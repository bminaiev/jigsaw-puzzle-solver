use std::{
    cmp::{max, min},
    collections::{hash_map::DefaultHasher, BTreeMap, BTreeSet, BinaryHeap, VecDeque},
    hash::{Hash, Hasher},
    i16::MAX,
    process,
};

use eframe::epaint::{Color32, ColorImage};
use itertools::Itertools;
use ndarray::Array4;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    border_matcher::is_picture_border,
    borders_graph::Graph,
    figure::BorderFigure,
    interactive_solutions_picker::InteractiveSolutionPicker,
    known_facts::{self, EdgeState, Fact, KnownFacts},
    known_positions::get_known_placement,
    parsed_puzzles::ParsedPuzzles,
    placement::{Placement, PotentialGroupLocation, Search3StateWithScore},
    point::{Point, PointF},
    positions_cache::PositionsCache,
    rects_fitter::get_bounding_box,
    search_states_cache::SearchStatesCache,
    surface_placer::{place_one_connected_component, rotate_component},
    topn::TopN,
    utils::{fmax, fmin, normalize_bounding_box, Side},
};

#[derive(Clone, Debug)]
struct Two {
    s0: Side,
    s1: Side,
    dist: f64,
}

#[derive(Clone, PartialEq, Debug)]
struct SearchState {
    all_edges: Vec<(Side, Side)>,
    vertices: BTreeSet<usize>,
    score: f64,
    av_edge_dist: f64,
}

impl SearchState {
    pub fn new(mut all_edges: Vec<(Side, Side)>, av_edge_dist: f64, bb: (i32, i32)) -> Self {
        all_edges.sort();
        all_edges.dedup();
        let mut vertices = BTreeSet::default();
        for (s1, s2) in all_edges.iter() {
            vertices.insert(s1.fig);
            vertices.insert(s2.fig);
        }
        let score = (all_edges.len() as f64) / (max(1, vertices.len()) as f64);
        Self {
            all_edges,
            vertices,
            score,
            av_edge_dist,
        }
    }
}

impl PartialOrd for SearchState {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for SearchState {}

impl Ord for SearchState {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score
            .total_cmp(&other.score)
            .then(self.av_edge_dist.total_cmp(&other.av_edge_dist).reverse())
            .then(self.all_edges.cmp(&other.all_edges))
    }
}

fn optimize_start_position(
    mut placement: Placement,
    n: usize,
    mut dist: impl FnMut(Side, Side) -> f64,
) -> Placement {
    let mut used = vec![false; n];
    for v in placement.get_all_used_figures() {
        used[v] = true;
    }
    loop {
        let mut changed = false;
        eprintln!("New loop!");
        for potential_location in placement.get_potential_locations(true) {
            let cur_figure = placement
                .get_figure_by_position(potential_location.pos)
                .unwrap();
            let mut cur_score = 0.0;
            for i in 0..4 {
                if let Some(another_side) = potential_location.neighbors[i] {
                    let cur_dist = dist(another_side, cur_figure[i]);
                    cur_score = fmax(cur_score, cur_dist);
                }
            }
            for fig in 0..n {
                if changed {
                    break;
                }
                if !used[fig] {
                    for offset in 0..4 {
                        let mut max_dist = 0.0;
                        let mut example_side_pair = None;
                        for i in 0..4 {
                            if let Some(existing_side) = potential_location.neighbors[i] {
                                let my_side = Side {
                                    fig,
                                    side: (i + offset) % 4,
                                };
                                let dist = dist(my_side, existing_side);
                                max_dist = fmax(max_dist, dist);
                                example_side_pair = Some((existing_side, my_side));
                            }
                        }
                        if max_dist < cur_score {
                            let old_figure = cur_figure[0].fig;
                            eprintln!(
                                "Replaced {old_figure} with {fig}: {cur_score} -> {max_dist}.Location = {:?}", potential_location.pos
                            );
                            changed = true;
                            placement.remove_figure(cur_figure[0].fig);
                            let (s1, s2) = example_side_pair.unwrap();
                            used[s2.fig] = true;
                            used[cur_figure[0].fig] = false;
                            if let Some(new_edges) = placement.join_sides(s1, s2) {
                                for &(s1, s2) in new_edges.iter() {
                                    eprintln!("New edge: {}", dist(s1, s2));
                                }
                            } else {
                                unreachable!();
                            }
                            break;
                        }
                    }
                }
            }
            if changed {
                break;
            }
        }
        if !changed {
            break;
        }
    }
    placement
}

#[derive(Clone, Debug)]
struct Edge {
    side1: Side,
    side2: Side,
    dist: f64,
}

#[derive(Clone, PartialEq)]
pub struct PlacedFigure {
    pub figure_id: usize,
    pub positions: Vec<PointF>,
}

#[derive(Clone, PartialEq)]
pub struct PotentialSolution {
    pub placed_figures: Vec<PlacedFigure>,
    pub placement_score: f64,
    pub text_offset: PointF,
    pub additional_text: String,
    pub debug_lines: Vec<Vec<PointF>>,
    pub bbox: (PointF, PointF),
    pub new_figures_used: Vec<usize>,
    pub new_edges_used: Vec<(Side, Side)>,
}

impl PartialOrd for PotentialSolution {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for PotentialSolution {}

impl Ord for PotentialSolution {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.placement_score.total_cmp(&other.placement_score)
    }
}

fn is_inside(border: &[PointF], p: PointF) -> bool {
    let mut cnt_intersections = 0;
    for (p1, p2) in border.iter().circular_tuple_windows() {
        let (p_min, p_max) = if p1.y < p2.y { (p1, p2) } else { (p2, p1) };
        if PointF::vect_mul(p_min, p_max, &p) <= 0.0 {
            if p.y >= p_min.y && p.y < p_max.y {
                cnt_intersections += 1;
            }
        }
    }
    cnt_intersections % 2 == 1
}

fn get_pixels_inside_figure(border: &[PointF]) -> Vec<Point> {
    let (p0, p1) = get_bounding_box(border);
    let mut res = vec![];
    for x in p0.x as usize..p1.x as usize {
        for y in p0.y as usize..p1.y as usize {
            if is_inside(
                border,
                PointF {
                    x: x as f64,
                    y: y as f64,
                },
            ) {
                res.push(Point { x, y });
            }
        }
    }
    res
}

impl PotentialSolution {
    pub fn gen_image(ps: &[Self]) -> ColorImage {
        eprintln!("Start generating solutions mask");
        let mut max_x = 0;
        let mut max_y = 0;
        for sol in ps.iter() {
            max_x = max(max_x, sol.bbox.1.x as usize + 2);
            max_y = max(max_y, sol.bbox.1.y as usize + 2);
        }

        let mut res = ColorImage::new([max_x, max_y], Color32::TRANSPARENT);

        for sol in ps.iter() {
            let to_color: Vec<_> = sol
                .placed_figures
                .par_iter()
                .map(|fig| {
                    if sol.new_figures_used.contains(&fig.figure_id) {
                        return vec![];
                    }
                    get_pixels_inside_figure(&fig.positions)
                })
                .flatten()
                .collect();
            for p in to_color.iter() {
                res[(p.x, p.y)] = Color32::LIGHT_GRAY;
            }
        }
        eprintln!("Image generated!");

        res
    }

    pub fn refresh(&self, known_facts: &KnownFacts) -> Option<Self> {
        let all_placed_vertices = known_facts.get_all_placed_vertices();
        let new_figures_used = self
            .new_figures_used
            .iter()
            .filter(|&fig| !all_placed_vertices.contains(fig))
            .cloned()
            .collect_vec();
        if new_figures_used.is_empty() {
            return None;
        }
        let mut new_edges_used = vec![];
        for &(s1, s2) in self.new_edges_used.iter() {
            match known_facts.get_edge_state(s1, s2) {
                EdgeState::Unknown => new_edges_used.push((s1, s2)),
                EdgeState::GoodEdge => {}
                EdgeState::WrongEdge => return None,
            }
        }
        if new_edges_used.is_empty() {
            return None;
        }
        let mut res = self.clone();
        res.new_figures_used = new_figures_used;
        res.new_edges_used = new_edges_used;
        Some(res)
    }
}

pub fn find_sides_by_known_edge(
    fig1: usize,
    fig2: usize,
    mut dist: impl FnMut(Side, Side) -> f64,
) -> (Side, Side) {
    let mut best_score = f64::MAX;
    let mut best_sides = None;
    for side1 in 0..4 {
        for side2 in 0..4 {
            let s1 = Side {
                fig: fig1,
                side: side1,
            };
            let s2 = Side {
                fig: fig2,
                side: side2,
            };
            let score = dist(s1, s2);
            if score < best_score {
                best_score = score;
                best_sides = Some((s1, s2));
            }
        }
    }
    best_sides.unwrap()
}

pub fn solve_graph(
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
    prev_state: Option<Graph>,
) -> Vec<PotentialSolution> {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    eprintln!("Hello there!");
    let n = graph.n;
    let base_dist = graph.gen_adj_matrix();
    eprintln!("nd array created!");

    let base_points_matrix = graph.get_base_points_matrix();

    let all_sides = (0..n)
        .cartesian_product(0..4)
        .map(|(fig, side)| Side { fig, side })
        .collect_vec();

    let base_dist = |s1: Side, s2: Side| -> f64 { base_dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    let sorted_by_dist = calc_sorted_by_dist(parsed_puzzles, base_dist);
    let mut multipliers = vec![vec![1.0; 4]; n];
    for fig in 0..n {
        for side in 0..4 {
            if sorted_by_dist[fig][side].is_empty() {
                continue;
            }
            let another_side = sorted_by_dist[fig][side][1];
            let cost = 1.0 / base_dist(Side { fig, side }, another_side);
            multipliers[fig][side] = cost;
        }
    }

    let mut dist = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 50.0);
    for fig1 in 0..n {
        for side1 in 0..4 {
            for fig2 in 0..n {
                for side2 in 0..4 {
                    dist[[fig1, side1, fig2, side2]] = base_dist(
                        Side {
                            fig: fig1,
                            side: side1,
                        },
                        Side {
                            fig: fig2,
                            side: side2,
                        },
                    ) * multipliers[fig1][side1]
                        * multipliers[fig2][side2];
                }
            }
        }
    }
    const OFFSET: f64 = 0.1;
    for fig in 0..n {
        for side in 0..4 {
            for (i, another_side) in sorted_by_dist[fig][side].iter().enumerate() {
                dist[[fig, side, another_side.fig, another_side.side]] += OFFSET * (i as f64);
                dist[[another_side.fig, another_side.side, fig, side]] += OFFSET * (i as f64);
            }
        }
    }
    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    let mut all_edges = vec![];
    for fig1 in 0..n {
        for side1 in 0..4 {
            for fig2 in fig1 + 1..n {
                for side2 in 0..4 {
                    let s1 = Side {
                        fig: fig1,
                        side: side1,
                    };
                    let s2 = Side {
                        fig: fig2,
                        side: side2,
                    };
                    all_edges.push(Edge {
                        side1: s1,
                        side2: s2,
                        dist: dist(s1, s2),
                    });
                }
            }
        }
    }
    all_edges.sort_by(|e1, e2| e1.dist.total_cmp(&e2.dist));
    for e in all_edges[..10].iter() {
        eprintln!("{:?}", e);
    }

    const MAX_DIST: f64 = 35.8;
    const MAX_DIST_MULIPLIERS: [f64; 5] = [0.0, 1.0, 1.1, 1.3, 2.0];

    let mut twos = vec![];
    for &s0 in &all_sides {
        for &s1 in &all_sides {
            let d0 = dist(s0, s1);
            if d0 > MAX_DIST * *MAX_DIST_MULIPLIERS.last().unwrap() || s0.fig >= s1.fig {
                continue;
            }
            twos.push(Two { s0, s1, dist: d0 });
        }
    }

    eprintln!("All edges len: {}", twos.len());
    twos.sort_by(|t1, t2| dist(t1.s0, t1.s1).total_cmp(&dist(t2.s0, t2.s1)));

    let gen_placement = |state: &SearchState| -> Placement {
        let mut placement = Placement::new();
        for &(s1, s2) in state.all_edges.iter() {
            placement.join_sides(s1, s2).unwrap();
        }
        placement
    };

    let gen_state = |new_placement: &Placement| -> SearchState {
        let new_state_edges = new_placement.get_all_neighbours();
        let sum_dists: f64 = new_state_edges.iter().map(|&(s1, s2)| dist(s1, s2)).sum();
        let bb = new_placement.get_bounding_box();
        let cnt_edges = new_state_edges.len();
        SearchState::new(new_state_edges, sum_dists / (cnt_edges as f64), bb)
    };

    const MAX_CNT: usize = 5000;
    let mut pq = vec![TopN::new(MAX_CNT); graph.n + 1];
    // pq[0].push(SearchState::new(vec![], 0.0));

    let start_vertices_num = {
        let mut start_state_edges = vec![];
        if let Some(prev_state) = prev_state {
            for e in prev_state.all_edges.iter() {
                eprintln!("WANT TO USE EDGE FROM PREV SOL: {:?}", e);
                start_state_edges.push((
                    Side {
                        fig: e.fig1,
                        side: e.side1,
                    },
                    Side {
                        fig: e.fig2,
                        side: e.side2,
                    },
                ))
            }
        } else {
            for e in graph.all_edges.iter() {
                if e.existing_edge {
                    eprintln!("WANT TO USE EDGE: {:?}", e);
                    start_state_edges.push((
                        Side {
                            fig: e.fig1,
                            side: e.side1,
                        },
                        Side {
                            fig: e.fig2,
                            side: e.side2,
                        },
                    ))
                }
            }
        }
        {
            const START_VERTEX: usize = 999;
            let placement = get_known_placement(graph);
            start_state_edges.extend(placement.get_all_neighbours_in_same_component(START_VERTEX));
        }
        if start_state_edges.is_empty() {
            start_state_edges.push((twos[0].s0, twos[0].s1));
        }
        let new_state = SearchState::new(start_state_edges, 0.0, (1, 1));
        let placement = gen_placement(&new_state);
        // let new_state = gen_state(&optimize_start_position(placement, n, dist));
        let new_state = gen_state(&placement);

        let start_vertices = new_state.vertices.len();
        eprintln!("start vertices: {start_vertices}");
        pq[start_vertices].insert(new_state);
        start_vertices
    };

    // let banned_figures = [
    //     190, 162, 568, 291, 132, 606, 745, 267, 179, 931, 989, 583, 751, 864, 934, 929, 822, 442,
    //     973, 427, 238, 947, 887, 783, 1035, 673, 890, 498, 324, 490,
    // ];

    let banned_figures: [usize; 0] = [];
    // let banned_figures = [716, 912, 928, 114, 656, 719];
    // let banned_figures = [
    //     864, 267, 783, 190, 275, 179, 427, 238, 890, 1035, 275, 583, 783, 673, 934, 929, 864, 989,
    //     931, 568, 190, 291, 162, 745, 132, 606,
    // ];

    let max_v = start_vertices_num + 1;
    for cnt_vertices in start_vertices_num..min(max_v, pq.len()) {
        let mut iter = 0;

        let cur_pq = std::mem::replace(&mut pq[cnt_vertices], TopN::new(0));

        for state in cur_pq.set.iter().rev() {
            let cur_placement = gen_placement(&state);
            let cur_bb = normalize_bounding_box(cur_placement.get_bounding_box());
            if iter < 30 {
                eprintln!(
                    "iter = {}. num vertices = {}, score = {}, av edge = {}, bbox : {:?}",
                    iter,
                    state.vertices.len(),
                    state.score,
                    state.av_edge_dist,
                    cur_bb
                );
            }
            iter += 1;
            let mut new_placement = cur_placement.clone();
            let mut can_use_figure = vec![true; n];
            for i in 0..n {
                if !parsed_puzzles.figures[i].is_good_puzzle() {
                    can_use_figure[i] = false;
                }
            }
            for &fig in banned_figures.iter() {
                can_use_figure[fig] = false;
            }
            for fig in cur_placement.get_all_used_figures() {
                can_use_figure[fig] = false;
            }
            for potential_location in cur_placement.get_potential_locations(false) {
                let mut ok_bb_increase = true;
                for fig in 0..n {
                    if !can_use_figure[fig] {
                        continue;
                    }
                    if !ok_bb_increase {
                        break;
                    }
                    for offset in 0..4 {
                        let mut max_dist = 0.0;
                        let mut cnt_connections = 0;
                        let mut example_side_pair = None;
                        for i in 0..4 {
                            if let Some(existing_side) = potential_location.neighbors[i] {
                                cnt_connections += 1;
                                let my_side = Side {
                                    fig,
                                    side: (i + offset) % 4,
                                };
                                let dist = dist(my_side, existing_side);
                                max_dist = fmax(max_dist, dist);
                                example_side_pair = Some((existing_side, my_side));
                            }
                        }
                        let max_dist_with_multiplier =
                            MAX_DIST * MAX_DIST_MULIPLIERS[cnt_connections];
                        if max_dist <= max_dist_with_multiplier {
                            let (s0, s1) = example_side_pair.unwrap();
                            if let Some(_new_edges) = new_placement.join_sides(s0, s1) {
                                let new_bounding_box =
                                    normalize_bounding_box(new_placement.get_bounding_box());
                                let mut ok = true;
                                if is_bad_bounding_box(new_bounding_box) {
                                    ok = false;
                                }
                                if new_bounding_box != cur_bb
                                    && !ok_to_increase_bb(cur_bb, cnt_vertices)
                                {
                                    ok = false;
                                }

                                if ok {
                                    let new_state = gen_state(&new_placement);
                                    // eprintln!(
                                    //     "With figure {fig} new score is {}",
                                    //     new_state.av_edge_dist
                                    // );
                                    pq[new_state.vertices.len()].insert(new_state);
                                } else {
                                    ok_bb_increase = false;
                                }
                                new_placement = cur_placement.clone();
                            }
                        }
                    }
                }
            }
        }
    }

    let mut solutions = vec![];
    let mut all = pq[max_v].set.iter().cloned().collect_vec();
    all.reverse();
    all.truncate(100);
    for state in all {
        eprintln!(
            "Try state... score = {}, cur_comp = {:?}",
            state.score, state.vertices
        );

        let cur_component = state.vertices.iter().cloned().collect_vec();
        let used_edges = state.all_edges.clone();
        let mut positions = vec![None; n];
        let placement_score = place_one_connected_component(
            parsed_puzzles,
            &cur_component,
            &used_edges,
            &mut positions,
            &base_points_matrix,
        )
        .unwrap_or(f64::MAX);
        eprintln!("placed score: {placement_score}");
        rotate_component(
            &cur_component,
            &mut positions,
            graph,
            parsed_puzzles,
            &[],
            &[],
        );

        let mut placed_figures = vec![];
        for (figure_id, positions) in positions.iter().enumerate() {
            if let Some(positions) = positions {
                placed_figures.push(PlacedFigure {
                    figure_id,
                    positions: positions.clone(),
                });
            }
        }

        solutions.push(PotentialSolution {
            placed_figures,
            placement_score,
            text_offset: PointF { x: 0.0, y: 0.0 },
            additional_text: "".to_owned(),
            debug_lines: vec![],
            bbox: (PointF::ZERO, PointF::ZERO),
            new_edges_used: vec![],
            new_figures_used: vec![],
        });
    }
    solutions.sort_by(|s1, s2| s1.placement_score.total_cmp(&s2.placement_score));
    solutions

    // let placement = gen_placement(&last_state.unwrap());
    // graph.get_subgraph(&placement)
}

fn is_bad_bounding_box((min, max): (i32, i32)) -> bool {
    min * 3 < max
}

fn ok_to_increase_bb((x, y): (i32, i32), cnt: usize) -> bool {
    let area = (x * y) as usize;
    let ok_from = area - area / 10;
    cnt >= ok_from
}

fn calc_sorted_by_dist(
    parsed_puzzles: &ParsedPuzzles,
    mut dist: impl FnMut(Side, Side) -> f64,
) -> Vec<Vec<Vec<Side>>> {
    let n = parsed_puzzles.figures.len();
    let mut sorted_by_dist = vec![vec![vec![]; 4]; n];
    for fig in 0..n {
        if !parsed_puzzles.figures[fig].is_good_puzzle() {
            continue;
        }

        for side in 0..4 {
            for fig2 in 0..n {
                if fig2 == fig || !parsed_puzzles.figures[fig2].is_good_puzzle() {
                    continue;
                }
                for side2 in 0..4 {
                    let another_side = Side {
                        fig: fig2,
                        side: side2,
                    };
                    if dist(Side { fig, side }, another_side) < 100500.0 {
                        sorted_by_dist[fig][side].push(another_side);
                    }
                }
            }
            let my_side = Side { fig, side };
            sorted_by_dist[fig][side]
                .sort_by(|s1, s2| dist(my_side, *s1).total_cmp(&dist(my_side, *s2)));
        }
    }
    sorted_by_dist
}

pub fn solve_graph_border(graph: &Graph, parsed_puzzles: &ParsedPuzzles) -> Graph {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    eprintln!("Hello there!");
    let n = graph.n;
    let dist = graph.gen_adj_matrix();
    eprintln!("nd array created!");

    let mut figures_on_border = parsed_puzzles.calc_figures_on_border();

    let sz = figures_on_border.len();
    eprintln!("Number of figures on border: {}", sz);

    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    let sorted_by_dist = calc_sorted_by_dist(parsed_puzzles, dist);

    const MAX_DIST: f64 = 4.5;

    let mut dist2 = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 50.0);
    for (iter, left_fig) in figures_on_border.iter().enumerate() {
        eprintln!("dist2 iter: {iter}/{}", figures_on_border.len());
        for right_fig in figures_on_border.iter() {
            if left_fig.figure_id == right_fig.figure_id {
                continue;
            }
            let s1 = left_fig.right_side;
            let s1_up = s1.pr();
            let s2 = right_fig.left_side;
            let s2_up = s2.ne();
            let mut cur_res = f64::MAX;

            let start_dist = dist(s1, s2);

            if start_dist > MAX_DIST {
                continue;
            }

            let list1 = &sorted_by_dist[s1.fig][s1_up.side];
            let list2 = &sorted_by_dist[s2.fig][s2_up.side];
            for sum_ix in 0..list1.len() + list2.len() - 1 {
                {
                    let ix1 = sum_ix / 2;
                    let ix2 = sum_ix - ix1;
                    let s3 = list1[ix1];
                    let s4 = list2[ix2];
                    if dist(s1_up, s3) > cur_res && dist(s2_up, s4) > cur_res {
                        // eprintln!(
                        //     "Break at idx: {sum_ix}. Dists: {} {}",
                        //     dist(s1_up, s3),
                        //     dist(s2_up, s4)
                        // );
                        break;
                    }
                }
                for ix1 in 0..min(sum_ix + 1, list1.len()) {
                    let ix2 = sum_ix - ix1;
                    if ix2 < list2.len() {
                        let s3 = list1[ix1];
                        let s4 = list2[ix2];
                        let d1_up = dist(s1_up, s3);
                        let d2_up = dist(s2_up, s4);
                        if d1_up > cur_res || d2_up > cur_res {
                            continue;
                        }
                        let check_res = fmax(start_dist, dist(s3.pr(), s4.ne()));
                        let check_res = fmax(fmax(d1_up, d2_up), check_res);
                        cur_res = fmin(cur_res, check_res);
                    }
                }
            }

            // if start_dist < 2.0 {
            //     eprintln!(
            //         "Calculated dist2: {:?}..{:?}: {start_dist} -> {cur_res}",
            //         left_fig, right_fig
            //     );
            // }
            dist2[[s1.fig, s1.side, s2.fig, s2.side]] = cur_res;
            dist2[[s2.fig, s2.side, s1.fig, s1.side]] = cur_res;
        }
    }

    let dist2 = |s1: Side, s2: Side| -> f64 { dist2[[s1.fig, s1.side, s2.fig, s2.side]] };

    let mut rnd = StdRng::seed_from_u64(123);

    let segm_cost = |border: &[BorderFigure],
                     mut i1: usize,
                     mut i2: usize,
                     mut i3: usize,
                     mut i4: usize|
     -> f64 {
        i1 %= sz;
        i2 %= sz;
        i3 %= sz;
        i4 %= sz;
        dist2(border[i1].right_side, border[i2].left_side)
            + dist2(border[i3].right_side, border[i4].left_side)
    };

    for iter in 0..100_000_000 {
        let start1 = rnd.gen_range(0..sz - 1);
        let len1 = rnd.gen_range(1..sz - 3);

        let offset = rnd.gen_range(1..sz - len1 - 2);
        let start2 = (start1 + len1 + offset) % sz;
        let len2 = rnd.gen_range(1..sz - len1 - offset);

        assert!(len1 + offset + len2 < sz);

        let start_cost = segm_cost(
            &figures_on_border,
            start1 + sz - 1,
            start1,
            start1 + len1 - 1,
            start1 + len1,
        ) + segm_cost(
            &figures_on_border,
            start2 + sz - 1,
            start2,
            start2 + len2 - 1,
            start2 + len2,
        );

        let end_cost = segm_cost(
            &figures_on_border,
            start1 + sz - 1,
            start2,
            start2 + len2 - 1,
            start1 + len1,
        ) + segm_cost(
            &figures_on_border,
            start2 + sz - 1,
            start1,
            start1 + len1 - 1,
            start2 + len2,
        );

        if end_cost < start_cost {
            eprintln!("Optimize (iter = {iter}): {start_cost} -> {end_cost}");

            figures_on_border.rotate_left(start1);
            let start2 = len1 + offset;

            let next_order = figures_on_border[start2..start2 + len2]
                .iter()
                .chain(figures_on_border[len1..start2].iter())
                .chain(figures_on_border[..len1].iter())
                .chain(figures_on_border[start2 + len2..].iter())
                .cloned()
                .collect_vec();
            assert_eq!(figures_on_border.len(), next_order.len());
            figures_on_border = next_order;
        }
    }

    let gen_placement = |order: &[BorderFigure]| -> Placement {
        let mut placement = Placement::new();
        for (f1, f2) in order.iter().circular_tuple_windows() {
            let d = dist(f1.right_side, f2.left_side);
            if d > MAX_DIST {
                continue;
            }
            eprintln!("Use score: {d}");
            if placement.join_sides(f1.right_side, f2.left_side).is_none() {
                eprintln!("Can't join sides");
            }
        }
        placement
    };

    let placement = gen_placement(&figures_on_border);
    graph.get_subgraph(&placement)
    //
}

fn gen_relative_dists(graph: &Graph, parsed_puzzles: &ParsedPuzzles) -> Array4<f64> {
    let n = graph.n;

    let base_dist = graph.gen_adj_matrix();
    let base_dist = |s1: Side, s2: Side| -> f64 { base_dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    let sorted_by_dist = calc_sorted_by_dist(parsed_puzzles, base_dist);
    let mut multipliers = vec![vec![1.0; 4]; n];
    for fig in 0..n {
        for side in 0..4 {
            if sorted_by_dist[fig][side].is_empty() {
                continue;
            }
            let another_side = sorted_by_dist[fig][side][1];
            let cost = 1.0 / base_dist(Side { fig, side }, another_side);
            multipliers[fig][side] = cost;
        }
    }

    let mut dist = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 50.0);
    for fig1 in 0..n {
        for side1 in 0..4 {
            for fig2 in 0..n {
                for side2 in 0..4 {
                    dist[[fig1, side1, fig2, side2]] = base_dist(
                        Side {
                            fig: fig1,
                            side: side1,
                        },
                        Side {
                            fig: fig2,
                            side: side2,
                        },
                    ) * multipliers[fig1][side1]
                        * multipliers[fig2][side2];
                }
            }
        }
    }
    const OFFSET: f64 = 0.1;
    for fig in 0..n {
        for side in 0..4 {
            for (i, another_side) in sorted_by_dist[fig][side].iter().enumerate() {
                dist[[fig, side, another_side.fig, another_side.side]] += OFFSET * (i as f64);
                dist[[another_side.fig, another_side.side, fig, side]] += OFFSET * (i as f64);
            }
        }
    }
    dist
}

fn find_best_next_rec(
    state: &PotentialGroupLocation,
    first_sides: &mut Vec<Side>,
    ways: &mut Vec<Search3StateWithScore>,
    dist: impl FnMut(Side, Side) -> f64 + Clone,
    to_check: &[Vec<Vec<Side>>],
    to_check_limit: usize,
    used: &[bool],
) {
    if state.locations.len() == first_sides.len() {
        ways.push(state.check_way(&first_sides, dist));
    } else {
        let cur_loc = first_sides.len();
        let mut good_sides = vec![];
        for offset in 0..4 {
            if let Some(other_side) = state.locations[cur_loc].neighbors[offset] {
                let to_check = &to_check[other_side.fig][other_side.side];
                let to_check = &to_check[..min(to_check_limit, to_check.len())];
                for my_side in to_check.iter() {
                    if used[my_side.fig] || first_sides.iter().any(|s| s.fig == my_side.fig) {
                        continue;
                    }
                    good_sides.push(my_side.ne_offset(4 - offset));
                }
            }
        }
        good_sides.sort();
        good_sides.dedup();
        for side in good_sides.into_iter() {
            first_sides.push(side);
            find_best_next_rec(
                state,
                first_sides,
                ways,
                dist.clone(),
                to_check,
                to_check_limit,
                used,
            );
            first_sides.pop();
        }
    }
}

fn find_best_next(
    state: &PotentialGroupLocation,
    limit_inside: usize,
    limit_res: usize,
    to_check: &[Vec<Vec<Side>>],
    used: &[bool],
    dist: impl FnMut(Side, Side) -> f64 + Clone,
) -> Vec<Search3StateWithScore> {
    let mut next = vec![];
    find_best_next_rec(
        state,
        &mut vec![],
        &mut next,
        dist,
        to_check,
        limit_inside,
        used,
    );
    next.sort();
    next.truncate(limit_res);
    while !next.is_empty() && next.last().unwrap().score > 10000.0 {
        // eprintln!("Remove last with score: {}", next.last().unwrap().score);
        next.pop();
    }
    next
}

const START_VERTEX: usize = 628;

pub fn solve_graph_add_by_3(
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
    _prev_state: Option<Graph>,
    known_facts: &mut KnownFacts,
) -> InteractiveSolutionPicker {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    eprintln!("Hello there!");
    let n = graph.n;
    eprintln!("nd array created!");
    let base_points_matrix = graph.get_base_points_matrix();
    let positions_cache = PositionsCache::load(&parsed_puzzles);

    let mut on_border = vec![false; n];
    for v in 0..n {
        if !parsed_puzzles.figures[v].is_good_puzzle() {
            continue;
        }
        for side in 0..4 {
            if is_picture_border(&parsed_puzzles.figures[v], side) {
                on_border[v] = true;
            }
        }
    }

    let all_sides = parsed_puzzles.gen_all_sides();

    let mut dist = gen_relative_dists(graph, parsed_puzzles);
    const MX: f64 = f64::MAX / 100000.0;
    for fact in known_facts.facts.iter() {
        if !fact.good_edge {
            dist[[
                fact.side1.fig,
                fact.side1.side,
                fact.side2.fig,
                fact.side2.side,
            ]] = MX;
            dist[[
                fact.side2.fig,
                fact.side2.side,
                fact.side1.fig,
                fact.side1.side,
            ]] = MX;
        }
    }
    let mut start_placement = Placement::new();
    for fact in known_facts.facts.iter() {
        if fact.good_edge {
            start_placement.join_sides(fact.side1, fact.side2).unwrap();
        }
    }
    for (s1, s2) in start_placement.get_all_neighbours() {
        for &s3 in all_sides.iter() {
            if s2 != s3 {
                dist[[s1.fig, s1.side, s3.fig, s3.side]] = MX;
            }
        }
    }
    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    let rot_positions;

    let states;
    let mut used = vec![false; graph.n];
    {
        let mut placement = Placement::new();
        for fact in known_facts.facts.iter() {
            if fact.good_edge {
                placement.join_sides(fact.side1, fact.side2).unwrap();
            }
        }
        for v in placement.get_all_used_figures() {
            used[v] = true;
        }
        let my_comp_placement = placement.get_only_one_component_placement(START_VERTEX);
        let edges = my_comp_placement.get_all_neighbours();
        rot_positions = get_correct_rotation_positions(
            &edges,
            parsed_puzzles,
            graph,
            &positions_cache,
            &base_points_matrix,
        );

        states = my_comp_placement.get_potential_group_locations(3);

        eprintln!("Start states: {}", states.len());
    }

    const CHECK_BEST: usize = 500;
    const TOO_BIG_COST: f64 = 100.0;
    let mut to_check = calc_sorted_by_dist(parsed_puzzles, dist);
    for fig in 0..to_check.len() {
        for side in 0..4 {
            to_check[fig][side].truncate(CHECK_BEST);
            while let Some(s2) = to_check[fig][side].last() {
                if dist(Side { fig, side }, *s2) > TOO_BIG_COST {
                    to_check[fig][side].pop();
                } else {
                    break;
                }
            }
        }
    }

    const LIMIT: usize = 250;
    const LIMIT_RES: usize = 50;
    let next_states: Vec<_> = states
        .par_iter()
        .map(|st| find_best_next(&st, LIMIT, LIMIT_RES, &to_check, &used, dist))
        .collect();
    let mut next_states = next_states.into_iter().flatten().collect_vec();
    next_states.sort();
    eprintln!("Generated {} start states.", next_states.len());
    let states_cache = SearchStatesCache::load();
    let next_states: Vec<_> = next_states
        .into_iter()
        .filter(|state| !states_cache.contains(state.get_hash(), 5.0))
        .collect_vec();
    eprintln!("After filtering: {} states", next_states.len());

    let more_res: Vec<_> = next_states
        .par_iter()
        .filter_map(|state| {
            let used_edges = state.all_edges();
            Some((
                state.clone(),
                gen_potential_solution(
                    &used_edges,
                    parsed_puzzles,
                    graph,
                    None,
                    known_facts,
                    &rot_positions,
                    &positions_cache,
                    &base_points_matrix,
                )?,
            ))
        })
        .collect();
    eprintln!("Generated All solutions.");
    for (a, b) in more_res.iter() {
        states_cache.insert(a.get_hash(), b.placement_score);
    }
    states_cache.save();
    InteractiveSolutionPicker::new(
        more_res,
        START_VERTEX,
        rot_positions,
        positions_cache,
        base_points_matrix,
        known_facts,
        parsed_puzzles,
        graph,
    )
}

fn get_correct_rotation_positions(
    edges: &[(Side, Side)],
    parsed_puzzles: &ParsedPuzzles,
    graph: &Graph,
    positions_cache: &PositionsCache,
    base_points_matrix: &Array4<[PointF; 2]>,
) -> Vec<Option<Vec<PointF>>> {
    let mut cur_component = edges
        .iter()
        .map(|(s1, s2)| [s1.fig, s2.fig])
        .flatten()
        .collect_vec();
    cur_component.sort();
    cur_component.dedup();

    let mut positions = vec![None; parsed_puzzles.figures.len()];
    let placement_score = positions_cache
        .place_one_connected_component_with_cache(
            parsed_puzzles,
            &cur_component,
            edges,
            &mut positions,
            base_points_matrix,
        )
        .unwrap_or(f64::MAX);
    eprintln!("placed score: {placement_score}");
    rotate_component(
        &cur_component,
        &mut positions,
        graph,
        parsed_puzzles,
        &[START_VERTEX],
        &[],
    );
    positions
}

pub fn gen_potential_solution(
    edges: &[(Side, Side)],
    parsed_puzzles: &ParsedPuzzles,
    graph: &Graph,
    set_fixed_score: Option<f64>,
    known_facts: &KnownFacts,
    rot_positions: &[Option<Vec<PointF>>],
    positions_cache: &PositionsCache,
    base_points_matrix: &Array4<[PointF; 2]>,
) -> Option<PotentialSolution> {
    let mut cur_component = edges
        .iter()
        .map(|(s1, s2)| [s1.fig, s2.fig])
        .flatten()
        .collect_vec();
    cur_component.sort();
    cur_component.dedup();

    let mut positions = vec![None; parsed_puzzles.figures.len()];
    let placement_score = positions_cache.place_one_connected_component_with_cache(
        parsed_puzzles,
        &cur_component,
        edges,
        &mut positions,
        base_points_matrix,
    )?;
    eprintln!("placed score: {placement_score}");
    rotate_component(
        &cur_component,
        &mut positions,
        graph,
        parsed_puzzles,
        &[],
        rot_positions,
    );

    let mut placed_figures = vec![];
    for (figure_id, positions) in positions.iter().enumerate() {
        if let Some(positions) = positions {
            placed_figures.push(PlacedFigure {
                figure_id,
                positions: positions.clone(),
            });
        }
    }

    let placed_vertices = known_facts.get_all_placed_vertices();
    let new_figures_used = cur_component
        .iter()
        .filter(|x| !placed_vertices.contains(x))
        .cloned()
        .collect_vec();
    let new_edges_used = edges
        .iter()
        .filter(|&&(s1, s2)| {
            let known = known_facts.facts.contains(&Fact::new(s1, s2, true));
            let from_old = placed_vertices.contains(&s1.fig) != placed_vertices.contains(&s2.fig);
            !known && from_old
        })
        .cloned()
        .collect_vec();

    Some(PotentialSolution {
        placed_figures,
        placement_score: set_fixed_score.unwrap_or(placement_score),
        text_offset: PointF { x: 0.0, y: 0.0 },
        additional_text: "".to_owned(),
        debug_lines: vec![],
        bbox: (PointF::ZERO, PointF::ZERO),
        new_figures_used,
        new_edges_used,
    })
}
