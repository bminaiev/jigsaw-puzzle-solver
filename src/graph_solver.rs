use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, BinaryHeap, VecDeque},
    i16::MAX,
    process,
};

use itertools::Itertools;
use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::{
    border_matcher::is_picture_border, borders_graph::Graph, parsed_puzzles::ParsedPuzzles,
    placement::Placement, topn::TopN, utils::Side,
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
        let max_bb_edge = max(bb.0, bb.1) as f64;
        let min_bb_edge = min(bb.0, bb.1) as f64;
        let bb_coef = min_bb_edge + 2.0 * max_bb_edge;
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

pub fn solve_graph(graph: &Graph, parsed_puzzles: &ParsedPuzzles) -> Graph {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    eprintln!("Hello there!");
    let n = graph.n;
    let dist = graph.gen_adj_matrix();
    eprintln!("nd array created!");

    let all_sides = (0..n)
        .cartesian_product(0..4)
        .map(|(fig, side)| Side { fig, side })
        .collect_vec();

    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    const MAX_DIST: f64 = 1.0;
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

    let gen_placement = |state: &SearchState| -> Placement {
        let mut placement = Placement::new();
        for &(s1, s2) in state.all_edges.iter() {
            placement.join_sides(s1, s2).unwrap();
        }
        placement
    };

    let gen_state =
        |new_state_edges: Vec<(Side, Side)>, new_placement: &Placement| -> SearchState {
            let sum_dists: f64 = new_state_edges.iter().map(|&(s1, s2)| dist(s1, s2)).sum();
            let bb = new_placement.get_bounding_box();
            let cnt_edges = new_state_edges.len();
            SearchState::new(new_state_edges, sum_dists / (cnt_edges as f64), bb)
        };

    const MAX_CNT: usize = 30;
    let mut pq = vec![TopN::new(MAX_CNT); graph.n + 1];
    // pq[0].push(SearchState::new(vec![], 0.0));

    let start_vertices_num = {
        let mut start_state_edges = vec![];
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
        let new_state = SearchState::new(start_state_edges, 0.0, (1, 1));
        let start_vertices = new_state.vertices.len();
        eprintln!("start vertices: {start_vertices}");
        pq[start_vertices].insert(new_state);
        start_vertices
    };

    let max_v = start_vertices_num + 50;
    let mut last_state = None;
    for cnt_vertices in start_vertices_num..min(max_v, pq.len()) {
        let mut iter = 0;
        let mut gg = 0;

        let cur_pq = std::mem::replace(&mut pq[cnt_vertices], TopN::new(0));

        if let Some(cur_best) = cur_pq.get_best() {
            last_state = Some(cur_best);
        }

        for state in cur_pq.set.iter().rev() {
            let cur_placement = gen_placement(&state);
            let bbox = cur_placement.get_bounding_box();
            if iter < 30 {
                eprintln!(
            "iter = {}. num vertices = {}, score = {}, av edge = {}, bbox : {:?}, gg = {gg}",
            iter,
            state.vertices.len(),
            state.score,
            state.av_edge_dist,
            bbox
        );
            }
            iter += 1;
            let mut new_placement = cur_placement.clone();
            for two in twos.iter() {
                let alr1 = state.vertices.contains(&two.s0.fig);
                let alr2 = state.vertices.contains(&two.s1.fig);
                if alr1 && alr2 {
                    continue;
                }
                if !alr1 && !alr2 && !state.vertices.is_empty() {
                    continue;
                }
                gg += 1;
                if let Some(new_edges) = new_placement.join_sides(two.s0, two.s1) {
                    let mut ok = true;
                    let mut new_state_edges = state.all_edges.clone();
                    let max_dist_with_multiplier = MAX_DIST * MAX_DIST_MULIPLIERS[new_edges.len()];
                    for &(s1, s2) in new_edges.iter() {
                        let d = dist(s1, s2);
                        if d > max_dist_with_multiplier {
                            ok = false;
                            break;
                        }
                        new_state_edges.push((s1, s2));
                    }

                    if ok {
                        let new_state = gen_state(new_state_edges, &new_placement);
                        pq[new_state.vertices.len()].insert(new_state);
                    }
                    new_placement = cur_placement.clone();
                }
            }
        }
    }

    let placement = gen_placement(&last_state.unwrap());
    graph.get_subgraph(&placement)
    //
}

#[derive(Clone, Copy)]
struct BorderFigure {
    figure_id: usize,
    left_side: Side,
    right_side: Side,
}

impl BorderFigure {
    pub fn new(figure_id: usize, left_side: usize, right_side: usize) -> Self {
        Self {
            figure_id,
            left_side: Side {
                fig: figure_id,
                side: (left_side % 4),
            },
            right_side: Side {
                fig: figure_id,
                side: (right_side % 4),
            },
        }
    }
}

pub fn solve_graph_border(graph: &Graph, parsed_puzzles: &ParsedPuzzles) -> Graph {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    eprintln!("Hello there!");
    let n = graph.n;
    let dist = graph.gen_adj_matrix();
    eprintln!("nd array created!");

    let mut figures_on_border = (0..n)
        .filter_map(|id| {
            let figure = &parsed_puzzles.figures[id];
            if !figure.is_good_puzzle() {
                return None;
            }
            let is_border = (0..4)
                .map(|border_id| is_picture_border(figure, border_id))
                .collect_vec();
            for i in 0..4 {
                if is_border[i] && is_border[(i + 1) % 4] {
                    return Some(BorderFigure::new(id, i + 2, i + 3));
                }
            }
            for i in 0..4 {
                if is_border[i] {
                    return Some(BorderFigure::new(id, i + 1, i + 3));
                }
            }
            None
        })
        .collect_vec();

    let sz = figures_on_border.len();
    eprintln!("Number of figures on border: {}", sz);

    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

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
        dist(border[i1].right_side, border[i2].left_side)
            + dist(border[i3].right_side, border[i4].left_side)
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

    const MAX_DIST: f64 = 1.5;

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
