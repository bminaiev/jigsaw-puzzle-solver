use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, BinaryHeap, VecDeque},
};

use itertools::Itertools;

use crate::{
    border_matcher::{local_optimize_coordinate_systems, match_borders, match_placed_borders},
    borders_graph::Graph,
    coordinate_system::CoordinateSystem,
    dsu::Dsu,
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
    placement::{self, Placement},
    point::{Point, PointF},
    rects_fitter::{self, RectsFitter},
    surface_placer::place_on_surface,
    topn::TopN,
    utils::{fmax, Side},
};

#[derive(Clone, Debug)]
struct Four {
    s0: Side,
    s1: Side,
    s2: Side,
    s3: Side,
    max_dist: f64,
}

impl Four {
    pub fn all(&self) -> [Side; 4] {
        [self.s0, self.s1, self.s2, self.s3]
    }

    pub fn neighbours(&self) -> [(Side, Side); 4] {
        [
            (self.s0, self.s1.ne2()),
            (self.s0.ne(), self.s2.pr()),
            (self.s1.ne(), self.s3.pr()),
            (self.s2, self.s3.ne2()),
        ]
    }
}

#[derive(Clone, Debug)]
struct Two {
    s0: Side,
    s1: Side,
    dist: f64,
}

fn fmax4(a: [f64; 4]) -> f64 {
    let mut res = a[0];
    for &x in a[1..].iter() {
        if x > res {
            res = x;
        }
    }
    res
}



pub fn solve_graph_simple(
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
) -> Vec<Option<Vec<PointF>>> {
    assert_eq!(graph.parsed_puzzles_hash, parsed_puzzles.calc_hash());

    let mut edges = graph.all_edges.clone();
    edges.sort_by(|e1, e2| e1.score.total_cmp(&e2.score));
    let mut placement = Placement::new();
    for e in edges.iter() {
        if e.score > 2.0 {
            break;
        }

        if placement
            .join_sides(
                Side {
                    fig: e.fig1,
                    side: e.side1,
                },
                Side {
                    fig: e.fig2,
                    side: e.side2,
                },
            )
            .is_some()
        {
            eprintln!("Join sides! {}", e.score);
        }
    }
    place_on_surface(graph, &placement, parsed_puzzles)
}

pub fn solve_graph_fours(
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
) -> Vec<Option<Vec<PointF>>> {
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

    const MAX_DIST: f64 = 3.5;
    const MAX_NEW_EDGE_DIST: f64 = 7.0;
    let mut fours = vec![];
    let mut twos = vec![];
    for &s0 in &all_sides {
        for &s1 in &all_sides {
            let d0 = dist(s0, s1.ne2());
            if d0 > MAX_DIST {
                continue;
            }
            twos.push(Two { s0, s1, dist: d0 });
            for &s2 in &all_sides {
                let d1 = dist(s0.ne(), s2.pr());
                if d1 > MAX_DIST {
                    continue;
                }
                for &s3 in &all_sides {
                    if s0.fig == s3.fig {
                        continue;
                    }
                    let d2 = dist(s1.ne(), s3.pr());
                    let d3 = dist(s2, s3.ne2());
                    if d2 > MAX_DIST || d3 > MAX_DIST {
                        continue;
                    }
                    fours.push(Four {
                        s0,
                        s1,
                        s2,
                        s3,
                        max_dist: fmax4([d0, d1, d2, d3]),
                    });
                }
            }
        }
    }
    fours.sort_by(|f1, f2| f1.max_dist.total_cmp(&f2.max_dist));
    twos.sort_by(|t1, t2| t1.dist.total_cmp(&t2.dist));
    eprintln!("total fours: {}, twos: {}", fours.len(), twos.len());

    let mut placement = Placement::new();

    let mut used = vec![false; graph.n];

    for (it, four) in fours.iter().enumerate() {
        if it % 1000 == 0 {
            eprintln!("fours: {}/{}", it, fours.len());
        }
        let mut cnt_used = 0;
        for s in four.all() {
            if used[s.fig] {
                cnt_used += 1;
            }
        }
        if cnt_used == 4 {
            continue;
        }
        let mut new_placement = placement.clone();
        let mut ok = true;
        for (s1, s2) in four.neighbours() {
            if let Some(edges) = new_placement.join_sides(s1, s2) {
                if edges
                    .iter()
                    .any(|&(s1, s2)| dist(s1, s2) > MAX_NEW_EDGE_DIST)
                {
                    ok = false;
                    break;
                }
            } else {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        placement = new_placement;
        for s in four.all() {
            used[s.fig] = true;
        }
    }
    eprintln!("finished fours..");
    const ADD_TWOS: bool = false;
    if ADD_TWOS {
        for two in twos.iter() {
            if placement.get_fig_comp_size(two.s0.fig) == 1
                || placement.get_fig_comp_size(two.s1.fig) == 1
            {
                continue;
            }
            let old_placement = placement.clone();
            if let Some(edges) = placement.join_sides(two.s0, two.s1) {
                if edges
                    .iter()
                    .any(|&(s1, s2)| dist(s1, s2) > MAX_NEW_EDGE_DIST)
                {
                    eprintln!("will not join, roll back!");
                    placement = old_placement;
                } else {
                    eprintln!("joined pair!");
                }
            }
        }
    }

    place_on_surface(graph, &placement, parsed_puzzles)
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
        let max_bb_edge = max(bb.0, bb.1);
        let min_bb_edge = min(bb.0, bb.1);
        let score = (all_edges.len() as f64)
            / (max(1, vertices.len()) as f64)
            / (max_bb_edge as f64)
            / (min_bb_edge as f64).powf(0.2);
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

pub fn solve_graph_queue_from_fours(
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
) -> Vec<Option<Vec<PointF>>> {
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

    const MAX_DIST: f64 = 3.5;

    let mut twos = vec![];
    for &s0 in &all_sides {
        for &s1 in &all_sides {
            // TODO: ne2???? REALLYYY?
            let d0 = dist(s0, s1.ne2());
            if d0 > MAX_DIST || s0.fig >= s1.fig {
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

    const MAX_CNT: usize = 2000;
    let mut pq = vec![TopN::new(MAX_CNT); graph.n + 1];
    // pq[0].push(SearchState::new(vec![], 0.0));

    for &s0 in &all_sides {
        for &s1 in &all_sides {
            let d0 = dist(s0, s1.ne2());
            if d0 > MAX_DIST {
                continue;
            }
            // TODO: why do I put twos twice?
            twos.push(Two { s0, s1, dist: d0 });
            for &s2 in &all_sides {
                let d1 = dist(s0.ne(), s2.pr());
                if d1 > MAX_DIST {
                    continue;
                }
                if s1.fig == s2.fig {
                    continue;
                }
                for &s3 in &all_sides {
                    if s0.fig == s3.fig {
                        continue;
                    }
                    let d2 = dist(s1.ne(), s3.pr());
                    let d3 = dist(s2, s3.ne2());
                    if d2 > MAX_DIST || d3 > MAX_DIST {
                        continue;
                    }

                    let four = Four {
                        s0,
                        s1,
                        s2,
                        s3,
                        max_dist: fmax4([d0, d1, d2, d3]),
                    };

                    let new_state_edges = four.neighbours().into_iter().collect_vec();
                    let max_edge_dist = new_state_edges
                        .iter()
                        .map(|&(s1, s2)| dist(s1, s2))
                        .max_by(f64::total_cmp)
                        .unwrap();
                    let new_state = SearchState::new(new_state_edges, max_edge_dist, (1, 1));
                    pq[new_state.vertices.len()].insert(new_state);
                }
            }
        }
    }

    const MAX_V: usize = 40;
    for cnt_vertices in 0..min(MAX_V, pq.len()) {
        let mut iter = 0;
        let mut gg = 0;

        let cur_pq = std::mem::replace(&mut pq[cnt_vertices], TopN::new(0));

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
                    for &(s1, s2) in new_edges.iter() {
                        let d = dist(s1, s2);
                        if d > MAX_DIST {
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

    let state = pq[MAX_V].set.iter().next().unwrap();
    let placement = gen_placement(state);
    place_on_surface(graph, &placement, parsed_puzzles)
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

    const MAX_DIST: f64 = 4.5;

    let mut twos = vec![];
    for &s0 in &all_sides {
        for &s1 in &all_sides {
            let d0 = dist(s0, s1);
            if d0 > MAX_DIST || s0.fig >= s1.fig {
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

    const MAX_CNT: usize = 1000;
    let mut pq = vec![TopN::new(MAX_CNT); graph.n + 1];
    // pq[0].push(SearchState::new(vec![], 0.0));

    {
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
        let max_edge_dist = start_state_edges
            .iter()
            .map(|&(s1, s2)| dist(s1, s2))
            .max_by(f64::total_cmp)
            .unwrap();
        let new_state = SearchState::new(start_state_edges, max_edge_dist, (1, 1));
        eprintln!("start vertices: {}", new_state.vertices.len());
        pq[new_state.vertices.len()].insert(new_state);
    }

    const MAX_V: usize = 25;
    for cnt_vertices in 0..min(MAX_V, pq.len()) {
        let mut iter = 0;
        let mut gg = 0;

        let cur_pq = std::mem::replace(&mut pq[cnt_vertices], TopN::new(0));

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
                    for &(s1, s2) in new_edges.iter() {
                        let d = dist(s1, s2);
                        if d > MAX_DIST {
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

    let state = pq[MAX_V].set.iter().next().unwrap();
    let placement = gen_placement(state);
    graph.get_subgraph(&placement)
    //
}
