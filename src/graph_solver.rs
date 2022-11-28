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
        let max_bb_edge = max(bb.0, bb.1);
        let min_bb_edge = min(bb.0, bb.1);
        let score = (all_edges.len() as f64)
            / (max(1, vertices.len()) as f64)
            / (max_bb_edge as f64).powf(2.0)
            / (min_bb_edge as f64).powf(1.0);
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

    const MAX_DIST: f64 = 1.5;
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

    const MAX_CNT: usize = 2000;
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
        let new_state = SearchState::new(start_state_edges, 0.0, (1, 1));
        eprintln!("start vertices: {}", new_state.vertices.len());
        pq[new_state.vertices.len()].insert(new_state);
    }

    const MAX_V: usize = 60;
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

    let state = pq[MAX_V].set.iter().next().unwrap();
    let placement = gen_placement(state);
    graph.get_subgraph(&placement)
    //
}
