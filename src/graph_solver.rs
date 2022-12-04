use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, BinaryHeap, VecDeque},
    i16::MAX,
    process,
};

use itertools::Itertools;
use ndarray::Array4;
use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::{
    border_matcher::is_picture_border,
    borders_graph::Graph,
    figure::BorderFigure,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    topn::TopN,
    utils::{fmax, fmin, Side},
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

    const MAX_DIST: f64 = 2.0;
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

    const MAX_CNT: usize = 300;
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

    let max_v = start_vertices_num + 30;
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
            let cur_bb = normalize_bounding_box(cur_placement.get_bounding_box());
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

                    let new_bounding_box = normalize_bounding_box(new_placement.get_bounding_box());
                    if is_bad_bounding_box(new_bounding_box) {
                        ok = false;
                    }
                    if new_bounding_box != cur_bb && !ok_to_increase_bb(cur_bb, cnt_vertices) {
                        ok = false;
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

fn normalize_bounding_box((x, y): (i32, i32)) -> (i32, i32) {
    (min(x, y), max(x, y))
}

fn is_bad_bounding_box((min, max): (i32, i32)) -> bool {
    min * 3 < max
}

fn ok_to_increase_bb((x, y): (i32, i32), cnt: usize) -> bool {
    let area = (x * y) as usize;
    let ok_from = area - area / 10;
    cnt >= ok_from
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
                    sorted_by_dist[fig][side].push(Side {
                        fig: fig2,
                        side: side2,
                    });
                }
            }
            let my_side = Side { fig, side };
            sorted_by_dist[fig][side]
                .sort_by(|s1, s2| dist(my_side, *s1).total_cmp(&dist(my_side, *s2)));
        }
    }

    const MAX_DIST: f64 = 3.5;

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
