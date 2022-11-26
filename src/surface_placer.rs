use std::collections::{BTreeMap, VecDeque};

use itertools::Itertools;

use crate::{
    border_matcher::{local_optimize_coordinate_systems, match_borders, match_placed_borders},
    borders_graph::Graph,
    coordinate_system::CoordinateSystem,
    dsu::Dsu,
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    point::PointF,
    rects_fitter::RectsFitter,
    utils::Side,
};

fn gen_basic_position(figure: &Figure) -> Vec<PointF> {
    figure.border.iter().map(|p| p.conv_f64()).collect_vec()
}

#[derive(Clone)]
struct MoveCS {
    from_cs: CoordinateSystem,
    to_cs: CoordinateSystem,
}

impl MoveCS {
    pub fn new(figure: &Figure, positions: &[PointF]) -> Self {
        let n = figure.border.len() / 2;
        assert!(figure.border[0] != figure.border[n]);
        assert_eq!(figure.border.len(), positions.len());
        Self {
            from_cs: CoordinateSystem::new(
                figure.border[0].conv_f64(),
                figure.border[n].conv_f64() - figure.border[0].conv_f64(),
            ),
            to_cs: CoordinateSystem::new(positions[0], positions[n] - positions[0]),
        }
    }

    pub fn conv_point(&self, p: PointF) -> PointF {
        self.to_cs.to_real(self.from_cs.create(p))
    }
}

#[derive(Debug)]
struct EdgesScores {
    all_scores: Vec<f64>,
}

impl EdgesScores {
    pub fn new() -> Self {
        Self { all_scores: vec![] }
    }

    pub fn add_score(&mut self, score: f64) {
        self.all_scores.push(score);
    }

    pub fn get_final_score(&self) -> f64 {
        self.all_scores.iter().sum()
    }
}

fn get_border(pos: &[PointF], parsed_puzzles: &ParsedPuzzles, side: Side) -> Vec<PointF> {
    let figure = &parsed_puzzles.figures[side.fig];
    let mut cur = figure.corner_positions[side.side];
    let end = figure.corner_positions[(side.side + 1) % figure.corner_positions.len()];
    let mut res = vec![];
    loop {
        res.push(pos[cur]);
        if cur == end {
            break;
        }
        cur = (cur + 1) % pos.len();
    }
    res
}

fn get_borders_dist(
    pos1: &[PointF],
    pos2: &[PointF],
    parsed_puzzles: &ParsedPuzzles,
    s1: Side,
    s2: Side,
) -> f64 {
    let border1 = get_border(pos1, parsed_puzzles, s1);
    let mut border2 = get_border(pos2, parsed_puzzles, s2);
    border2.reverse();
    match_placed_borders(&border1, &border2)
}

fn local_optimize_positions(
    all_edges: &[(Side, Side)],
    positions: &mut [Option<Vec<PointF>>],
    parsed_puzzles: &ParsedPuzzles,
    cur_component: &[usize],
) {
    let local_id_mapping: BTreeMap<usize, usize> = cur_component
        .iter()
        .enumerate()
        .map(|(local_id, &global_id)| (global_id, local_id))
        .collect();

    let all_edges = all_edges
        .iter()
        .cloned()
        .filter(|(s1, s2)| {
            positions[s1.fig].is_some()
                && positions[s2.fig].is_some()
                && local_id_mapping.contains_key(&s1.fig)
                && local_id_mapping.contains_key(&s2.fig)
                && s1.fig < s2.fig
        })
        .collect_vec();

    for glob_iter in 0..1 {
        eprintln!("global iter: {glob_iter}");

        let move_cs = cur_component
            .iter()
            .map(|&v| MoveCS::new(&parsed_puzzles.figures[v], &positions[v].as_ref().unwrap()))
            .collect_vec();

        let calc_new_positions = |to_cs: &[CoordinateSystem]| {
            (0..to_cs.len())
                .map(|local_id| {
                    parsed_puzzles.figures[cur_component[local_id]]
                        .border
                        .iter()
                        .map(|&p| {
                            to_cs[local_id].to_real(move_cs[local_id].from_cs.create(p.conv_f64()))
                        })
                        .collect_vec()
                })
                .collect_vec()
        };

        let calc_score = |to_cs: &[CoordinateSystem]| {
            let my_pos = calc_new_positions(to_cs);
            let mut res = EdgesScores::new();
            for &(s1, s2) in all_edges.iter() {
                let dist = get_borders_dist(
                    &my_pos[local_id_mapping[&s1.fig]],
                    &my_pos[local_id_mapping[&s2.fig]],
                    parsed_puzzles,
                    s1,
                    s2,
                );
                res.add_score(dist);
            }
            res
        };

        let to_cs = move_cs
            .iter()
            .map(|move_cs| move_cs.to_cs.clone())
            .collect_vec();
        let start_score = calc_score(&to_cs);
        let to_cs =
            local_optimize_coordinate_systems(&to_cs, |to_cs| calc_score(to_cs).get_final_score());
        let new_score = calc_score(&to_cs);
        if new_score.get_final_score() < start_score.get_final_score() {
            let new_pos = calc_new_positions(&to_cs);
            for local_id in 0..cur_component.len() {
                positions[cur_component[local_id]] = Some(new_pos[local_id].clone());
            }
        } else {
            break;
        }
    }
}

pub fn place_on_surface(
    graph: &Graph,
    placement: &Placement,
    parsed_puzzles: &ParsedPuzzles,
) -> Vec<Option<Vec<PointF>>> {
    assert_eq!(graph.n, parsed_puzzles.figures.len());
    eprintln!("Start placing on the surface!");
    let dists = graph.gen_adj_matrix();
    let mut all_edges = placement.get_all_neighbours();
    let dist = |s1: Side, s2: Side| dists[[s1.fig, s1.side, s2.fig, s2.side]];
    all_edges.sort_by(|&(s1, s2), &(s3, s4)| dist(s1, s2).total_cmp(&dist(s3, s4)));
    let mut dsu = Dsu::new(graph.n);
    let mut tree_edges = vec![vec![]; graph.n];
    for &(s1, s2) in all_edges.iter() {
        if dsu.unite(s1.fig, s2.fig) {
            tree_edges[s1.fig].push((s1, s2));
            tree_edges[s2.fig].push((s2, s1));
        }
    }
    let mut rects_fitter = RectsFitter::new();
    let mut positions: Vec<Option<Vec<PointF>>> = vec![None; graph.n];

    const DO_LOCAL_OPTIMIZATIONS_ON_EVERY_MOVE: bool = false;

    // TODO: smarter logic for root choosing
    for root in 0..positions.len() {
        if positions[root].is_some() || tree_edges[root].is_empty() {
            continue;
        }
        eprintln!("Start root: {}", root);
        // TODO: choose good starting position
        positions[root] = Some(gen_basic_position(&parsed_puzzles.figures[root]));
        let mut queue = VecDeque::new();
        queue.push_back(root);
        let mut cur_component = vec![root];
        while let Some(v) = queue.pop_front() {
            for (s1, s2) in tree_edges[v].iter() {
                if positions[s2.fig].is_some() {
                    continue;
                }
                cur_component.push(s2.fig);
                queue.push_back(s2.fig);

                let match_res = match_borders(
                    &parsed_puzzles.figures[s1.fig],
                    s1.side,
                    &parsed_puzzles.figures[s2.fig],
                    s2.side,
                    s1.fig,
                    s2.fig,
                )
                .unwrap();

                let move_cs = MoveCS::new(
                    &parsed_puzzles.figures[s1.fig],
                    &positions[s1.fig].as_ref().unwrap(),
                );

                eprintln!("score = {}", match_res.score);

                positions[s2.fig] = Some(
                    match_res
                        .rhs
                        .iter()
                        .map(|&p| move_cs.conv_point(p - match_res.shifted))
                        .collect(),
                );
                if DO_LOCAL_OPTIMIZATIONS_ON_EVERY_MOVE {
                    local_optimize_positions(
                        &all_edges,
                        &mut positions,
                        parsed_puzzles,
                        &cur_component,
                    );
                }
            }
        }
        if !DO_LOCAL_OPTIMIZATIONS_ON_EVERY_MOVE {
            // local_optimize_positions(&all_edges, &mut positions, parsed_puzzles, &cur_component);
        }
        {
            let new_pts = cur_component
                .iter()
                .map(|&id| positions[id].as_ref().unwrap())
                .cloned()
                .flatten()
                .collect_vec();
            let shift = rects_fitter.add_points(&new_pts);
            for &fig in cur_component.iter() {
                for p in positions[fig].as_mut().unwrap() {
                    *p = *p + shift
                }
            }
        }
    }
    positions
}
