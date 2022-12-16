use std::{
    cmp::min,
    collections::{BTreeMap, BTreeSet, VecDeque},
    f64::consts::PI,
};

use eframe::egui::plot::Corner;
use itertools::Itertools;

use crate::{
    border_matcher::{
        local_optimize_coordinate_systems, match_borders, match_placed_borders, BorderAndNeighbors,
        MatchResult,
    },
    borders_graph::Graph,
    coordinate_system::CoordinateSystem,
    dsu::Dsu,
    figure::Figure,
    graph_solver::PotentialSolution,
    parsed_puzzles::ParsedPuzzles,
    placement::{Placement, PotentialLocation},
    point::{Point, PointF},
    rects_fitter::{get_bounding_box, RectsFitter},
    utils::{dedup_edges, fmax, fmin, Side},
};

fn gen_basic_position(figure: &Figure) -> Vec<PointF> {
    figure.border.iter().map(|p| p.conv_f64()).collect_vec()
}

#[derive(Clone)]
struct MoveCS {
    from_cs: CoordinateSystem,
    to_cs: CoordinateSystem,
}

fn from_cs_from_figure(figure: &Figure) -> CoordinateSystem {
    let (p1, p2) = figure.get_cs_points();
    assert!(p1 != p2);
    CoordinateSystem::new(p1.conv_f64(), p2.conv_f64() - p1.conv_f64())
}

impl MoveCS {
    pub fn from_from_to(from_cs: CoordinateSystem, to_cs: CoordinateSystem) -> Self {
        Self { from_cs, to_cs }
    }

    pub fn new(figure: &Figure, positions: &[PointF]) -> Self {
        let n = figure.border.len() / 2;
        let from_cs = from_cs_from_figure(figure);
        assert!(figure.border[0] != figure.border[n]);
        assert_eq!(figure.border.len(), positions.len());
        Self {
            from_cs,
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
        if self.all_scores.is_empty() {
            return 0.0;
        }
        self.all_scores.iter().sum::<f64>() / self.all_scores.len() as f64
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

pub fn get_border_and_neighbors(
    pos: &[PointF],
    parsed_puzzles: &ParsedPuzzles,
    side: Side,
) -> BorderAndNeighbors {
    let mut prev = get_border(pos, parsed_puzzles, side.pr());
    {
        let sz = min(BorderAndNeighbors::NEIGHBORS, prev.len());
        prev.rotate_right(sz);
        prev.truncate(sz);
        prev.reverse();
    }
    let mut next = get_border(pos, parsed_puzzles, side.ne());
    next.truncate(BorderAndNeighbors::NEIGHBORS);

    BorderAndNeighbors {
        border: get_border(pos, parsed_puzzles, side),
        prev,
        next,
    }
}

fn get_borders_dist(
    pos1: &[PointF],
    pos2: &[PointF],
    parsed_puzzles: &ParsedPuzzles,
    s1: Side,
    s2: Side,
) -> f64 {
    let border1 = get_border_and_neighbors(pos1, parsed_puzzles, s1);
    let border2 = get_border_and_neighbors(pos2, parsed_puzzles, s2);
    let border2 = border2.reverse();
    match_placed_borders(&border1, &border2)
}

fn local_optimize_positions(
    all_edges: &[(Side, Side)],
    positions: &mut [Option<Vec<PointF>>],
    parsed_puzzles: &ParsedPuzzles,
    cur_component: &[usize],
) -> f64 {
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
        })
        .collect_vec();

    let mut final_score = 0.0;

    for glob_iter in 0..1 {
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
        final_score = new_score.get_final_score();
        if new_score.get_final_score() < start_score.get_final_score() {
            let new_pos = calc_new_positions(&to_cs);
            for local_id in 0..cur_component.len() {
                positions[cur_component[local_id]] = Some(new_pos[local_id].clone());
            }
        } else {
            break;
        }
    }
    final_score
}

pub fn rotate_component(
    component: &[usize],
    positions: &mut Vec<Option<Vec<PointF>>>,
    graph: &Graph,
    parsed_puzzles: &ParsedPuzzles,
    base: &[usize],
    rot_positions: &[Option<Vec<PointF>>],
) {
    let mut all_points: Vec<PointF> = vec![];
    for &c in component.iter() {
        all_points.extend(&positions[c].as_ref().unwrap().clone());
    }
    let mut probably_correct_dir = graph.get_puzzles_with_probably_correct_directions();
    for &v in base.iter() {
        probably_correct_dir[v] = true;
    }
    for i in 0..rot_positions.len() {
        if rot_positions[i].is_some() {
            probably_correct_dir[i] = true;
        }
    }
    {
        if !component.iter().any(|&c| probably_correct_dir[c]) {
            probably_correct_dir[component[0]] = true;
        }
    }
    let calc_score = |angle: f64| -> f64 {
        let mut res = 0.0;

        for &c in component.iter() {
            if !probably_correct_dir[c] {
                continue;
            }
            let (p1, p2) = parsed_puzzles.figures[c].get_cs_points();
            let (mut p1, mut p2) = (p1.conv_f64(), p2.conv_f64());
            let (i1, i2) = parsed_puzzles.figures[c].get_cs_points_indexes();

            if rot_positions.len() > c && rot_positions[c].is_some() {
                p1 = rot_positions[c].as_ref().unwrap()[i1];
                p2 = rot_positions[c].as_ref().unwrap()[i2];
            }

            let p3 = positions[c].as_ref().unwrap()[i1].rotate(angle);
            let p4 = positions[c].as_ref().unwrap()[i2].rotate(angle);
            let dir1 = (p2 - p1).norm();
            let dir2 = (p4 - p3).norm();
            res += dir1.x * dir2.x + dir1.y * dir2.y;
        }
        res
    };
    let mut best_score = -f64::MAX;
    let mut best_angle = 0.0;
    const STEPS: usize = 1000;
    for angle in (0..STEPS).map(|step| (step as f64) * PI * 2.0 / (STEPS as f64)) {
        let score = calc_score(angle);
        if score > best_score {
            best_score = score;
            best_angle = angle;
        }
    }
    for &c in component.iter() {
        for p in positions[c].as_mut().unwrap().iter_mut() {
            *p = p.rotate(best_angle);
        }
    }
}

pub fn place_one_connected_component(
    parsed_puzzles: &ParsedPuzzles,
    component: &[usize],
    used_edges: &[(Side, Side)],
    positions: &mut Vec<Option<Vec<PointF>>>,
) -> f64 {
    let mut matched_borders = BTreeMap::new();

    let used_edges = dedup_edges(used_edges);
    let mut used_edges_two_sides = vec![];
    for &(s1, s2) in used_edges.iter() {
        used_edges_two_sides.push((s1, s2));
        used_edges_two_sides.push((s2, s1));
    }

    for &(s1, s2) in used_edges_two_sides.iter() {
        let match_res = match_borders(parsed_puzzles, s1, s2).unwrap();

        matched_borders.insert((s1, s2), match_res);
    }

    let root = component[0];
    positions[root] = Some(gen_basic_position(&parsed_puzzles.figures[root]));

    let from_cs = parsed_puzzles
        .figures
        .iter()
        .map(|f| from_cs_from_figure(f))
        .collect_vec();
    let mut to_cs = vec![None; positions.len()];
    to_cs[root] = Some(from_cs[root].clone());

    let predict_point_based_on_edge =
        |from_v: usize, to_cs: &[Option<CoordinateSystem>], p: PointF| {
            let move_cs = MoveCS::from_from_to(
                from_cs[from_v].clone(),
                to_cs[from_v].as_ref().unwrap().clone(),
            );
            move_cs.conv_point(p)
        };

    let predict_to_cs_based_on_edge =
        |from_v: usize,
         to_cs: &[Option<CoordinateSystem>],
         to_v: usize,
         match_result: &MatchResult| {
            let (i1, i2) = parsed_puzzles.figures[to_v].get_cs_points_indexes();
            let p1_placed = predict_point_based_on_edge(from_v, to_cs, match_result.rhs[i1]);
            let p2_placed = predict_point_based_on_edge(from_v, to_cs, match_result.rhs[i2]);
            CoordinateSystem::new(p1_placed, p2_placed - p1_placed)
        };

    loop {
        let mut changed = false;

        for &(s1, s2) in used_edges_two_sides.iter() {
            if to_cs[s1.fig].is_some() && to_cs[s2.fig].is_none() {
                changed = true;

                let match_res = &matched_borders[&(s1, s2)];

                to_cs[s2.fig] = Some(predict_to_cs_based_on_edge(
                    s1.fig, &to_cs, s2.fig, match_res,
                ));
            }
        }

        if !changed {
            break;
        }
    }

    {
        let mut new_start = vec![PointF::ZERO; positions.len()];
        let mut new_x_dir = vec![PointF::ZERO; positions.len()];
        let mut cnt_edges = vec![0; positions.len()];
        for _iter in 0..1000 {
            for &c in component.iter() {
                new_start[c] = PointF::ZERO;
                new_x_dir[c] = PointF::ZERO;
                cnt_edges[c] = 0;
            }
            for &(s1, s2) in used_edges_two_sides.iter() {
                let match_res = &matched_borders[&(s1, s2)];

                let new_cs = predict_to_cs_based_on_edge(s1.fig, &to_cs, s2.fig, match_res);

                new_start[s2.fig] = new_start[s2.fig] + new_cs.start;
                new_x_dir[s2.fig] = new_x_dir[s2.fig] + new_cs.x_dir;
                cnt_edges[s2.fig] += 1;
            }
            for &c in component.iter() {
                let cs = CoordinateSystem::new(
                    new_start[c] / (cnt_edges[c] as f64),
                    new_x_dir[c] / (cnt_edges[c] as f64),
                );
                to_cs[c] = Some(cs);
            }
        }
    }

    for &v in component.iter() {
        assert!(to_cs[v].is_some(), "no position for {v}");
        let move_cs = MoveCS::from_from_to(from_cs[v].clone(), to_cs[v].as_ref().unwrap().clone());
        positions[v] = Some(
            parsed_puzzles.figures[v]
                .border
                .iter()
                .map(|&p| move_cs.conv_point(p.conv_f64()))
                .collect(),
        );
    }

    local_optimize_positions(&used_edges, positions, parsed_puzzles, component)
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
    for &(s1, s2) in all_edges.iter() {
        dsu.unite(s1.fig, s2.fig);
    }
    let mut rects_fitter = RectsFitter::new();
    let mut positions: Vec<Option<Vec<PointF>>> = vec![None; graph.n];

    for cur_component in dsu.get_components().into_iter() {
        if cur_component.len() <= 1 {
            continue;
        }
        let comp_set: BTreeSet<_> = cur_component.iter().collect();
        let used_edges = all_edges
            .iter()
            .filter(|(s1, s2)| comp_set.contains(&s1.fig) && comp_set.contains(&s2.fig))
            .cloned()
            .collect_vec();

        place_one_connected_component(parsed_puzzles, &cur_component, &used_edges, &mut positions);

        eprintln!("Rotate component!");
        rotate_component(
            &cur_component,
            &mut positions,
            graph,
            parsed_puzzles,
            &[],
            &[],
        );

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

pub fn put_solutions_on_surface(solutions: &mut [PotentialSolution]) {
    let mut rects_fitter = RectsFitter::new();
    for i in 0..solutions.len() {
        let new_pts: Vec<PointF> = solutions[i]
            .placed_figures
            .iter()
            .flat_map(|pf| pf.positions.iter())
            .cloned()
            .collect();
        let bbox = get_bounding_box(&new_pts);
        let shift = rects_fitter.add_points(&new_pts);
        for fig in solutions[i].placed_figures.iter_mut() {
            for p in fig.positions.iter_mut() {
                *p = *p + shift;
            }
        }
        for line in solutions[i].debug_lines.iter_mut() {
            for p in line.iter_mut() {
                *p = *p + shift;
            }
        }
        solutions[i].text_offset = bbox.0 + shift;
        solutions[i].bbox = (bbox.0 + shift, bbox.1 + shift)
    }
}
