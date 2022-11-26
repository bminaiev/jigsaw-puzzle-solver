use std::collections::{BTreeMap, BTreeSet, VecDeque};

use eframe::egui::plot::Corner;
use itertools::Itertools;

use crate::{
    border_matcher::{
        local_optimize_coordinate_systems, match_borders, match_placed_borders, MatchResult,
    },
    borders_graph::Graph,
    coordinate_system::CoordinateSystem,
    dsu::Dsu,
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    point::{Point, PointF},
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

fn place_one_connected_component(
    parsed_puzzles: &ParsedPuzzles,
    component: &[usize],
    used_edges: &[(Side, Side)],
    positions: &mut Vec<Option<Vec<PointF>>>,
) {
    let mut matched_borders = BTreeMap::new();

    for &(s1, s2) in used_edges.iter() {
        let match_res = match_borders(
            &parsed_puzzles.figures[s1.fig],
            s1.side,
            &parsed_puzzles.figures[s2.fig],
            s2.side,
            s1.fig,
            s2.fig,
        )
        .unwrap();

        matched_borders.insert((s1, s2), match_res);
    }

    eprintln!("Built all borders default matchers.");

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

        for &(s1, s2) in used_edges.iter() {
            if to_cs[s1.fig].is_some() && to_cs[s2.fig].is_none() {
                changed = true;

                let match_res = &matched_borders[&(s1, s2)];

                eprintln!("score = {}", match_res.score);

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
            for &(s1, s2) in used_edges.iter() {
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

    eprintln!("Do local optimizations!");
    local_optimize_positions(used_edges, positions, parsed_puzzles, component)
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
        eprintln!("unite {} {}", s1.fig, s2.fig);
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
