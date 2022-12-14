use std::collections::BTreeMap;

use itertools::Itertools;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    border_matcher::{match_borders, match_side_borders_v2},
    borders_graph::Graph,
    graph_solver::{find_sides_by_known_edge, PlacedFigure, PotentialSolution},
    known_positions::get_known_placement,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    point::PointF,
    surface_placer::{get_border_and_neighbors, place_one_connected_component, rotate_component},
    utils::{dedup_edges, Side},
};

#[derive(Default)]
struct ScoringFunctionQuality {
    sum_log: f64,
    cnt: usize,
}

impl ScoringFunctionQuality {
    pub fn add_score(&mut self, position: usize) {
        assert!(position >= 1);
        self.sum_log += (position as f64).log2();
        self.cnt += 1;
    }

    pub fn get_quality(&self) -> f64 {
        if self.cnt == 0 {
            0.0
        } else {
            self.sum_log / (self.cnt as f64)
        }
    }
}

pub fn optimize_edge_scores(
    parsed_puzzles: &ParsedPuzzles,
    graph: &Graph,
    calc_new_scores: bool,
) -> Vec<PotentialSolution> {
    let dist = graph.gen_adj_matrix();
    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };
    let placement = get_known_placement(graph);
    let known_sides = dedup_edges(&placement.get_all_neighbours());

    let all_sides: Vec<Side> = parsed_puzzles.gen_all_sides();

    let mut places = BTreeMap::new();

    let mut existing_scoring_fun = ScoringFunctionQuality::default();
    let mut new_scoring_fun = ScoringFunctionQuality::default();
    for &(s1, s2) in known_sides.iter() {
        let cur_dist = dist(s1, s2);
        for &stay_edge in [s1, s2].iter() {
            let better = all_sides
                .iter()
                .filter(|another_side| dist(**another_side, stay_edge) <= cur_dist)
                .count();
            let other_side = if stay_edge == s1 { s2 } else { s1 };
            places.insert((stay_edge, other_side), (better, None));
            existing_scoring_fun.add_score(better);
            if calc_new_scores {
                // eprintln!("Checking {:?} {:?}. Start dist: {cur_dist}", s1, s2);
                // eprintln!("THIS IS FAIL.");
                // eprintln!("For {:?}, on position {better}", stay_edge);
                let new_my_score = match_borders(parsed_puzzles, stay_edge, other_side)
                    .unwrap()
                    .score;
                let new_better = all_sides
                    .par_iter()
                    .filter(|another_side| {
                        match_borders(parsed_puzzles, stay_edge, **another_side)
                            .map(|mr| mr.score)
                            .unwrap_or(f64::MAX)
                            <= new_my_score
                    })
                    .count();
                new_scoring_fun.add_score(new_better);
                eprintln!("({}, {}): {better} -> {new_better}", s1.fig, s2.fig);
                places.insert((stay_edge, other_side), (better, Some(new_better)));
            }
        }
    }
    eprintln!(
        "Scoring functions quality: {} {}",
        existing_scoring_fun.get_quality(),
        new_scoring_fun.get_quality()
    );

    let mut solutions = vec![];

    for &(s1, s2) in known_sides.iter() {
        let mut positions = vec![None; graph.n];
        let cur_component = [s1.fig, s2.fig];
        let used_edges = [(s1, s2)];
        let placement_score = place_one_connected_component(
            parsed_puzzles,
            &cur_component,
            &used_edges,
            &mut positions,
        );
        rotate_component(&cur_component, &mut positions, graph, parsed_puzzles);

        let mut placed_figures = vec![];
        for (figure_id, positions) in positions.iter().enumerate() {
            if let Some(positions) = positions {
                placed_figures.push(PlacedFigure {
                    figure_id,
                    positions: positions.clone(),
                });
            }
        }

        let (better, new_better) = places[&(s1, s2)];
        let (better2, new_better2) = places[&(s2, s1)];
        let mut additional_text = if let Some(new_better) = new_better {
            format!(
                " (place {better}/{better2} -> {new_better}/{})",
                new_better2.unwrap()
            )
        } else {
            format!(" (place {better}/{better2})")
        };

        let mut debug_lines = vec![];
        {
            let pos1 = positions[s1.fig].as_ref().unwrap();
            let pos2 = positions[s2.fig].as_ref().unwrap();

            let border1 = get_border_and_neighbors(pos1, parsed_puzzles, s1);
            let border2 = get_border_and_neighbors(pos2, parsed_puzzles, s2);
            let border2 = border2.reverse();

            {
                let res1 = match_side_borders_v2(&border1.prev, &border2.prev, true);
                debug_lines.push(res1.1);
                let res2 = match_side_borders_v2(&border1.next, &border2.next, true);
                debug_lines.push(res2.1);

                additional_text += &format!(". borders: {:.3}", res1.0 + res2.0);
            }
        }

        solutions.push(PotentialSolution {
            placed_figures,
            placement_score,
            text_offset: PointF { x: 0.0, y: 0.0 },
            additional_text,
            debug_lines,
        });
    }
    solutions.sort_by(|s1, s2| s1.placement_score.total_cmp(&s2.placement_score));
    solutions
}
