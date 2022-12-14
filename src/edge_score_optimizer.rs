use std::collections::BTreeMap;

use itertools::Itertools;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    border_matcher::{match_borders, match_side_borders_v2},
    borders_graph::Graph,
    graph_solver::{find_sides_by_known_edge, PlacedFigure, PotentialSolution},
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    point::PointF,
    surface_placer::{get_border_and_neighbors, place_one_connected_component, rotate_component},
    utils::{dedup_edges, Side},
};

pub fn optimize_edge_scores(
    parsed_puzzles: &ParsedPuzzles,
    graph: &Graph,
    calc_new_scores: bool,
) -> Vec<PotentialSolution> {
    let known_edges = vec![
        (471, 974),
        (417, 471),
        (713, 974),
        (871, 974),
        (713, 756),
        // (756, 436),<- TODO: fix correct edges
        // (436, 815),
        (999, 926),
        (999, 574),
        (926, 637),
        (574, 963),
        (439, 570),
        (570, 548),
        (926, 421),
        (683, 421),
    ];

    let dist = graph.gen_adj_matrix();
    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };
    let mut known_sides = known_edges
        .iter()
        .map(|&(fig1, fig2)| find_sides_by_known_edge(fig1, fig2, dist))
        .collect_vec();

    known_sides.push((Side { fig: 756, side: 3 }, Side { fig: 436, side: 3 }));
    known_sides.push((Side { fig: 815, side: 3 }, Side { fig: 436, side: 1 }));

    let mut placement = Placement::new();
    for &(s1, s2) in known_sides.iter() {
        placement.join_sides(s1, s2).unwrap();
    }
    let known_sides = dedup_edges(&placement.get_all_neighbours());

    let all_sides: Vec<Side> = parsed_puzzles.gen_all_sides();

    let mut places = BTreeMap::new();

    for &(s1, s2) in known_sides.iter() {
        let cur_dist = dist(s1, s2);
        for &stay_edge in [s1, s2].iter() {
            let better = all_sides
                .iter()
                .filter(|another_side| dist(**another_side, stay_edge) <= cur_dist)
                .count();
            let other_side = if stay_edge == s1 { s2 } else { s1 };
            places.insert((stay_edge, other_side), (better, None));
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
                eprintln!("({}, {}): {better} -> {new_better}", s1.fig, s2.fig);
                places.insert((stay_edge, other_side), (better, Some(new_better)));
            }
        }
    }

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
