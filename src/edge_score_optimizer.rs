use itertools::Itertools;

use crate::{
    border_matcher::match_borders, borders_graph::Graph, graph_solver::find_sides_by_known_edge,
    parsed_puzzles::ParsedPuzzles, utils::Side,
};

pub fn optimize_edge_scores(parsed_puzzles: &ParsedPuzzles, graph: &Graph) {
    let known_edges = vec![
        (471, 974),
        (417, 471),
        (713, 974),
        (871, 974),
        (713, 756),
        (756, 436),
        (436, 815),
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
    let known_sides = known_edges
        .iter()
        .map(|&(fig1, fig2)| find_sides_by_known_edge(fig1, fig2, dist))
        .collect_vec();
    let all_sides: Vec<Side> = parsed_puzzles.gen_all_sides();
    for &(s1, s2) in known_sides.iter() {
        let cur_dist = dist(s1, s2);
        for &stay_edge in [s1, s2].iter() {
            let better = all_sides
                .iter()
                .filter(|another_side| dist(**another_side, stay_edge) <= cur_dist)
                .count();
            {
                // eprintln!("Checking {:?} {:?}. Start dist: {cur_dist}", s1, s2);
                // eprintln!("THIS IS FAIL.");
                // eprintln!("For {:?}, on position {better}", stay_edge);
                let other_side = if stay_edge == s1 { s2 } else { s1 };
                let new_my_score = match_borders(parsed_puzzles, stay_edge, other_side)
                    .unwrap()
                    .score;
                let new_better = all_sides
                    .iter()
                    .filter(|another_side| {
                        match_borders(parsed_puzzles, stay_edge, **another_side)
                            .map(|mr| mr.score)
                            .unwrap_or(f64::MAX)
                            <= new_my_score
                    })
                    .count();
                eprintln!("({}, {}): {better} -> {new_better}", s1.fig, s2.fig);
            }
        }
    }
}
