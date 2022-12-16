use std::cmp::{max, min};

use crate::{
    borders_graph::Graph,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    utils::{normalize_bounding_box, Side},
};

fn placement_by_field(field: &[Vec<i32>], mut dist: impl FnMut(Side, Side) -> f64) -> Placement {
    let cnt_known = field.iter().flatten().filter(|&id| *id >= 0).count();
    let mut all_potential_edges = vec![];
    for r in 0..field.len() {
        for c in 0..field[r].len() {
            for dr in [0, 1] {
                let dc = 1 - dr;
                if r + dr < field.len() && c + dc < field[r].len() {
                    let id1 = field[r][c];
                    let id2 = field[r + dr][c + dc];
                    if id1 >= 0 && id2 >= 0 {
                        let id1 = id1 as usize;
                        let id2 = id2 as usize;
                        for side1 in 0..4 {
                            for side2 in 0..4 {
                                let s1 = Side {
                                    fig: id1,
                                    side: side1,
                                };
                                let s2 = Side {
                                    fig: id2,
                                    side: side2,
                                };
                                all_potential_edges.push((s1, s2));
                            }
                        }
                    }
                }
            }
        }
    }
    all_potential_edges.sort_by(|&(s1, s2), &(s3, s4)| dist(s1, s2).total_cmp(&dist(s3, s4)));
    let mut placement = Placement::new();
    for &(s1, s2) in all_potential_edges.iter() {
        let placement_copy = placement.clone();
        if placement.join_sides(s1, s2).is_some() {
            if !placement.get_all_neighbours().iter().all(|&(s1, s2)| {
                all_potential_edges.contains(&(s1, s2)) || all_potential_edges.contains(&(s2, s1))
            }) {
                placement = placement_copy;
            } else {
                eprintln!("JOINED {:?} {:?}", s1, s2);
            }
        }
    }
    assert_eq!(placement.get_cnt_figures(), cnt_known);
    let (b_min, b_max) = normalize_bounding_box(placement.get_bounding_box());
    let (exp_min, exp_max) = (
        min(field.len(), field[0].len()),
        max(field.len(), field[0].len()),
    );
    assert_eq!(b_min as usize, exp_min);
    assert_eq!(b_max as usize, exp_max);
    placement
}

pub fn get_known_placement(graph: &Graph) -> Placement {
    let mut total_placement = Placement::new();
    let dist = graph.gen_adj_matrix();
    let dist = |s1: Side, s2: Side| dist[[s1.fig, s1.side, s2.fig, s2.side]];
    {
        let field = [
            vec![/*686, */ 717, 102, 365, 401, 969, 343, 637, 926, 999],
            vec![/*648, */ 542, 1006, 577, 979, 676, 142, 683, 421, 574],
            vec![1052, 620, 737, 832, 292, 215, 827, 424, 963],
            vec![-1, -1, -1, -1, -1, -1, 882, 448, 995],
            vec![-1, -1, -1, -1, -1, -1, 980, 508, 628],
            // vec![-1, -1, -1, -1, -1, -1, -1, 794, 157, 623],
        ];
        let placement = placement_by_field(&field, &dist);
        assert!(total_placement.join_with(&placement));
    }
    {
        let field = [vec![439, 570, 548]];
        let placement = placement_by_field(&field, &dist);
        assert!(total_placement.join_with(&placement));
    }
    {
        let field = [vec![458, 777]];
        let placement = placement_by_field(&field, &dist);
        assert!(total_placement.join_with(&placement));
    }
    {
        let field = [
            vec![-1, 258, 383, 337, 871, 417],
            vec![815, 436, 756, 713, 974, 471],
        ];
        let placement = placement_by_field(&field, &dist);
        assert!(total_placement.join_with(&placement));
    }
    total_placement
}
