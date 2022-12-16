use std::collections::BTreeSet;

use crate::{
    border_matcher::{match_borders, match_borders_without_move},
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
    placement::{self, Placement},
    utils::Side,
};

use itertools::Itertools;
use ndarray::Array4;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Edge {
    pub fig1: usize,
    pub fig2: usize,
    pub side1: usize,
    pub side2: usize,
    pub score: f64,
    pub existing_edge: bool,
}

impl Edge {
    pub fn sides(&self) -> (Side, Side) {
        (
            Side {
                fig: self.fig1,
                side: self.side1,
            },
            Side {
                fig: self.fig2,
                side: self.side2,
            },
        )
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Graph {
    pub n: usize,
    pub all_edges: Vec<Edge>,
    pub parsed_puzzles_hash: u64,
}
impl Graph {
    pub fn get_subgraph(&self, placement: &Placement) -> Self {
        let all_sides: BTreeSet<_> = placement.get_all_neighbours().into_iter().collect();
        let all_edges = self
            .all_edges
            .iter()
            .filter(|e| all_sides.contains(&e.sides()))
            .cloned()
            .collect_vec();
        Self {
            n: self.n,
            parsed_puzzles_hash: self.parsed_puzzles_hash,
            all_edges,
        }
    }

    pub fn new(parsed_puzzles: &ParsedPuzzles, only_borders: bool) -> Self {
        let mut all_edges = vec![];
        let figures = &parsed_puzzles.figures;

        let mut ok_sides = BTreeSet::new();
        if only_borders {
            let figures = parsed_puzzles.calc_figures_on_border();
            for fig in figures.iter() {
                ok_sides.insert(fig.left_side);
                ok_sides.insert(fig.right_side);
            }
        } else {
            for fig in 0..figures.len() {
                for side in 0..4 {
                    ok_sides.insert(Side { fig, side });
                }
            }
        }

        for fig1 in 0..figures.len() {
            eprintln!("{}/{}", fig1, figures.len());
            let more_edges: Vec<_> = (fig1 + 1..figures.len())
                .into_par_iter()
                .map(|fig2| {
                    if !figures[fig1].is_good_puzzle() || !figures[fig2].is_good_puzzle() {
                        return vec![];
                    }
                    let mut new_edges = vec![];
                    for side1 in 0..4 {
                        for side2 in 0..4 {
                            let s1 = Side {
                                fig: fig1,
                                side: side1,
                            };
                            let s2 = Side {
                                fig: fig2,
                                side: side2,
                            };
                            if !ok_sides.contains(&s1) || !ok_sides.contains(&s2) {
                                continue;
                            }
                            let existing_edge = match_borders_without_move(
                                &figures[fig1],
                                side1,
                                &figures[fig2],
                                side2,
                                fig1,
                                fig2,
                            )
                            .is_some();
                            if let Some(res) = match_borders(parsed_puzzles, s1, s2) {
                                let score = res.score;
                                new_edges.push(Edge {
                                    fig1,
                                    fig2,
                                    side1,
                                    side2,
                                    score,
                                    existing_edge,
                                });
                                if existing_edge {
                                    eprintln!("Add existing edge: {fig1} {fig2}");
                                }
                            }
                        }
                    }
                    new_edges
                })
                .flatten()
                .collect();
            all_edges.extend(more_edges);
        }
        Graph {
            n: parsed_puzzles.figures.len(),
            all_edges,
            parsed_puzzles_hash: parsed_puzzles.calc_hash(),
        }
    }

    pub fn get_puzzles_with_probably_correct_directions(&self) -> Vec<bool> {
        let mut probably_correct_puzzle_direction = vec![false; self.n];
        for e in self.all_edges.iter() {
            if e.existing_edge {
                probably_correct_puzzle_direction[e.fig1] = true;
                probably_correct_puzzle_direction[e.fig2] = true;
            }
        }
        probably_correct_puzzle_direction
    }

    pub fn gen_adj_matrix(&self) -> Array4<f64> {
        let n = self.n;
        let mut dist = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 50.0);
        for edge in self.all_edges.iter() {
            dist[[edge.fig1, edge.side1, edge.fig2, edge.side2]] = edge.score;
            dist[[edge.fig2, edge.side2, edge.fig1, edge.side1]] = edge.score;
        }
        dist
    }

    pub(crate) fn num_vertices(&self) -> usize {
        let mut seen = BTreeSet::new();
        for e in self.all_edges.iter() {
            seen.insert(e.fig1);
            seen.insert(e.fig2);
        }
        seen.len()
    }
}
