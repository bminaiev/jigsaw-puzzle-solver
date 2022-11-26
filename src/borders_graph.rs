use crate::{
    border_matcher::{match_borders, match_borders_without_move},
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
};

use ndarray::Array4;
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

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
    pub n: usize,
    pub all_edges: Vec<Edge>,
    pub parsed_puzzles_hash: u64,
}
impl Graph {
    pub fn new(parsed_puzzles: &ParsedPuzzles) -> Self {
        let mut all_edges = vec![];
        let figures = &parsed_puzzles.figures;

        for fig1 in 0..figures.len() {
            eprintln!("{}/{}", fig1, figures.len());
            for fig2 in fig1 + 1..figures.len() {
                if !figures[fig1].is_good_puzzle() || !figures[fig2].is_good_puzzle() {
                    continue;
                }
                for side1 in 0..4 {
                    for side2 in 0..4 {
                        let existing_edge = match_borders_without_move(
                            &figures[fig1],
                            side1,
                            &figures[fig2],
                            side2,
                            fig1,
                            fig2,
                        )
                        .is_some();
                        if let Some(res) =
                            match_borders(&figures[fig1], side1, &figures[fig2], side2, fig1, fig2)
                        {
                            let score = res.score;
                            all_edges.push(Edge {
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
            }
        }
        Graph {
            n: parsed_puzzles.figures.len(),
            all_edges,
            parsed_puzzles_hash: parsed_puzzles.calc_hash(),
        }
    }

    pub fn gen_adj_matrix(&self) -> Array4<f64> {
        let n = self.n;
        let mut dist = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 10.0);
        for edge in self.all_edges.iter() {
            dist[[edge.fig1, edge.side1, edge.fig2, edge.side2]] = edge.score;
            dist[[edge.fig2, edge.side2, edge.fig1, edge.side1]] = edge.score;
        }
        dist
    }
}
