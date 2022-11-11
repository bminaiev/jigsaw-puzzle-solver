use crate::{border_matcher::match_borders, parsed_puzzles::ParsedPuzzles};

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Edge {
    pub fig1: usize,
    pub fig2: usize,
    pub side1: usize,
    pub side2: usize,
    pub score: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
    pub n: usize,
    pub all_edges: Vec<Edge>,
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
                            });
                        }
                    }
                }
            }
        }
        Graph {
            n: parsed_puzzles.figures.len(),
            all_edges,
        }
    }
}
