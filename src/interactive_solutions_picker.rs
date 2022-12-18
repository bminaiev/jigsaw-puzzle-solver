use std::time::Instant;

use eframe::epaint::ColorImage;
use egui_extras::RetainedImage;
use itertools::Itertools;
use ndarray::Array4;

use crate::{
    borders_graph::Graph,
    graph_solver::{gen_potential_solution, PotentialSolution},
    known_facts::KnownFacts,
    parsed_puzzles::ParsedPuzzles,
    placement::{Placement, PotentialGroupLocation, Search3StateWithScore},
    point::PointF,
    positions_cache::PositionsCache,
    surface_placer::put_solutions_on_surface,
};

pub struct InteractiveSolutionPicker {
    pub solutions_to_show: Vec<PotentialSolution>,
    pub mask_image: RetainedImage,
    all_solutions: Vec<(Search3StateWithScore, PotentialSolution)>,
    start_vertex: usize,
    rot_positions: Vec<Option<Vec<PointF>>>,
    positions_cache: PositionsCache,
    base_points_matrix: Array4<[PointF; 2]>,
}

impl InteractiveSolutionPicker {
    pub fn new(
        all_solutions: Vec<(Search3StateWithScore, PotentialSolution)>,
        start_vertex: usize,
        rot_positions: Vec<Option<Vec<PointF>>>,
        positions_cache: PositionsCache,
        base_points_matrix: Array4<[PointF; 2]>,
        known_facts: &KnownFacts,
        parsed_puzzles: &ParsedPuzzles,
        graph: &Graph,
    ) -> Self {
        let mut res = Self {
            solutions_to_show: vec![],
            mask_image: RetainedImage::from_color_image("", ColorImage::default()),
            all_solutions,
            start_vertex,
            rot_positions,
            positions_cache,
            base_points_matrix,
        };
        res.refresh(known_facts, parsed_puzzles, graph);
        res
    }

    pub fn refresh(
        &mut self,
        known_facts: &KnownFacts,
        parsed_puzzles: &ParsedPuzzles,
        graph: &Graph,
    ) {
        eprintln!("Start interactive solution picker refresh!");
        let start = Instant::now();
        self.all_solutions = self
            .all_solutions
            .iter()
            .filter_map(|(state, sol)| {
                let new_sol = sol.refresh(known_facts)?;
                Some((state.clone(), new_sol.clone()))
            })
            .collect_vec();
        self.solutions_to_show.clear();
        self.solutions_to_show.push(self.gen_basic_known_solution(
            known_facts,
            parsed_puzzles,
            graph,
        ));
        self.all_solutions.sort_by(|(a, b), (c, d)| {
            a.cmp(&c)
                .then(b.placement_score.total_cmp(&d.placement_score))
        });

        {
            let mut i = 0;
            while i != self.all_solutions.len() {
                let mut j = i;
                while j != self.all_solutions.len()
                    && self.all_solutions[i].0.get_key() == self.all_solutions[j].0.get_key()
                {
                    j += 1;
                }
                let mut r = self.all_solutions[i].1.clone();
                if i + 1 != j {
                    let rat = self.all_solutions[i + 1].1.placement_score
                        / self.all_solutions[i].1.placement_score;
                    r.placement_score /= rat.powf(2.0);
                    r.additional_text += &format!(
                        ". next = {:.3}",
                        self.all_solutions[i + 1].1.placement_score
                    );
                }
                self.solutions_to_show.push(r);
                i = j;
            }
        }
        self.solutions_to_show.sort();
        self.solutions_to_show.truncate(25);
        put_solutions_on_surface(&mut self.solutions_to_show);

        let mask_image = RetainedImage::from_color_image(
            "solutions mask",
            PotentialSolution::gen_image(&self.solutions_to_show),
        );
        self.mask_image = mask_image;
        eprintln!("Finished refresh in {:?}", start.elapsed());
    }

    fn gen_basic_known_solution(
        &self,
        known_facts: &KnownFacts,
        parsed_puzzles: &ParsedPuzzles,
        graph: &Graph,
    ) -> PotentialSolution {
        let mut placement = Placement::new();
        for fact in known_facts.facts.iter() {
            if fact.good_edge {
                placement.join_sides(fact.side1, fact.side2).unwrap();
            }
        }
        let edges = placement.get_all_neighbours_in_same_component(self.start_vertex);
        gen_potential_solution(
            &edges,
            parsed_puzzles,
            graph,
            Some(0.0),
            known_facts,
            &self.rot_positions,
            &self.positions_cache,
            &self.base_points_matrix,
        )
        .unwrap()
    }
}
