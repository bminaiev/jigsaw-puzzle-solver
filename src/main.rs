#![feature(slice_group_by)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use std::fs;

use eframe::{egui, epaint::pos2};

use crate::{
    borders_graph::Graph,
    crop::crop,
    edge_score_optimizer::optimize_edge_scores,
    graph_solver::{solve_graph, solve_graph_border, PotentialSolution},
    my_widget::MyWidget,
    parsed_puzzles::ParsedPuzzles,
    placement::Placement,
    point::PointF,
    surface_placer::{place_on_surface, put_solutions_on_surface},
    utils::load_image_from_path,
};

mod average_color;
mod border_matcher;
mod borders_graph;
mod coordinate_system;
mod crop;
mod dsu;
mod edge_score_optimizer;
mod figure;
mod graph_solver;
mod known_positions;
mod matcher_tests;
mod my_widget;
mod parsed_puzzles;
mod placement;
mod point;
mod rects_fitter;
mod surface_placer;
mod topn;
mod utils;

const BEFORE_CROP_PATH: &str = "img/prod2/24.jpg";
const PATH: &str = "img/prod2/crop2_full.jpg";
const GRAPH_PATH: &str = "graph_with_start.json";
const GRAPH_SOLUTION_PATH: &str = "graph_solution.json";
const LOAD_EXISTING_SOLUTION: bool = false;

const PUZZLE_PIXEL_WHITE_THRESHOLD: usize = 460;

// TODO: nicer type
fn main_ui(
    solutions: Vec<PotentialSolution>,
    path: &str,
    show_parsed: bool,
    show_image: bool,
    show_matched_borders: bool,
    crop_enabled: bool,
) {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1400.0, 1100.0)),
        ..Default::default()
    };
    let app_created = Box::new(MyApp::new(
        solutions,
        path,
        show_parsed,
        show_image,
        show_matched_borders,
        crop_enabled,
    ));
    eframe::run_native("jigsaw solver", options, Box::new(|_| app_created));
}

fn main_build_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph = Graph::new(&parsed_puzzles, false);
    fs::write(GRAPH_PATH, serde_json::to_string(&graph).unwrap()).unwrap();
}

fn main_load_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph: Graph = serde_json::from_str(&fs::read_to_string(GRAPH_PATH).unwrap()).unwrap();
    eprintln!("graph loaded! n = {}", graph.n);

    let mut solution_graph = {
        let mut prev_state: Option<Graph> = if LOAD_EXISTING_SOLUTION {
            Some(serde_json::from_str(&fs::read_to_string(GRAPH_SOLUTION_PATH).unwrap()).unwrap())
        } else {
            None
        };
        solve_graph(&graph, &parsed_puzzles, prev_state.clone())
    };
    // fs::write(
    //     GRAPH_SOLUTION_PATH,
    //     serde_json::to_string(&solution_graph).unwrap(),
    // )
    // .unwrap();
    // let placement = Placement::from_full_graph(&solution_graph);
    // let positions = place_on_surface(&graph, &placement, &parsed_puzzles);
    put_solutions_on_surface(&mut solution_graph);
    eprintln!("positions generated!");
    main_ui(solution_graph, PATH, true, false, true, false);
}

fn main_before_crop() {
    main_ui(vec![], BEFORE_CROP_PATH, false, true, false, true);
}

fn main_check_parsing() {
    main_ui(vec![], PATH, true, true, true, false);
}

fn main_check_crop() {
    let pts = [
        pos2(1576.5, 1072.4),
        pos2(2631.0, 1469.2),
        pos2(1682.6, 3007.4),
        pos2(364.4, 2241.2),
    ];
    crop(BEFORE_CROP_PATH, &pts);
}

fn main_optimize_edge_scoring() {
    let color_image = load_image_from_path(PATH).unwrap();
    let graph: Graph = serde_json::from_str(&fs::read_to_string(GRAPH_PATH).unwrap()).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);

    let mut solutions = optimize_edge_scores(&parsed_puzzles, &graph, true);
    put_solutions_on_surface(&mut solutions);
    main_ui(solutions, PATH, true, false, true, false);
}

fn main() {
    // main_before_crop();
    // main_check_parsing();
    main_build_graph();
    // main_load_graph();
    // main_optimize_edge_scoring();
}

struct MyApp {
    my_widget: MyWidget,
}

impl MyApp {
    fn new(
        solutions: Vec<PotentialSolution>,
        path: &str,
        show_parsed: bool,
        show_image: bool,
        show_matched_borders: bool,
        crop_enabled: bool,
    ) -> Self {
        Self {
            my_widget: MyWidget::new(
                path,
                solutions,
                show_parsed,
                show_image,
                show_matched_borders,
                crop_enabled,
            ),
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            self.my_widget.ui(ui);
        });
    }
}
