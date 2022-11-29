#![feature(slice_group_by)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use std::fs;

use eframe::{egui, epaint::pos2};

use crate::{
    borders_graph::Graph, crop::crop, graph_solver::solve_graph, my_widget::MyWidget,
    parsed_puzzles::ParsedPuzzles, placement::Placement, point::PointF,
    surface_placer::place_on_surface, utils::load_image_from_path,
};

mod border_matcher;
mod borders_graph;
mod coordinate_system;
mod crop;
mod dsu;
mod figure;
mod graph_solver;
mod matcher_tests;
mod my_widget;
mod parsed_puzzles;
mod placement;
mod point;
mod rects_fitter;
mod surface_placer;
mod topn;
mod utils;

const BEFORE_CROP_PATH: &str = "img/prod/22.jpg";
const PATH: &str = "img/crop.jpg";
const GRAPH_PATH: &str = "graph_with_start.json";
const GRAPH_SOLUTION_PATH: &str = "graph_solution.json";
const LOAD_EXISTING_SOLUTION: bool = false;

const PUZZLE_PIXEL_WHITE_THRESHOLD: usize = 460;

// TODO: nicer type
fn main_ui(
    positions: Vec<Option<Vec<PointF>>>,
    path: &str,
    show_parsed: bool,
    show_image: bool,
    show_matched_borders: bool,
) {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1400.0, 1100.0)),
        ..Default::default()
    };
    let app_created = Box::new(MyApp::new(
        positions,
        path,
        show_parsed,
        show_image,
        show_matched_borders,
    ));
    eframe::run_native("jigsaw solver", options, Box::new(|_| app_created));
}

fn main_build_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph = Graph::new(&parsed_puzzles);
    fs::write(GRAPH_PATH, serde_json::to_string(&graph).unwrap()).unwrap();
}

fn main_load_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph: Graph = serde_json::from_str(&fs::read_to_string(GRAPH_PATH).unwrap()).unwrap();
    eprintln!("graph loaded! n = {}", graph.n);
    let solution_graph = if LOAD_EXISTING_SOLUTION {
        serde_json::from_str(&fs::read_to_string(GRAPH_SOLUTION_PATH).unwrap()).unwrap()
    } else {
        let res = solve_graph(&graph, &parsed_puzzles);
        fs::write(GRAPH_SOLUTION_PATH, serde_json::to_string(&res).unwrap()).unwrap();
        res
    };
    let placement = Placement::from_full_graph(&solution_graph);
    let positions = place_on_surface(&graph, &placement, &parsed_puzzles);
    eprintln!("positions generated!");
    main_ui(positions, PATH, true, false, true);
}

fn main_before_crop() {
    main_ui(vec![], BEFORE_CROP_PATH, false, true, false);
}

fn main_check_parsing() {
    main_ui(vec![], PATH, true, false, true);
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

fn main() {
    main_before_crop();
    // main_check_parsing();
    // main_build_graph();
    // main_load_graph();
}

struct MyApp {
    my_widget: MyWidget,
}

impl MyApp {
    fn new(
        positions: Vec<Option<Vec<PointF>>>,
        path: &str,
        show_parsed: bool,
        show_image: bool,
        show_matched_borders: bool,
    ) -> Self {
        Self {
            my_widget: MyWidget::new(
                path,
                positions,
                show_parsed,
                show_image,
                show_matched_borders,
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
