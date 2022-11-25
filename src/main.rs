#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use std::fs;

use eframe::egui;

use crate::{
    borders_graph::Graph, graph_solver::solve_graph, my_widget::MyWidget,
    parsed_puzzles::ParsedPuzzles, point::PointF, utils::load_image_from_path,
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
mod topn;
mod utils;

const PATH: &str = "img/crop_only_white.jpg";

// TODO: nicer type
fn main_ui(positions: Vec<Option<Vec<PointF>>>) {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1400.0, 1100.0)),
        ..Default::default()
    };
    eframe::run_native(
        "jigsaw solver",
        options,
        Box::new(|_cc| Box::new(MyApp::new(positions, false))),
    );
}

fn main_build_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph = Graph::new(&parsed_puzzles);
    fs::write("graph.json", serde_json::to_string(&graph).unwrap()).unwrap();
}

fn main_load_graph() {
    let color_image = load_image_from_path(PATH).unwrap();
    let parsed_puzzles = ParsedPuzzles::new(&color_image);
    let graph: Graph = serde_json::from_str(&fs::read_to_string("graph.json").unwrap()).unwrap();
    eprintln!("graph loaded! n = {}", graph.n);
    let positions = solve_graph(&graph, &parsed_puzzles);
    eprintln!("positions generated!");
    main_ui(positions);
}

fn main() {
    // main_build_graph();
    main_load_graph();
}

struct MyApp {
    my_widget: MyWidget,
}

impl MyApp {
    fn new(positions: Vec<Option<Vec<PointF>>>, show_parsed: bool) -> Self {
        Self {
            my_widget: MyWidget::new(PATH, positions, show_parsed),
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
