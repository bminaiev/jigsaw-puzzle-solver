//! A good way of displaying an SVG image in egui.
//!
//! Requires the dependency `egui_extras` with the `svg` feature.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use eframe::{egui, epaint::pos2};

use crate::my_widget::MyWidget;

mod crop;
mod dsu;
mod my_widget;

fn main() {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1400.0, 1100.0)),
        ..Default::default()
    };
    eframe::run_native(
        "svg example",
        options,
        Box::new(|_cc| Box::new(MyApp::default())),
    );
}

struct MyApp {
    my_widget: MyWidget,
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            my_widget: MyWidget::new("img/crop.jpg"),
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
