#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use eframe::egui;

use crate::my_widget::MyWidget;

mod border_matcher;
mod crop;
mod dsu;
mod figure;
mod my_widget;
mod point;

fn main() {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1400.0, 1100.0)),
        ..Default::default()
    };
    eframe::run_native(
        "jigsaw solver",
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
