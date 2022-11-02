//! A good way of displaying an SVG image in egui.
//!
//! Requires the dependency `egui_extras` with the `svg` feature.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use eframe::{
    egui::{self, Image, ScrollArea, Widget},
    epaint::{pos2, Color32, Mesh, Rect, Rounding, Shape},
};

use crate::my_widget::MyWidget;

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
    svg_image: egui_extras::RetainedImage,
    my_widget: MyWidget,
}

fn load_image_from_path(path: &str) -> Result<egui::ColorImage, image::ImageError> {
    let image = image::io::Reader::open(path)?.decode()?;
    let size = [image.width() as _, image.height() as _];
    let image_buffer = image.to_rgba8();
    let pixels = image_buffer.as_flat_samples();
    Ok(egui::ColorImage::from_rgba_unmultiplied(
        size,
        pixels.as_slice(),
    ))
}

impl Default for MyApp {
    fn default() -> Self {
        let color_image = load_image_from_path("img/4.jpg").unwrap();
        let color_image2 = load_image_from_path("img/4.jpg").unwrap();
        let svg_image = egui_extras::RetainedImage::from_color_image("test", color_image);
        let svg_image2 = egui_extras::RetainedImage::from_color_image("test", color_image2);
        Self {
            svg_image: svg_image,
            my_widget: MyWidget::new(svg_image2),
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
