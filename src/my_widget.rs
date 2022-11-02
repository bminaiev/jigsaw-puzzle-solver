use eframe::{
    egui::{Key, Sense},
    epaint::{pos2, vec2, CircleShape, Color32, ColorImage, Mesh, Pos2, Rect, Shape, Stroke, Vec2},
};
use egui_extras::RetainedImage;

use crate::crop::crop;

pub struct MyWidget {
    offset: Vec2,
    zoom_log: f32,
    image: RetainedImage,
    frame: Vec<Pos2>,
    image_path: String,
}

const ZOOME_DELTA_COEF: f32 = 500.0;

fn load_image_from_path(path: &str) -> Result<ColorImage, image::ImageError> {
    let image = image::io::Reader::open(path)?.decode()?;
    let size = [image.width() as _, image.height() as _];
    let image_buffer = image.to_rgba8();
    let pixels = image_buffer.as_flat_samples();
    Ok(ColorImage::from_rgba_unmultiplied(size, pixels.as_slice()))
}

impl MyWidget {
    pub fn new(path: &str) -> Self {
        let color_image2 = load_image_from_path(path).unwrap();
        let image = egui_extras::RetainedImage::from_color_image("test", color_image2);
        let img_size = image.size_vec2();
        Self {
            offset: vec2(0.0, 0.0),
            zoom_log: -1.0,
            frame: vec![
                pos2(0.0, 0.0),
                pos2(img_size.x, 0.0),
                pos2(img_size.x, img_size.y),
                pos2(0.0, img_size.y),
            ],
            image,
            image_path: path.to_owned(),
        }
    }

    fn get_zoom(&self) -> f32 {
        self.zoom_log.exp()
    }

    pub fn convert_to_screen(&self, p: Pos2) -> Pos2 {
        let zoom = self.get_zoom();
        pos2((p.x + self.offset.x) * zoom, (p.y + self.offset.y) * zoom)
    }

    pub fn convert_from_screen(&self, p: Pos2) -> Pos2 {
        let zoom = self.get_zoom();
        pos2(p.x / zoom, p.y / zoom) - self.offset
    }

    pub fn change_zoom(&mut self, ui: &eframe::egui::Ui) {
        let delta = ui.input().scroll_delta.y;
        if delta == 0.0 {
            return;
        }
        if let Some(mouse) = ui.input().pointer.hover_pos() {
            let real_mouse_pos = self.convert_from_screen(mouse);
            self.zoom_log += delta / ZOOME_DELTA_COEF;
            let zoom = self.get_zoom();
            // (real.x + offset.x) * zoom = mouse.x
            self.offset = pos2(mouse.x / zoom, mouse.y / zoom) - real_mouse_pos;
        }
    }

    fn find_closest_object(&self, pts: &[Pos2], ui: &eframe::egui::Ui) -> Option<usize> {
        let mut options = vec![];

        if let Some(mouse) = ui.input().pointer.hover_pos() {
            for (_i, p) in pts.iter().enumerate() {
                let screen_p = self.convert_to_screen(p.clone());
                let dist_to_mouse = screen_p.distance(mouse);
                // TODO: change to const
                if dist_to_mouse < 30.0 {
                    options.push((_i, dist_to_mouse));
                }
            }
        }

        options.sort_by(|c1, c2| c1.1.partial_cmp(&c2.1).unwrap());

        options.get(0).map(|(x, _)| *x)
    }

    pub fn ui(&mut self, ui: &mut eframe::egui::Ui) -> eframe::egui::Response {
        let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::drag());
        // {
        //     ui.painter()
        //         .add(Shape::rect_filled(rect, Rounding::none(), Color32::YELLOW));
        // }
        {
            self.change_zoom(ui);
        }

        let img_size = self.image.size_vec2();
        {
            let min = self.convert_to_screen(pos2(0.0, 0.0));
            let max = self.convert_to_screen(pos2(img_size.x, img_size.y));
            let rect = Rect::from_min_max(min, max);
            let texture_id = self.image.texture_id(&ui.ctx());
            let mut mesh = Mesh::with_texture(texture_id);

            let uv = Rect::from_min_max(pos2(0.0, 0.0), pos2(1.0, 1.0));
            mesh.add_rect_with_uv(rect, uv, Color32::WHITE);
            ui.painter().add(Shape::mesh(mesh));
        }

        for pos in self.frame.iter() {
            ui.painter().add(Shape::Circle(CircleShape::filled(
                self.convert_to_screen(pos.clone()),
                10.0,
                Color32::RED,
            )));
        }

        for i in 0..self.frame.len() {
            let p1 = self.convert_to_screen(self.frame[i]);
            let p2 = self.convert_to_screen(self.frame[(i + 1) % self.frame.len()]);
            ui.painter().add(Shape::line_segment(
                [p1, p2],
                Stroke::new(2.0, Color32::RED),
            ));
        }

        let drag_delta = response.drag_delta() / self.get_zoom();
        if let Some(id) = self.find_closest_object(&self.frame, ui) {
            let screen_p = self.convert_to_screen(self.frame[id]);
            ui.painter().add(Shape::circle_stroke(
                screen_p,
                13.0,
                Stroke::new(2.0, Color32::GREEN),
            ));

            self.frame[id] += drag_delta;
        } else {
            self.offset += drag_delta;
        }

        let frame = self.frame.clone();
        if ui.input().key_pressed(Key::Enter) {
            dbg!("CROP!", &frame);
            crop(&self.image_path, &frame);
            dbg!("CROPED!");
        }

        response
    }
}
