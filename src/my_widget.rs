use std::f32::consts::PI;

use eframe::{
    egui::{Key, Sense},
    emath::Align2,
    epaint::{
        pos2, vec2, CircleShape, Color32, ColorImage, FontFamily, FontId, Mesh, Pos2, Rect,
        Rounding, Shape, Stroke, Vec2,
    },
};
use egui_extras::RetainedImage;
use image::ImageBuffer;
use rand::Rng;

use crate::{
    border_matcher::match_borders,
    crop::crop,
    dsu::Dsu,
    figure::Figure,
    point::{Point, PointF},
};

use itertools::Itertools;

pub struct MyWidget {
    offset: Vec2,
    zoom_log: f32,
    image: RetainedImage,
    frame: Vec<Pos2>,
    image_path: String,
    mask_image: RetainedImage,
    parsed_puzzles: ParsedPuzzles,
}

const ZOOM_DELTA_COEF: f32 = 500.0;

fn load_image_from_path(path: &str) -> Result<ColorImage, image::ImageError> {
    let image = image::io::Reader::open(path)?.decode()?;
    let size = [image.width() as _, image.height() as _];
    let image_buffer = image.to_rgba8();
    let pixels = image_buffer.as_flat_samples();
    Ok(ColorImage::from_rgba_unmultiplied(size, pixels.as_slice()))
}

fn is_puzzle_color(color: Color32) -> bool {
    (color.r() as usize) + (color.g() as usize) + (color.b() as usize) >= 630
}

fn save_color_image(color_image: &ColorImage, path: &str) {
    let mut new_img = ImageBuffer::new(color_image.size[0] as u32, color_image.size[1] as u32);
    for x in 0..new_img.width() {
        for y in 0..new_img.height() {
            let color = color_image[(x as usize, y as usize)];
            new_img[(x, y)] = image::Rgba([color.r(), color.g(), color.b(), color.a()]);
        }
    }
    new_img.save(path).unwrap();
}

fn calc_dist_to_not_puzzle(
    width: usize,
    height: usize,
    mut id: impl FnMut(usize, usize) -> usize,
    is_puzzle: &[bool],
) -> Vec<usize> {
    let mut dist_to_not_puzzle = vec![std::usize::MAX / 2; width * height];

    fn update_min(with: usize, place: &mut usize) {
        if *place > with {
            *place = with;
        }
    }

    for x in 0..width {
        for y in 0..height {
            if is_puzzle[id(x, y)] {
                if x > 0 {
                    update_min(
                        1 + dist_to_not_puzzle[id(x - 1, y)],
                        &mut dist_to_not_puzzle[id(x, y)],
                    );
                }
                if y > 0 {
                    update_min(
                        1 + dist_to_not_puzzle[id(x, y - 1)],
                        &mut dist_to_not_puzzle[id(x, y)],
                    );
                }
            } else {
                dist_to_not_puzzle[id(x, y)] = 0;
            }
        }
    }
    for x in (0..width).rev() {
        for y in (0..height).rev() {
            if is_puzzle[id(x, y)] {
                if x + 1 < width {
                    update_min(
                        1 + dist_to_not_puzzle[id(x + 1, y)],
                        &mut dist_to_not_puzzle[id(x, y)],
                    );
                }
                if y + 1 < height {
                    update_min(
                        1 + dist_to_not_puzzle[id(x, y + 1)],
                        &mut dist_to_not_puzzle[id(x, y)],
                    );
                }
            } else {
                dist_to_not_puzzle[id(x, y)] = 0;
            }
        }
    }
    dist_to_not_puzzle
}

struct ParsedPuzzles {
    width: usize,
    height: usize,
    figures: Vec<Figure>,
}

impl ParsedPuzzles {
    pub fn new(color_image: &ColorImage) -> Self {
        let width = color_image.size[0];
        let height = color_image.size[1];

        let id = |x: usize, y: usize| -> usize { x + y * width };
        let mut is_puzzle = vec![false; width * height];

        for x in 0..width {
            for y in 0..height {
                if is_puzzle_color(color_image[(x, y)]) {
                    is_puzzle[id(x, y)] = true;
                }
            }
        }

        let dist_to_not_puzzle = calc_dist_to_not_puzzle(width, height, id, &is_puzzle);

        let mut dsu = Dsu::new(width * height);
        for x in 0..width {
            for y in 0..height {
                if is_puzzle[id(x, y)] {
                    if x + 1 < width && is_puzzle[id(x + 1, y)] {
                        dsu.unite(id(x, y), id(x + 1, y));
                    }
                    if y + 1 < height && is_puzzle[id(x, y + 1)] {
                        dsu.unite(id(x, y), id(x, y + 1));
                    }
                }
            }
        }

        let mut figures_borders = vec![vec![]; width * height];
        let mut figures_inside = vec![vec![]; width * height];
        for x in 0..width {
            for y in 0..height {
                if is_puzzle[id(x, y)] {
                    if dist_to_not_puzzle[id(x, y)] == 1 {
                        figures_borders[dsu.get(id(x, y))].push(Point { x, y });
                    }
                    figures_inside[dsu.get(id(x, y))].push(Point { x, y });
                }
            }
        }

        let mut res_figures = vec![];

        for i in 0..figures_borders.len() {
            if figures_borders[i].len() < 20 {
                continue;
            }
            res_figures.push(Figure::new(&figures_inside[i], &figures_borders[i]));
        }

        res_figures.sort_by_key(|f| (f.center.y, f.center.x));

        let f1 = &res_figures[243];

        for f2_id in 0..res_figures.len() {
            if res_figures[f2_id].good_border && res_figures[f2_id].corner_positions.len() == 4 {
                for j in 0..4 {
                    let result = match_borders(f1, 2, &res_figures[f2_id], j);
                    if result.score < 50.0 {
                        eprintln!("{}x{} -> {}", f2_id, j, result.score);
                    }
                }
            }
        }

        Self {
            width,
            height,
            figures: res_figures,
        }
    }

    pub fn gen_image(&self) -> ColorImage {
        let mut res = ColorImage::new([self.width, self.height], Color32::WHITE);

        let mut rng = rand::thread_rng();
        for figure in self.figures.iter() {
            let color = Color32::from_rgb(rng.gen(), rng.gen(), rng.gen());
            for p in figure.all_pts.iter() {
                res[(p.x, p.y)] = color;
            }
            let border_color = if figure.good_border {
                Color32::BLACK
            } else {
                Color32::RED
            };
            for p in figure.border.iter() {
                res[(p.x, p.y)] = border_color;
            }
            // res[(figure.center.x, figure.center.y)] = Color32::BLUE;
        }

        res
    }
}

impl MyWidget {
    pub fn new(path: &str) -> Self {
        let color_image = load_image_from_path(path).unwrap();

        let parsed_puzzles = ParsedPuzzles::new(&color_image);
        let mask_image = {
            let color = parsed_puzzles.gen_image();
            save_color_image(&color, "img/puzzle.jpg");
            RetainedImage::from_color_image("mask", color)
        };

        let image = RetainedImage::from_color_image("test", color_image);
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
            mask_image,
            parsed_puzzles,
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
            self.zoom_log += delta / ZOOM_DELTA_COEF;
            let zoom = self.get_zoom();
            // (real.x + offset.x) * zoom = mouse.x
            self.offset = pos2(mouse.x / zoom, mouse.y / zoom) - real_mouse_pos;
        }
    }

    fn find_closest_object(
        &self,
        pts: &[Pos2],
        ui: &eframe::egui::Ui,
        max_dist: f32,
    ) -> Option<usize> {
        let mut options = vec![];

        if let Some(mouse) = ui.input().pointer.hover_pos() {
            for (_i, p) in pts.iter().enumerate() {
                let screen_p = self.convert_to_screen(p.clone());
                let dist_to_mouse = screen_p.distance(mouse);
                // TODO: change to const
                if dist_to_mouse < max_dist {
                    options.push((_i, dist_to_mouse));
                }
            }
        }

        options.sort_by(|c1, c2| c1.1.partial_cmp(&c2.1).unwrap());

        options.get(0).map(|(x, _)| *x)
    }

    fn which_figure_hovered(&self, ui: &eframe::egui::Ui) -> Option<usize> {
        let centers: Vec<_> = self
            .parsed_puzzles
            .figures
            .iter()
            .map(|f| f.center.pos2())
            .collect();
        let closest_figure_id = self.find_closest_object(&centers, ui, f32::MAX)?;
        if let Some(mouse) = ui.input().pointer.hover_pos() {
            for p in self.parsed_puzzles.figures[closest_figure_id]
                .all_pts
                .iter()
            {
                if mouse.distance(self.convert_to_screen(p.pos2())) <= 25.0 {
                    return Some(closest_figure_id);
                }
            }
        }
        None
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
            if false {
                ui.painter().add(Shape::mesh(mesh));
            }

            if true {
                let texture_id2 = self.mask_image.texture_id(&ui.ctx());
                let mut mesh2 = Mesh::with_texture(texture_id2);
                mesh2.add_rect_with_uv(rect, uv, Color32::WHITE);
                ui.painter().add(Shape::mesh(mesh2));
            }
        }

        if false {
            for pos in self.frame.iter() {
                ui.painter().add(Shape::Circle(CircleShape::filled(
                    self.convert_to_screen(pos.clone()),
                    10.0,
                    Color32::RED,
                )));
            }
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
        if let Some(id) = self.find_closest_object(&self.frame, ui, 30.0) {
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

        let show_border = |figure: &Figure| {
            for i in 0..figure.border.len() {
                let p1 = self.convert_to_screen(figure.border[i].pos2());

                let color = Color32::GREEN;
                ui.painter().circle_filled(p1, 3.0, color);
                let p2 =
                    self.convert_to_screen(figure.border[(i + 1) % figure.border.len()].pos2());
                ui.painter()
                    .add(Shape::line_segment([p1, p2], Stroke::new(2.0, color)));
            }
        };

        if let Some(figure_id) = self.which_figure_hovered(ui) {
            let figure = &self.parsed_puzzles.figures[figure_id];
            show_border(figure);
        }

        for &figure_id in [243, 278].iter() {
            show_border(&self.parsed_puzzles.figures[figure_id]);
        }

        {
            let min = self.convert_to_screen(pos2(5000.0, 0.0));
            let max = self.convert_to_screen(pos2(8000.0, 4000.0));
            ui.painter().rect_filled(
                Rect::from_min_max(min, max),
                Rounding::default(),
                Color32::WHITE,
            );

            let f1 = &self.parsed_puzzles.figures[243];
            let f2 = &self.parsed_puzzles.figures[278];
            for i in 0..4 {
                for j in 0..4 {
                    let result = match_borders(f1, i, f2, j);

                    let conv_point = |p: PointF| -> Pos2 {
                        pos2(
                            p.x as f32 + 5000.0 + (i as f32) * 200.0,
                            p.y as f32 + (j as f32) * 200.0,
                        )
                    };

                    let draw = |pts: &[PointF], color: Color32| {
                        for seg in pts.windows(2) {
                            let p1 = conv_point(seg[0]);
                            let p2 = conv_point(seg[1]);
                            ui.painter().line_segment(
                                [self.convert_to_screen(p1), self.convert_to_screen(p2)],
                                Stroke::new(2.0, color),
                            );
                            ui.painter()
                                .circle_filled(self.convert_to_screen(p1), 3.0, color);
                        }
                    };

                    draw(&result.lhs, Color32::BLUE);
                    draw(&result.rhs, Color32::RED);
                    // eprintln!("{}x{} -> {}", i, j, result.score);
                }
            }
        }

        for (figure_id, figure) in self.parsed_puzzles.figures.iter().enumerate() {
            if figure.good_border {
                for &pos in figure.corner_positions.iter() {
                    let p = self.convert_to_screen(figure.border[pos].pos2());
                    ui.painter().add(Shape::circle_filled(p, 4.0, Color32::RED));
                }

                let center = self.convert_to_screen(figure.center.pos2());
                ui.painter().text(
                    center,
                    Align2::CENTER_CENTER,
                    figure_id.to_string(),
                    FontId::new(20.0, FontFamily::Monospace),
                    Color32::BLACK,
                );
            }
        }

        let frame = self.frame.clone();
        if ui.input().key_pressed(Key::Enter) {
            // dbg!("CROP!", &frame);
            // crop(&self.image_path, &frame);
            // dbg!("CROPED!");
        }

        response
    }
}
