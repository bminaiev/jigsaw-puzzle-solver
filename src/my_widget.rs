use std::f32::consts::PI;

use eframe::{
    egui::{Key, Sense},
    emath::Align2,
    epaint::{
        pos2, vec2, CircleShape, Color32, ColorImage, FontFamily, FontId, Mesh, Pos2, Rect, Shape,
        Stroke, Vec2,
    },
};
use egui_extras::RetainedImage;
use image::ImageBuffer;
use rand::Rng;

use crate::{crop::crop, dsu::Dsu};

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

#[derive(Clone, Copy)]
struct Point {
    x: usize,
    y: usize,
}

impl Point {
    pub fn dist2(&self, other: &Point) -> usize {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    pub fn pos2(&self) -> Pos2 {
        pos2(self.x as f32, self.y as f32)
    }

    pub fn angle_to(&self, other: &Self) -> f32 {
        let dy = (other.y as f32) - (self.y as f32);
        let dx = (other.x as f32) - (self.x as f32);
        dy.atan2(dx)
    }
}

fn find_cycle(points: &[Point]) -> Vec<Point> {
    let mut used = vec![false; points.len()];
    used[0] = true;
    let mut path = vec![points[0]];
    loop {
        let last = *path.last().unwrap();
        let mut next = None;
        for i in 0..used.len() {
            if !used[i] {
                let ndist = last.dist2(&points[i]);
                if next.is_none() || last.dist2(&points[next.unwrap()]) > ndist {
                    next = Some(i);
                }
            }
        }
        if let Some(next) = next {
            path.push(points[next]);
            used[next] = true;
        } else {
            break;
        }
    }
    optimize_cycle(&mut path);
    let area = calc_signed_area(&path);
    if area > 0.0 {
        // make counter-clock-wise
        path.reverse();
    }
    path
}

fn cycle_rev(a: &mut [Point], mut from: usize, mut len: usize) {
    while len > 1 {
        let end = (from + len - 1) % a.len();
        a.swap(from, end);
        from = (from + 1) % a.len();
        len -= 2;
    }
}

fn optimize_cycle(cycle: &mut Vec<Point>) {
    // reverse a part of the cycle to minimize sum of dist^2 between neighbours
    loop {
        let mut changed = false;
        for i in 0..cycle.len() {
            for len in 2..(cycle.len() + 1) / 2 {
                let first = cycle[i];
                let prev = cycle[(i + cycle.len() - 1) % cycle.len()];
                let last = cycle[(i + len - 1) % cycle.len()];
                let next = cycle[(i + len) % cycle.len()];

                let cur_dist = prev.dist2(&first) + last.dist2(&next);
                let potential_dist = prev.dist2(&last) + first.dist2(&next);
                if potential_dist < cur_dist {
                    changed = true;
                    cycle_rev(cycle, i, len);
                }
            }
        }

        if !changed {
            break;
        }
    }
}

fn find_center(points: &[Point]) -> Point {
    let sum_x: usize = points.iter().map(|p| p.x).sum();
    let sum_y: usize = points.iter().map(|p| p.y).sum();
    let cnt = points.len();
    Point {
        x: sum_x / cnt,
        y: sum_y / cnt,
    }
}

fn calc_signed_area(pts: &[Point]) -> f32 {
    let mut res = 0.0;
    for i in 0..pts.len() {
        let p1 = pts[i];
        let p2 = pts[(i + 1) % pts.len()];
        res += (p1.x as f32) * (p2.y as f32) - (p1.y as f32) * (p2.x as f32);
    }
    res
}

fn fmax(x: f32, y: f32) -> f32 {
    if x > y {
        x
    } else {
        y
    }
}

fn fmin(x: f32, y: f32) -> f32 {
    if x < y {
        x
    } else {
        y
    }
}

fn filter_4_corners(all: &[Point], idxs: &[usize]) -> Vec<usize> {
    let mut best_score = f32::MAX;
    let n = idxs.len();
    let mut best = vec![];

    fn norm_angle(mut a: f32) -> f32 {
        while a < 0.0 {
            a += PI * 2.0;
        }
        while a >= PI * 2.0 {
            a -= PI * 2.0;
        }
        fmin(a, PI * 2.0 - a)
    }

    fn find_angle(p1: Point, p2: Point, p3: Point) -> f32 {
        let a1 = p1.angle_to(&p2);
        let a2 = p2.angle_to(&p3);
        norm_angle(a1 - a2)
    }

    fn angle_score(a: f32) -> f32 {
        (PI / 2.0 - a).abs()
    }

    let calc_score = |pts: &[Point]| -> f32 {
        let sum_angles: f32 = (0..4)
            .map(|i| angle_score(find_angle(pts[i], pts[(i + 1) % 4], pts[(i + 2) % 4])))
            .sum();

        let dists: Vec<_> = (0..4).map(|i| pts[i].dist2(&pts[(i + 1) % 4])).collect();
        let max_dist = *dists.iter().max().unwrap();
        let min_dist = *dists.iter().min().unwrap();
        let dist_coef = (max_dist as f32) / (min_dist as f32);

        sum_angles * dist_coef
    };

    for i1 in 0..n {
        for i2 in i1 + 1..n {
            for i3 in i2 + 1..n {
                for i4 in i3 + 1..n {
                    let pts = [all[idxs[i1]], all[idxs[i2]], all[idxs[i3]], all[idxs[i4]]];
                    let score = calc_score(&pts);
                    if score < best_score {
                        best_score = score;
                        best = vec![idxs[i1], idxs[i2], idxs[i3], idxs[i4]];
                    }
                }
            }
        }
    }
    best
}

fn find_corners_positions(points: &[Point]) -> Vec<usize> {
    let center = find_center(points);
    let dists: Vec<_> = points.iter().map(|p| p.dist2(&center)).collect();
    const CHECK_LEN: usize = 5;
    let mut corners = vec![];
    for i in 0..dists.len() {
        let (_, max_pos) = (0..CHECK_LEN * 2)
            .map(|shift| (dists[(i + shift) % dists.len()], shift))
            .max()
            .unwrap();
        if max_pos == CHECK_LEN {
            corners.push((i + max_pos) % dists.len());
        }
    }
    filter_4_corners(&points, &corners)
}

#[derive(Clone)]
struct Figure {
    all_pts: Vec<Point>,
    border: Vec<Point>,
    good_border: bool,
    center: Point,
}

impl Figure {
    pub fn new(all_pts: &[Point], border: &[Point]) -> Self {
        let cycle = find_cycle(border);
        let max_dist2 = cycle
            .iter()
            .circular_tuple_windows()
            .map(|(w0, w1)| w0.dist2(&w1))
            .max()
            .unwrap();
        let center = find_center(&cycle);

        Self {
            all_pts: all_pts.iter().cloned().collect(),
            border: cycle,
            good_border: max_dist2 <= 10,
            center,
        }
    }
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

        for (figure_id, figure) in self.parsed_puzzles.figures.iter().enumerate() {
            if figure.good_border {
                let cornre_pos = find_corners_positions(&figure.border);
                for pos in cornre_pos.into_iter() {
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
