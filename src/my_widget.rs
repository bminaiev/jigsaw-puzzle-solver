use eframe::{
    egui::{Key, Sense},
    epaint::{pos2, vec2, CircleShape, Color32, ColorImage, Mesh, Pos2, Rect, Shape, Stroke, Vec2},
};
use egui_extras::RetainedImage;
use image::ImageBuffer;
use rand::Rng;

use crate::{crop::crop, dsu::Dsu};

pub struct MyWidget {
    offset: Vec2,
    zoom_log: f32,
    image: RetainedImage,
    frame: Vec<Pos2>,
    image_path: String,
    mask_image: RetainedImage,
}

const ZOOME_DELTA_COEF: f32 = 500.0;

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
            let new_dist = last.dist2(&points[next]);
            if new_dist > 2 {
                let mut best_sum_dist = new_dist * 2;
                let mut best_insert_pos = path.len();
                for pos in 0..path.len() - 1 {
                    let cur_dist =
                        path[pos].dist2(&points[next]) + path[pos + 1].dist2(&points[next]);
                    if cur_dist < best_sum_dist {
                        best_sum_dist = cur_dist;
                        best_insert_pos = pos + 1;
                    }
                }
                path.insert(best_insert_pos, points[next]);
                // eprintln!(
                //     "inserted in pos {}/{}, dist = {}",
                //     best_insert_pos,
                //     path.len(),
                //     best_sum_dist
                // );
            } else {
                path.push(points[next]);
            }
            used[next] = true;
        } else {
            break;
        }
    }
    // TODO: optimize
    path
}

fn optimize_cycle(cycle: &mut Vec<Point>) {
    loop {
        let mut changed = false;
        for i in 0..cycle.len() {
            let prev = cycle[(i + cycle.len() - 1) % cycle.len()];
            let cur = cycle[i];
            let next = cycle[(i + 1) % cycle.len()];
            let d_prev = cur.dist2(&prev);
            let d_next = cur.dist2(&next);
            let d_prev_next = prev.dist2(&next);
            if d_prev + d_next > d_prev_next {
                let mut bi = i;
                let mut bsum = d_prev + d_next;
                for j in 0..cycle.len() {
                    if j + 1 == i {
                        continue;
                    }
                    let p1 = cycle[j];
                    let p2 = cycle[(j + 1) % cycle.len()];
                    let check_dist = cur.dist2(&p1) + cur.dist2(&p2);
                    if check_dist < bsum {
                        bsum = check_dist;
                        bi = j + 1;
                    }
                }
                if bi != i {
                    changed = true;
                    cycle.remove(i);
                    if bi > i {
                        bi -= 1;
                    }
                    cycle.insert(bi, cur);
                }
            }
        }

        if !changed {
            break;
        }
    }
}

fn gen_mask(color_image: &ColorImage) -> ColorImage {
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

    let mut figures = vec![vec![]; width * height];
    for x in 0..width {
        for y in 0..height {
            if is_puzzle[id(x, y)] && dist_to_not_puzzle[id(x, y)] == 1 {
                figures[dsu.get(id(x, y))].push(Point { x, y });
            }
        }
    }

    let mut use_color = vec![Color32::WHITE; width * height];
    let mut rng = rand::thread_rng();
    for i in 0..use_color.len() {
        use_color[i] = Color32::from_rgb(rng.gen(), rng.gen(), rng.gen());
    }

    let mut res = ColorImage::new(color_image.size, Color32::WHITE);
    for x in 0..width {
        for y in 0..height {
            if is_puzzle[id(x, y)] {
                if dist_to_not_puzzle[id(x, y)] == 1 {
                    res[(x, y)] = Color32::BLACK;
                } else {
                    res[(x, y)] = use_color[dsu.get(id(x, y))];
                }
            }
        }
    }
    for i in 0..figures.len() {
        if figures[i].len() < 20 {
            continue;
        }
        let mut cycle = find_cycle(&figures[i]);
        let max_dist2 = cycle.windows(2).map(|w| w[0].dist2(&w[1])).max().unwrap();
        // optimize_cycle(&mut cycle);
        let max_dist2_after_opt = cycle.windows(2).map(|w| w[0].dist2(&w[1])).max().unwrap();
        eprintln!(
            "figure len: {}, max dist2: {}, after_opt : {}",
            figures[i].len(),
            max_dist2,
            max_dist2_after_opt
        );
        if max_dist2 > 10 {
            for j in 0..cycle.len() {
                res[(cycle[j].x, cycle[j].y)] =
                    mid_color(j, cycle.len(), Color32::RED, Color32::RED);
            }
        }
    }
    res
}

fn mid_color(pos: usize, len: usize, start: Color32, end: Color32) -> Color32 {
    let calc =
        |l: u8, r: u8| -> u8 { (((l as usize) * (len - pos - 1) + (r as usize) * pos) / len) as _ };
    let r = calc(start.r(), end.r());
    let g = calc(start.g(), end.g());
    let b = calc(start.b(), end.b());
    Color32::from_rgb(r, g, b)
}

impl MyWidget {
    pub fn new(path: &str) -> Self {
        let color_image2 = load_image_from_path(path).unwrap();
        let mask_image = {
            let color = gen_mask(&color_image2);
            save_color_image(&color, "img/puzzle.jpg");
            RetainedImage::from_color_image("mask", color)
        };

        let image = RetainedImage::from_color_image("test", color_image2);
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
