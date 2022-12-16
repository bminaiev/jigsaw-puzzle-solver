use std::cmp::min;

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
use rand::{rngs::ThreadRng, Rng};

use crate::{
    border_matcher::{match_borders, MatchResult},
    crop::crop,
    dsu::Dsu,
    figure::Figure,
    graph_solver::PotentialSolution,
    known_facts::{self, EdgeState, Fact, KnownFacts},
    parsed_puzzles::ParsedPuzzles,
    point::{Point, PointF},
    utils::{load_image_from_path, save_color_image, Side},
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
    matched_borders: Vec<Vec<MatchResult>>,
    solutions: Vec<PotentialSolution>,
    show_parsed: bool,
    show_image: bool,
    fig_colors: Vec<Color32>,
    show_matched_borders: bool,
    piece_picker: String,
    crop_enabled: bool,
    selected_solution: Option<usize>,
    new_edges: Vec<EdgeState>,
    known_facts: KnownFacts,
}

const ZOOM_DELTA_COEF: f32 = 500.0;

fn gen_good_color(rng: &mut ThreadRng) -> Color32 {
    let parts: [u8; 3] = [rng.gen(), rng.gen(), rng.gen()];
    let mut parts_sorted = parts.clone();
    parts_sorted.sort();
    if parts_sorted[2] - parts_sorted[0] < 100 {
        // avoid too grey colors
        return gen_good_color(rng);
    }
    Color32::from_rgb(parts[0], parts[1], parts[2])
}

impl MyWidget {
    pub fn new(
        path: &str,
        solutions: Vec<PotentialSolution>,
        show_parsed: bool,
        show_image: bool,
        show_matched_borders: bool,
        crop_enabled: bool,
        known_facts: KnownFacts,
    ) -> Self {
        let color_image = load_image_from_path(path).unwrap();

        let parsed_puzzles = ParsedPuzzles::new(&color_image);
        let mask_image = if !show_parsed {
            RetainedImage::from_color_image("test", color_image.clone())
        } else {
            let color = parsed_puzzles.gen_image();
            save_color_image(&color, "img/puzzle.jpg");
            RetainedImage::from_color_image("mask", color)
        };

        let image = RetainedImage::from_color_image("test", color_image.clone());
        let img_size = image.size_vec2();
        let mut rng = rand::thread_rng();
        let fig_colors = (0..parsed_puzzles.figures.len())
            .map(|_| gen_good_color(&mut rng))
            .collect_vec();
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
            matched_borders: vec![vec![]; 4],
            solutions,
            show_parsed,
            fig_colors,
            show_image,
            show_matched_borders,
            piece_picker: String::new(),
            crop_enabled,
            selected_solution: None,
            new_edges: vec![],
            known_facts,
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
                if mouse.distance(self.convert_to_screen(p.pos2())) <= 12.0 {
                    return Some(closest_figure_id);
                }
            }
        }
        None
    }

    fn show_border(&self, figure: &Figure, ui: &mut eframe::egui::Ui) {
        for i in 0..figure.border.len() {
            let p1 = self.convert_to_screen(figure.border[i].pos2());

            let color = Color32::GREEN;
            ui.painter().circle_filled(p1, 1.5, color);
            let p2 = self.convert_to_screen(figure.border[(i + 1) % figure.border.len()].pos2());
            ui.painter()
                .add(Shape::line_segment([p1, p2], Stroke::new(1.0, color)));
        }
    }

    fn show_white_background(&mut self, ui: &mut eframe::egui::Ui) {
        let img_size = self.image.size_vec2();
        let min = self.convert_to_screen(pos2(0.0, 0.0));
        let max = self.convert_to_screen(pos2(img_size.x, img_size.y));
        let rect = Rect::from_min_max(min, max);
        ui.painter()
            .rect_filled(rect, Rounding::default(), Color32::WHITE);
    }

    fn show_parsed(&mut self, ui: &mut eframe::egui::Ui) {
        {
            let img_size = self.image.size_vec2();
            let min = self.convert_to_screen(pos2(0.0, 0.0));
            let max = self.convert_to_screen(pos2(img_size.x, img_size.y));
            let rect = Rect::from_min_max(min, max);
            let texture_id = self.image.texture_id(&ui.ctx());

            let uv = Rect::from_min_max(pos2(0.0, 0.0), pos2(1.0, 1.0));

            if self.show_image {
                let mut mesh = Mesh::with_texture(texture_id);

                mesh.add_rect_with_uv(rect, uv, Color32::WHITE);
                ui.painter().add(Shape::mesh(mesh));
            } else {
            }

            {
                let texture_id2 = self.mask_image.texture_id(&ui.ctx());
                let mut mesh2 = Mesh::with_texture(texture_id2);
                mesh2.add_rect_with_uv(rect, uv, Color32::WHITE);
                ui.painter().add(Shape::mesh(mesh2));
            }
        }

        if !self.show_parsed {
            return;
        }

        if let Some(figure_id) = self.which_figure_hovered(ui) {
            let figure = &self.parsed_puzzles.figures[figure_id];
            self.show_border(figure, ui);
            if ui.input().pointer.primary_clicked() && self.show_matched_borders {
                for i in 0..4 {
                    self.matched_borders[i].clear();
                    if self.parsed_puzzles.figures[figure_id].is_good_puzzle() {
                        for other_figure_id in 0..self.parsed_puzzles.figures.len() {
                            if other_figure_id == figure_id {
                                continue;
                            }
                            let figure = &self.parsed_puzzles.figures[figure_id];
                            let other_figure = &self.parsed_puzzles.figures[other_figure_id];
                            if other_figure.is_good_puzzle() {
                                for j in 0..4 {
                                    if let Some(result) = match_borders(
                                        &self.parsed_puzzles,
                                        Side {
                                            fig: figure_id,
                                            side: i,
                                        },
                                        Side {
                                            fig: other_figure_id,
                                            side: j,
                                        },
                                    ) {
                                        self.matched_borders[i].push(result);
                                    }
                                }
                            }
                        }
                        self.matched_borders[i].sort_by(|a, b| a.score.total_cmp(&b.score));
                    }
                }
            }
        }

        if let Ok(id) = self.piece_picker.parse::<usize>() {
            if id < self.parsed_puzzles.figures.len() {
                self.show_border(&self.parsed_puzzles.figures[id], ui);
            }
        }

        if let Some(sol_id) = self.selected_solution {
            for &new_fig in self.solutions[sol_id].new_figures_used.iter() {
                self.show_border(&self.parsed_puzzles.figures[new_fig], ui);
            }
        }
    }

    fn show_best_matched_borders(&mut self, ui: &mut eframe::egui::Ui) {
        if !self.show_matched_borders {
            return;
        }

        let font_id = FontId::new(10.0, FontFamily::Monospace);
        let x_offset = self.image.size()[0] as f32 + 100.0;
        const MAX_OPTIONS: usize = 20;
        const ONE_OPTION_SIZE: f32 = 250.0;

        {
            let min = self.convert_to_screen(pos2(x_offset, 0.0));
            let max = self.convert_to_screen(pos2(
                x_offset + 1000.0,
                (MAX_OPTIONS as f32) * ONE_OPTION_SIZE,
            ));
            ui.painter().rect_filled(
                Rect::from_min_max(min, max),
                Rounding::default(),
                Color32::WHITE,
            );
        }

        for i in 0..4 {
            for j in 0..min(MAX_OPTIONS, self.matched_borders[i].len()) {
                let result = &self.matched_borders[i][j];

                let offset = result.get_offset();

                let conv_point = |p: PointF| -> Pos2 {
                    let p = p + offset;
                    pos2(
                        p.x as f32 + x_offset + (i as f32) * ONE_OPTION_SIZE + 20.0,
                        p.y as f32 + (j as f32) * ONE_OPTION_SIZE + 20.0,
                    )
                };

                let draw = |pts: &[PointF], color: Color32| {
                    for (p1, p2) in pts.iter().circular_tuple_windows() {
                        let p1 = conv_point(*p1);
                        let p2 = conv_point(*p2);
                        ui.painter().line_segment(
                            [self.convert_to_screen(p1), self.convert_to_screen(p2)],
                            Stroke::new(1.0, color),
                        );
                        ui.painter()
                            .circle_filled(self.convert_to_screen(p1), 3.0, color);
                    }
                };

                draw(&result.lhs, Color32::BLUE);
                draw(&result.rhs, Color32::RED);
                ui.painter().text(
                    self.convert_to_screen(conv_point(result.lhs_center)),
                    Align2::CENTER_CENTER,
                    result.lhs_id.to_string(),
                    font_id.clone(),
                    Color32::BLUE,
                );
                ui.painter().text(
                    self.convert_to_screen(conv_point(result.rhs_center)),
                    Align2::CENTER_CENTER,
                    result.rhs_id.to_string(),
                    font_id.clone(),
                    Color32::RED,
                );
                ui.painter().text(
                    self.convert_to_screen(conv_point(PointF::ZERO - offset)),
                    Align2::LEFT_TOP,
                    format!("score = {:.3}", result.score),
                    font_id.clone(),
                    Color32::BLACK,
                );
                // eprintln!("{}x{} -> {}", i, j, result.score);
            }
        }
    }

    fn show_solution(&mut self, ui: &mut eframe::egui::Ui) {
        if self.solutions.is_empty() {
            return;
        }

        let offset = vec2(0.0, self.image.size()[1] as f32 + 100.0);

        {
            let img_size = self.image.size_vec2();
            let min = self.convert_to_screen(pos2(0.0, 0.0) + offset);
            let max = self.convert_to_screen(pos2(img_size.x, img_size.y) + offset);
            let rect = Rect::from_min_max(min, max);
            ui.painter()
                .rect_filled(rect, Rounding::default(), Color32::WHITE);
        }

        for (sol_id, sol) in self.solutions.iter().enumerate() {
            for fig in sol.placed_figures.iter() {
                let color = self.fig_colors[fig.figure_id];
                for (p1, p2) in fig.positions.iter().circular_tuple_windows() {
                    let p1 = self.convert_to_screen(p1.pos2() + offset);
                    let p2 = self.convert_to_screen(p2.pos2() + offset);

                    ui.painter().line_segment([p1, p2], Stroke::new(1.5, color))
                }

                for &pos in self.parsed_puzzles.figures[fig.figure_id]
                    .corner_positions
                    .iter()
                {
                    let p = self.convert_to_screen(fig.positions[pos].pos2() + offset);
                    ui.painter().add(Shape::circle_filled(p, 2.0, color));
                }

                let center = self.convert_to_screen(calc_center(&fig.positions).pos2() + offset);
                ui.painter().text(
                    center,
                    Align2::CENTER_CENTER,
                    fig.figure_id.to_string(),
                    FontId::new(10.0, FontFamily::Monospace),
                    Color32::BLACK,
                );
            }
            for debug_line in sol.debug_lines.iter() {
                for w in debug_line.windows(2) {
                    let p1 = self.convert_to_screen(w[0].pos2() + offset);
                    let p2 = self.convert_to_screen(w[1].pos2() + offset);
                    ui.painter()
                        .line_segment([p1, p2], Stroke::new(1.5, Color32::BLACK));
                }
            }
            {
                let p0 = self.convert_to_screen(sol.bbox.0.pos2() + offset);
                let p1 = self.convert_to_screen(sol.bbox.1.pos2() + offset);
                let is_mouse_inside = if let Some(mouse) = ui.input().pointer.hover_pos() {
                    mouse.x >= p0.x && mouse.x <= p1.x && mouse.y >= p0.y && mouse.y <= p1.y
                } else {
                    false
                };
                if is_mouse_inside
                    && ui.input().pointer.primary_clicked()
                    && self.selected_solution != Some(sol_id)
                {
                    self.selected_solution = Some(sol_id);
                    self.new_edges.clear();
                    for &(s1, s2) in sol.new_edges_used.iter() {
                        self.new_edges.push(self.known_facts.get_edge_state(s1, s2));
                    }
                }
                let selected = self.selected_solution == Some(sol_id);
                if is_mouse_inside || selected {
                    let color = if selected {
                        Color32::BLUE
                    } else {
                        Color32::BLACK
                    };
                    ui.painter().rect_stroke(
                        Rect::from_min_max(p0, p1),
                        Rounding::default(),
                        Stroke::new(3.0, color),
                    )
                }
            }
            ui.painter().text(
                self.convert_to_screen(sol.text_offset.pos2() + offset + vec2(10.0, 10.0)),
                Align2::LEFT_TOP,
                format!("{:.3}{}", sol.placement_score, sol.additional_text),
                FontId::new(10.0, FontFamily::Monospace),
                Color32::BLACK,
            );
        }
    }

    fn show_elements(&mut self, ui: &mut eframe::egui::Ui) {
        ui.text_edit_singleline(&mut self.piece_picker);
        if let Some(sol_id) = self.selected_solution {
            for i in 0..self.new_edges.len() {
                ui.horizontal(|ui| {
                    let (s1, s2) = self.solutions[sol_id].new_edges_used[i];
                    ui.label(format!("{} - {}", s1.fig, s2.fig));
                    if ui
                        .radio_value(&mut self.new_edges[i], EdgeState::Unknown, "Unknown")
                        .clicked()
                    {
                        self.known_facts.remove_fact(&Fact::new(s1, s2, false));
                        self.known_facts.remove_fact(&Fact::new(s1, s2, true));
                    }
                    if ui
                        .radio_value(&mut self.new_edges[i], EdgeState::GoodEdge, "Good")
                        .clicked()
                    {
                        self.known_facts.add_fact(&Fact::new(s1, s2, true));
                    }
                    if ui
                        .radio_value(&mut self.new_edges[i], EdgeState::WrongEdge, "Wrong")
                        .clicked()
                    {
                        self.known_facts.add_fact(&Fact::new(s1, s2, false));
                    }
                });
            }
        }
    }

    pub fn ui(&mut self, ui: &mut eframe::egui::Ui) -> eframe::egui::Response {
        self.show_elements(ui);

        let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::drag());
        {
            self.change_zoom(ui);
        }

        if !self.show_image {
            self.show_white_background(ui);
        }
        self.show_parsed(ui);

        if self.show_image {
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

        self.show_best_matched_borders(ui);

        if self.show_parsed {
            for (figure_id, figure) in self.parsed_puzzles.figures.iter().enumerate() {
                if figure.good_border {
                    for &pos in figure.corner_positions.iter() {
                        let p = self.convert_to_screen(figure.border[pos].pos2());
                        ui.painter().add(Shape::circle_filled(p, 2.0, Color32::RED));
                    }

                    let center = self.convert_to_screen(figure.center.pos2());
                    let color = Color32::BLACK;
                    ui.painter().text(
                        center,
                        Align2::CENTER_CENTER,
                        figure_id.to_string(),
                        FontId::new(10.0, FontFamily::Monospace),
                        color,
                    );
                }
            }
        }

        self.show_solution(ui);

        let frame = self.frame.clone();
        if self.crop_enabled {
            if ui.input().key_pressed(Key::Enter) {
                dbg!("CROP!", &frame);
                crop(&self.image_path, &frame);
                dbg!("CROPED!");
                std::process::exit(0);
            }
        }

        response
    }
}

fn calc_center(positions: &[PointF]) -> PointF {
    let mut res = PointF::ZERO;
    for p in positions.iter() {
        res = res + *p;
    }
    res / (positions.len() as f64)
}
