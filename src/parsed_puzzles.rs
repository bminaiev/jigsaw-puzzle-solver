use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
};

use eframe::epaint::{Color32, ColorImage};
use rand::Rng;

use crate::{dsu::Dsu, figure::Figure, point::Point};

#[derive(Hash)]
pub struct ParsedPuzzles {
    pub width: usize,
    pub height: usize,
    pub figures: Vec<Figure>,
}

fn is_puzzle_color(color: Color32) -> bool {
    (color.r() as usize) + (color.g() as usize) + (color.b() as usize) >= 630
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

        let mut dsu_figures = vec![vec![]; width * height];
        for x in 0..width {
            for y in 0..height {
                if is_puzzle[id(x, y)] {
                    dsu_figures[dsu.get(id(x, y))].push(Point { x, y });
                }
            }
        }

        let mut res_figures = vec![];

        for i in 0..dsu_figures.len() {
            if let Some(figure) = Figure::new(&dsu_figures[i]) {
                res_figures.push(figure);
            }
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

    pub fn calc_hash(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }
}
