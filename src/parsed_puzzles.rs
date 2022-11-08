use eframe::epaint::{Color32, ColorImage};
use rand::Rng;

use crate::{dsu::Dsu, figure::Figure, point::Point};

pub struct ParsedPuzzles {
    pub width: usize,
    pub height: usize,
    pub figures: Vec<Figure>,
}

fn is_puzzle_color(color: Color32) -> bool {
    (color.r() as usize) + (color.g() as usize) + (color.b() as usize) >= 630
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
