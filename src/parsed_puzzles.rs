use std::{
    cmp::min,
    collections::{hash_map::DefaultHasher, BTreeMap},
    hash::{Hash, Hasher},
    ops::Range,
    process,
};

use eframe::epaint::{Color32, ColorImage};
use itertools::Itertools;
use rand::Rng;

use crate::{
    average_color::{self, AverareColor},
    border_matcher::is_picture_border,
    dsu::Dsu,
    figure::{BorderFigure, Figure},
    point::Point,
    utils::{save_color_image, Side},
    PUZZLE_PIXEL_WHITE_THRESHOLD,
};

#[derive(Hash)]
pub struct ParsedPuzzles {
    pub width: usize,
    pub height: usize,
    pub figures: Vec<Figure>,
}

fn get_color_score(color: Color32) -> usize {
    (color.r() as usize) + (color.g() as usize) + (color.b() as usize)
}

fn is_puzzle_color(color: Color32) -> bool {
    get_color_score(color) >= PUZZLE_PIXEL_WHITE_THRESHOLD
}

fn split_into_figures(
    pts: &[Point],
    figure_size_limit: usize,
    color_image: &ColorImage,
) -> Vec<Figure> {
    if pts.len() <= figure_size_limit {
        if let Some(figure) = Figure::new(pts) {
            return vec![figure];
        } else {
            return vec![];
        }
    }
    let mut id_by_p = BTreeMap::new();
    for (id, p) in pts.iter().enumerate() {
        id_by_p.insert(p.clone(), id);
    }
    let mut dsu = Dsu::new(pts.len());
    let color_score = pts
        .iter()
        .map(|p| {
            let mut res = get_color_score(color_image[(p.x, p.y)]);
            for p2 in p.neighbours() {
                if p2.x < color_image.size[0] && p2.y < color_image.size[1] {
                    res = min(res, get_color_score(color_image[(p2.x, p2.y)]));
                }
            }
            res
        })
        .collect_vec();
    let mut queue = (0..pts.len()).collect_vec();
    queue.sort_by_key(|id| color_score[*id]);
    queue.reverse();
    for &id in queue.iter() {
        let p = pts[id];
        for p2 in p.neighbours() {
            if let Some(&id2) = id_by_p.get(&p2) {
                let dsu1 = dsu.get(id);
                let dsu2 = dsu.get(id2);
                if dsu1 != dsu2 && dsu.get_size(dsu1) + dsu.get_size(dsu2) < figure_size_limit {
                    dsu.unite(dsu1, dsu2);
                }
            }
        }
    }
    let mut res = vec![];
    for fig in dsu.get_components() {
        assert!(fig.len() <= figure_size_limit);
        let fig_pts = fig.iter().map(|&id| pts[id]).collect_vec();
        if let Some(fig) = Figure::new(&fig_pts) {
            res.push(fig);
        }
    }
    return res;
}

fn get_sum(color: Color32) -> i32 {
    (color.r() as i32) + (color.g() as i32) + (color.b() as i32)
}

fn is_puzzle_pixel(
    color_image: &ColorImage,
    average_color: &AverareColor,
    r: usize,
    offset: i32,
    x: usize,
    y: usize,
) -> bool {
    let c1 = color_image[(x, y)];
    let c_av = average_color.get_averare_color(x, y, r);

    get_sum(c1) > get_sum(c_av) + offset
}

fn save_debug(color_image: &ColorImage, average_color: &AverareColor, r: usize, offset: i32) {
    let mut res = ColorImage::new(color_image.size, Color32::WHITE);

    for x in 0..color_image.size[0] {
        for y in 0..color_image.size[1] {
            let c1 = color_image[(x, y)];
            let c_av = average_color.get_averare_color(x, y, r);

            if get_sum(c1) > get_sum(c_av) + offset {
                res[(x, y)] = c1;
            } else {
                res[(x, y)] = Color32::BLACK;
            }
        }
    }

    save_color_image(&res, &format!("img/debug/debug-r={r},o={offset}.jpg"));
}

fn good_size_range(figures: &[Figure]) -> Range<usize> {
    let mut sizes = vec![];
    for figure in figures.iter() {
        if figure.is_good_puzzle() {
            sizes.push(figure.all_pts.len());
        }
    }
    sizes.sort();
    let med = sizes[sizes.len() / 2];
    med * 2 / 3..med * 4 / 3
}

impl ParsedPuzzles {
    pub fn new(color_image: &ColorImage) -> Self {
        let average_color = AverareColor::new(color_image);
        const R: usize = 200;
        const OFFSET: i32 = 200;

        // for &r in [25, 50, 100, 200].iter() {
        //     for &offset in [20, 50, 100, 200].iter() {
        //         save_debug(color_image, &average_color, r, offset);
        //     }
        // }

        // process::exit(0);

        let width = color_image.size[0];
        let height = color_image.size[1];

        let id = |x: usize, y: usize| -> usize { x + y * width };
        let mut is_puzzle = vec![false; width * height];

        for x in 0..width {
            for y in 0..height {
                if is_puzzle_pixel(color_image, &average_color, R, OFFSET, x, y) {
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

        eprintln!("Start parsing figures...");
        let mut figure_sizes = vec![];
        for fig in dsu_figures.iter() {
            if fig.len() >= 333 {
                figure_sizes.push(fig.len());
            }
        }
        figure_sizes.sort();
        let med_figure_size = if figure_sizes.is_empty() {
            0
        } else {
            figure_sizes[figure_sizes.len() / 2]
        };
        let figure_size_limit = med_figure_size * 3 / 2;
        eprintln!("figure size limit: {}", figure_size_limit);

        for i in 0..dsu_figures.len() {
            res_figures.extend(split_into_figures(
                &dsu_figures[i],
                figure_size_limit,
                color_image,
            ));
        }

        eprintln!("Found {} figures", res_figures.len());
        let good_sizes = good_size_range(&res_figures);
        for fig in res_figures.iter_mut() {
            if !good_sizes.contains(&fig.all_pts.len()) {
                fig.good_size = false;
            }
        }
        let good_figures = res_figures.iter().filter(|f| f.is_good_puzzle()).count();
        eprintln!("Good figures: {good_figures}");

        res_figures.sort_by_key(|f| (f.center.y, f.center.x));

        Self {
            width,
            height,
            figures: res_figures,
        }
    }

    pub fn gen_image(&self) -> ColorImage {
        let mut res = ColorImage::new([self.width, self.height], Color32::TRANSPARENT);

        let mut rng = rand::thread_rng();
        for figure in self.figures.iter() {
            let color = Color32::from_rgb(rng.gen(), rng.gen(), rng.gen());
            for p in figure.all_pts.iter() {
                if p.x < self.width && p.y < self.height {
                    // res[(p.x, p.y)] = color;
                }
            }
            let border_color = if figure.good_border {
                Color32::BLUE
            } else {
                Color32::RED
            };
            for p in figure.border.iter() {
                if p.x < self.width && p.y < self.height {
                    res[(p.x, p.y)] = border_color;
                }
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

    pub fn calc_figures_on_border(&self) -> Vec<BorderFigure> {
        (0..self.figures.len())
            .filter_map(|id| {
                let figure = &self.figures[id];
                if !figure.is_good_puzzle() {
                    return None;
                }
                let is_border = (0..4)
                    .map(|border_id| is_picture_border(figure, border_id))
                    .collect_vec();
                for i in 0..4 {
                    if is_border[i] && is_border[(i + 1) % 4] {
                        return Some(BorderFigure::new(id, i + 2, i + 3));
                    }
                }
                for i in 0..4 {
                    if is_border[i] {
                        return Some(BorderFigure::new(id, i + 1, i + 3));
                    }
                }
                None
            })
            .collect_vec()
    }

    pub fn gen_all_sides(&self) -> Vec<Side> {
        let mut res = vec![];
        for fig in 0..self.figures.len() {
            if self.figures[fig].is_good_puzzle() {
                for side in 0..4 {
                    res.push(Side { fig, side });
                }
            }
        }
        res
    }
}
