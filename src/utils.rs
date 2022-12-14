use std::cmp::{max, min};

use eframe::epaint::ColorImage;
use image::ImageBuffer;
use itertools::Itertools;

pub fn load_image_from_path(path: &str) -> Result<ColorImage, image::ImageError> {
    let image = image::io::Reader::open(path)?.decode()?;
    let size = [image.width() as _, image.height() as _];
    let image_buffer = image.to_rgba8();
    let pixels = image_buffer.as_flat_samples();
    Ok(ColorImage::from_rgba_unmultiplied(size, pixels.as_slice()))
}

// TODO: move to a separate file
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Side {
    pub fig: usize,
    pub side: usize,
}

fn norm_side(x: usize) -> usize {
    if x >= 4 {
        x - 4
    } else {
        x
    }
}

impl Side {
    pub fn ne(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 1),
        }
    }

    pub fn ne2(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 2),
        }
    }

    pub fn pr(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 3),
        }
    }
}

pub fn fmax(x: f64, y: f64) -> f64 {
    if x > y {
        x
    } else {
        y
    }
}

pub fn fmin(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

pub fn save_color_image(color_image: &ColorImage, path: &str) {
    let mut new_img = ImageBuffer::new(color_image.size[0] as u32, color_image.size[1] as u32);
    for x in 0..new_img.width() {
        for y in 0..new_img.height() {
            let color = color_image[(x as usize, y as usize)];
            new_img[(x, y)] = image::Rgba([color.r(), color.g(), color.b(), color.a()]);
        }
    }
    new_img.save(path).unwrap();
}

pub fn dedup_edges(used_edges: &[(Side, Side)]) -> Vec<(Side, Side)> {
    let mut used_edges = used_edges
        .iter()
        .map(|&(s1, s2)| if s1 < s2 { (s1, s2) } else { (s2, s1) })
        .collect_vec();
    used_edges.sort();
    used_edges.dedup();
    used_edges
}

pub fn gauss(matrix: &mut [Vec<f64>]) -> Vec<f64> {
    for c in 0..matrix[0].len() - 1 {
        let mut best_r = c;
        for r in c..matrix.len() {
            if matrix[r][c].abs() > matrix[best_r][c].abs() {
                best_r = r;
            }
        }
        matrix.swap(c, best_r);
        if matrix[c][c] == 0.0 {
            continue;
        }
        let mul = 1.0 / matrix[c][c];
        for x in matrix[c].iter_mut() {
            *x *= mul;
        }
        for r2 in 0..matrix.len() {
            if r2 == c {
                continue;
            }
            let mul = matrix[r2][c];
            for c2 in 0..matrix[r2].len() {
                matrix[r2][c2] -= mul * matrix[c][c2];
            }
        }
    }
    let mut res = vec![0.0; matrix.len()];
    for i in 0..res.len() {
        res[i] = -matrix[i].last().unwrap();
    }
    res
}

pub fn normalize_bounding_box((x, y): (i32, i32)) -> (i32, i32) {
    (min(x, y), max(x, y))
}
