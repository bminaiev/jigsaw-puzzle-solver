use std::cmp::min;

use eframe::epaint::{Color32, ColorImage};

pub struct AverareColor {
    width: usize,
    height: usize,
    pref_sum: Vec<[i64; 3]>,
}

impl AverareColor {
    fn id(&self, x: usize, y: usize) -> usize {
        x + (self.width + 1) * y
    }

    pub fn new(image: &ColorImage) -> Self {
        let size = image.size;
        let pref_sum = vec![[0; 3]; (size[0] + 1) * (size[1] + 1)];
        let mut res = Self {
            width: size[0],
            height: size[1],
            pref_sum,
        };
        for x in 0..size[0] {
            for y in 0..size[1] {
                let color = image[(x, y)];
                let id = res.id(x + 1, y + 1);
                res.pref_sum[id][0] += color.r() as i64;
                res.pref_sum[id][1] += color.g() as i64;
                res.pref_sum[id][2] += color.b() as i64;
                for idx in 0..3 {
                    assert!(res.pref_sum[id][idx] <= 255);
                    res.pref_sum[id][idx] -= res.pref_sum[res.id(x, y)][idx];
                    res.pref_sum[id][idx] += res.pref_sum[res.id(x + 1, y)][idx];
                    res.pref_sum[id][idx] += res.pref_sum[res.id(x, y + 1)][idx];
                }
            }
        }
        res
    }

    pub fn get_averare_color(&self, x: usize, y: usize, radius: usize) -> Color32 {
        let x_start = if x > radius { x - radius } else { 0 };
        let y_start = if y > radius { y - radius } else { 0 };
        let x_end = min(self.width, x + radius + 1);
        let y_end = min(self.height, y + radius + 1);
        let id00 = self.id(x_start, y_start);
        let id01 = self.id(x_start, y_end);
        let id10 = self.id(x_end, y_start);
        let id11 = self.id(x_end, y_end);
        let mut sum = [0; 3];
        let total_size = ((x_end - x_start) * (y_end - y_start)) as i64;
        for idx in 0..3 {
            sum[idx] += self.pref_sum[id11][idx];
            sum[idx] += self.pref_sum[id00][idx];
            sum[idx] -= self.pref_sum[id01][idx];
            sum[idx] -= self.pref_sum[id10][idx];
            sum[idx] /= total_size;
            assert!(
                sum[idx] <= 255,
                "total_size: {total_size}, av: {}",
                sum[idx]
            );
        }
        Color32::from_rgb(sum[0] as u8, sum[1] as u8, sum[2] as u8)
    }
}
