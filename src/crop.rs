use eframe::epaint::{pos2, vec2, Pos2};
use image::{GenericImageView, ImageBuffer};

use crate::utils::gauss;

#[derive(Clone, Copy)]
struct Line {
    a: f32,
    b: f32,
    c: f32,
}

impl Line {
    pub fn new(p1: Pos2, p2: Pos2) -> Self {
        let a = p2.y - p1.y;
        let b = p1.x - p2.x;
        let c = -(p1.x * a + p1.y * b);
        Self { a, b, c }
    }

    fn dist(&self, p: Pos2) -> f32 {
        (self.a * p.x + self.b * p.y + self.c) / (self.a * self.a + self.b * self.b).sqrt()
    }
}

fn intersection(l1: Line, l2: Line) -> Pos2 {
    let denom = l1.a * l2.b - l2.a * l1.b;
    let x = l2.c * l1.b - l1.c * l2.b;
    let y = -(l2.c * l1.a - l1.c * l2.a);
    pos2(x / denom, y / denom)
}

fn vec_mul(p1: Pos2, p2: Pos2, p3: Pos2) -> f32 {
    (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
}

fn fmin(a: f32, b: f32) -> f32 {
    if a < b {
        a
    } else {
        b
    }
}

fn fmax(a: f32, b: f32) -> f32 {
    if a > b {
        a
    } else {
        b
    }
}

fn inside(mid: Pos2, a: Pos2, b: Pos2) -> bool {
    mid.x >= fmin(a.x, b.x)
        && mid.x <= fmax(a.x, b.x)
        && mid.y >= fmin(a.y, b.y)
        && mid.y <= fmax(a.y, b.y)
}

fn cut(pts: &[Pos2], fr: Pos2, to: Pos2) -> Vec<Pos2> {
    let mut res = vec![];
    for i in 0..pts.len() {
        let p1 = pts[i];
        let p2 = pts[(i + 1) % pts.len()];
        let v_cur = vec_mul(fr, to, p1);
        let v_next = vec_mul(fr, to, p2);
        if v_cur >= 0.0 {
            res.push(p1);
        }
        if v_cur != 0.0 && v_next != 0.0 && v_cur.signum() != v_next.signum() {
            let new_p = intersection(Line::new(p1, p2), Line::new(fr, to));
            // assert!(inside(new_p, p1, p2));
            res.push(new_p);
        }
    }

    res
}

fn area(pts: &[Pos2]) -> f32 {
    let mut res = 0.0;
    for i in 0..pts.len() {
        let p1 = pts[i];
        let p2 = pts[(i + 1) % pts.len()];
        res += vec_mul(pts[0], p1, p2);
    }
    res.abs()
}

fn calc_intersection_area(mut pts: Vec<Pos2>, x: f32, y: f32, debug: bool) -> f32 {
    let p1 = pos2(x, y);
    let p2 = pos2(x + 1.0, y);
    let p3 = pos2(x + 1.0, y + 1.0);
    let p4 = pos2(x, y + 1.0);
    for w in [p1, p2, p3, p4, p1].windows(2) {
        pts = cut(&pts, w[0], w[1]);
        if debug {
            dbg!(&pts);
        }
    }
    area(&pts)
}

struct TransofmationMatrix {
    a: Vec<Vec<f64>>,
}

impl TransofmationMatrix {
    pub fn from_gauss_result(gauss_res: &[f64]) -> Self {
        let mut a = vec![vec![1.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let id = i * 3 + j;
                if id < 8 {
                    a[i][j] = gauss_res[id];
                }
            }
        }
        Self { a }
    }

    fn transform_point(&self, p: Pos2) -> Pos2 {
        let mult = [p.x as f64, p.y as f64, 1.0];
        let mut res = [0.0; 3];
        for row in 0..3 {
            for col in 0..3 {
                res[row] += self.a[row][col] * mult[col];
            }
        }
        pos2((res[0] / res[2]) as f32, (res[1] / res[2]) as f32)
    }
}

fn find_transformation_matrix(frame: &[Pos2], new_width: f32) -> TransofmationMatrix {
    let sz = 13;
    let mut a = vec![vec![vec![0.0; sz]; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let id = i * 3 + j;
            if id == 8 {
                a[i][j][sz - 1] = 1.0;
            } else {
                a[i][j][id] = 1.0;
            }
        }
    }
    let target_points = [
        pos2(0.0, 0.0),
        pos2(new_width, 0.0),
        pos2(new_width, new_width),
        pos2(0.0, new_width),
    ];
    let mut matrix = vec![];
    for (pt_id, (my, target)) in frame.iter().zip(target_points.iter()).enumerate() {
        for row in 0..3 {
            let mut new_row = vec![0.0; sz];
            let b = [my.x as f64, my.y as f64, 1.0];

            for col in 0..3 {
                for t in 0..sz {
                    new_row[t] += a[row][col][t] * b[col];
                }
            }

            let k_coef = if row == 0 {
                target.x as f64
            } else if row == 1 {
                target.y as f64
            } else {
                1.0
            };
            new_row[8 + pt_id] -= k_coef;
            matrix.push(new_row);
        }
    }
    let gauss_res = gauss(&mut matrix);
    let tranformation_matrix = TransofmationMatrix::from_gauss_result(&gauss_res);
    for (my, target) in frame.iter().zip(target_points.iter()) {
        let new_target = tranformation_matrix.transform_point(*my);
        eprintln!("target = {:?}, got = {:?}", target, new_target);
    }
    tranformation_matrix
}

pub fn crop(path: &str, frame: &[Pos2]) {
    let image = image::open(path).unwrap();
    eprintln!("crop points: {:?}", frame);

    const NEW_WIDTH: usize = 2000;

    let matrix = find_transformation_matrix(frame, NEW_WIDTH as f32);

    let w = image.width();
    let h = image.height();

    let mut new_img_comps = vec![vec![[0.0; 3]; NEW_WIDTH]; NEW_WIDTH];
    let mut new_img_area = vec![vec![0.0; NEW_WIDTH]; NEW_WIDTH];

    for x in 0..w {
        for y in 0..h {
            let pixel = image.get_pixel(x, y);

            let p = pos2(x as f32, y as f32);
            let p1 = matrix.transform_point(p);
            let p2 = matrix.transform_point(p + vec2(1.0, 0.0));
            let p3 = matrix.transform_point(p + vec2(1.0, 1.0));
            let p4 = matrix.transform_point(p + vec2(0.0, 1.0));

            let pts = vec![p1, p2, p3, p4];
            let calc = |f: fn(&Pos2) -> usize| -> (usize, usize) {
                let mut coords: Vec<_> = pts.iter().map(f).collect();
                coords.sort();
                (coords[0], *coords.last().unwrap() + 1)
            };
            let (xmin, xmax) = calc(|p| p.x as usize);
            let (ymin, ymax) = calc(|p| p.y as usize);
            let mut sum_area = 0.0;
            let expected_sum_area = area(&pts);
            for x in xmin..=xmax {
                for y in ymin..=ymax {
                    let area = calc_intersection_area(pts.clone(), x as f32, y as f32, false);
                    sum_area += area;
                    if x >= NEW_WIDTH as usize || y >= NEW_WIDTH as usize {
                        continue;
                    }
                    for i in 0..3 {
                        new_img_comps[x][y][i] += area * (pixel.0[i] as f32);
                    }
                    new_img_area[x][y] += area;
                }
            }
            if (sum_area - expected_sum_area).abs() > 0.1 {
                // dbg!(sum_area, expected_sum_area, xmin, xmax, ymin, ymax);
                // assert!(false);
            }
        }
    }

    let mut new_img = ImageBuffer::new(NEW_WIDTH as u32, NEW_WIDTH as u32);
    for x in 0..NEW_WIDTH {
        for y in 0..NEW_WIDTH {
            let mut res_color = new_img_comps[x][y].clone();
            let sum_area = new_img_area[x][y];
            assert!(sum_area > 0.001);
            for i in 0..3 {
                res_color[i] /= sum_area;
                if res_color[i] > 255.0 {
                    res_color[i] = 255.0;
                }
                if !(res_color[i] <= 256.0) {
                    dbg!(res_color[i]);
                    assert!(false);
                }
            }
            let pixel = image::Rgb([res_color[0] as u8, res_color[1] as u8, res_color[2] as u8]);
            new_img.put_pixel(x as u32, y as u32, pixel);
        }
    }
    new_img.save("img/crop.jpg").unwrap();
}
