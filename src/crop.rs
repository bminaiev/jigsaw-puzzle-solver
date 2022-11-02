use eframe::epaint::{pos2, Pos2};
use image::{GenericImageView, ImageBuffer};

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

pub fn crop(path: &str, frame: &[Pos2]) {
    let image = image::open(path).unwrap();
    let w = image.width() as f32;
    let h = image.height() as f32;

    // let frame = vec![pos2(100.0, 100.0), pos2(w, 0.0), pos2(w, h), pos2(0.0, h)];

    let new_width = ((frame[0].distance(frame[1]) + frame[2].distance(frame[3])) / 2.0) as usize;
    let new_height = ((frame[1].distance(frame[2]) + frame[3].distance(frame[0])) / 2.0) as usize;

    let pt = |from: Pos2, to: Pos2, pos: usize, max_pos: usize| -> Pos2 {
        let left = pos as f32;
        let total = max_pos as f32;
        let right = total - left;
        pos2(
            (from.x * right + to.x * left) / total,
            (from.y * right + to.y * left) / total,
        )
    };
    dbg!(new_height, new_width);
    let mut new_img = ImageBuffer::new(new_width as u32, new_height as u32);
    for x in 0..new_width {
        for y in 0..new_height {
            let offsets = [x, y, new_width - x - 1, new_height - y - 1];
            let max_len = [new_width, new_height, new_width, new_height];
            let mut points = Vec::with_capacity(offsets.len() * 2);
            for i in 0..offsets.len() {
                let fr = frame[i];
                let to = frame[(i + 1) % frame.len()];
                let p1 = pt(fr, to, offsets[i], max_len[i]);
                let p2 = pt(fr, to, offsets[i] + 1, max_len[i]);
                points.push(p1);
                points.push(p2);
            }
            let l1 = Line::new(points[0], points[5]);
            let l2 = Line::new(points[1], points[4]);
            let l3 = Line::new(points[2], points[7]);
            let l4 = Line::new(points[3], points[6]);
            let p1 = intersection(l1, l3);
            let p2 = intersection(l2, l3);
            let p3 = intersection(l2, l4);
            let p4 = intersection(l1, l4);
            let pts = vec![p1, p2, p3, p4];
            let calc = |f: fn(&Pos2) -> usize| -> (usize, usize) {
                let mut coords: Vec<_> = pts.iter().map(f).collect();
                coords.sort();
                (coords[0], *coords.last().unwrap() + 1)
            };
            let (xmin, xmax) = calc(|p| p.x as usize);
            let (ymin, ymax) = calc(|p| p.y as usize);
            let mut sum_area = 0.0;
            let mut res_color = [0.0, 0.0, 0.0];
            for x in xmin..=xmax {
                for y in ymin..=ymax {
                    let area = calc_intersection_area(pts.clone(), x as f32, y as f32, false);
                    sum_area += area;
                    if x >= image.width() as usize || y >= image.height() as usize {
                        continue;
                    }
                    let pixel = image.get_pixel(x as u32, y as u32);
                    for i in 0..3 {
                        res_color[i] += (pixel.0[i] as f32) * area;
                    }
                }
            }
            if (sum_area == 0.0) {
                dbg!(sum_area);
                dbg!(x);
                dbg!(y);

                for x in xmin..xmax {
                    for y in ymin..ymax {
                        let area = calc_intersection_area(pts.clone(), x as f32, y as f32, true);
                        dbg!(x, y, area);
                    }
                }
            }
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
