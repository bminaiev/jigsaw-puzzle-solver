use crate::{point::PointF, utils::fmax};

#[derive(Clone)]
struct Rect {
    start: PointF,
    end: PointF,
}

impl Rect {
    pub fn intersects(&self, other: &Self) -> bool {
        if self.start.x >= other.end.x || self.end.x <= other.start.x {
            return false;
        }
        if other.start.x >= self.end.x || other.end.x <= self.start.x {
            return false;
        }
        if self.start.y >= other.end.y || self.end.y <= other.start.y {
            return false;
        }
        if other.start.y >= self.end.y || other.end.y <= self.start.y {
            return false;
        }
        true
    }
}

pub struct RectsFitter {
    rects: Vec<Rect>,
}

const OFFSET: f64 = 20.0;

impl RectsFitter {
    pub fn new() -> Self {
        Self {
            rects: vec![Rect {
                start: PointF::ZERO,
                end: PointF::ZERO,
            }],
        }
    }

    // returns shift, which we need to apply to points
    pub fn add_points(&mut self, pts: &[PointF]) -> PointF {
        let min_max = |f: fn(&PointF) -> f64| -> (f64, f64) {
            (
                pts.iter().map(f).min_by(f64::total_cmp).unwrap(),
                pts.iter().map(f).max_by(f64::total_cmp).unwrap(),
            )
        };
        let (min_x, max_x) = min_max(|p| p.x);
        let (min_y, max_y) = min_max(|p| p.y);
        let cur_start = PointF {
            x: min_x - OFFSET,
            y: min_y - OFFSET,
        };
        let cur_end = PointF {
            x: max_x + OFFSET,
            y: max_y + OFFSET,
        };
        let mut options = vec![];
        for x in self.rects.iter().map(|r| r.end.x) {
            for y in self.rects.iter().map(|r| r.end.y) {
                options.push(PointF { x, y });
            }
        }
        options.sort_by(|p1, p2| fmax(p1.x, p1.y).total_cmp(&fmax(p2.x, p2.y)));
        for &start in options.iter() {
            let mut ok = true;
            let rect_to_check = Rect {
                start,
                end: start + (cur_end - cur_start),
            };
            for r in self.rects.iter() {
                if r.intersects(&rect_to_check) {
                    ok = false;
                    break;
                }
            }
            if ok {
                self.rects.push(rect_to_check.clone());
                return rect_to_check.start - cur_start;
            }
        }
        unreachable!("Should have found position");
    }
}
