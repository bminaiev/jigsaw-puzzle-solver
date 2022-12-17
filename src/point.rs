use eframe::epaint::{pos2, Pos2};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Hash, PartialEq, PartialOrd, Ord, Eq, Debug)]
pub struct Point {
    pub x: usize,
    pub y: usize,
}

impl Point {
    pub fn dist2(&self, other: &Point) -> usize {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    pub fn pos2(&self) -> Pos2 {
        pos2(self.x as f32, self.y as f32)
    }

    pub fn angle_to(&self, other: &Self) -> f32 {
        let dy = (other.y as f32) - (self.y as f32);
        let dx = (other.x as f32) - (self.x as f32);
        dy.atan2(dx)
    }

    pub fn conv_f64(&self) -> PointF {
        PointF {
            x: self.x as f64,
            y: self.y as f64,
        }
    }

    pub fn neighbours(&self) -> Vec<Self> {
        let mut res = vec![
            Point {
                x: self.x + 1,
                y: self.y,
            },
            Point {
                x: self.x,
                y: self.y + 1,
            },
        ];
        if self.x > 0 {
            res.push(Point {
                x: self.x - 1,
                y: self.y,
            });
        }
        if self.y > 0 {
            res.push(Point {
                x: self.x,
                y: self.y - 1,
            });
        }
        res
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PointF {
    pub x: f64,
    pub y: f64,
}

impl PointF {
    pub const ZERO: Self = PointF { x: 0.0, y: 0.0 };

    pub fn dist2(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    pub fn rotate_ccw90(&self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
        }
    }

    pub fn len2(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    pub fn len(&self) -> f64 {
        self.len2().sqrt()
    }

    pub fn norm(&self) -> Self {
        let len = self.len();
        Self {
            x: self.x / len,
            y: self.y / len,
        }
    }

    pub fn scal_mul(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    pub fn pos2(&self) -> Pos2 {
        pos2(self.x as f32, self.y as f32)
    }

    pub fn rotate(&self, angle: f64) -> Self {
        Self {
            x: angle.cos() * self.x + angle.sin() * self.y,
            y: -angle.sin() * self.x + angle.cos() * self.y,
        }
    }

    pub fn vect_mul(p1: &Self, p2: &Self, p3: &Self) -> f64 {
        (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    }
}

impl std::ops::Add for PointF {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl std::ops::Div<f64> for PointF {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl std::ops::Mul<f64> for PointF {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl std::ops::Sub for PointF {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

pub fn find_center(pts: &[PointF]) -> PointF {
    let mut res = PointF::ZERO;
    for p in pts.iter() {
        res = res + *p;
    }
    res / (pts.len() as f64)
}
