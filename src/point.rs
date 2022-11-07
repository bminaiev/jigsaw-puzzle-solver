use eframe::epaint::{pos2, Pos2};

#[derive(Clone, Copy)]
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
}

#[derive(Clone, Copy)]
pub struct PointF {
    pub x: f64,
    pub y: f64,
}

impl PointF {
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
