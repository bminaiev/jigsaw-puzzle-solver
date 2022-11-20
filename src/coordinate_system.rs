use crate::point::PointF;

#[derive(Clone)]
pub struct CoordinateSystem {
    pub start: PointF,
    pub x_dir: PointF,
    pub y_dir: PointF,
}

impl CoordinateSystem {
    pub fn new(start: PointF, x_dir: PointF) -> Self {
        assert!(x_dir.len2() != 0.0);
        let x_dir = x_dir.norm();
        let y_dir = x_dir.rotate_ccw90();
        Self {
            start,
            x_dir,
            y_dir,
        }
    }

    pub fn create(&self, p: PointF) -> PointF {
        let p = p - self.start;
        PointF {
            x: self.x_dir.scal_mul(&p),
            y: self.y_dir.scal_mul(&p),
        }
    }

    pub fn to_real(&self, p: PointF) -> PointF {
        self.start + self.x_dir * p.x + self.y_dir * p.y
    }
}
