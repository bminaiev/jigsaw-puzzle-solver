use itertools::Itertools;

use crate::{figure::Figure, point::PointF};

// smaller -> better
pub fn match_placed_borders(lhs: &[PointF], rhs: &[PointF]) -> f64 {
    // TODO: smarter logic
    let expected_dist2 =
        (lhs[0].dist2(&rhs[0]) + lhs.last().unwrap().dist2(&rhs.last().unwrap())) / 2.0;

    // For each point from each side, we find the closest point from the other side, and check
    // how much it differs from `expected_dist2`.
    //
    // Then we calculate the average of all points.

    let score_one_side = |lhs: &[PointF], rhs: &[PointF]| -> f64 {
        lhs.iter()
            .map(|p| {
                let min_dist = rhs
                    .iter()
                    .map(|other| p.dist2(other))
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                (min_dist - expected_dist2).abs()
            })
            .sum::<f64>()
            / (lhs.len() as f64)
    };

    score_one_side(lhs, rhs) + score_one_side(rhs, lhs)
}

#[derive(Clone)]
pub struct MatchResult {
    pub score: f64,
    pub lhs: Vec<PointF>,
    pub rhs: Vec<PointF>,
    pub lhs_id: usize,
    pub rhs_id: usize,
    pub lhs_center: PointF,
    pub rhs_center: PointF,
    // TODO: position
}

impl MatchResult {
    pub fn new(
        score: f64,
        lhs: Vec<PointF>,
        rhs: Vec<PointF>,
        lhs_id: usize,
        rhs_id: usize,
    ) -> Self {
        let min_x = lhs
            .iter()
            .chain(rhs.iter())
            .map(|p| p.x)
            .min_by(|a, b| a.total_cmp(b))
            .unwrap();
        let min_y = lhs
            .iter()
            .chain(rhs.iter())
            .map(|p| p.y)
            .min_by(|a, b| a.total_cmp(b))
            .unwrap();
        let move_pts = |pts: Vec<PointF>| -> Vec<PointF> {
            pts.iter()
                .map(|p| *p - PointF { x: min_x, y: min_y })
                .collect()
        };
        let find_center = |pts: &[PointF]| -> PointF {
            let mut res = PointF::ZERO;
            for p in pts.iter() {
                res = res + *p;
            }
            res / (pts.len() as f64)
        };
        let lhs = move_pts(lhs);
        let rhs = move_pts(rhs);
        Self {
            lhs_center: find_center(&lhs),
            rhs_center: find_center(&rhs),
            score,
            lhs,
            rhs,
            lhs_id,
            rhs_id,
        }
    }
}

fn get_figure_border(figure: &Figure, border_id: usize) -> Vec<PointF> {
    let mut res = vec![];
    let mut cur = figure.corner_positions[border_id];
    let to = figure.corner_positions[(border_id + 1) % figure.corner_positions.len()];
    loop {
        let p = figure.border[cur];
        res.push(p.conv_f64());
        if cur == to {
            break;
        }
        cur = (cur + 1) % figure.border.len();
    }
    res
}

struct CoordinateSystem {
    start: PointF,
    x_dir: PointF,
    y_dir: PointF,
}

impl CoordinateSystem {
    pub fn new(start: PointF, x_dir: PointF) -> Self {
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

pub fn match_borders(
    lhs_figure: &Figure,
    lhs_border_id: usize,
    rhs_figure: &Figure,
    rhs_border_id: usize,
    lhs_id: usize,
    rhs_id: usize,
) -> MatchResult {
    let lhs = get_figure_border(&lhs_figure, lhs_border_id);
    let mut rhs = get_figure_border(&rhs_figure, rhs_border_id);
    rhs.reverse();

    let mid = (lhs[0] + *lhs.last().unwrap()) / 2.0;
    let x_dir = *lhs.last().unwrap() - lhs[0];
    let dir = x_dir.rotate_ccw90().norm();

    // TODO: 5.0 is not a bullet-proof solution.
    let start_mid = mid; //+ (dir * 5.0);

    let to_cs = CoordinateSystem::new(start_mid, x_dir);

    let rhs_mid = (rhs[0] + *rhs.last().unwrap()) / 2.0;
    let from_cs = CoordinateSystem::new(rhs_mid, *rhs.last().unwrap() - rhs[0]);

    let conv_point = |p: PointF| -> PointF {
        let in_system = from_cs.create(p);
        to_cs.to_real(in_system)
    };

    let rhs_moved = rhs.into_iter().map(conv_point).collect_vec();

    MatchResult::new(
        match_placed_borders(&lhs, &rhs_moved),
        lhs_figure.border.iter().map(|p| p.conv_f64()).collect_vec(),
        rhs_figure
            .border
            .iter()
            .map(|p| conv_point(p.conv_f64()))
            .collect_vec(),
        lhs_id,
        rhs_id,
    )
}
