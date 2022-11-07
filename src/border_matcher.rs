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

pub struct MatchResult {
    pub score: f64,
    pub lhs: Vec<PointF>,
    pub rhs: Vec<PointF>,
    // TODO: position
}

fn get_figure_border(figure: &Figure, border_id: usize) -> Vec<PointF> {
    let mut res = vec![];
    let mut cur = figure.corner_positions[border_id];
    let to = figure.corner_positions[(border_id + 1) % figure.corner_positions.len()];
    loop {
        let p = figure.border[cur];
        res.push(PointF {
            x: p.x as f64,
            y: p.y as f64,
        });
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

    MatchResult {
        score: match_placed_borders(&lhs, &rhs_moved),
        lhs,
        rhs: rhs_moved,
    }
}
