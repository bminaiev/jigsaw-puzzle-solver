use std::cmp::min;

use itertools::Itertools;

use crate::{
    coordinate_system::CoordinateSystem,
    figure::Figure,
    point::{find_center, PointF},
};

// smaller -> better
pub fn match_placed_borders(lhs: &[PointF], rhs: &[PointF]) -> f64 {
    // TODO: smarter logic

    const CHECK_NEXT: usize = 5;
    let score_one_side = |lhs: &[PointF], rhs: &[PointF]| -> f64 {
        let mut iter = 0;
        let mut sum_dists = 0.0;
        for p in lhs.iter() {
            let dists = rhs[iter..min(iter + CHECK_NEXT, rhs.len())]
                .iter()
                .map(|rhs_p| p.dist2(rhs_p))
                .collect_vec();
            let mut min_shift = 0;
            for i in 0..dists.len() {
                if dists[i] < dists[min_shift] {
                    min_shift = i;
                }
            }
            iter += min_shift;
            sum_dists += dists[min_shift];
        }
        sum_dists / (lhs.len() as f64)
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
    pub shifted: PointF,
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
        let shifted = PointF {
            x: -min_x,
            y: -min_y,
        };
        let move_pts =
            |pts: Vec<PointF>| -> Vec<PointF> { pts.iter().map(|p| *p + shifted).collect() };

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
            shifted,
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

fn estimate_coordinate_system_by_border(border: &[PointF]) -> Option<CoordinateSystem> {
    const OFFSET: usize = 5;
    if border.len() <= OFFSET * 2 + 3 {
        return None;
    }
    let mid = border.len() / 2;
    let p1 = find_center(&border[OFFSET..mid]);
    let p2 = find_center(&border[mid..border.len() - OFFSET]);
    if p1 == p2 {
        return None;
    }
    Some(CoordinateSystem::new(p1, p2 - p1))
}

pub fn local_optimize_coordinate_system(
    start_cs: CoordinateSystem,
    mut scorer: impl FnMut(&CoordinateSystem) -> f64,
) -> CoordinateSystem {
    let start_score = scorer(&start_cs);
    let mut last_score = start_score;

    // TODO: think about constants
    let mut start_coord_step = 10.0;
    let mut dir_step = 0.1;
    const MIN_EPS: f64 = 1e-2;
    const MULT: f64 = 0.3;
    let moves = vec![
        PointF { x: 1.0, y: 0.0 },
        PointF { x: -1.0, y: 0.0 },
        PointF { x: 0.0, y: 1.0 },
        PointF { x: 0.0, y: -1.0 },
    ];
    let mut cs = start_cs;
    while start_coord_step > MIN_EPS || dir_step > MIN_EPS {
        {
            let mut changed = false;
            for mv in moves.iter() {
                let ncs = CoordinateSystem::new(cs.start + *mv * start_coord_step, cs.x_dir);
                let new_score = scorer(&ncs);
                if new_score < last_score {
                    changed = true;
                    last_score = new_score;
                    cs = ncs;
                }
            }
            if !changed {
                start_coord_step *= MULT;
            }
        }
        {
            let mut changed = false;
            for mv in moves.iter() {
                let ncs = CoordinateSystem::new(cs.start, cs.x_dir + *mv * dir_step);
                let new_score = scorer(&ncs);
                if new_score < last_score {
                    changed = true;
                    last_score = new_score;
                    cs = ncs;
                }
            }
            if !changed {
                dir_step *= MULT;
            }
        }
    }
    cs
}

// TODO: use `Side` type
pub fn match_borders(
    lhs_figure: &Figure,
    lhs_border_id: usize,
    rhs_figure: &Figure,
    rhs_border_id: usize,
    lhs_id: usize,
    rhs_id: usize,
) -> Option<MatchResult> {
    let lhs = get_figure_border(&lhs_figure, lhs_border_id);
    let mut rhs = get_figure_border(&rhs_figure, rhs_border_id);
    rhs.reverse();

    let to_cs = estimate_coordinate_system_by_border(&lhs)?;
    let from_cs_estimation = estimate_coordinate_system_by_border(&rhs)?;

    let conv_point =
        |from_cs: &CoordinateSystem, p: PointF| -> PointF { to_cs.to_real(from_cs.create(p)) };

    let move_rhs = |from_cs: &CoordinateSystem| -> Vec<PointF> {
        rhs.iter().map(|p| conv_point(from_cs, *p)).collect_vec()
    };

    let from_cs = if match_placed_borders(&lhs, &move_rhs(&from_cs_estimation)) > 100.0 {
        from_cs_estimation
    } else {
        local_optimize_coordinate_system(from_cs_estimation, |from_cs| {
            match_placed_borders(&lhs, &move_rhs(from_cs))
        })
    };

    let res = MatchResult::new(
        match_placed_borders(&lhs, &move_rhs(&from_cs)),
        lhs_figure.border.iter().map(|p| p.conv_f64()).collect_vec(),
        rhs_figure
            .border
            .iter()
            .map(|p| conv_point(&from_cs, p.conv_f64()))
            .collect_vec(),
        lhs_id,
        rhs_id,
    );
    Some(res)
}
