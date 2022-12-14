use std::cmp::min;

use itertools::Itertools;

use crate::{
    coordinate_system::CoordinateSystem,
    figure::Figure,
    parsed_puzzles::ParsedPuzzles,
    point::{find_center, PointF},
    utils::{fmax, fmin, gauss, Side},
};

pub struct BorderAndNeighbors {
    pub border: Vec<PointF>,
    pub prev: Vec<PointF>,
    pub next: Vec<PointF>,
}

impl BorderAndNeighbors {
    pub const NEIGHBORS: usize = 10;

    pub fn reverse(&self) -> Self {
        let mut border = self.border.clone();
        border.reverse();
        Self {
            border,
            prev: self.next.clone(),
            next: self.prev.clone(),
        }
    }
}

// smaller -> better
pub fn match_placed_borders_center(lhs: &[PointF], rhs: &[PointF]) -> f64 {
    // TODO: smarter logic

    const CHECK_NEXT: usize = 5;
    const OK_MIN_DIFF: f64 = 0.2;
    let score_one_side = |lhs: &[PointF], rhs: &[PointF]| -> f64 {
        let socket_detector = SocketDetector::new(lhs);

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
            let coef = if socket_detector.is_socket_point(p) {
                1.0
            } else {
                4.0
            };

            sum_dists += coef * fmax(0.0, dists[min_shift] - OK_MIN_DIFF);
        }
        sum_dists / (lhs.len() as f64)
    };

    score_one_side(lhs, rhs) + score_one_side(rhs, lhs)
}

#[derive(Clone, Copy)]
struct Line {
    a: f64,
    b: f64,
    c: f64,
}

impl Line {
    pub fn new(p1: &PointF, p2: &PointF) -> Self {
        let a = p2.y - p1.y;
        let b = p1.x - p2.x;
        let c = -(p1.x * a + p1.y * b);
        Self { a, b, c }
    }

    fn dist(&self, p: PointF) -> f64 {
        (self.a * p.x + self.b * p.y + self.c).abs() / (self.a * self.a + self.b * self.b).sqrt()
    }
}

fn match_side_borders(lhs: &[PointF], rhs: &[PointF]) -> f64 {
    let mut res = f64::MAX / 100.0;
    for p1 in lhs.iter() {
        for p2 in rhs.iter() {
            if p1 != p2 {
                let line = Line::new(p1, p2);
                let mut sum_dist = 0.0;
                for &np in lhs.iter() {
                    sum_dist += line.dist(np);
                }
                for &np in rhs.iter() {
                    sum_dist += line.dist(np);
                }
                res = fmin(res, sum_dist);
            }
        }
    }
    let cnt = (lhs.len() + rhs.len()) as f64;
    res / cnt
}

fn local_optimize(coefs: &mut [f64], mut f: impl FnMut(&[f64]) -> f64, max_iters: usize) {
    let mut change = 1.0;
    const MULT: f64 = 0.5;
    const MIN_DELTA: f64 = 1e-5;
    let mut last_score = f(coefs);
    for iters in 0..max_iters {
        if change < MIN_DELTA {
            // eprintln!("finished after {iters}");
            return;
        }
        let mut changed = false;
        for pos in 0..coefs.len() {
            for sign in [-1.0, 1.0].iter() {
                coefs[pos] += sign * change;
                let new_score = f(coefs);
                if new_score < last_score {
                    last_score = new_score;
                    changed = true;
                } else {
                    coefs[pos] -= sign * change;
                }
            }
        }
        if !changed {
            change *= MULT;
        }
    }
}

// sum (ax^2 + bx + c - y) * (ax^2 + bx + c - y) -> min
// z0 * a^2 + z1 * b^2 + z2 * c^2 + z3 * a * b + ... -> min
// df/da = 0
// 2a * z0 + z3 * b + z4 * c + z6 = 0
fn find_polynomial_coefs(pts: &[PointF]) -> Vec<f64> {
    let mut z = [0.0; 9];
    for p in pts.iter() {
        let x = p.x;
        let y = p.y;
        let x2 = x * x;
        z[0] += x2 * x2;
        z[1] += x * x;
        z[2] += 1.0;
        z[3] += 2.0 * x2 * x;
        z[4] += 2.0 * x2;
        z[5] += 2.0 * x;
        z[6] += -2.0 * y * x2;
        z[7] += -2.0 * y * x;
        z[8] += -2.0 * y;
    }
    let mut equations = vec![];
    {
        // f/da = 0
        equations.push(vec![2.0 * z[0], z[3], z[4], z[6]]);
        // f/db = 0
        equations.push(vec![z[3], 2.0 * z[1], z[5], z[7]]);
        // f/dc = 0
        equations.push(vec![z[4], z[5], 2.0 * z[2], z[8]]);
    }
    let mut res = gauss(&mut equations);
    res.reverse();
    res
}

fn interpolate_with_polynomial_score(pts: &[PointF], debug: bool) -> (f64, Vec<PointF>) {
    let first = pts[0];
    let last = *pts.last().unwrap();
    let dir = (last - first).norm();
    let rev_dir = dir.rotate_ccw90();
    let normalized = pts
        .iter()
        .map(|&p| PointF {
            x: dir.scal_mul(&(p - first)),
            y: rev_dir.scal_mul(&(p - first)),
        })
        .collect_vec();
    // let mut coefs = vec![0.0; 3];

    const OK_MIN_DIFF: f64 = 0.3;
    const MULTIPLIER: f64 = 50.0;

    let calc_score = |coefs: &[f64]| -> f64 {
        normalized
            .iter()
            .map(|p| {
                let poly_y = coefs[0] + coefs[1] * p.x + coefs[2] * p.x * p.x;
                let diff_y = (poly_y - p.y).abs();
                let diff_y = fmax(0.0, diff_y - OK_MIN_DIFF);
                diff_y * diff_y
            })
            .sum::<f64>()
            * MULTIPLIER
            / normalized.len() as f64
    };
    // local_optimize(&mut coefs, calc_score, 1000);
    // eprintln!("COEFS: {:?}", coefs);
    let another_coefs = find_polynomial_coefs(&normalized);
    // eprintln!("FOUND {:?}", another_coefs);
    let coefs = another_coefs;
    let debug_points = if !debug {
        vec![]
    } else {
        normalized
            .iter()
            .map(|p| {
                let poly_y = coefs[0] + coefs[1] * p.x + coefs[2] * p.x * p.x;
                first + dir * p.x + rev_dir * poly_y
            })
            .collect_vec()
    };
    (calc_score(&coefs), debug_points)
}

pub fn match_side_borders_v2(lhs: &[PointF], rhs: &[PointF], debug: bool) -> (f64, Vec<PointF>) {
    const REM_FIRST: usize = 2;
    if lhs.len() <= REM_FIRST || rhs.len() <= REM_FIRST {
        return (f64::MAX / 100.0, vec![]);
    }
    let lhs = &lhs[REM_FIRST..];
    let rhs = &rhs[REM_FIRST..];

    let sum_pts = lhs.iter().rev().chain(rhs.iter()).cloned().collect_vec();
    interpolate_with_polynomial_score(&sum_pts, debug)
}

pub fn match_placed_borders(lhs: &BorderAndNeighbors, rhs: &BorderAndNeighbors) -> f64 {
    let center = match_placed_borders_center(&lhs.border, &rhs.border);
    let prev = match_side_borders_v2(&lhs.prev, &rhs.prev, false).0;
    let next = match_side_borders_v2(&lhs.next, &rhs.next, false).0;
    center + (prev + next)
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
        // let move_pts = |pts: Vec<PointF>| -> Vec<PointF> { pts.iter().map(|p| *p).collect() };

        // let lhs = move_pts(lhs);
        // let rhs = move_pts(rhs);
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

    pub fn get_offset(&self) -> PointF {
        let min_x = self
            .lhs
            .iter()
            .chain(self.rhs.iter())
            .map(|p| p.x)
            .min_by(|a, b| a.total_cmp(b))
            .unwrap();
        let min_y = self
            .lhs
            .iter()
            .chain(self.rhs.iter())
            .map(|p| p.y)
            .min_by(|a, b| a.total_cmp(b))
            .unwrap();
        PointF {
            x: -min_x,
            y: -min_y,
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

pub fn local_optimize_coordinate_systems(
    start_cs: &[CoordinateSystem],
    mut scorer: impl FnMut(&[CoordinateSystem]) -> f64,
) -> Vec<CoordinateSystem> {
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
    let mut cs = start_cs.to_vec();
    let mut changed_steps = 0;
    const MAX_CHANGED_STEPS: usize = 5;
    while start_coord_step > MIN_EPS || dir_step > MIN_EPS {
        changed_steps += 1;
        if changed_steps > MAX_CHANGED_STEPS {
            changed_steps = 0;
        }
        {
            let mut changed = false;
            for cs_id in 0..cs.len() {
                for mv in moves.iter() {
                    let ncs = CoordinateSystem::new(
                        cs[cs_id].start + *mv * start_coord_step,
                        cs[cs_id].x_dir,
                    );
                    let prev_cs = cs[cs_id].clone();
                    cs[cs_id] = ncs;
                    let new_score = scorer(&cs);
                    if new_score < last_score {
                        changed = true;
                        last_score = new_score;
                    } else {
                        cs[cs_id] = prev_cs;
                    }
                }
            }
            if !changed || changed_steps == MAX_CHANGED_STEPS {
                start_coord_step *= MULT;
            }
        }
        {
            let mut changed = false;
            for cs_id in 0..cs.len() {
                for mv in moves.iter() {
                    let ncs =
                        CoordinateSystem::new(cs[cs_id].start, cs[cs_id].x_dir + *mv * dir_step);
                    let prev_cs = cs[cs_id].clone();
                    cs[cs_id] = ncs;
                    let new_score = scorer(&cs);
                    if new_score < last_score {
                        changed = true;
                        last_score = new_score;
                    } else {
                        cs[cs_id] = prev_cs;
                    }
                }
            }
            if !changed || changed_steps == MAX_CHANGED_STEPS {
                dir_step *= MULT;
            }
        }
    }
    cs
}

struct SocketDetector {
    start: PointF,
    dir: PointF,
    len: f64,
}

impl SocketDetector {
    pub fn new(pts: &[PointF]) -> Self {
        let dir = *pts.last().unwrap() - pts[0];
        let len = dir.len();
        let dir = dir.norm();
        Self {
            dir,
            len,
            start: pts[0],
        }
    }

    pub fn is_socket_point(&self, p: &PointF) -> bool {
        ((p.x - self.start.x) * self.dir.y - self.dir.x * (p.y - self.start.y)).abs()
            > 0.1 * self.len
    }
}

fn is_picture_border_impl(pts: &[PointF]) -> bool {
    let socket_detector = SocketDetector::new(pts);
    !pts.iter().any(|p| socket_detector.is_socket_point(p))
}

pub fn is_picture_border(figure: &Figure, border_id: usize) -> bool {
    let pts = get_figure_border(figure, border_id);
    is_picture_border_impl(&pts)
}

fn get_figure_border_and_neighbors(figure: &Figure, border_id: usize) -> BorderAndNeighbors {
    let mut prev = get_figure_border(figure, (border_id + 3) % 4);
    {
        let sz = min(BorderAndNeighbors::NEIGHBORS, prev.len());
        prev.rotate_right(sz);
        prev.truncate(sz);
        prev.reverse();
    }
    let mut next = get_figure_border(figure, (border_id + 1) % 4);
    next.truncate(BorderAndNeighbors::NEIGHBORS);

    BorderAndNeighbors {
        border: get_figure_border(figure, border_id),
        prev,
        next,
    }
}

pub fn match_borders(
    parsed_puzzles: &ParsedPuzzles,
    side1: Side,
    side2: Side,
) -> Option<MatchResult> {
    let lhs_figure = &parsed_puzzles.figures[side1.fig];
    let rhs_figure = &parsed_puzzles.figures[side2.fig];

    let lhs = get_figure_border_and_neighbors(&lhs_figure, side1.side);
    let rhs = get_figure_border_and_neighbors(&rhs_figure, side2.side);
    let rhs = rhs.reverse();

    if is_picture_border_impl(&lhs.border) || is_picture_border_impl(&rhs.border) {
        return None;
    }

    if is_picture_border_impl(&get_figure_border(lhs_figure, (side1.side + 1) % 4))
        != is_picture_border_impl(&get_figure_border(rhs_figure, (side2.side + 3) % 4))
    {
        return None;
    }

    if is_picture_border_impl(&get_figure_border(lhs_figure, (side1.side + 3) % 4))
        != is_picture_border_impl(&get_figure_border(rhs_figure, (side2.side + 1) % 4))
    {
        return None;
    }

    let to_cs = estimate_coordinate_system_by_border(&lhs.border)?;
    let from_cs_estimation = estimate_coordinate_system_by_border(&rhs.border)?;

    let conv_point =
        |from_cs: &CoordinateSystem, p: PointF| -> PointF { to_cs.to_real(from_cs.create(p)) };

    let move_rhs = |from_cs: &CoordinateSystem| -> BorderAndNeighbors {
        let conv_v = |v: &[PointF]| -> Vec<PointF> {
            v.iter().map(|p| conv_point(from_cs, *p)).collect_vec()
        };
        BorderAndNeighbors {
            border: conv_v(&rhs.border),
            prev: conv_v(&rhs.prev),
            next: conv_v(&rhs.next),
        }
    };

    let from_cs = if match_placed_borders(&lhs, &move_rhs(&from_cs_estimation)) > 100.0 {
        from_cs_estimation
    } else {
        local_optimize_coordinate_systems(&[from_cs_estimation], |from_cs| {
            match_placed_borders(&lhs, &move_rhs(&from_cs[0]))
        })[0]
            .clone()
    };

    let res = MatchResult::new(
        match_placed_borders(&lhs, &move_rhs(&from_cs)),
        lhs_figure.border.iter().map(|p| p.conv_f64()).collect_vec(),
        rhs_figure
            .border
            .iter()
            .map(|p| conv_point(&from_cs, p.conv_f64()))
            .collect_vec(),
        side1.fig,
        side2.fig,
    );
    Some(res)
}

pub fn match_borders_without_move(
    lhs_figure: &Figure,
    lhs_border_id: usize,
    rhs_figure: &Figure,
    rhs_border_id: usize,
    lhs_id: usize,
    rhs_id: usize,
) -> Option<MatchResult> {
    let lhs = get_figure_border_and_neighbors(&lhs_figure, lhs_border_id);
    let mut rhs = get_figure_border_and_neighbors(&rhs_figure, rhs_border_id);
    let rhs = rhs.reverse();

    let score = match_placed_borders(&lhs, &rhs);

    if score > 30.0 {
        return None;
    }

    let res = MatchResult::new(
        score,
        lhs_figure.border.iter().map(|p| p.conv_f64()).collect_vec(),
        rhs_figure.border.iter().map(|p| p.conv_f64()).collect_vec(),
        lhs_id,
        rhs_id,
    );
    Some(res)
}
