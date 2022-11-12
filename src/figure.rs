use std::{
    collections::{BTreeMap, VecDeque},
    f32::consts::PI,
};

use itertools::Itertools;

use crate::point::Point;

#[derive(Clone, Hash)]
pub struct Figure {
    pub all_pts: Vec<Point>,
    pub border: Vec<Point>,
    pub good_border: bool,
    pub center: Point,
    pub corner_positions: Vec<usize>,
}

fn find_cycle(points: &[Point]) -> Vec<Point> {
    let mut used = vec![false; points.len()];
    used[0] = true;
    let mut path = vec![points[0]];
    loop {
        let last = *path.last().unwrap();
        let mut next = None;
        for i in 0..used.len() {
            if !used[i] {
                let ndist = last.dist2(&points[i]);
                if next.is_none() || last.dist2(&points[next.unwrap()]) > ndist {
                    next = Some(i);
                }
            }
        }
        if let Some(next) = next {
            path.push(points[next]);
            used[next] = true;
        } else {
            break;
        }
    }
    optimize_cycle(&mut path);
    let area = calc_signed_area(&path);
    if area > 0.0 {
        // make clock-wise
        path.reverse();
    }
    path
}

fn cycle_rev(a: &mut [Point], mut from: usize, mut len: usize) {
    while len > 1 {
        let end = (from + len - 1) % a.len();
        a.swap(from, end);
        from = (from + 1) % a.len();
        len -= 2;
    }
}

fn optimize_cycle(cycle: &mut Vec<Point>) {
    // reverse a part of the cycle to minimize sum of dist^2 between neighbours
    loop {
        let mut changed = false;
        for i in 0..cycle.len() {
            for len in 2..(cycle.len() + 1) / 2 {
                let first = cycle[i];
                let prev = cycle[(i + cycle.len() - 1) % cycle.len()];
                let last = cycle[(i + len - 1) % cycle.len()];
                let next = cycle[(i + len) % cycle.len()];

                let cur_dist = prev.dist2(&first) + last.dist2(&next);
                let potential_dist = prev.dist2(&last) + first.dist2(&next);
                if potential_dist < cur_dist {
                    changed = true;
                    cycle_rev(cycle, i, len);
                }
            }
        }

        if !changed {
            break;
        }
    }
}

fn find_center(points: &[Point]) -> Point {
    let sum_x: usize = points.iter().map(|p| p.x).sum();
    let sum_y: usize = points.iter().map(|p| p.y).sum();
    let cnt = points.len();
    Point {
        x: sum_x / cnt,
        y: sum_y / cnt,
    }
}

fn calc_signed_area(pts: &[Point]) -> f32 {
    let mut res = 0.0;
    for i in 0..pts.len() {
        let p1 = pts[i];
        let p2 = pts[(i + 1) % pts.len()];
        res += (p1.x as f32) * (p2.y as f32) - (p1.y as f32) * (p2.x as f32);
    }
    res
}

fn fmax(x: f32, y: f32) -> f32 {
    if x > y {
        x
    } else {
        y
    }
}

fn fmin(x: f32, y: f32) -> f32 {
    if x < y {
        x
    } else {
        y
    }
}

fn filter_4_corners(all: &[Point], idxs: &[usize]) -> Vec<usize> {
    let mut best_score = f32::MAX;
    let n = idxs.len();
    let mut best = vec![];

    fn norm_angle(mut a: f32) -> f32 {
        while a < 0.0 {
            a += PI * 2.0;
        }
        while a >= PI * 2.0 {
            a -= PI * 2.0;
        }
        fmin(a, PI * 2.0 - a)
    }

    fn find_angle(p1: Point, p2: Point, p3: Point) -> f32 {
        let a1 = p1.angle_to(&p2);
        let a2 = p2.angle_to(&p3);
        norm_angle(a1 - a2)
    }

    fn angle_score(a: f32) -> f32 {
        (PI / 2.0 - a).abs()
    }

    let calc_score = |pts: &[Point]| -> f32 {
        let sum_angles: f32 = (0..4)
            .map(|i| angle_score(find_angle(pts[i], pts[(i + 1) % 4], pts[(i + 2) % 4])))
            .sum();

        let dists: Vec<_> = (0..4).map(|i| pts[i].dist2(&pts[(i + 1) % 4])).collect();
        let max_dist = *dists.iter().max().unwrap();
        let min_dist = *dists.iter().min().unwrap();
        let dist_coef = (max_dist as f32) / (min_dist as f32);

        sum_angles * dist_coef
    };

    for i1 in 0..n {
        for i2 in i1 + 1..n {
            for i3 in i2 + 1..n {
                for i4 in i3 + 1..n {
                    let pts = [all[idxs[i1]], all[idxs[i2]], all[idxs[i3]], all[idxs[i4]]];
                    let score = calc_score(&pts);
                    if score < best_score {
                        best_score = score;
                        best = vec![idxs[i1], idxs[i2], idxs[i3], idxs[i4]];
                    }
                }
            }
        }
    }
    best
}

fn find_corners_positions(points: &[Point]) -> Vec<usize> {
    let center = find_center(points);
    let dists: Vec<_> = points.iter().map(|p| p.dist2(&center)).collect();
    const CHECK_LEN: usize = 5;
    let mut corners = vec![];
    for i in 0..dists.len() {
        let (_, max_pos) = (0..CHECK_LEN * 2)
            .map(|shift| (dists[(i + shift) % dists.len()], shift))
            .max()
            .unwrap();
        if max_pos == CHECK_LEN {
            corners.push((i + max_pos) % dists.len());
        }
    }
    filter_4_corners(&points, &corners)
}

fn bfs(pts: &[Point], max_dist: usize) -> BTreeMap<Point, usize> {
    let mut res = BTreeMap::new();
    let mut queue = VecDeque::new();
    for p in pts.iter() {
        res.insert(*p, 0);
        queue.push_back(*p);
    }
    while let Some(p) = queue.pop_front() {
        let cur_dist = res[&p];
        if cur_dist == max_dist {
            continue;
        }
        for next in p.neighbours() {
            if !res.contains_key(&next) {
                res.insert(next, cur_dist + 1);
                queue.push_back(next);
            }
        }
    }
    res
}

impl Figure {
    pub fn new(pts: &[Point]) -> Option<Self> {
        const MAX_DIST: usize = 2;
        let dist_by_pts = bfs(pts, MAX_DIST);
        let border = dist_by_pts
            .iter()
            .filter_map(|(&k, &v)| if v == MAX_DIST { Some(k) } else { None })
            .collect_vec();

        if border.len() <= 50 {
            return None;
        }

        let cycle = find_cycle(&border);
        let max_dist2 = cycle
            .iter()
            .circular_tuple_windows()
            .map(|(w0, w1)| w0.dist2(&w1))
            .max()
            .unwrap();
        let center = find_center(&cycle);

        let res = Self {
            all_pts: dist_by_pts.keys().cloned().collect_vec(),
            corner_positions: find_corners_positions(&cycle),
            border: cycle,
            good_border: max_dist2 <= 10,
            center,
        };
        Some(res)
    }

    pub fn is_good_puzzle(&self) -> bool {
        self.good_border && self.corner_positions.len() == 4
    }
}
