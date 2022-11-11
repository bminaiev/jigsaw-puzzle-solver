use itertools::Itertools;
use ndarray::Array4;

use crate::{
    borders_graph::Graph,
    placement::{self, Placement},
    utils::Side,
};

#[derive(Clone, Debug)]
struct Four {
    s0: Side,
    s1: Side,
    s2: Side,
    s3: Side,
    max_dist: f64,
}

impl Four {
    pub fn all(&self) -> [Side; 4] {
        [self.s0, self.s1, self.s2, self.s3]
    }
}

fn fmax4(a: [f64; 4]) -> f64 {
    let mut res = a[0];
    for &x in a[1..].iter() {
        if x > res {
            res = x;
        }
    }
    res
}

pub fn solve_graph(graph: &Graph) {
    let mut edges = graph.all_edges.clone();
    edges.sort_by(|e1, e2| e1.score.total_cmp(&e2.score));
    let mut placement = Placement::new(graph.n);
    for e in edges.iter() {
        if placement
            .join_sides(
                Side {
                    fig: e.fig1,
                    side: e.side1,
                },
                Side {
                    fig: e.fig2,
                    side: e.side2,
                },
            )
            .is_some()
        {
            eprintln!("Join sides! {}", e.score);
        }
    }
}

pub fn solve_graph2(graph: &Graph) {
    eprintln!("Hello there!");
    let n = graph.n;
    let mut dist = Array4::<f64>::from_elem((n, 4, n, 4), f64::MAX / 10.0);
    for edge in graph.all_edges.iter() {
        dist[[edge.fig1, edge.side1, edge.fig2, edge.side2]] = edge.score;
        dist[[edge.fig2, edge.side2, edge.fig1, edge.side1]] = edge.score;
    }
    eprintln!("nd array created!");

    let all_sides = (0..n)
        .cartesian_product(0..4)
        .map(|(fig, side)| Side { fig, side })
        .collect_vec();

    let dist = |s1: Side, s2: Side| -> f64 { dist[[s1.fig, s1.side, s2.fig, s2.side]] };

    const MAX_DIST: f64 = 3.0;
    let mut fours = vec![];
    for &s0 in &all_sides {
        for &s1 in &all_sides {
            let d0 = dist(s0, s1.ne2());
            if d0 > MAX_DIST {
                continue;
            }
            for &s2 in &all_sides {
                let d1 = dist(s0.ne(), s2.pr());
                if d1 > MAX_DIST {
                    continue;
                }
                for &s3 in &all_sides {
                    let d2 = dist(s1.ne(), s3.pr());
                    let d3 = dist(s2, s3.ne2());
                    if d2 > MAX_DIST || d3 > MAX_DIST {
                        continue;
                    }
                    fours.push(Four {
                        s0,
                        s1,
                        s2,
                        s3,
                        max_dist: fmax4([d0, d1, d2, d3]),
                    });
                }
            }
        }
    }
    fours.sort_by(|f1, f2| f1.max_dist.total_cmp(&f2.max_dist));

    let mut rot_by_figure = vec![None; n];

    for four in fours.iter() {
        let mut ok = true;
        for s in four.all() {
            if let Some(rot) = rot_by_figure[s.fig] {
                if rot != s.side {
                    ok = false;
                }
            }
        }
        if ok {
            for s in four.all() {
                rot_by_figure[s.fig] = Some(s.side);
            }
            eprintln!("{:?}", four);
        }
    }
    eprintln!("finished..");
}
