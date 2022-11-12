use std::collections::{HashMap, HashSet};

use itertools::Itertools;

use crate::utils::Side;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Pos {
    x: i32,
    y: i32,
}

impl std::ops::Sub for Pos {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl std::ops::Add for Pos {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Pos {
    pub fn rotate(&self) -> Self {
        Self {
            x: self.y,
            y: -self.x,
        }
    }
}

#[derive(Clone)]
pub struct Placement {
    positions: Vec<[Pos; 4]>,
    components: Vec<Vec<usize>>,
}

impl Placement {
    pub fn new(n: usize) -> Self {
        let default_pos = [
            Pos { x: 0, y: 0 },
            Pos { x: 0, y: 1 },
            Pos { x: 1, y: 1 },
            Pos { x: 1, y: 0 },
        ];
        let components = (0..n).map(|i| vec![i]).collect_vec();
        Self {
            positions: vec![default_pos; n],
            components,
        }
    }

    fn get_comp_id(&self, fig: usize) -> usize {
        for (comp_id, comp) in self.components.iter().enumerate() {
            if comp.contains(&fig) {
                return comp_id;
            }
        }
        unreachable!("Figure doesn't exist");
    }

    fn get_side_positions(&self, s: Side) -> [Pos; 2] {
        let all_pos = self.positions[s.fig];
        [all_pos[s.side], all_pos[(s.side + 1) % 4]]
    }

    fn rotate_comp(&mut self, comp_id: usize) {
        for id in self.components[comp_id].clone().into_iter() {
            for p in self.positions[id].iter_mut() {
                *p = p.rotate();
            }
        }
    }

    fn shift_comp(&mut self, comp_id: usize, shift: Pos) {
        for id in self.components[comp_id].clone().into_iter() {
            for p in self.positions[id].iter_mut() {
                *p = *p + shift;
            }
        }
    }

    fn get_top_left_corner(&self, id: usize) -> Pos {
        *self.positions[id].iter().min().unwrap()
    }

    fn can_joins_comps(&self, comp1: usize, comp2: usize) -> bool {
        let mut used = HashSet::new();
        for &id in self.components[comp1].iter() {
            used.insert(self.get_top_left_corner(id));
        }
        for &id in self.components[comp2].iter() {
            if used.contains(&self.get_top_left_corner(id)) {
                return false;
            }
        }
        true
    }

    fn gen_edges_between_comps(&self, comp1: usize, comp2: usize) -> Vec<(Side, Side)> {
        let mut edges = vec![];
        let mut sides = HashMap::new();
        for &fig in self.components[comp1].iter() {
            for side in 0..4 {
                let p1 = self.positions[fig][side];
                let p2 = self.positions[fig][(side + 1) % 4];
                sides.insert((p1, p2), Side { fig, side });
            }
        }

        for &fig in self.components[comp2].iter() {
            for side in 0..4 {
                let p1 = self.positions[fig][side];
                let p2 = self.positions[fig][(side + 1) % 4];

                if let Some(existing_side) = sides.get(&(p2, p1)) {
                    edges.push((*existing_side, Side { fig, side }));
                }
            }
        }

        edges
    }

    fn join_comps(&mut self, comp1: usize, comp2: usize) {
        let old_comp2 = std::mem::replace(&mut self.components[comp2], vec![]);
        self.components[comp1].extend(old_comp2);
    }

    // None if can't union.
    // Some list of new edges enforced if united
    pub fn join_sides(&mut self, s1: Side, s2: Side) -> Option<Vec<(Side, Side)>> {
        let comp1 = self.get_comp_id(s1.fig);
        let comp2 = self.get_comp_id(s2.fig);
        if comp1 == comp2 {
            let pos1 = self.get_side_positions(s1);
            let pos2 = self.get_side_positions(s2);
            if pos1[0] == pos2[1] && pos1[1] == pos2[0] {
                return Some(vec![(s1, s2)]);
            } else {
                return None;
            }
        } else {
            loop {
                let pos1 = self.get_side_positions(s1);
                let pos2 = self.get_side_positions(s2);
                let delta1 = pos1[1] - pos1[0];
                let delta2 = pos2[0] - pos2[1];
                if delta1 == delta2 {
                    let shift = pos1[0] - pos2[1];
                    self.shift_comp(self.get_comp_id(s2.fig), shift);
                    let comp1 = self.get_comp_id(s1.fig);
                    let comp2 = self.get_comp_id(s2.fig);
                    if !self.can_joins_comps(comp1, comp2) {
                        return None;
                    }
                    let new_edges = self.gen_edges_between_comps(
                        self.get_comp_id(s1.fig),
                        self.get_comp_id(s2.fig),
                    );
                    self.join_comps(self.get_comp_id(s1.fig), self.get_comp_id(s2.fig));
                    return Some(new_edges);
                } else {
                    self.rotate_comp(self.get_comp_id(s2.fig));
                }
            }
        }
    }

    pub fn get_all_neighbours(&self) -> Vec<(Side, Side)> {
        let mut res = vec![];
        for comp_id in 0..self.components.len() {
            res.extend(self.gen_edges_between_comps(comp_id, comp_id));
        }
        res
    }
}
