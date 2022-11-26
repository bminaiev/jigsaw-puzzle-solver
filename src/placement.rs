use itertools::Itertools;

use crate::{borders_graph::Graph, utils::Side};

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
struct FigureInfo {
    figure_id: usize,
    positions: [Pos; 4],
    comp_id: usize,
}

#[derive(Clone)]
pub struct Placement {
    figures: Vec<FigureInfo>,
}

const DEFAULT_POS: [Pos; 4] = [
    Pos { x: 0, y: 0 },
    Pos { x: 0, y: 1 },
    Pos { x: 1, y: 1 },
    Pos { x: 1, y: 0 },
];

impl Placement {
    pub fn new() -> Self {
        Self { figures: vec![] }
    }

    pub fn from_full_graph(graph: &Graph) -> Self {
        let mut res = Self::new();
        for edge in graph.all_edges.iter() {
            let (s1, s2) = edge.sides();
            res.join_sides(s1, s2).unwrap();
        }
        res
    }

    pub fn get_bounding_box(&self) -> (i32, i32) {
        let all_pos = self
            .figures
            .iter()
            .map(|f| f.positions.iter())
            .flatten()
            .collect_vec();
        let max_x = all_pos.iter().map(|p| p.x).max().unwrap();
        let min_x = all_pos.iter().map(|p| p.x).min().unwrap();
        let max_y = all_pos.iter().map(|p| p.y).max().unwrap();
        let min_y = all_pos.iter().map(|p| p.y).min().unwrap();
        (max_x - min_x, max_y - min_y)
    }

    fn get_fig_index(&self, fig: usize) -> Option<usize> {
        self.figures
            .binary_search_by(|check| check.figure_id.cmp(&fig))
            .ok()
    }

    fn get_comp_id(&self, fig: usize) -> usize {
        self.figures[self.get_fig_index(fig).unwrap()].comp_id
    }

    pub fn get_fig_comp_size(&self, fig: usize) -> usize {
        let comp_id = self.get_comp_id(fig);
        self.figures.iter().filter(|f| f.comp_id == comp_id).count()
    }

    fn get_side_positions(&self, s: Side) -> [Pos; 2] {
        let all_pos = self.figures[self.get_fig_index(s.fig).unwrap()].positions;
        [all_pos[s.side], all_pos[(s.side + 1) % 4]]
    }

    fn rotate_comp(&mut self, comp_id: usize) {
        for f in self.figures.iter_mut() {
            if f.comp_id == comp_id {
                for p in f.positions.iter_mut() {
                    *p = p.rotate();
                }
            }
        }
    }

    fn rotate_idx(&mut self, idx: usize) {
        for p in self.figures[idx].positions.iter_mut() {
            *p = p.rotate();
        }
    }

    fn shift_comp(&mut self, comp_id: usize, shift: Pos) {
        for f in self.figures.iter_mut() {
            if f.comp_id == comp_id {
                for p in f.positions.iter_mut() {
                    *p = *p + shift;
                }
            }
        }
    }

    fn get_top_left_corner_by_idx(&self, idx: usize) -> Pos {
        *self.figures[idx].positions.iter().min().unwrap()
    }

    fn can_joins_comps(&self, comp1: usize, comp2: usize) -> bool {
        let mut corners = vec![];
        for (idx, f) in self.figures.iter().enumerate() {
            if f.comp_id == comp1 || f.comp_id == comp2 {
                let corner = self.get_top_left_corner_by_idx(idx);
                corners.push(corner);
            }
        }
        corners.sort();
        let sz = corners.len();
        corners.dedup();
        corners.len() == sz
    }

    fn can_joins_comps_one_fig(&self, comp1: usize, idx: usize) -> bool {
        let search_corner = self.get_top_left_corner_by_idx(idx);
        for (idx, f) in self.figures.iter().enumerate() {
            if f.comp_id == comp1 {
                let corner = self.get_top_left_corner_by_idx(idx);
                if corner == search_corner {
                    return false;
                }
            }
        }
        true
    }

    fn gen_edges_between_comps(&self, comp1: usize, comp2: usize) -> Vec<(Side, Side)> {
        let mut edges = vec![];
        let mut all_sides = vec![];
        for fig in self.figures.iter() {
            if fig.comp_id == comp1 {
                for side in 0..4 {
                    let p1 = fig.positions[side];
                    let p2 = fig.positions[(side + 1) % 4];
                    all_sides.push((
                        p1,
                        p2,
                        Side {
                            fig: fig.figure_id,
                            side,
                        },
                    ));
                }
            }
        }
        for fig in self.figures.iter() {
            if fig.comp_id == comp2 {
                for side in 0..4 {
                    let p1 = fig.positions[side];
                    let p2 = fig.positions[(side + 1) % 4];

                    all_sides.push((
                        p2,
                        p1,
                        Side {
                            fig: fig.figure_id,
                            side,
                        },
                    ));
                }
            }
        }

        all_sides.sort();
        for w in all_sides.windows(2) {
            let e1 = w[0];
            let e2 = w[1];
            if e1.0 == e2.0 && e1.1 == e2.1 {
                edges.push((e1.2, e2.2));
            }
        }

        edges
    }

    fn gen_edges_between_comps2(&self, comp1: usize, s_idx: usize) -> Vec<(Side, Side)> {
        let mut edges = vec![];
        let mut all_sides = vec![];
        {
            let fig = &self.figures[s_idx];
            for side in 0..4 {
                let p1 = fig.positions[side];
                let p2 = fig.positions[(side + 1) % 4];
                all_sides.push((
                    p1,
                    p2,
                    Side {
                        fig: fig.figure_id,
                        side,
                    },
                ));
            }
        }
        for fig in self.figures.iter() {
            if fig.comp_id == comp1 {
                for side in 0..4 {
                    let p1 = fig.positions[side];
                    let p2 = fig.positions[(side + 1) % 4];

                    for &(prev_p1, prev_p2, s2) in all_sides.iter() {
                        if prev_p1 == p2 && prev_p2 == p1 {
                            edges.push((
                                Side {
                                    fig: fig.figure_id,
                                    side,
                                },
                                s2,
                            ));
                        }
                    }
                }
            }
        }

        edges
    }

    fn join_comps(&mut self, comp1: usize, comp2: usize) {
        for f in self.figures.iter_mut() {
            if f.comp_id == comp1 {
                f.comp_id = comp2;
            }
        }
    }

    fn force_exist(&mut self, fig_id: usize) {
        if self.get_fig_index(fig_id).is_some() {
            return;
        }
        self.figures.push(FigureInfo {
            figure_id: fig_id,
            positions: DEFAULT_POS.clone(),
            comp_id: self.figures.len(),
        });
        self.figures.sort_by_key(|f| f.figure_id);
    }

    fn remove_single_comps(&mut self) {
        self.figures.sort_by_key(|f| f.comp_id);
        for idx in (0..self.figures.len()).rev() {
            let mut ok = false;
            if idx + 1 < self.figures.len()
                && self.figures[idx + 1].comp_id == self.figures[idx].comp_id
            {
                ok = true;
            }
            if idx > 0 && self.figures[idx - 1].comp_id == self.figures[idx].comp_id {
                ok = true;
            }
            if !ok {
                self.figures.remove(idx);
            }
        }
        self.figures.sort_by_key(|f| f.figure_id);
    }

    fn join_sides_new_fig(&mut self, s1: Side, s2: Side) -> Option<Vec<(Side, Side)>> {
        self.force_exist(s1.fig);
        self.force_exist(s2.fig);
        let s2_idx = self.get_fig_index(s2.fig).unwrap();

        loop {
            let pos1 = self.get_side_positions(s1);
            let pos2 = self.get_side_positions(s2);
            let delta1 = pos1[1] - pos1[0];
            let delta2 = pos2[0] - pos2[1];
            if delta1 == delta2 {
                let shift = pos1[0] - pos2[1];
                self.shift_comp(self.get_comp_id(s2.fig), shift);
                let comp1 = self.get_comp_id(s1.fig);
                if !self.can_joins_comps_one_fig(comp1, s2_idx) {
                    self.figures.remove(s2_idx);
                    return None;
                }
                let new_edges = self.gen_edges_between_comps2(self.get_comp_id(s1.fig), s2_idx);
                self.join_comps(self.get_comp_id(s1.fig), self.get_comp_id(s2.fig));
                return Some(new_edges);
            } else {
                self.rotate_idx(s2_idx);
            }
        }
    }

    // None if can't union.
    // Some list of new edges enforced if united
    pub fn join_sides(&mut self, s1: Side, s2: Side) -> Option<Vec<(Side, Side)>> {
        if self.get_fig_index(s1.fig).is_none() {
            return self.join_sides_new_fig(s2, s1);
        }
        if self.get_fig_index(s2.fig).is_none() {
            return self.join_sides_new_fig(s1, s2);
        }
        self.remove_single_comps();
        self.force_exist(s1.fig);
        self.force_exist(s2.fig);
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

    fn get_all_comp_ids(&self) -> Vec<usize> {
        let mut res = self.figures.iter().map(|f| f.comp_id).collect_vec();
        res.sort();
        res.dedup();
        res
    }

    pub fn get_all_neighbours(&self) -> Vec<(Side, Side)> {
        let mut res = vec![];
        for comp_id in self.get_all_comp_ids() {
            res.extend(self.gen_edges_between_comps(comp_id, comp_id));
        }
        res
    }
}
