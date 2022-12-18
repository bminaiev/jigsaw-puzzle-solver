use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
};

use eframe::epaint::util::hash;
use itertools::Itertools;

use crate::{borders_graph::Graph, utils::Side};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Pos {
    x: i32,
    y: i32,
}

impl Pos {
    pub fn shift(&self, dx: i32, dy: i32) -> Self {
        Self {
            x: self.x + dx,
            y: self.y + dy,
        }
    }
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

fn calc_bbox(all_pos: &[Pos]) -> (Pos, Pos) {
    if all_pos.is_empty() {
        return (Pos { x: 0, y: 0 }, Pos { x: 0, y: 0 });
    }
    let max_x = all_pos.iter().map(|p| p.x).max().unwrap();
    let min_x = all_pos.iter().map(|p| p.x).min().unwrap();
    let max_y = all_pos.iter().map(|p| p.y).max().unwrap();
    let min_y = all_pos.iter().map(|p| p.y).min().unwrap();
    (Pos { x: min_x, y: min_y }, Pos { x: max_x, y: max_y })
}

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

    pub fn get_bounding_box_points(&self) -> (Pos, Pos) {
        let all_pos = self
            .figures
            .iter()
            .map(|f| f.positions.iter())
            .flatten()
            .cloned()
            .collect_vec();
        calc_bbox(&all_pos)
    }

    pub fn get_bounding_box(&self) -> (i32, i32) {
        let (p0, p1) = self.get_bounding_box_points();
        (p1.x - p0.x, p1.y - p0.y)
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
                    eprintln!(
                        "Already exist thing in {:?}, figure_id = {}",
                        corner, f.figure_id
                    );
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
                        0,
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
                        1,
                    ));
                }
            }
        }

        all_sides.sort();
        // Sorry, I am just too lazy to write it nicely.
        for w in all_sides.group_by(|s1, s2| (s1.0, s1.1) == (s2.0, s2.1)) {
            for e1 in w.iter() {
                for e2 in w.iter() {
                    if e1.3 == 0 && e2.3 == 1 {
                        edges.push((e1.2, e2.2));
                    }
                }
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
        let max_comp_id = self.figures.iter().map(|f| f.comp_id).max().unwrap_or(0);
        self.figures.push(FigureInfo {
            figure_id: fig_id,
            positions: DEFAULT_POS.clone(),
            comp_id: max_comp_id + 1,
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

    // each edge is returned two times: (s1, s2) and (s2, s1)
    pub fn get_all_neighbours(&self) -> Vec<(Side, Side)> {
        let mut res = vec![];
        for comp_id in self.get_all_comp_ids() {
            res.extend(self.gen_edges_between_comps(comp_id, comp_id));
        }
        res
    }

    pub fn get_all_neighbours_in_same_component(&self, v: usize) -> Vec<(Side, Side)> {
        let comp_id = self.get_comp_id(v);
        self.gen_edges_between_comps(comp_id, comp_id)
    }

    pub fn get_only_one_component_placement(&self, vertex: usize) -> Self {
        let edges = self.get_all_neighbours_in_same_component(vertex);
        let mut res = Placement::new();
        for &(s1, s2) in edges.iter() {
            res.join_sides(s1, s2).unwrap();
        }
        res
    }

    fn get_corner_by_side(p1: Pos, p2: Pos) -> (Pos, usize) {
        for i in 0..4 {
            let check_p1 = DEFAULT_POS[i];
            let check_p2 = DEFAULT_POS[(i + 1) % 4];
            if check_p1 - check_p2 == p1 - p2 {
                return (p1 - check_p1, i);
            }
        }
        unreachable!();
    }

    fn get_all_corners(&self) -> Vec<Pos> {
        (0..self.figures.len())
            .map(|idx| self.get_top_left_corner_by_idx(idx))
            .collect_vec()
    }

    pub fn get_potential_locations(&self, only_existing_locations: bool) -> Vec<PotentialLocation> {
        let mut potential_locations = vec![];
        for fig in self.figures.iter() {
            for i in 0..4 {
                let p2 = fig.positions[i];
                let p1 = fig.positions[(i + 1) % 4];
                let (corner, side_id) = Self::get_corner_by_side(p1, p2);
                potential_locations.push(PotentialLocationSide {
                    corner,
                    side_id,
                    other_side: Side {
                        fig: fig.figure_id,
                        side: i,
                    },
                });
            }
        }
        potential_locations.sort();

        let mut existing_corners = self.get_all_corners();
        existing_corners.sort();

        let mut existing_iter = 0;

        let mut i = 0;
        let mut res = vec![];
        while i != potential_locations.len() {
            while existing_iter != existing_corners.len()
                && existing_corners[existing_iter] < potential_locations[i].corner
            {
                existing_iter += 1;
            }

            let mut j = i;
            while j != potential_locations.len()
                && potential_locations[j].corner == potential_locations[i].corner
            {
                j += 1;
            }
            let mut neighbors = [None; 4];
            while i != j {
                let pl = potential_locations[i];
                i += 1;
                neighbors[pl.side_id] = Some(pl.other_side);
            }
            let existing_location = existing_iter != existing_corners.len()
                && potential_locations[i - 1].corner == existing_corners[existing_iter];
            if existing_location == only_existing_locations {
                res.push(PotentialLocation {
                    neighbors,
                    pos: potential_locations[i - 1].corner,
                });
            }

            i = j;
        }
        res
    }

    pub fn get_potential_group_locations(&self, max_cnt: u32) -> Vec<PotentialGroupLocation> {
        for i in 1..self.figures.len() {
            assert_eq!(self.figures[i].comp_id, self.figures[0].comp_id);
        }
        let all_neighbours = self.get_all_neighbours();
        let mut res = vec![];
        let locations = self.get_potential_locations(false);
        let (p0, p1) = self.get_bounding_box_points();
        const WINDOW_SIZE: i32 = 3;
        for x0 in p0.x - 1..p1.x {
            for y0 in p0.y - 1..p1.y {
                let mut inside = vec![];
                for loc in locations.iter() {
                    if loc.pos.x >= x0 && loc.pos.x < x0 + WINDOW_SIZE {
                        if loc.pos.y >= y0 && loc.pos.y < y0 + WINDOW_SIZE {
                            inside.push(loc.clone());
                        }
                    }
                }
                for mask in 0u32..(1 << inside.len()) {
                    if mask.count_ones() > 0 && mask.count_ones() <= max_cnt {
                        let mut used_inside = vec![];
                        for i in 0..inside.len() {
                            if ((1 << i) & mask) != 0 {
                                used_inside.push(inside[i].clone());
                            }
                        }
                        let (p0, p1) = calc_bbox(&used_inside.iter().map(|l| l.pos).collect_vec());
                        let p0 = p0.shift(-1, -1);
                        let p1 = p1.shift(1, 1);
                        let existing_figures = (0..self.figures.len())
                            .filter_map(|idx| {
                                let corner = self.get_top_left_corner_by_idx(idx);
                                if corner.x >= p0.x
                                    && corner.x <= p1.x
                                    && corner.y >= p0.y
                                    && corner.y <= p1.y
                                {
                                    Some(self.figures[idx].figure_id)
                                } else {
                                    None
                                }
                            })
                            .collect_vec();
                        let existing_edges = all_neighbours
                            .iter()
                            .filter(|(s1, s2)| {
                                existing_figures.contains(&s1.fig)
                                    && existing_figures.contains(&s2.fig)
                                    && s1.fig < s2.fig
                            })
                            .cloned()
                            .collect_vec();
                        let new_figures_edges = calc_new_figures_edges(&used_inside);
                        if new_figures_edges.len() < used_inside.len() - 1 {
                            // checking that everything is connected.
                            // otherwise this check will not work.
                            assert!(max_cnt <= 3);
                            continue;
                        }
                        // we only care about sides nearby, so we can hash the state
                        for loc in used_inside.iter_mut() {
                            loc.pos = Pos { x: 0, y: 0 };
                        }
                        let potential_group = PotentialGroupLocation {
                            locations: used_inside,
                            existing_edges,
                            new_figures_edges,
                        };
                        if potential_group.each_new_figure_has_at_least_two_edges() {
                            res.push(potential_group);
                        }
                    }
                }
            }
        }
        res.sort();
        res.dedup();
        res
    }

    pub fn get_figure_by_position(&self, top_left: Pos) -> Option<[Side; 4]> {
        for (idx, fig) in self.figures.iter().enumerate() {
            if self.get_top_left_corner_by_idx(idx) == top_left {
                let offset = fig
                    .positions
                    .iter()
                    .position(|pos| *pos == top_left)
                    .unwrap();
                let mut res = [Side::default(); 4];
                for i in 0..4 {
                    res[i] = Side {
                        fig: fig.figure_id,
                        side: (offset + i) % 4,
                    };
                }
                return Some(res);
            }
        }
        None
    }

    pub fn get_all_used_figures(&self) -> Vec<usize> {
        self.figures.iter().map(|f| f.figure_id).collect_vec()
    }

    pub fn get_cnt_figures(&self) -> usize {
        self.figures.len()
    }

    pub fn remove_figure(&mut self, figure: usize) {
        for i in 0..self.figures.len() {
            if self.figures[i].figure_id == figure {
                self.figures.remove(i);
                return;
            }
        }
    }

    pub fn join_with(&mut self, another: &Self) -> bool {
        for (s1, s2) in another.get_all_neighbours() {
            if self.join_sides(s1, s2).is_none() {
                return false;
            }
        }
        true
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PotentialLocationSide {
    corner: Pos,
    side_id: usize,
    other_side: Side,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PotentialLocation {
    pub neighbors: [Option<Side>; 4],
    pub pos: Pos,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NewFigureEdge {
    location_ids: [usize; 2],
    sides: [usize; 2],
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone)]
pub struct PotentialGroupLocation {
    pub locations: Vec<PotentialLocation>,
    pub existing_edges: Vec<(Side, Side)>,
    pub new_figures_edges: Vec<NewFigureEdge>,
}

fn calc_new_figures_edges(locations: &[PotentialLocation]) -> Vec<NewFigureEdge> {
    let mut res = vec![];
    for i in 0..locations.len() {
        for j in i + 1..locations.len() {
            for s1 in 0..4 {
                for s2 in 0..4 {
                    let p1 = locations[i].pos + DEFAULT_POS[s1];
                    let p2 = locations[i].pos + DEFAULT_POS[(s1 + 1) % 4];
                    let p3 = locations[j].pos + DEFAULT_POS[s2];
                    let p4 = locations[j].pos + DEFAULT_POS[(s2 + 1) % 4];
                    if p1 == p4 && p2 == p3 {
                        res.push(NewFigureEdge {
                            location_ids: [i, j],
                            sides: [s1, s2],
                        })
                    }
                }
            }
        }
    }
    res
}

impl PotentialGroupLocation {
    pub fn each_new_figure_has_at_least_two_edges(&self) -> bool {
        for i in 0..self.locations.len() {
            let mut cnt = 0;
            for j in 0..4 {
                if self.locations[i].neighbors[j].is_some() {
                    cnt += 1;
                }
            }
            for edge in self.new_figures_edges.iter() {
                if edge.location_ids.contains(&i) {
                    cnt += 1;
                }
            }
            if cnt < 2 {
                return false;
            }
        }
        true
    }

    pub fn get_way_edges(&self, first_sides: &[Side]) -> Vec<(Side, Side)> {
        // TODO: check two-way.
        let mut res = vec![];
        res.extend(self.existing_edges.clone());
        for i in 0..first_sides.len() {
            for offset in 0..4 {
                if let Some(s) = self.locations[i].neighbors[offset] {
                    res.push((s, first_sides[i].ne_offset(offset)))
                }
            }
        }
        for edge in self.new_figures_edges.iter() {
            let side1 = first_sides[edge.location_ids[0]].ne_offset(edge.sides[0]);
            let side2 = first_sides[edge.location_ids[1]].ne_offset(edge.sides[1]);
            res.push((side1, side2));
        }
        res
    }

    pub fn check_way(
        &self,
        first_sides: &[Side],
        mut dist: impl FnMut(Side, Side) -> f64,
    ) -> Search3StateWithScore {
        let edges = self.get_way_edges(first_sides);
        let score = edges.iter().map(|&(s1, s2)| dist(s1, s2)).sum::<f64>() / (edges.len() as f64);
        Search3StateWithScore {
            state: self.clone(),
            first_sides: first_sides.to_vec(),
            score,
        }
    }
}

#[derive(Clone, PartialEq)]
pub struct Search3StateWithScore {
    pub state: PotentialGroupLocation,
    pub first_sides: Vec<Side>,
    pub score: f64,
}

impl Eq for Search3StateWithScore {}

impl PartialOrd for Search3StateWithScore {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Search3StateWithScore {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.total_cmp(&other.score)
    }
}

impl Search3StateWithScore {
    pub fn get_hash(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.state.hash(&mut hasher);
        self.first_sides.hash(&mut hasher);
        hasher.finish()
    }

    pub fn all_edges(&self) -> Vec<(Side, Side)> {
        self.state.get_way_edges(&self.first_sides)
    }

    pub fn get_key(&self) -> &PotentialGroupLocation {
        &self.state
    }
}
