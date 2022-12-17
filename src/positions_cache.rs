use std::{
    fs,
    sync::{Arc, Mutex},
};

use ndarray::Array4;
use serde::{Deserialize, Serialize};

use crate::{
    parsed_puzzles::ParsedPuzzles, point::PointF, surface_placer::place_one_connected_component,
    utils::Side,
};

#[derive(Serialize, Deserialize)]
pub struct PositionsCacheData {
    pts: Vec<Option<Vec<PointF>>>,
}

pub struct PositionsCache {
    data: Arc<Mutex<PositionsCacheData>>,
}

const DEFAULT_PATH: &str = "positions_cache.json";

const CACHE_GRAPH_SIZE: usize = 15;

impl PositionsCache {
    pub fn load(parsed_puzzles: &ParsedPuzzles) -> Self {
        let data = if let Ok(content) = fs::read_to_string(DEFAULT_PATH) {
            serde_json::from_str(&content).unwrap()
        } else {
            PositionsCacheData {
                pts: vec![None; parsed_puzzles.figures.len()],
            }
        };
        Self {
            data: Arc::new(Mutex::new(data)),
        }
    }

    fn save(data: &PositionsCacheData) {
        fs::write(DEFAULT_PATH, serde_json::to_string(data).unwrap()).unwrap();
    }

    pub fn place_one_connected_component_with_cache(
        &self,
        parsed_puzzles: &ParsedPuzzles,
        component: &[usize],
        used_edges: &[(Side, Side)],
        positions: &mut Vec<Option<Vec<PointF>>>,
        base_points_matrix: &Array4<[PointF; 2]>,
    ) -> Option<f64> {
        if component.len() < CACHE_GRAPH_SIZE {
            return place_one_connected_component(
                parsed_puzzles,
                component,
                used_edges,
                positions,
                base_points_matrix,
            );
        }
        let mut data = self.data.lock().unwrap();
        eprintln!("Will use cache!");
        for &fig in component.iter() {
            positions[fig] = data.pts[fig].clone();
        }
        let res = place_one_connected_component(
            parsed_puzzles,
            component,
            used_edges,
            positions,
            base_points_matrix,
        );
        for &fig in component.iter() {
            data.pts[fig] = positions[fig].clone();
        }
        Self::save(&data);
        res
    }
}
