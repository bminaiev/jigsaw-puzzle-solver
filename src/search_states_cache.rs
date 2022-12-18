use std::{fs, sync::Arc};

use eframe::epaint::{ahash::HashMap, mutex::Mutex};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct SearchStateCacheData {
    hashes: HashMap<u64, f64>,
}

pub struct SearchStatesCache {
    data: Arc<Mutex<SearchStateCacheData>>,
}

const DEFAULT_PATH: &str = "states_cache.json";

impl SearchStatesCache {
    pub fn load() -> Self {
        let data = if let Ok(content) = fs::read_to_string(DEFAULT_PATH) {
            serde_json::from_str(&content).unwrap()
        } else {
            SearchStateCacheData {
                hashes: HashMap::default(),
            }
        };
        Self {
            data: Arc::new(Mutex::new(data)),
        }
    }

    fn save_data(data: &SearchStateCacheData) {
        fs::write(DEFAULT_PATH, serde_json::to_string(data).unwrap()).unwrap();
    }

    pub fn save(&self) {
        eprintln!(
            "Saved bad states cache. {} entries",
            self.data.lock().hashes.len()
        );
        Self::save_data(&self.data.lock());
    }

    pub fn contains(&self, hash: u64, bound: f64) -> bool {
        *self.data.lock().hashes.get(&hash).unwrap_or(&0.0) > bound
    }

    pub fn insert(&self, hash: u64, value: f64) {
        self.data.lock().hashes.insert(hash, value);
    }
}
