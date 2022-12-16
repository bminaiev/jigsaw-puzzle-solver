use std::fs;

use serde::{Deserialize, Serialize};

use crate::utils::Side;

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug)]
pub struct Fact {
    pub side1: Side,
    pub side2: Side,
    pub good_edge: bool,
}

impl Fact {
    pub fn new(side1: Side, side2: Side, good_edge: bool) -> Self {
        if side1 < side2 {
            Self {
                side1,
                side2,
                good_edge,
            }
        } else {
            Self {
                side1: side2,
                side2: side1,
                good_edge,
            }
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct KnownFacts {
    pub facts: Vec<Fact>,
}

const DEFAULT_PATH: &str = "facts.json";

impl KnownFacts {
    pub fn load() -> Self {
        if let Ok(content) = fs::read_to_string(DEFAULT_PATH) {
            serde_json::from_str(&content).unwrap()
        } else {
            Self { facts: vec![] }
        }
    }

    pub fn save(&self) {
        fs::write(DEFAULT_PATH, serde_json::to_string(self).unwrap()).unwrap();
    }

    pub fn add_fact(&mut self, fact: &Fact) {
        if self.facts.contains(&fact) {
            return;
        }
        eprintln!("ADD NEW FACT: {:?}", fact);
        for i in 0..self.facts.len() {
            if self.facts[i].side1 == fact.side1 && self.facts[i].side2 == fact.side2 {
                self.facts.remove(i);
                break;
            }
        }
        self.facts.push(fact.clone());
        self.save();
    }

    pub fn get_all_placed_vertices(&self) -> Vec<usize> {
        let mut res = vec![];
        for fact in self.facts.iter() {
            res.push(fact.side1.fig);
            res.push(fact.side2.fig);
        }
        res.sort();
        res.dedup();
        res
    }

    pub fn remove_fact(&mut self, fact: &Fact) {
        if let Some(idx) = self.facts.iter().position(|f| f == fact) {
            self.facts.remove(idx);
            self.save();
        }
    }

    pub fn get_edge_state(&self, side1: Side, side2: Side) -> EdgeState {
        if self.facts.contains(&Fact::new(side1, side2, true)) {
            return EdgeState::GoodEdge;
        }
        if self.facts.contains(&Fact::new(side1, side2, false)) {
            return EdgeState::WrongEdge;
        }
        EdgeState::Unknown
    }
}

#[derive(Clone, PartialEq, Eq)]
pub enum EdgeState {
    Unknown,
    GoodEdge,
    WrongEdge,
}
