use std::{cmp::Ordering, collections::BTreeSet};

#[derive(Clone)]
pub struct TopN<T: Clone + Ord> {
    pub set: BTreeSet<T>,
    max_n: usize,
}

impl<T: Clone + Ord> TopN<T> {
    pub fn new(max_n: usize) -> Self {
        Self {
            set: BTreeSet::default(),
            max_n,
        }
    }

    pub fn insert(&mut self, elem: T) {
        if self.set.len() == self.max_n {
            if elem.cmp(self.set.iter().next().unwrap()) != Ordering::Greater {
                return;
            }
        }
        self.set.insert(elem);
        if self.set.len() > self.max_n {
            self.set.remove(&self.set.iter().next().unwrap().clone());
        }
    }

    pub(crate) fn get_best(&self) -> Option<T> {
        self.set.iter().next_back().cloned()
    }
}
