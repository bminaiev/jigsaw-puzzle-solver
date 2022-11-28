use itertools::Itertools;

pub struct Dsu {
    parent: Vec<usize>,
    sizes: Vec<usize>,
}

impl Dsu {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            sizes: vec![1; n],
        }
    }

    pub fn get(&mut self, v: usize) -> usize {
        if self.parent[v] == v {
            return v;
        }
        self.parent[v] = self.get(self.parent[v]);
        self.parent[v]
    }

    pub fn unite(&mut self, mut x: usize, mut y: usize) -> bool {
        x = self.get(x);
        y = self.get(y);
        if x == y {
            return false;
        }
        self.parent[x] = y;
        self.sizes[y] += self.sizes[x];
        true
    }

    pub fn get_components(&mut self) -> Vec<Vec<usize>> {
        let mut res = vec![vec![]; self.parent.len()];
        for v in 0..self.parent.len() {
            res[self.get(v)].push(v);
        }
        res.into_iter().filter(|v| !v.is_empty()).collect_vec()
    }

    pub fn get_size(&mut self, v: usize) -> usize {
        let v = self.get(v);
        self.sizes[v]
    }
}
