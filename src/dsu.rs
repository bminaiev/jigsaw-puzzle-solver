pub struct Dsu {
    parent: Vec<usize>,
}

impl Dsu {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    pub fn get(&mut self, v: usize) -> usize {
        if self.parent[v] == v {
            return v;
        }
        self.parent[v] = self.get(self.parent[v]);
        self.parent[v]
    }

    pub fn unite(&mut self, mut x: usize, mut y: usize) {
        x = self.get(x);
        y = self.get(y);
        self.parent[x] = y;
    }
}
