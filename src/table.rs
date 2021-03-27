use scoped_threadpool::Pool;
use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct Table<T> {
    pub array: Vec<T>,
}

impl<T> Table<T> {
    pub fn new(array: Vec<T>) -> Self {
        Self { array }
    }

    pub fn foreach(&mut self, f: impl Fn(usize, &T, &Vec<T>) -> T) {
        let mut res = Vec::with_capacity(self.array.len());
        for (index, elem) in self.array.iter().enumerate() {
            res.push(f(index, elem, &self.array));
        }
        self.array = res;
    }

    pub fn get(&self, index: usize) -> Option<&T> {
        self.array.get(index)
    }

    pub fn iter<'a>(&'a self) -> std::slice::Iter<'a, T> {
        self.array.iter()
    }
}

const GROUP_BY: usize = 50;

impl<T: Send + Clone + Sync> Table<T> {
    pub fn foreach_threaded(
        &mut self,
        n_threads: u32,
        f: impl (Fn(usize, &T, &Vec<T>) -> T) + Send + Copy,
    ) {
        let mut pool = Pool::new(n_threads);

        let mut res = Vec::with_capacity(self.array.len());
        for _ in 0..self.array.len() {
            res.push(None);
        }
        let res = Mutex::new(res);
        pool.scoped(|scope| {
            let res = &res;
            let array = &self.array;
            let mut next_pool = Vec::with_capacity(GROUP_BY);
            for (index, elem) in self.array.iter().enumerate() {
                next_pool.push((index, elem));
                if index % GROUP_BY == GROUP_BY - 1 {
                    scope.execute(move || {
                        let mut results = Vec::new();
                        for (index, elem) in next_pool {
                            results.push((index, f(index, elem, array)));
                        }
                        if let Ok(ref mut guard) = res.lock() {
                            for (i, r) in results {
                                guard[i] = Some(r);
                            }
                        }
                    });
                    next_pool = Vec::with_capacity(GROUP_BY);
                }
            }
            scope.execute(move || {
                let mut results = Vec::new();
                for (index, elem) in next_pool {
                    results.push((index, f(index, elem, array)));
                }
                if let Ok(ref mut guard) = res.lock() {
                    for (i, r) in results {
                        guard[i] = Some(r);
                    }
                }
            });
        });
        self.array = res
            .into_inner()
            .unwrap()
            .into_iter()
            .map(|x| x.unwrap())
            .collect();
    }
}
