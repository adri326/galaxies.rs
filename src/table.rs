use scoped_threadpool::Pool;
use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct Table<T> {
    pub array: Vec<T>,
}

impl<T> Table<T> {
    pub fn new(array: Vec<T>) -> Self {
        Self {
            array
        }
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

impl<T: Send + Clone + Sync> Table<T> {
    pub fn foreach_threaded(&mut self, n_threads: u32, f: impl (Fn(usize, &T, &Vec<T>) -> T) + Send + Copy) {
        let mut pool = Pool::new(n_threads);

        let res = Mutex::new(vec![None; self.array.len()]);
        pool.scoped(|scope| {
            let res = &res;
            let array = &self.array;
            for (index, elem) in self.array.iter().enumerate() {
                scope.execute(move || {
                    let value = f(index, elem, array);
                    *res.lock().unwrap().get_mut(index).unwrap() = Some(value);
                });
            }
        });
        self.array = res.into_inner().unwrap().into_iter().map(|x| x.unwrap()).collect();
    }
}
