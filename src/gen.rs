use super::table::Table;
use super::Prec;
use rand::distributions::{Distribution, Uniform};

/// Rejection uniform sampling on a sphere of radius `radius`
pub fn uniform(radius: Prec, number: usize) -> Table<(Prec, Prec, Prec)> {
    let mut res = Vec::with_capacity(number);
    let mut rng = rand::thread_rng();

    let radius_squared = radius * radius;
    let range = Uniform::new_inclusive(-radius, radius);
    for _ in 0..number {
        res.push(loop {
            let x = range.sample(&mut rng);
            let y = range.sample(&mut rng);
            let z = range.sample(&mut rng);

            if x * x + y * y + z * z <= radius_squared {
                break (x, y, z);
            }
        });
    }

    Table::new(res)
}

pub fn null(number: usize) -> Table<(Prec, Prec, Prec)> {
    Table::new(vec![(0.0, 0.0, 0.0); number])
}

pub fn scalar_uniform(from: Prec, to: Prec, number: usize) -> Table<Prec> {
    let mut res = Vec::with_capacity(number);
    let mut rng = rand::thread_rng();

    let range = Uniform::new_inclusive(from, to);
    for _ in 0..number {
        res.push(range.sample(&mut rng));
    }

    Table::new(res)
}
