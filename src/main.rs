extern crate minifb;
extern crate rand;
extern crate raqote;
extern crate scoped_threadpool;
extern crate noise;

pub mod gen;
pub mod table;
pub mod voxel;
pub mod draw;

use table::Table;
use voxel::VoxelGrid;
use std::time::Instant;
use std::sync::{mpsc, Mutex};
use std::thread;
use noise::{NoiseFn, Perlin};
use rand::Rng;

pub type Prec = f32;

const N_ELEMS: usize = 100000;
const N_THREADS: u32 = 15;
const GRID_SCALE: Prec = 2000.0;
const SUB_GRID_SCALE: Prec = 150.0;
const PRECISION_RADIUS: Prec = 100.0;
const PRECISION_RADIUS_SQUARED: Prec = PRECISION_RADIUS * PRECISION_RADIUS;
const SUB_GRID_RADIUS: Prec = 700.0;
const SUB_GRID_RADIUS_SQUARED: Prec = SUB_GRID_RADIUS * SUB_GRID_RADIUS;
const G: Prec = 20.0;
const DT: Prec = 0.01;
const COLLISION_DIST: Prec = 200.0;
const GONE_RADIUS: Prec = 2.0 * RADIUS;
const DISTANCE_EPSILON: Prec = 0.01;
const SMOOTHING_SQUARED: Prec = 100.0;

const STEPS_PER_FRAME: usize = 9;
const ACCEL_MID_FREQ: usize = 3;
const RENDER_EVERY_N: usize = 2;

const RADIUS: Prec = 7500.0;
const RANDOM_SPEED: Prec = 5.0;
const BLACK_HOLE_MASS: Prec = 3000000.0;
const GALAXY_WEIGHT_FACTOR: Prec = 1.1;
const Z_SCALE: Prec = 0.3;

// const NORM_MULT: Prec = 0.001;
const WIDTH: usize = 1920;
const HEIGHT: usize = 1080;
const SCALE: Prec = 6000.0;
const MIN_X: Prec = -SCALE * WIDTH as Prec / HEIGHT as Prec;
const MAX_X: Prec = SCALE * WIDTH as Prec / HEIGHT as Prec;
const MIN_Y: Prec = -SCALE;
const MAX_Y: Prec = SCALE;

const OUTPUT_WINDOW: bool = true;
const OUTPUT_FILE: bool = true;

const NOISE_POS_SCALE: f64 = 2000.0;
const NOISE_POS_AMP: Prec = 1000.0;

const NOISE_MASS_SCALE: f64 = 200.0;
const NOISE_MASS_AMP: Prec = 2.0;

const MASSIVE_STAR_RARITY: u32 = 100;
const MASSIVE_STAR_MASS: Prec = 15.0;

fn main() {
    let (tx, rx) = mpsc::channel();
    let perlin = Perlin::new();
    let mut position = gen::uniform(RADIUS, N_ELEMS);
    position.foreach(|_, p, _| (
        p.0 + NOISE_POS_AMP * perlin.get([p.0 as f64 / NOISE_POS_SCALE, p.1 as f64 / NOISE_POS_SCALE, p.2 as f64 / NOISE_POS_SCALE, 0.0]) as Prec,
        p.1 + NOISE_POS_AMP * perlin.get([p.0 as f64 / NOISE_POS_SCALE, p.1 as f64 / NOISE_POS_SCALE, p.2 as f64 / NOISE_POS_SCALE, 2.0]) as Prec,
        (p.2 + NOISE_POS_AMP * perlin.get([p.0 as f64 / NOISE_POS_SCALE, p.1 as f64 / NOISE_POS_SCALE, p.2 as f64 / NOISE_POS_SCALE, -2.0]) as Prec) * Z_SCALE
    ));

    let mut mass = gen::scalar_uniform(0.5, 1.5, N_ELEMS);
    mass.foreach(|i, m, _| {
        let mut rng = rand::thread_rng();
        let mut m = *m;
        let p = position.array[i];
        if rng.gen_ratio(1, MASSIVE_STAR_RARITY) {
            m += MASSIVE_STAR_MASS;
        }
        m + NOISE_MASS_AMP * perlin.get([p.0 as f64 / NOISE_MASS_SCALE, p.1 as f64 / NOISE_MASS_SCALE, p.2 as f64 / NOISE_MASS_SCALE, 6.0]) as Prec
    });

    let mut speed = gen::uniform(RANDOM_SPEED, N_ELEMS);
    speed.foreach(|i, s, _| {
        let p = position.get(i).unwrap();
        let n = Prec::sqrt(p.0 * p.0 + p.1 * p.1);
        let alpha = Prec::atan2(p.1 / n, p.0 / n);
        let distance = Prec::sqrt(distance_squared(*p, (0.0, 0.0, 0.0)));
        let distance_sqrt = Prec::sqrt(distance);
        (
            s.0 - Prec::sqrt(BLACK_HOLE_MASS * G + distance / RADIUS * N_ELEMS as Prec * GALAXY_WEIGHT_FACTOR) * Prec::sin(alpha) / distance_sqrt,
            s.1 + Prec::sqrt(BLACK_HOLE_MASS * G + distance / RADIUS * N_ELEMS as Prec * GALAXY_WEIGHT_FACTOR) * Prec::cos(alpha) / distance_sqrt,
            s.2,
        )
    });
    let mut accel_close = gen::null(N_ELEMS);
    let mut accel_mid = gen::null(N_ELEMS);
    let mut accel_far = gen::null(N_ELEMS);

    let mut jerk_mid = gen::null(N_ELEMS);
    let mut jerk_far = gen::null(N_ELEMS);

    position.array.push((0.0, 0.0, 0.0));
    mass.array.push(BLACK_HOLE_MASS);
    speed.array.push((0.0, 0.0, 0.0));
    accel_close.array.push((0.0, 0.0, 0.0));
    accel_mid.array.push((0.0, 0.0, 0.0));
    jerk_mid.array.push((0.0, 0.0, 0.0));
    accel_far.array.push((0.0, 0.0, 0.0));
    jerk_far.array.push((0.0, 0.0, 0.0));

    let mut previous_position = position.array.clone();

    let mut voxels = VoxelGrid::new(GRID_SCALE, SUB_GRID_SCALE, GONE_RADIUS);
    voxels.populate(&position);
    voxels.update(&position, &speed, &mass);

    let mut step = 0;

    thread::spawn(move || draw::draw(rx));

    loop {
        let start = Instant::now();
        step += 1;
        // Advance N step

        // Half-step for the position before calculating the forces for the furthest stars
        position.foreach(|i, prev, _| {
            add(*prev, mul(speed.array[i], DT / 2.0))
        });

        accel_far.foreach_threaded(N_THREADS, |i, _, _| {
            let my_position = position.array[i];
            let mut res = (0.0, 0.0, 0.0);

            if mass.array[i] == 0.0 {
                return (0.0, 0.0, 0.0)
            }

            for voxel in voxels.iter() {
                let voxel_distance_squared = distance_squared(my_position, voxel.position);
                if voxel_distance_squared > SUB_GRID_RADIUS_SQUARED && voxel.total_mass > 0.0 {
                    let distance_squared = distance_squared(my_position, voxel.center_of_mass);
                    let norm = voxel.total_mass / distance_squared; // F/m
                    if distance_squared < DISTANCE_EPSILON {
                        continue
                    }
                    if norm.is_nan() || norm.is_infinite() {
                        panic!("{} {:?} {:?} {:?} {:?}", norm, voxel.total_mass, distance_squared, voxel.center_of_mass, my_position);
                    }
                    res = add(res, normalize_difference(my_position, voxel.position, norm));
                }
            }

            mul(res, G)
        });

        jerk_far.foreach_threaded(N_THREADS, |i, _, _| {
            let my_position = position.array[i];
            let my_speed = speed.array[i];
            let mut res = (0.0, 0.0, 0.0);

            if mass.array[i] == 0.0 {
                return (0.0, 0.0, 0.0)
            }

            for voxel in voxels.iter() {
                let voxel_distance_squared = distance_squared(my_position, voxel.position);
                if voxel_distance_squared > SUB_GRID_RADIUS_SQUARED && voxel.total_mass > 0.0 {
                    let distance_squared = distance_squared(my_position, voxel.center_of_mass);
                    if distance_squared < DISTANCE_EPSILON {
                        continue
                    }
                    let distance_cubed = distance_squared * distance_squared.sqrt();
                    let speed_diff = sub(my_speed, voxel.average_speed);
                    let pos_diff = sub(my_position, voxel.center_of_mass);

                    res = add(res, mul(speed_diff, 1.0 / distance_cubed));

                    let speed_pos = prod(speed_diff, pos_diff);
                    res = add(res, normalize_difference(my_position, voxel.position, 3.0 * speed_pos / distance_squared / distance_squared));

                    if
                        res.0.is_nan() || res.0.is_infinite() || res.0.abs() > 1000000.0
                        || res.1.is_nan() || res.1.is_infinite() || res.1.abs() > 1000000.0
                        || res.2.is_nan() || res.2.is_infinite() || res.2.abs() > 1000000.0
                    {
                        panic!(
                            "{:?} {:?} {:?} {:?} {:?} {:?} {:?}",
                            res,
                            voxel.total_mass,
                            distance_squared,
                            voxel.center_of_mass,
                            my_position,
                            voxel.average_speed,
                            my_speed,
                        );
                    }
                }
            }

            mul(res, -G)
        });

        for n_step in 0..STEPS_PER_FRAME {
            // Half-step for xₜ, except on the first loop, as this half-step was already done while calculating the forces for the furthest stars
            if n_step > 0 {
                position.foreach(|i, prev, _| {
                    add(*prev, mul(speed.array[i], DT / 2.0))
                });
            }

            if n_step % ACCEL_MID_FREQ == 0 {
                accel_mid.foreach_threaded(N_THREADS, |i, _, _| {
                    let my_position = position.array[i];
                    let my_coords = voxel::get_coords(GRID_SCALE, my_position);
                    let mut res = (0.0, 0.0, 0.0);

                    if mass.array[i] == 0.0 {
                        return (0.0, 0.0, 0.0)
                    }

                    for voxel in voxels.iter() {
                        if my_coords == voxel.coords || distance_squared(my_position, voxel.position) <= SUB_GRID_RADIUS_SQUARED {
                            let my_coords = voxel::get_coords(SUB_GRID_SCALE, my_position);
                            for voxel in voxel.children.iter() {
                                if voxel.total_mass > 0.0 && distance_squared(my_position, voxel.position) > PRECISION_RADIUS_SQUARED && my_coords != voxel.coords {
                                    let distance_squared = distance_squared(my_position, add(voxel.center_of_mass, mul(voxel.average_speed, DT * n_step as Prec)));
                                    if distance_squared < DISTANCE_EPSILON {
                                        continue
                                    }
                                    let norm = voxel.total_mass / distance_squared; // F/m
                                    if norm.is_nan() || norm.is_infinite() {
                                        panic!("{} {:?} {:?} {:?} {:?}", norm, voxel.total_mass, distance_squared, voxel.center_of_mass, my_position);
                                    }
                                    res = add(res, normalize_difference(my_position, voxel.position, norm));
                                }
                            }
                        }
                    }

                    mul(res, G)
                });

                // TODO: move this up before the x-half step?
                jerk_mid.foreach_threaded(N_THREADS, |i, _, _| {
                    let my_position = position.array[i];
                    let my_speed = speed.array[i];
                    let my_coords = voxel::get_coords(GRID_SCALE, my_position);
                    let mut res = (0.0, 0.0, 0.0);

                    if mass.array[i] == 0.0 {
                        return (0.0, 0.0, 0.0)
                    }

                    for voxel in voxels.iter() {
                        if my_coords == voxel.coords || distance_squared(my_position, voxel.position) <= SUB_GRID_RADIUS_SQUARED {
                            let my_coords = voxel::get_coords(SUB_GRID_SCALE, my_position);
                            for voxel in voxel.children.iter() {
                                if voxel.total_mass > 0.0 && distance_squared(my_position, voxel.position) > PRECISION_RADIUS_SQUARED && my_coords != voxel.coords {
                                    let distance_squared = distance_squared(my_position, voxel.center_of_mass);
                                    if distance_squared < DISTANCE_EPSILON {
                                        continue
                                    }
                                    let distance_cubed = distance_squared * distance_squared.sqrt();
                                    let speed_diff = sub(my_speed, voxel.average_speed);
                                    let pos_diff = sub(my_position, voxel.center_of_mass);

                                    res = add(res, mul(speed_diff, 1.0 / distance_cubed));

                                    let speed_pos = prod(speed_diff, pos_diff);
                                    res = add(res, normalize_difference(my_position, voxel.position, 3.0 * speed_pos / distance_squared / distance_squared));

                                    if
                                        res.0.is_nan() || res.0.is_infinite()
                                        || res.1.is_nan() || res.1.is_infinite()
                                        || res.2.is_nan() || res.2.is_infinite()
                                    {
                                        panic!(
                                            "{:?} {:?} {:?} {:?} {:?} {:?} {:?}",
                                            res,
                                            voxel.total_mass,
                                            distance_squared,
                                            voxel.center_of_mass,
                                            my_position,
                                            voxel.average_speed,
                                            my_speed,
                                        );
                                    }
                                }
                            }
                        }
                    }

                    mul(res, -G)
                });
            }

            let collisions = Mutex::new(Vec::new());
            {
                let collisions = &collisions;
                accel_close.foreach_threaded(N_THREADS, |i, _, _| {
                    let my_position = position.array[i];
                    let my_coords = voxel::get_coords(GRID_SCALE, my_position);
                    let mut res = (0.0, 0.0, 0.0);

                    if mass.array[i] == 0.0 {
                        return (0.0, 0.0, 0.0)
                    }

                    for voxel in voxels.iter() {
                        if my_coords == voxel.coords || distance_squared(my_position, voxel.position) <= SUB_GRID_RADIUS_SQUARED {
                            let my_coords = voxel::get_coords(SUB_GRID_SCALE, my_position);
                            for voxel in voxel.children.iter() {
                                if my_coords == voxel.coords || distance_squared(my_position, voxel.position) <= PRECISION_RADIUS_SQUARED {
                                    for &index in voxel.elements.iter() {
                                        let their_position = position.array[index];
                                        let mut distance_squared = distance_squared(my_position, their_position) + SMOOTHING_SQUARED;
                                        if mass.array[index] == 0.0 {
                                            continue
                                        }
                                        if distance_squared < COLLISION_DIST * mass.array[i].max(mass.array[index]).sqrt() {
                                            if mass.array[i] > mass.array[index] {
                                                collisions.lock().unwrap().push((i, index));
                                            }
                                            return (0.0, 0.0, 0.0)
                                        }
                                        if distance_squared < DISTANCE_EPSILON {
                                            continue
                                        }
                                        let their_mass = mass.array[index];
                                        let norm = G * their_mass / distance_squared; // F/m
                                        if norm.is_nan() || norm.is_infinite() {
                                            panic!("{} {:?} {:?} {:?} {:?}", norm, their_mass, distance_squared, their_position, my_position);
                                        }
                                        res = add(res, normalize_difference(my_position, their_position, norm));
                                    }
                                }
                            }
                        }
                    }

                    res
                });
            }

            for (swallowing, swallowed) in collisions.into_inner().unwrap() {
                // println!("{} and {} collided!", swallowing, swallowed);
                let mass_a = mass.array[swallowing];
                let mass_b = mass.array[swallowed];
                let coeff_a = mass_a / (mass_a + mass_b);
                let coeff_b = mass_b / (mass_a + mass_b);
                accel_far.array[swallowing] = add(mul(accel_far.array[swallowing], coeff_a), mul(accel_far.array[swallowed], coeff_b));
                accel_mid.array[swallowing] = add(mul(accel_mid.array[swallowing], coeff_a), mul(accel_mid.array[swallowed], coeff_b));
                accel_close.array[swallowing] = add(mul(accel_close.array[swallowing], coeff_a), mul(accel_close.array[swallowed], coeff_b));
                jerk_far.array[swallowing] = add(mul(jerk_far.array[swallowing], coeff_a), mul(jerk_far.array[swallowed], coeff_b));
                jerk_mid.array[swallowing] = add(mul(jerk_mid.array[swallowing], coeff_a), mul(jerk_mid.array[swallowed], coeff_b));
                speed.array[swallowing] = add(mul(speed.array[swallowing], coeff_a), mul(speed.array[swallowed], coeff_b));
                position.array[swallowing] = add(mul(position.array[swallowing], coeff_a), mul(position.array[swallowed], coeff_b));
                mass.array[swallowing] += mass.array[swallowed];
                mass.array[swallowed] = 0.0;
            }

            speed.foreach_threaded(N_THREADS, |i, prev, _| {
                if mass.array[i] == 0.0 {
                    return (0.0, 0.0, 0.0)
                }

                let mut res = accel_close.array[i];
                res = add(res, accel_mid.array[i]);
                res = add(res, mul(jerk_mid.array[i], (n_step % ACCEL_MID_FREQ) as Prec * DT));
                res = add(res, accel_far.array[i]);
                res = add(res, mul(jerk_far.array[i], n_step as Prec * DT));

                add(*prev, mul(res, DT))
            });

            position.foreach(|i, prev, _| {
                add(*prev, mul(speed.array[i], DT / 2.0))
            });

            voxels.update(&position, &speed, &mass);
        }

        voxels.relocate(&position);
        voxels.gc();
        voxels.update(&position, &speed, &mass);
        println!("{:?}", start.elapsed());
        if step % RENDER_EVERY_N == 0 {
            tx.send((position.clone(), previous_position, mass.iter().map(|&x| x > 0.0).collect())).unwrap();
            previous_position = position.array.clone();
        }
    }
}

#[inline]
fn distance_squared(a: (Prec, Prec, Prec), b: (Prec, Prec, Prec)) -> Prec {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let dz = b.2 - a.2;
    dx * dx + dy * dy + dz * dz
}

#[inline]
fn normalize_difference(
    a: (Prec, Prec, Prec),
    b: (Prec, Prec, Prec),
    unit: Prec,
) -> (Prec, Prec, Prec) {
    if unit == 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let dz = b.2 - a.2;
    let d2 = dx * dx + dy * dy + dz * dz;

    if d2 == 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let d = unit / Prec::sqrt(d2);

    (dx * d, dy * d, dz * d)
}

#[inline]
fn add(a: (Prec, Prec, Prec), b: (Prec, Prec, Prec)) -> (Prec, Prec, Prec) {
    (a.0 + b.0, a.1 + b.1, a.2 + b.2)
}

#[inline]
fn mul(a: (Prec, Prec, Prec), by: Prec) -> (Prec, Prec, Prec) {
    (a.0 * by, a.1 * by, a.2 * by)
}

#[inline]
fn norm(a: (Prec, Prec, Prec)) -> Prec {
    Prec::sqrt(a.0 * a.0 + a.1 * a.1 + a.2 * a.2)
}

fn sigma(x: Prec) -> Prec {
    1.0 / (1.0 + Prec::exp(-x))
}

#[inline]
fn sub(a: (Prec, Prec, Prec), b: (Prec, Prec, Prec)) -> (Prec, Prec, Prec) {
    (
        a.0 - b.0,
        a.1 - b.1,
        a.2 - b.2,
    )
}

#[inline]
fn prod(a: (Prec, Prec, Prec), b: (Prec, Prec, Prec)) -> Prec {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}
