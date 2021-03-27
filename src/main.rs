extern crate minifb;
extern crate rand;
extern crate raqote;
extern crate scoped_threadpool;

pub mod gen;
pub mod table;
pub mod voxel;

use minifb::{MouseMode, Scale, ScaleMode, Window, WindowOptions};
use raqote::*;
use table::Table;
use voxel::VoxelGrid;
use std::time::Instant;

pub type Prec = f32;

const N_ELEMS: usize = 100000;
const N_THREADS: u32 = 16;
const GRID_SCALE: Prec = 100.0;
const PRECISION_RADIUS: Prec = 150.0;
const PRECISION_RADIUS_SQUARED: Prec = PRECISION_RADIUS * PRECISION_RADIUS;
const G: Prec = 40.0;
const DT: Prec = 0.01;
const BOUNCE_DIST: Prec = 5.0;
const BOUNCE_COEFF: Prec = 0.2;

const STEPS_PER_FRAME: usize = 10;

const RADIUS: Prec = 1250.0;
const RANDOM_SPEED: Prec = 5.0;
const BLACK_HOLE_MASS: Prec = 500000.0;
const Z_SCALE: Prec = 0.3;

// const NORM_MULT: Prec = 0.001;
const WIDTH: usize = 1920;
const HEIGHT: usize = 1080;
const SCALE: Prec = 1000.0;
const MIN_X: Prec = -SCALE * WIDTH as Prec / HEIGHT as Prec;
const MAX_X: Prec = SCALE * WIDTH as Prec / HEIGHT as Prec;
const MIN_Y: Prec = -SCALE;
const MAX_Y: Prec = SCALE;

fn main() {
    let mut position = gen::uniform(RADIUS, N_ELEMS);
    position.foreach(|_, p, _| (p.0, p.1, p.2 * Z_SCALE));
    let mut mass = gen::scalar_uniform(0.5, 1.5, N_ELEMS);
    let mut speed = gen::uniform(RANDOM_SPEED, N_ELEMS);
    speed.foreach(|i, s, _| {
        let p = position.get(i).unwrap();
        let n = Prec::sqrt(p.0 * p.0 + p.1 * p.1);
        let alpha = Prec::atan2(p.1 / n, p.0 / n);
        let distance = Prec::sqrt(distance_squared(*p, (0.0, 0.0, 0.0)));
        let distance_sqrt = Prec::sqrt(distance);
        (
            s.0 - Prec::sqrt(BLACK_HOLE_MASS * G + distance / RADIUS * N_ELEMS as Prec) * Prec::sin(alpha) / distance_sqrt,
            s.1 + Prec::sqrt(BLACK_HOLE_MASS * G + distance / RADIUS * N_ELEMS as Prec) * Prec::cos(alpha) / distance_sqrt,
            s.2,
        )
    });
    let mut accel_close = gen::null(N_ELEMS);
    let mut accel_far = gen::null(N_ELEMS);

    position.array.push((0.0, 0.0, 0.0));
    mass.array.push(BLACK_HOLE_MASS);
    speed.array.push((0.0, 0.0, 0.0));
    accel_close.array.push((0.0, 0.0, 0.0));
    accel_far.array.push((0.0, 0.0, 0.0));

    let mut voxels = VoxelGrid::new(GRID_SCALE);
    voxels.populate(&position);
    voxels.update(&position, &speed, &mass);

    let mut window = Window::new(
        "Galaxy goes brr",
        WIDTH,
        HEIGHT,
        WindowOptions {
            ..WindowOptions::default()
        },
    )
    .unwrap();

    let mut dt = DrawTarget::new(WIDTH as i32, HEIGHT as i32);

    {
        let mut pb = PathBuilder::new();

        pb.rect(0.0, 0.0, WIDTH as f32, HEIGHT as f32);
        dt.fill(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(
                0xff, 0x00, 0x00, 0x00,
            )),
            &DrawOptions::new(),
        );
    }

    let mut step = 0;

    loop {
        let start = Instant::now();
        step += 1;
        // Advance N step

        accel_far.foreach_threaded(N_THREADS, |i, _, _| {
            let my_position = *position.get(i).unwrap();
            let mut res = (0.0, 0.0, 0.0);

            for voxel in voxels.iter() {
                let voxel_distance_squared = distance_squared(my_position, voxel.position);
                if voxel_distance_squared > PRECISION_RADIUS_SQUARED {
                    let distance_squared = distance_squared(my_position, voxel.center_of_mass);
                    let norm = G * voxel.total_mass / distance_squared; // F/m
                    res = add(res, normalize_difference(my_position, voxel.position, norm));
                }
            }

            res
        });

        for _ in 0..STEPS_PER_FRAME {
            accel_close.foreach_threaded(N_THREADS, |i, _, _| {
                let my_position = *position.get(i).unwrap();
                let my_coords = voxel::get_coords(GRID_SCALE, my_position);
                let mut res = (0.0, 0.0, 0.0);

                // Old code:
                // for (index, &position) in position.iter().enumerate() {
                //     if index == i {
                //         continue
                //     }
                //     let mut distance_squared = distance_squared(my_position, position);
                //     if distance_squared < BOUNCE_DIST {
                //         distance_squared *= -BOUNCE_COEFF;
                //     }
                //     let their_mass = mass.get(index).unwrap();
                //     let norm = G * their_mass / distance_squared; // F/m
                //     res = add(res, normalize_difference(my_position, position, norm))
                // }

                for voxel in voxels.iter() {
                    if my_coords == voxel.coords || distance_squared(my_position, voxel.position) <= PRECISION_RADIUS_SQUARED {
                        for &index in voxel.elements.iter() {
                            let their_position = *position.get(index).unwrap();
                            let mut distance_squared = distance_squared(my_position, their_position);
                            if distance_squared < BOUNCE_DIST {
                                distance_squared *= -BOUNCE_COEFF;
                            }
                            let their_mass = mass.get(index).unwrap();
                            let norm = G * their_mass / distance_squared; // F/m
                            res = add(res, normalize_difference(my_position, their_position, norm));
                        }
                    }
                }

                res
            });

            speed.foreach(|i, prev, _| add(*prev, mul(add(*accel_close.get(i).unwrap(), *accel_far.get(i).unwrap()), DT)));
            position.foreach(|i, prev, _| add(*prev, mul(*speed.get(i).unwrap(), DT)));
            // voxels.update(&position, &speed, &mass);
        }

        voxels.relocate(&position);
        voxels.gc();
        voxels.update(&position, &speed, &mass);
        // println!("{:?}", start.elapsed());

        // Start drawing
        let mut pb = PathBuilder::new();

        pb.rect(0.0, 0.0, WIDTH as f32, HEIGHT as f32);
        dt.fill(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(0x20, 0x00, 0x00, 0x00)),
            &DrawOptions::new()
        );

        let mut pb = PathBuilder::new();

        let mut n = 0;
        for p in position.iter() {
            // let light = sigma(norm(*speed) * NORM_MULT);
            // println!("{:?}", (light * 255.0) as u8);
            let c = position.array.last().unwrap();
            let delta_x = p.0 - c.0;
            let delta_y = p.1 - c.1;
            let delta_z = p.2 - c.2;

            let x = (delta_x - MIN_X) / (MAX_X - MIN_X) * WIDTH as Prec;
            let y = ((delta_y * 0.333 + delta_z * 0.666) - MIN_Y) / (MAX_Y - MIN_Y) * HEIGHT as Prec;

            if x >= 0.0 && x <= WIDTH as Prec && y >= 0.0 && y <= HEIGHT as Prec {
                pb.rect(Prec::round(x) as f32, Prec::round(y) as f32, 1.0, 1.0);
                n += 1;
                if n % 100 == 0 {
                    dt.fill(
                        &pb.finish(),
                        &Source::Solid(SolidSource::from_unpremultiplied_argb(0xd0, 0xff, 0xff, 0xff)),
                        &DrawOptions::new()
                    );
                    pb = PathBuilder::new();
                }
            }
        }

        dt.fill(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(0xd0, 0xff, 0xff, 0xff)),
            &DrawOptions::new()
        );

        // Update image
        window.update_with_buffer(dt.get_data(), WIDTH, HEIGHT).unwrap();
        dt.write_png(format!("./out/{}.png", step)).unwrap();
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
