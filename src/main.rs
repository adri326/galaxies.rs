extern crate scoped_threadpool;
extern crate rand;
extern crate raqote;
extern crate minifb;

pub mod table;
pub mod gen;

use table::Table;
use raqote::*;
use minifb::{MouseMode, Window, WindowOptions, ScaleMode, Scale};

pub type Prec = f32;

const N_ELEMS: usize = 25000;
const N_THREADS: u32 = 16;
const G: Prec = 20.0;
const DT: Prec = 0.02;
const BOUNCE_DIST: Prec = 10.0;
const BOUNCE_COEFF: Prec = 0.025;

const STEPS_PER_FRAME: usize = 10;

const RADIUS: Prec = 1000.0;
const RANDOM_SPEED: Prec = 5.0;
const TANGENT_SPEED: Prec = BLACK_HOLE_MASS * G;
const BLACK_HOLE_MASS: Prec = 100000.0;
const Z_SCALE: Prec = 0.3;

const NORM_MULT: Prec = 0.001;
const WIDTH: usize = 1920;
const HEIGHT: usize = 1080;
const MIN_X: Prec = -1000.0 * WIDTH as Prec / HEIGHT as Prec;
const MAX_X: Prec = 1000.0 * WIDTH as Prec / HEIGHT as Prec;
const MIN_Y: Prec = -1000.0;
const MAX_Y: Prec = 1000.0;

fn main() {
    let mut position = gen::uniform(RADIUS, N_ELEMS);
    position.foreach(|_, p, _| (p.0, p.1, p.2 * Z_SCALE));
    let mut mass = gen::scalar_uniform(0.5, 1.0, N_ELEMS);
    let mut speed = gen::uniform(RANDOM_SPEED, N_ELEMS);
    speed.foreach(|i, s, _| {
        let p = position.get(i).unwrap();
        let n = Prec::sqrt(p.0 * p.0 + p.1 * p.1);
        let alpha = Prec::atan2(p.1 / n, p.0 / n);
        let distance_sqrt = Prec::sqrt(Prec::sqrt(distance_squared(*p, (0.0, 0.0, 0.0))));
        (
            s.0 - Prec::sqrt(TANGENT_SPEED) * Prec::sin(alpha) / distance_sqrt,
            s.1 + Prec::sqrt(TANGENT_SPEED) * Prec::cos(alpha) / distance_sqrt,
            s.2
        )
    });
    let mut accel = gen::null(N_ELEMS);

    position.array.push((0.0, 0.0, 0.0));
    mass.array.push(BLACK_HOLE_MASS);
    speed.array.push((0.0, 0.0, 0.0));
    accel.array.push((0.0, 0.0, 0.0));

    let mut window = Window::new("Galaxy goes brr", WIDTH, HEIGHT, WindowOptions {
        ..WindowOptions::default()
    }).unwrap();


    let mut dt = DrawTarget::new(WIDTH as i32, HEIGHT as i32);

    {
        let mut pb = PathBuilder::new();

        pb.rect(0.0, 0.0, WIDTH as f32, HEIGHT as f32);
        dt.fill(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(0xff, 0x00, 0x00, 0x00)),
            &DrawOptions::new()
        );
    }

    let mut step = 0;

    loop {
        step += 1;
        // Advance N step
        for _ in 0..STEPS_PER_FRAME {
            accel.foreach_threaded(N_THREADS, |i, _, _| {
                let my_position = position.get(i).unwrap();
                let mut res = (0.0, 0.0, 0.0);
                for (index, position) in position.iter().enumerate() {
                    if index == i {
                        continue
                    }
                    let mut distance_squared = distance_squared(*my_position, *position);
                    if distance_squared < BOUNCE_DIST {
                        distance_squared *= -BOUNCE_COEFF;
                    }
                    let their_mass = mass.get(index).unwrap();
                    let norm = G * their_mass / distance_squared; // F/m
                    res = add(res, normalize_difference(*my_position, *position, norm))
                }
                res
            });
            speed.foreach(|i, prev, _| add(*prev, mul(*accel.get(i).unwrap(), DT)));
            position.foreach(|i, prev, _| add(*prev, mul(*speed.get(i).unwrap(), DT)));
        }

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
        dt.write_png(format!("./out/{}.png", step));
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
fn normalize_difference(a: (Prec, Prec, Prec), b: (Prec, Prec, Prec), unit: Prec) -> (Prec, Prec, Prec) {
    if unit == 0.0 {
        return (0.0, 0.0, 0.0)
    }

    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let dz = b.2 - a.2;
    let d2 = dx * dx + dy * dy + dz * dz;

    if d2 == 0.0 {
        return (0.0, 0.0, 0.0)
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
