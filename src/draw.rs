use super::*;
use minifb::{MouseMode, Scale, ScaleMode, Window, WindowOptions};
use raqote::*;

pub fn draw(rx: std::sync::mpsc::Receiver<(Table<(Prec, Prec, Prec)>, Vec<(Prec, Prec, Prec)>)>) {
    let mut window = if OUTPUT_WINDOW {
        Some(Window::new(
            "Galaxy goes brr",
            WIDTH,
            HEIGHT,
            WindowOptions {
                ..WindowOptions::default()
            },
        ).unwrap())
    } else {
        None
    };

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
        let (position, previous_position) = rx.recv().unwrap();
        step += 1;
        // Start drawing
        let mut pb = PathBuilder::new();

        pb.rect(0.0, 0.0, WIDTH as f32, HEIGHT as f32);
        dt.fill(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(0x50, 0x00, 0x00, 0x00)),
            &DrawOptions::new()
        );

        let mut pb = PathBuilder::new();

        let mut n = 0;
        for (i, p) in position.iter().enumerate() {
            // let light = sigma(norm(*speed) * NORM_MULT);
            // println!("{:?}", (light * 255.0) as u8);
            let c = position.array.last().unwrap();
            let prev_c = previous_position.last().unwrap();
            let prev_p = previous_position[i];
            let delta_x = p.0 - c.0;
            let delta_y = p.1 - c.1;
            let delta_z = p.2 - c.2;

            let prev_delta_x = prev_p.0 - prev_c.0;
            let prev_delta_y = prev_p.1 - prev_c.1;
            let prev_delta_z = prev_p.2 - prev_c.2;

            let x = (delta_x - MIN_X) / (MAX_X - MIN_X) * WIDTH as Prec;
            let y = ((delta_y * 0.333 + delta_z * 0.666) - MIN_Y) / (MAX_Y - MIN_Y) * HEIGHT as Prec;

            let prev_x = (prev_delta_x - MIN_X) / (MAX_X - MIN_X) * WIDTH as Prec;
            let prev_y = ((prev_delta_y * 0.333 + prev_delta_z * 0.666) - MIN_Y) / (MAX_Y - MIN_Y) * HEIGHT as Prec;

            if
                x >= 0.0 && x < WIDTH as Prec && y >= 0.0 && y < HEIGHT as Prec
                && prev_x >= 0.0 && prev_x < WIDTH as Prec && prev_y >= 0.0 && prev_y < HEIGHT as Prec
            {
                pb.move_to(Prec::round(prev_x) as f32, Prec::round(prev_y) as f32);
                pb.line_to(Prec::round(x) as f32, Prec::round(y) as f32);
                n += 1;
                if n % 1000 == 0 {
                    dt.stroke(
                        &pb.finish(),
                        &Source::Solid(SolidSource::from_unpremultiplied_argb(0x40, 0xff, 0xff, 0xff)),
                        &StrokeStyle {
                            width: 1.0,
                            cap: LineCap::Square,
                            join: LineJoin::Bevel,
                            miter_limit: 0.0,
                            dash_array: vec![],
                            dash_offset: 0.0
                        },
                        &DrawOptions::new()
                    );
                    pb = PathBuilder::new();
                }
            }
        }

        dt.stroke(
            &pb.finish(),
            &Source::Solid(SolidSource::from_unpremultiplied_argb(0x40, 0xff, 0xff, 0xff)),
            &StrokeStyle {
                width: 1.0,
                cap: LineCap::Square,
                join: LineJoin::Bevel,
                miter_limit: 0.0,
                dash_array: vec![],
                dash_offset: 0.0
            },
            &DrawOptions::new()
        );

        // Update image
        if let Some(ref mut window) = window {
            window.update_with_buffer(dt.get_data(), WIDTH, HEIGHT).unwrap();
        }
        if OUTPUT_FILE {
            dt.write_png(format!("./out/{}.png", step)).unwrap();
        }
    }
}
