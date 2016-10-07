#![allow(dead_code, unused_variables, unused)]

extern crate sdl2;
extern crate time;

use sdl2::pixels::PixelFormatEnum;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::render::{Renderer, Texture};

use time::PreciseTime;

use std::ops::{Add, Sub, Mul, MulAssign};

const SCREEN_WIDTH: u32 = 800;
const SCREEN_HEIGHT: u32 = 600;
const UPSAMPLE_FACTOR: u32 = 1;

const FRACTALS: [fn(usize, usize, FractalBounds) -> FractalParams; 2] = [mandelbrot, julia1];
const FRACTAL_BOUNDS: [FractalBounds; 2] = [
        FractalBounds {
            lo: ComplexNumber { r: -2.5, j: -1.33 },
            hi: ComplexNumber {r: 1.0, j: 1.33 } },
        FractalBounds {
            lo: ComplexNumber { r: -1.333333, j: -1.333333 },
            hi: ComplexNumber { r: 1.333333, j: 1.333333 } }];
const START_FRACTAL: usize = 0;

const PI: FloatPrecision = 3.14159265359;

const PALETTE_HEIGHT: u32 = 10;

type FloatPrecision = f64;

#[derive(Clone, Copy, PartialEq, Debug)]
struct ComplexNumber {
    r: FloatPrecision,
    j: FloatPrecision,
}

impl ComplexNumber {
    fn square_module(&self) -> FloatPrecision {
        self.r * self.r + self.j * self.j
    }
}

impl Add for ComplexNumber {
    type Output = ComplexNumber;
    fn add(self, rhs: ComplexNumber) -> Self::Output {
        ComplexNumber {r: self.r + rhs.r, j: self.j + rhs.j }
    }
}

impl Sub for ComplexNumber {
    type Output = ComplexNumber;
    fn sub(self, rhs: ComplexNumber) -> Self::Output {
        ComplexNumber {r: self.r - rhs.r, j: self.j - rhs.j }
    }
}

impl Mul for ComplexNumber {
    type Output = ComplexNumber;
    fn mul(self, rhs: ComplexNumber) -> Self::Output {
        ComplexNumber {r: self.r * rhs.r - self.j * rhs.j,
                       j: self.j * rhs.r + self.r * rhs.j}
    }
}

impl MulAssign for ComplexNumber {
    fn mul_assign(&mut self, rhs: ComplexNumber) {
        let new_r = self.r * rhs.r - self.j * rhs.j;
        let new_j = self.j * rhs.r + self.r * rhs.j;
        self.r = new_r;
        self.j = new_j;
    }
}

#[derive(Copy, Clone, Debug)]
struct FractalBounds {
    lo: ComplexNumber,
    hi: ComplexNumber,
}

#[derive(Copy, Clone, Debug)]
struct FractalParams {
    z0: ComplexNumber,
    c: ComplexNumber,
    p: u32,
    max_value: FloatPrecision,
    max_iter: u64,
}

#[derive(Copy, Clone, Debug)]
struct Color(FloatPrecision, FloatPrecision, FloatPrecision);

fn powi(z: ComplexNumber, p: u32) -> ComplexNumber {
    let mut res = z;
    let mut pp: u32 = p;
    while pp > 1 {
        res *= z;
        pp -= 1;
    }
    res
}

fn f(z: ComplexNumber, c: ComplexNumber, p: u32) -> ComplexNumber {
    powi(z, p) + c
}

#[test]
fn test_powi() {
    let c = ComplexNumber { r: 1.0, j: 3.0 };
    let a = c * c * c;
    let b = powi(c, 3);
    assert_eq!(a, b);
}

fn to_fract_coordinates(screen_x: usize, screen_y: usize, fractal_bounds: FractalBounds)
  -> ComplexNumber {
    let fractal_width = fractal_bounds.hi.r - fractal_bounds.lo.r;
    let fractal_height = fractal_bounds.hi.j - fractal_bounds.lo.j;
    let screen_midpix_x = (screen_x as FloatPrecision) + 0.5;
    let screen_midpix_y = (screen_y as FloatPrecision) + 0.5;
    let r = (screen_midpix_x * fractal_width / (SCREEN_WIDTH as FloatPrecision)) + fractal_bounds.lo.r;
    let j = (screen_midpix_y * fractal_height / (SCREEN_HEIGHT as FloatPrecision)) + fractal_bounds.lo.j;
    ComplexNumber { r: r, j: j }
}

fn escape_time(z0: ComplexNumber, c: ComplexNumber, p: u32, max_value: FloatPrecision, max_iter: u64) -> FloatPrecision {
    let mut n = 0;
    let mut z = z0;
    while z.square_module() < max_value && n < max_iter {
        z = f(z, c, p);

        // z = new_z;
        n += 1;
    }

    z = f(z, c, p);
    z = f(z, c, p);
    let modulus = z.square_module().sqrt();
    let mu: FloatPrecision = (n as FloatPrecision) - modulus.ln().log(2.0);
    return mu;
}

fn index(u: FloatPrecision, max_iter: u64) -> FloatPrecision {
    2.0 * (clamp(u / 100.0) - 0.15)
    // 2.0 * clamp(u / max_iter as FloatPrecision)
}

fn palette(t: FloatPrecision) -> Color {
    let r: FloatPrecision = 0.5 + 0.5 * (2.0 * PI * (t + 0.04)).cos();
    let g: FloatPrecision = 0.5 + 0.5 * (2.0 * PI * (t + 0.10)).cos();
    let b: FloatPrecision = 0.5 + 0.5 * (2.0 * PI * (t + 0.20)).cos();
    return Color(r, g, b);
}

fn mandelbrot(screen_x: usize, screen_y: usize, bounds: FractalBounds) -> FractalParams {
    let z0 = ComplexNumber { r: 0.0, j: 0.0 };
    let c = to_fract_coordinates(screen_x, screen_y, bounds);
    return FractalParams { z0: z0, c: c, p: 2, max_value: 10.0, max_iter: 128 };
}

fn julia1(screen_x: usize, screen_y: usize, bounds: FractalBounds) -> FractalParams {
    let z0 = to_fract_coordinates(screen_x, screen_y, bounds);
    let c = ComplexNumber { r: 0.5, j: 0.25 };
    return FractalParams { z0: z0, c: c, p: 2, max_value: 4.0, max_iter: 256 };
}

fn clamp(x: FloatPrecision) -> FloatPrecision {
    if x <= 0.0 {
        return 0.0;
    } else if x >= 1.0 {
       return 1.0;
    } else {
        return x;
    }
}

fn render_fractal(renderer: &mut Renderer, texture: &mut Texture, current_fractal: usize, current_bounds: FractalBounds) {
    let start = PreciseTime::now();

    texture.with_lock(None, |buffer: &mut [u8], pitch: usize| {
        for y in 0..((SCREEN_HEIGHT - PALETTE_HEIGHT) as usize) {
            for x in 0..(SCREEN_WIDTH as usize) {
                let FractalParams{z0, c, p, max_value, max_iter} =
                    FRACTALS[current_fractal](x, y, current_bounds);

                let u = escape_time(z0, c, p, max_value, max_iter);
                let t = index(u, max_iter);
                let Color(r, g, b) = palette(t);
                let offset = y * pitch + x * 3;
                buffer[offset + 0] = (r * 255.0) as u8;
                buffer[offset + 1] = (g * 255.0) as u8;
                buffer[offset + 2] = (b * 255.0) as u8;
            }
        }

        for y in ((SCREEN_HEIGHT - PALETTE_HEIGHT) as usize)..(SCREEN_HEIGHT as usize) {
            for x in 0..(SCREEN_WIDTH as usize) {
                let t = x as FloatPrecision / SCREEN_WIDTH as FloatPrecision;
                let Color(r, g, b) = palette(t);
                let offset = y * pitch + x * 3;
                buffer[offset + 0] = (r * 255.0) as u8;
                buffer[offset + 1] = (g * 255.0) as u8;
                buffer[offset + 2] = (b * 255.0) as u8;
            }
        }
    }).unwrap();
    renderer.clear();
    renderer.copy(&texture, None, None);
    renderer.present();

    let render_time = start.to(PreciseTime::now()).num_milliseconds();
    println!("Frame time: {}", render_time);
}

pub fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("rust-sdl2 demo: Video", SCREEN_WIDTH * UPSAMPLE_FACTOR, SCREEN_HEIGHT * UPSAMPLE_FACTOR)
    .position_centered()
    .opengl()
    .build()
    .unwrap();

    let mut renderer = window.renderer().build().unwrap();

    let mut texture = renderer.create_texture_streaming(
        PixelFormatEnum::RGB24, SCREEN_WIDTH, SCREEN_HEIGHT).unwrap();

    let t = PreciseTime::now();
    let mut current_fractal: usize = START_FRACTAL;
    let mut current_bounds = FRACTAL_BOUNDS[current_fractal];
    let mut something_changed: bool = false;

    render_fractal(&mut renderer, &mut texture, current_fractal, current_bounds);

    let mut event_pump = sdl_context.event_pump().unwrap();


    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                Event::KeyDown { keycode: Some(Keycode::Down), .. } => {
                    something_changed = true;
                    current_fractal = (current_fractal + 1) % FRACTALS.len();
                    current_bounds = FRACTAL_BOUNDS[current_fractal];
                },
                Event::KeyDown { keycode: Some(Keycode::R), ..}  => {
                    something_changed = true;
                    current_bounds = FRACTAL_BOUNDS[current_fractal];
                }
                Event::MouseButtonDown { x, y, mouse_btn, .. } => {
                    println!("Click {} {}", x, y);

                    let mut fractal_width = current_bounds.hi.r - current_bounds.lo.r;
                    let mut fractal_height = current_bounds.hi.j - current_bounds.lo.j;

                    let x_proportion = x as FloatPrecision / SCREEN_WIDTH as FloatPrecision;
                    let y_proportion = y as FloatPrecision / SCREEN_HEIGHT as FloatPrecision;
                    // println!("{} {} ", x_proportion, y_proportion);
                    let new_center = current_bounds.lo +
                        ComplexNumber { r: x_proportion * fractal_width, j: y_proportion * fractal_height };

                    const ZOOM_FACTOR: FloatPrecision = 10.0;
                    let mut diagonal = ComplexNumber { r: fractal_width / 2.0, j: fractal_height / 2.0 };
                    if mouse_btn == sdl2::mouse::Mouse::Left {
                        diagonal.r /= ZOOM_FACTOR;
                        diagonal.j /= ZOOM_FACTOR;
                    } else if mouse_btn == sdl2::mouse::Mouse::Right {
                        diagonal.r *= ZOOM_FACTOR;
                        diagonal.j *= ZOOM_FACTOR;
                    }
                    current_bounds.lo = new_center - diagonal;
                    current_bounds.hi = new_center + diagonal;
                    something_changed = true;
                    // println!("New center {:?}", diagonal);
                }
                _ => {}
            }
        }
        if something_changed {
            render_fractal(&mut renderer, &mut texture, current_fractal, current_bounds);
            something_changed = false;
        }

    }
}
