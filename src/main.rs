use clap::Parser;
use std::fmt;
use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::path::PathBuf;

#[derive(Debug, Copy, Clone)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

type Point = Vec3;
type Color = Vec3;

impl Vec3 {
    const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn length_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    fn dot(self, rhs: Vec3) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    fn _cross(self, rhs: Vec3) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    fn unit_vector(self) -> Self {
        self / self.length()
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Vec3) {
        *self = *self + rhs;
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Vec3 {
        self + (-rhs)
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Vec3) {
        *self = *self - rhs;
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Vec3 {
        Vec3 {
            x: rhs * self.x,
            y: rhs * self.y,
            z: rhs * self.z,
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        rhs * self
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Vec3 {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, rhs: f64) {
        *self = *self / rhs;
    }
}

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} {} {}", self.x, self.y, self.z)
    }
}

fn write_color(mut output: impl Write, mut color: Color) -> anyhow::Result<()> {
    color *= 255.999;
    writeln!(
        output,
        "{} {} {}",
        color.x as u64, color.y as u64, color.z as u64,
    )?;
    Ok(())
}

#[derive(Copy, Clone, Debug)]
struct Ray {
    orig: Point,
    dir: Vec3,
}

impl Ray {
    fn at(self, t: f64) -> Point {
        self.orig + t * self.dir
    }
}

#[derive(Copy, Clone, Debug)]
struct HitRecord {
    point: Point,
    t: f64,
    normal: Vec3,
    front_face: bool,
}

impl HitRecord {
    fn new(point: Point, t: f64, r: Ray, outward_normal: Vec3) -> Self {
        debug_assert_eq!(
            outward_normal.length() as f32, // approximate
            1.0,
            "should be a unit vector",
        );
        let front_face = r.dir.dot(outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        Self {
            point,
            t,
            normal,
            front_face,
        }
    }
}

trait Hittable {
    fn hit(&self, r: Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

struct Sphere {
    center: Point,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, r: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = self.center - r.orig;
        let a = r.dir.length_squared();
        let h = r.dir.dot(oc);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = h * h - a * c;
        if discriminant < 0.0 {
            return None;
        }
        let discriminant_sqrt = discriminant.sqrt();
        let mut root = (h - discriminant_sqrt) / a;
        if !(t_min..=t_max).contains(&root) {
            root = (h + discriminant_sqrt) / a;
        }
        if !(t_min..=t_max).contains(&root) {
            return None;
        }
        let t = root;
        let point = r.at(t);
        let outward_normal = (point - self.center) / self.radius;
        Some(HitRecord::new(point, t, r, outward_normal))
    }
}

impl Hittable for Vec<Box<dyn Hittable>> {
    fn hit(&self, r: Ray, t_min: f64, mut t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = None;
        for hittable in self {
            if let Some(hit) = hittable.hit(r, t_min, t_max) {
                assert!(hit.t <= t_max);
                t_max = hit.t;
                closest_so_far = Some(hit);
            }
        }
        closest_so_far
    }
}

fn ray_color(r: Ray, world: &Vec<Box<dyn Hittable>>) -> Color {
    if let Some(hit) = world.hit(r, 0.0, f64::MAX) {
        return 0.5 * (hit.normal + Color::new(1.0, 1.0, 1.0));
    }
    // If we didn't hit anything, paint the blue sky background.
    let a = 0.5 * (r.dir.unit_vector().y + 1.0);
    (1.0 - a) * Color::new(1.0, 1.0, 1.0) + a * Color::new(0.5, 0.7, 1.0)
}

#[derive(Parser)]
struct Args {
    path: Option<PathBuf>,
}

const IMAGE_WIDTH: u64 = 400;
const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_HEIGHT: u64 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u64;
const FOCAL_LENGTH: f64 = 1.0;
const VIEWPORT_HEIGHT: f64 = 2.0;
const VIEWPORT_WIDTH: f64 = VIEWPORT_HEIGHT * IMAGE_WIDTH as f64 / IMAGE_HEIGHT as f64;
const CAMERA_CENTER: Vec3 = Point::new(0.0, 0.0, 0.0);

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let mut outfile: Box<dyn Write>;
    let bar;
    if let Some(path) = args.path {
        outfile = Box::new(BufWriter::new(File::create(path)?));
        bar = None;
    } else {
        outfile = Box::new(std::io::stdout());
        bar = Some(indicatif::ProgressBar::new(IMAGE_HEIGHT));
    };

    let mut world: Vec<Box<dyn Hittable>> = Vec::new();
    world.push(Box::new(Sphere {
        center: Point::new(0.0, 0.0, -1.0),
        radius: 0.5,
    }));
    world.push(Box::new(Sphere {
        center: Point::new(0.0, -100.5, -1.0),
        radius: 100.0,
    }));

    let viewport_u = Vec3::new(VIEWPORT_WIDTH, 0.0, 0.0);
    let viewport_v = Vec3::new(0.0, -VIEWPORT_HEIGHT, 0.0);
    let pixel_delta_u = viewport_u / IMAGE_WIDTH as f64;
    let pixel_delta_v = viewport_v / IMAGE_HEIGHT as f64;
    let viewport_upper_left =
        CAMERA_CENTER - Vec3::new(0.0, 0.0, FOCAL_LENGTH) - viewport_u / 2.0 - viewport_v / 2.0;
    let pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    write!(outfile, "P3\n{IMAGE_WIDTH} {IMAGE_HEIGHT}\n255\n")?;

    for j in 0..IMAGE_HEIGHT {
        for i in 0..IMAGE_WIDTH {
            let pixel_center =
                pixel00_loc + (i as f64 * pixel_delta_u) + (j as f64 * pixel_delta_v);
            let ray_direction = pixel_center - CAMERA_CENTER;
            let r = Ray {
                orig: CAMERA_CENTER,
                dir: ray_direction,
            };
            let color = ray_color(r, &world);
            write_color(&mut outfile, color)?;
        }
        if let Some(bar) = &bar {
            bar.inc(1);
        }
    }

    outfile.flush()?;

    Ok(())
}
