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

    const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
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

#[derive(Debug, Clone, Copy)]
struct Interval {
    min: f64,
    max: f64,
}

impl Interval {
    fn new(min: f64, max: f64) -> Self {
        Self { min, max }
    }

    fn _contains(self, value: f64) -> bool {
        self.min <= value && value <= self.max
    }

    fn surrounds(self, value: f64) -> bool {
        self.min < value && value < self.max
    }

    fn clamp(self, value: f64) -> f64 {
        value.clamp(self.min, self.max)
    }
}

#[derive(Copy, Clone, Debug)]
struct HitRecord {
    point: Point,
    t: f64,
    normal: Vec3,
    _front_face: bool,
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
            _front_face: front_face,
        }
    }
}

trait Hittable {
    fn hit(&self, r: Ray, t_range: Interval) -> Option<HitRecord>;
}

struct Sphere {
    center: Point,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, r: Ray, t_range: Interval) -> Option<HitRecord> {
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
        if !t_range.surrounds(root) {
            root = (h + discriminant_sqrt) / a;
        }
        if !t_range.surrounds(root) {
            return None;
        }
        let t = root;
        let point = r.at(t);
        let outward_normal = (point - self.center) / self.radius;
        Some(HitRecord::new(point, t, r, outward_normal))
    }
}

type World = Vec<Box<dyn Hittable>>;

impl Hittable for World {
    fn hit(&self, r: Ray, mut t_range: Interval) -> Option<HitRecord> {
        let mut closest_so_far = None;
        for hittable in self {
            if let Some(hit) = hittable.hit(r, t_range) {
                assert!(hit.t <= t_range.max);
                t_range.max = hit.t;
                closest_so_far = Some(hit);
            }
        }
        closest_so_far
    }
}

struct Camera {
    pub aspect_ratio: f64,
    pub image_width: u64,
}

impl Camera {
    fn image_height(&self) -> u64 {
        (self.image_width as f64 / self.aspect_ratio) as u64
    }

    fn render(
        &self,
        world: &World,
        mut outfile: impl Write,
        progress_bar: Option<&indicatif::ProgressBar>,
    ) -> anyhow::Result<()> {
        let focal_length: f64 = 1.0;
        let viewport_height = 2.0;
        let viewport_width = viewport_height * self.image_width as f64 / self.image_height() as f64;
        let camera_center = Point::new(0.0, 0.0, 0.0);
        let viewport_u = Vec3::new(viewport_width, 0.0, 0.0);
        let viewport_v = Vec3::new(0.0, -viewport_height, 0.0);
        let pixel_delta_u = viewport_u / self.image_width as f64;
        let pixel_delta_v = viewport_v / self.image_height() as f64;
        let viewport_upper_left =
            camera_center - Vec3::new(0.0, 0.0, focal_length) - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        write!(
            outfile,
            "P3\n{} {}\n255\n",
            self.image_width,
            self.image_height(),
        )?;

        for j in 0..self.image_height() {
            for i in 0..self.image_width {
                let samples_per_pixel = 100;
                let mut color = Color::zero();
                for _ in 0..samples_per_pixel {
                    let x_random_offset = rand::random::<f64>() - 0.5;
                    let y_random_offset = rand::random::<f64>() - 0.5;
                    let pixel_center = pixel00_loc
                        + ((i as f64 + x_random_offset) * pixel_delta_u)
                        + ((j as f64 + y_random_offset) * pixel_delta_v);
                    let ray_direction = pixel_center - camera_center;
                    let r = Ray {
                        orig: camera_center,
                        dir: ray_direction,
                    };
                    let max_bounces = 10;
                    color += ray_color(r, &world, max_bounces);
                }
                color /= samples_per_pixel as f64;
                write_color(&mut outfile, color)?;
            }
            if let Some(bar) = &progress_bar {
                bar.inc(1);
            }
        }

        outfile.flush()?;
        Ok(())
    }
}

fn random_unit_vector() -> Vec3 {
    loop {
        let p = Point {
            x: rand::random_range(-1.0..1.0),
            y: rand::random_range(-1.0..1.0),
            z: rand::random_range(-1.0..1.0),
        };
        let lensq = p.length_squared();
        if 0.0 < lensq && lensq <= 1.0 {
            return p / lensq;
        }
    }
}

fn _random_on_hemisphere(normal: Vec3) -> Vec3 {
    let on_unit_sphere = random_unit_vector();
    if on_unit_sphere.dot(normal) > 0.0 {
        on_unit_sphere
    } else {
        -on_unit_sphere
    }
}

fn ray_color(r: Ray, world: &World, remaining_bounces: u64) -> Color {
    if remaining_bounces == 0 {
        return Color::zero();
    }
    // 0.001 here fixes "shadow acne".
    if let Some(hit) = world.hit(r, Interval::new(0.001, f64::MAX)) {
        let orig = hit.point;
        let dir = hit.normal + random_unit_vector(); // Lambertian distribution
        return 0.5 * ray_color(Ray { orig, dir }, world, remaining_bounces - 1);
    }
    // If we didn't hit anything, paint the blue sky background.
    let a = 0.5 * (r.dir.unit_vector().y + 1.0);
    (1.0 - a) * Color::new(1.0, 1.0, 1.0) + a * Color::new(0.5, 0.7, 1.0)
}

fn linear_to_gamma(linear: f64) -> f64 {
    if linear > 0.0 { linear.sqrt() } else { 0.0 }
}

fn write_color(mut output: impl Write, color: Color) -> anyhow::Result<()> {
    let intensity = Interval {
        min: 0.0,
        max: 0.999,
    };
    let red = (256.0 * intensity.clamp(linear_to_gamma(color.x))) as u64;
    let blue = (256.0 * intensity.clamp(linear_to_gamma(color.y))) as u64;
    let green = (256.0 * intensity.clamp(linear_to_gamma(color.z))) as u64;
    writeln!(output, "{red} {blue} {green}")?;
    Ok(())
}

#[derive(Parser)]
struct Args {
    path: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let image_width = 400;
    let aspect_ratio = 16.0 / 9.0;
    let camera = Camera {
        image_width,
        aspect_ratio,
    };

    let mut outfile: Box<dyn Write>;
    let progress_bar;
    if let Some(path) = args.path {
        outfile = Box::new(BufWriter::new(File::create(path)?));
        progress_bar = None;
    } else {
        outfile = Box::new(std::io::stdout());
        progress_bar = Some(indicatif::ProgressBar::new(camera.image_height()));
    };

    let mut world: World = Vec::new();
    world.push(Box::new(Sphere {
        center: Point::new(0.0, 0.0, -1.0),
        radius: 0.5,
    }));
    world.push(Box::new(Sphere {
        center: Point::new(0.0, -100.5, -1.0),
        radius: 100.0,
    }));

    camera.render(&world, &mut outfile, progress_bar.as_ref())?;

    outfile.flush()?;

    Ok(())
}
