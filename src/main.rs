use clap::Parser;
use rand::Rng;
use rand::SeedableRng;
use rayon::prelude::*;
use riddance::Registry;
use std::fmt;
use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

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

    fn cross(self, rhs: Vec3) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    fn unit_vector(self) -> Self {
        self / self.length()
    }

    fn near_zero(self) -> bool {
        let s = 1e-8;
        self.x.abs() < s && self.y.abs() < s && self.z.abs() < s
    }

    fn reflect(self, normal: Self) -> Self {
        self - 2.0 * self.dot(normal) * normal
    }

    fn refract(self, normal: Self, etai_over_etat: f64) -> Self {
        let cos_theta = self.dot(-normal).min(1.0);
        let r_out_perp = etai_over_etat * (self + cos_theta * normal);
        let r_out_parallel = -(1.0 - r_out_perp.length_squared()).abs().sqrt() * normal;
        return r_out_perp + r_out_parallel;
    }

    fn elementwise_mul(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
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
    from_outside: bool,
    material: Material,
    attenuation: Color,
    glow: Color,
}

impl HitRecord {
    fn new(
        point: Point,
        t: f64,
        r: Ray,
        outward_normal: Vec3,
        material: Material,
        attenuation: Color,
        glow: Color,
    ) -> Self {
        debug_assert_eq!(
            outward_normal.length() as f32, // approximate
            1.0,
            "should be a unit vector",
        );
        let from_outside = r.dir.dot(outward_normal) < 0.0;
        let normal = if from_outside {
            outward_normal
        } else {
            -outward_normal
        };
        Self {
            point,
            t,
            normal,
            from_outside,
            material,
            attenuation,
            glow,
        }
    }
}

struct Sphere {
    center: Point,
    radius: f64,
    attenuation: Color,
    material: Material,
    glow: Color,
}

impl Sphere {
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
        Some(HitRecord::new(
            point,
            t,
            r,
            outward_normal,
            self.material,
            self.attenuation,
            self.glow,
        ))
    }
}

enum Hittable {
    Sphere(Sphere),
}

impl Hittable {
    fn hit(&self, r: Ray, t_range: Interval) -> Option<HitRecord> {
        match self {
            Hittable::Sphere(sphere) => sphere.hit(r, t_range),
        }
    }
}

#[derive(Copy, Clone, Debug)]
enum Material {
    Lambertian,
    Metal { fuzz: f64 },
    Dielectric { refraction_index: f64 },
}

impl Material {
    fn scatter(self, ray: Ray, hit: HitRecord) -> Option<Ray> {
        let orig = hit.point;
        match self {
            Material::Lambertian => {
                let mut dir = hit.normal + random_unit_vector();
                if dir.near_zero() {
                    dir = hit.normal;
                }
                Some(Ray { orig, dir })
            }
            Material::Metal { fuzz } => {
                let mut dir = ray.dir.reflect(hit.normal);
                dir = dir.unit_vector() + (fuzz * random_unit_vector());
                if dir.dot(hit.normal) > 0.0 {
                    Some(Ray { orig, dir })
                } else {
                    // Fuzziness pointed the reflection inwards.
                    None
                }
            }
            Material::Dielectric { refraction_index } => {
                let ri = if hit.from_outside {
                    1.0 / refraction_index
                } else {
                    refraction_index
                };
                let unit_direction = ray.dir.unit_vector();
                let cos_theta = unit_direction.dot(-hit.normal).min(1.0);
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
                let cannot_refract = ri * sin_theta > 1.0;
                let orig = hit.point;
                let dir = if cannot_refract || reflectance(cos_theta, ri) > rand::random() {
                    unit_direction.reflect(hit.normal)
                } else {
                    unit_direction.refract(hit.normal, ri)
                };
                Some(Ray { orig, dir })
            }
        }
    }
}

fn reflectance(cos: f64, ri: f64) -> f64 {
    // Use Schlick's approximation for reflectance.
    let mut r0 = (1.0 - ri) / (1.0 + ri);
    r0 *= r0;
    r0 + (1.0 - r0) * (1.0 - cos).powi(5)
}

struct World {
    hittables: Registry<Hittable>,
}

impl World {
    fn new() -> Self {
        Self {
            hittables: Registry::new(),
        }
    }
}

impl World {
    fn hit(&self, r: Ray, mut t_range: Interval) -> Option<HitRecord> {
        let mut closest_so_far = None;
        for hittable in self.hittables.values() {
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
    pub vfov: f64,
    pub lookfrom: Point,
    pub lookat: Point,
    pub vup: Vec3,
    pub defocus_angle: f64,
    pub focus_distance: f64,
    pub samples_per_pixel: u64,
    pub max_bounces: u64,
}

impl Camera {
    fn image_height(&self) -> u64 {
        (self.image_width as f64 / self.aspect_ratio) as u64
    }

    fn render(
        &self,
        world: &World,
        mut outfile: impl Write,
        progress_counter: &AtomicU64,
        progress_bar: Option<&indicatif::ProgressBar>,
    ) -> anyhow::Result<()> {
        let camera_center = self.lookfrom;
        let theta = self.vfov.to_radians();
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h * self.focus_distance;
        let viewport_width = viewport_height * self.image_width as f64 / self.image_height() as f64;

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        let w = (self.lookfrom - self.lookat).unit_vector();
        let u = self.vup.cross(w).unit_vector();
        let v = w.cross(u);

        let viewport_u = viewport_width * u;
        let viewport_v = viewport_height * -v;
        let pixel_delta_u = viewport_u / self.image_width as f64;
        let pixel_delta_v = viewport_v / self.image_height() as f64;
        let viewport_upper_left =
            camera_center - self.focus_distance * w - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        let defocus_radius = self.focus_distance * (self.defocus_angle / 2.0).to_radians().tan();
        let defocus_disk_u = u * defocus_radius;
        let defocus_disk_v = v * defocus_radius;

        write!(
            outfile,
            "P3\n{} {}\n255\n",
            self.image_width,
            self.image_height(),
        )?;

        let mut indexes = Vec::new();
        for j in 0..self.image_height() {
            for i in 0..self.image_width {
                indexes.push((i, j));
            }
        }
        let pixels: Vec<Color> = indexes
            .into_par_iter()
            .map(|(i, j)| {
                let mut color = Color::zero();
                for _ in 0..self.samples_per_pixel {
                    // Construct a camera ray originating from the defocus disk and directed at a
                    // randomly sampled point around the pixel location i, j.
                    let x_random_offset = rand::random::<f64>() - 0.5;
                    let y_random_offset = rand::random::<f64>() - 0.5;
                    let pixel_center = pixel00_loc
                        + ((i as f64 + x_random_offset) * pixel_delta_u)
                        + ((j as f64 + y_random_offset) * pixel_delta_v);
                    let ray_origin = if self.defocus_angle <= 0.0 {
                        camera_center
                    } else {
                        let p = random_in_unit_disk();
                        camera_center + p.x * defocus_disk_u + p.y * defocus_disk_v
                    };
                    let ray_direction = pixel_center - ray_origin;
                    let r = Ray {
                        orig: ray_origin,
                        dir: ray_direction,
                    };
                    color += ray_color(r, &world, self.max_bounces);
                }
                color /= self.samples_per_pixel as f64;
                let progress = progress_counter.fetch_add(1, Relaxed);
                if progress % PROGRESS_FRACTION == 0 {
                    if let Some(bar) = &progress_bar {
                        bar.inc(1);
                    }
                }
                color
            })
            .collect();

        for pixel in pixels {
            write_color(&mut outfile, pixel)?;
        }
        outfile.flush()?;
        Ok(())
    }
}

fn random_unit_vector() -> Vec3 {
    loop {
        let p = Vec3 {
            x: rand::random_range(-1.0..1.0),
            y: rand::random_range(-1.0..1.0),
            z: rand::random_range(-1.0..1.0),
        };
        let lensq = p.length_squared();
        if 1e-160 < lensq && lensq <= 1.0 {
            return p / lensq.sqrt();
        }
    }
}

fn random_in_unit_disk() -> Vec3 {
    loop {
        let p = Vec3 {
            x: rand::random_range(-1.0..1.0),
            y: rand::random_range(-1.0..1.0),
            z: 0.0,
        };
        let lensq = p.length_squared();
        if 1e-160 < lensq && lensq <= 1.0 {
            return p / lensq.sqrt();
        }
    }
}

fn ray_color(r: Ray, world: &World, remaining_bounces: u64) -> Color {
    if remaining_bounces == 0 {
        return Color::zero();
    }
    // 0.001 here fixes "shadow acne".
    if let Some(hit) = world.hit(r, Interval::new(0.001, f64::MAX)) {
        let mut color = if let Some(scatter) = hit.material.scatter(r, hit) {
            ray_color(scatter, world, remaining_bounces - 1)
        } else {
            Color::zero()
        };
        color = hit.attenuation.elementwise_mul(color);
        if hit.from_outside {
            color.x = hit.glow.x + color.x * (1.0 - hit.glow.x);
            color.y = hit.glow.y + color.y * (1.0 - hit.glow.y);
            color.z = hit.glow.z + color.z * (1.0 - hit.glow.z);
        }
        return color;
    }
    // If we didn't hit anything, the sky is black.
    Color::new(0.0, 0.0, 0.0)
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

// Number of pixels rendered in between updates to the progress bar.
const PROGRESS_FRACTION: u64 = 1024;

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let camera = Camera {
        image_width: 2000,
        aspect_ratio: 16.0 / 9.0,
        vfov: 20.0,
        lookfrom: Point::new(13.0, 2.0, -3.0),
        lookat: Point::zero(),
        vup: Vec3::new(0.0, 1.0, 0.0),
        defocus_angle: 0.0,
        focus_distance: 10.0,
        samples_per_pixel: 500,
        max_bounces: 50,
    };

    let mut outfile: Box<dyn Write>;
    let progress_counter = AtomicU64::new(0);
    let progress_bar;
    if let Some(path) = args.path {
        outfile = Box::new(BufWriter::new(File::create(path)?));
        progress_bar = None;
    } else {
        outfile = Box::new(std::io::stdout());
        let bar = indicatif::ProgressBar::new(
            camera.image_width * camera.image_height() / PROGRESS_FRACTION,
        );
        bar.set_style(
            indicatif::style::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {wide_bar}{percent:>3}%")?,
        );
        progress_bar = Some(bar);
    };

    let mut world: World = World::new();

    // the ground
    _ = world.hittables.insert(Hittable::Sphere(Sphere {
        center: Point::new(0.0, -1000.0, 0.0),
        radius: 1000.0,
        attenuation: Color::new(0.5, 0.5, 0.5),
        material: Material::Lambertian,
        glow: Color::zero(),
    }));

    // lots of little spheres
    let mut sphere_rng = rand_chacha::ChaCha8Rng::seed_from_u64(6);
    for a in -11..11 {
        for b in -11..11 {
            let center = Point::new(
                a as f64 + 0.9 * sphere_rng.random::<f64>(),
                0.2,
                b as f64 + 0.9 * sphere_rng.random::<f64>(),
            );
            if (center - Point::new(4.0, 0.2, 0.0)).length() > 0.9 {
                match sphere_rng.random() {
                    0.0..0.8 => {
                        // diffuse
                        _ = world.hittables.insert(Hittable::Sphere(Sphere {
                            center,
                            radius: 0.2,
                            attenuation: Color {
                                x: sphere_rng.random(),
                                y: sphere_rng.random(),
                                z: sphere_rng.random(),
                            },
                            material: Material::Lambertian,
                            glow: Color::zero(),
                        }));
                    }
                    0.8..0.95 => {
                        // metal
                        _ = world.hittables.insert(Hittable::Sphere(Sphere {
                            center,
                            radius: 0.2,
                            attenuation: Color {
                                x: sphere_rng.random_range(0.5..1.0),
                                y: sphere_rng.random_range(0.5..1.0),
                                z: sphere_rng.random_range(0.5..1.0),
                            },
                            material: Material::Metal {
                                fuzz: sphere_rng.random_range(0.0..0.5),
                            },
                            glow: Color::zero(),
                        }));
                    }
                    _ => {
                        // glass
                        let glow_red = sphere_rng.random();
                        let glow_green = sphere_rng.random_range(0.0..glow_red);
                        let glow_blue = sphere_rng.random_range(0.0..glow_green);
                        _ = world.hittables.insert(Hittable::Sphere(Sphere {
                            center,
                            radius: 0.2,
                            attenuation: Color::new(1.0, 1.0, 1.0),
                            material: Material::Dielectric {
                                refraction_index: 1.5,
                            },
                            glow: Color::new(glow_red, glow_green, glow_blue),
                        }));
                    }
                }
            }
        }
    }

    // four big spheres
    _ = world.hittables.insert(Hittable::Sphere(Sphere {
        center: Point::new(0.0, 1.0, 0.0),
        radius: 1.0,
        attenuation: Color::new(1.0, 1.0, 1.0),
        material: Material::Dielectric {
            refraction_index: 1.5,
        },
        glow: Color::zero(),
    }));
    _ = world.hittables.insert(Hittable::Sphere(Sphere {
        center: Point::new(-4.0, 1.0, 0.0),
        radius: 1.0,
        attenuation: Color::new(0.7, 0.3, 0.1),
        material: Material::Lambertian,
        // glow: Color::zero(),
        glow: Color::new(0.7, 0.3, 0.1),
    }));
    _ = world.hittables.insert(Hittable::Sphere(Sphere {
        center: Point::new(4.0, 1.0, 0.0),
        radius: 1.0,
        attenuation: Color::new(0.7, 0.6, 0.5),
        material: Material::Metal { fuzz: 0.0 },
        glow: Color::zero(),
    }));
    _ = world.hittables.insert(Hittable::Sphere(Sphere {
        center: Point::new(7.0, 4.0, -6.0),
        radius: 4.0,
        attenuation: Color::new(0.7, 0.6, 0.5),
        material: Material::Metal { fuzz: 0.0 },
        glow: Color::new(0.5, 0.5, 0.5),
    }));

    camera.render(
        &world,
        &mut outfile,
        &progress_counter,
        progress_bar.as_ref(),
    )?;

    outfile.flush()?;

    Ok(())
}
