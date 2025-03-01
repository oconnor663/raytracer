const IMAGE_WIDTH: u64 = 256;
const IMAGE_HEIGHT: u64 = 256;

fn main() {
    print!("P3\n{IMAGE_WIDTH} {IMAGE_HEIGHT}\n255\n");

    for j in 0..IMAGE_HEIGHT {
        for i in 0..IMAGE_WIDTH {
            let r = i as f64 / (IMAGE_WIDTH - 1) as f64;
            let g = j as f64 / (IMAGE_HEIGHT - 1) as f64;
            let b = 0.0;

            let ir = (255.999 * r) as u64;
            let ig = (255.999 * g) as u64;
            let ib = (255.999 * b) as u64;

            println!("{ir} {ig} {ib}");
        }
    }
}
