use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
struct Args {
    path: Option<PathBuf>,
}

const IMAGE_WIDTH: u64 = 256;
const IMAGE_HEIGHT: u64 = 256;

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

    write!(outfile, "P3\n{IMAGE_WIDTH} {IMAGE_HEIGHT}\n255\n")?;

    for j in 0..IMAGE_HEIGHT {
        for i in 0..IMAGE_WIDTH {
            let r = i as f64 / (IMAGE_WIDTH - 1) as f64;
            let g = j as f64 / (IMAGE_HEIGHT - 1) as f64;
            let b = 0.0;

            let ir = (255.999 * r) as u64;
            let ig = (255.999 * g) as u64;
            let ib = (255.999 * b) as u64;

            writeln!(outfile, "{ir} {ig} {ib}")?;
        }
        if let Some(bar) = &bar {
            bar.inc(1);
        }
    }

    outfile.flush()?;

    Ok(())
}
