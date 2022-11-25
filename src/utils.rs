use eframe::epaint::ColorImage;

pub fn load_image_from_path(path: &str) -> Result<ColorImage, image::ImageError> {
    let image = image::io::Reader::open(path)?.decode()?;
    let size = [image.width() as _, image.height() as _];
    let image_buffer = image.to_rgba8();
    let pixels = image_buffer.as_flat_samples();
    Ok(ColorImage::from_rgba_unmultiplied(size, pixels.as_slice()))
}

// TODO: move to a separate file
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Side {
    pub fig: usize,
    pub side: usize,
}

fn norm_side(x: usize) -> usize {
    if x >= 4 {
        x - 4
    } else {
        x
    }
}

impl Side {
    pub fn ne(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 1),
        }
    }

    pub fn ne2(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 2),
        }
    }

    pub fn pr(&self) -> Self {
        Self {
            fig: self.fig,
            side: norm_side(self.side + 3),
        }
    }
}

pub fn fmax(x: f64, y: f64) -> f64 {
    if x > y {
        x
    } else {
        y
    }
}
