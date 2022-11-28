use eframe::epaint::ColorImage;
use image::ImageBuffer;

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

pub fn fmin(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

pub fn save_color_image(color_image: &ColorImage, path: &str) {
    let mut new_img = ImageBuffer::new(color_image.size[0] as u32, color_image.size[1] as u32);
    for x in 0..new_img.width() {
        for y in 0..new_img.height() {
            let color = color_image[(x as usize, y as usize)];
            new_img[(x, y)] = image::Rgba([color.r(), color.g(), color.b(), color.a()]);
        }
    }
    new_img.save(path).unwrap();
}
