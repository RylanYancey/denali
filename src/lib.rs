
#![allow(unused_parens)]

mod gen;
use gen::*;

pub struct SimplexNoise {
    octaves      : u8,
    x_frequency  : f32,
    y_frequency  : f32,
    lacunarity   : f32,
    persistence  : f32,

    max: f32,
    min: f32,

    perm: [u8; 512],
}

impl SimplexNoise {
    pub fn new(
        octaves: u8, x_frequency: f32, y_frequency: f32,
        lacunarity: f32, persistence: f32, max: f32, min: f32, seed: u128
    ) -> Self {
        Self { octaves, x_frequency, y_frequency, lacunarity, 
               persistence, max, min, perm: get_perm(seed) }
    }

    pub fn generate2D (&self, x: f32, y: f32) -> f32 {
        let mut output: f32 = 0.0;
        let mut denom : f32 = 0.0;
    
        let mut xfreq = self.x_frequency;
        let mut yfreq = self.y_frequency;
        let mut amp = 1.0;
    
        for i in 0..self.octaves {
            output += amp * simplex2d(x * xfreq, y * yfreq, &self.perm);
            denom += amp;

            xfreq += self.lacunarity;
            yfreq += self.lacunarity;

            amp *= self.persistence;
        }
    
        return (((output / denom) + 1.0) * (self.max - self.min)) / 2.0 + self.min;
    }

    pub fn generate3D (&self, x: f32, y: f32, z: f32) -> f32 {
        let mut output: f32 = 0.0;
        let mut denom : f32 = 0.0;

        let mut xfreq = self.x_frequency;
        let mut yfreq = self.y_frequency;
        let mut amp = 1.0;

        for i in 0..self.octaves {
            output += amp * simplex2d(x * xfreq, y * yfreq, &self.perm);
            denom += amp;

            xfreq += self.lacunarity;
            yfreq += self.lacunarity;

            amp *= self.persistence;
        }

        return (((output / denom) + 1.0) * (self.max - self.min)) / 2.0 + self.min;
    }
}
