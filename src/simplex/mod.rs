
/// Generate module for Simplex Noise
/// Generates raw simplex noise
/// Not intended for public use, but it's here if you need it.
pub mod gen;
use gen::*;

/// Interface for working with Simplex Noise and Fractal Brownian Motion. \
/// Can be used for both 2D and 3D noise values. \
/// # Examples
/// 
/// Include denali and crate a Simplex object. 
/// ```
/// use denali::*;
/// 
/// let noise = Simplex::new(
///     3, // octaves
///     0.01, // x_freq
///     0.01, // y_freq
///     0.01, // z_freq
///     2.5, // lacunarity
///     0.5, // persistence
///     255.0, // max
///     0.0, // min
///     67893402, // Seed
/// );
/// ```
/// Denali can generate single noise values in 2D or 3D:
/// ```rust,ignore
/// let n: f32 = noise.generate2D(5.0, 10.0);
/// 
/// let m: f32 = noise.generate3D(5.0, 10.0, 8.0);
/// ```
/// Or generate noisemaps 
/// ```rust,ignore
/// // 10 x 10 array
/// let mut map: [f32; 100] = [0.0; 100];
/// noise.generate_noisemap2D(5.0, 5.0, 5.0, &mut map, 10);
/// 
/// // 10 x 10 x 10 array
/// let mut map: [f32; 1000] = [0.0; 1000];
/// noise.generate_noisemap3D(5.0, 5.0, 5.0, &mut map, 10, 10);
/// ```
/// ## More Info
/// Simplex implements Send and Sync.\
/// It also derives Clone and Copy.\
/// it also implements PartialEq, which compares the seeds of two SimplexNoise objects.
#[derive(Clone, Copy)]
pub struct Simplex {
    /// The number of waves to combine together.\
    /// As octaves increases, level of detail generally increases.\
    /// Octaves has a profound impact on lacunarity and persistence.\
    /// It is best practices for octaves to stay between 1 and 8.
    pub octaves      : u8,

    /// the starting x frequency.\
    /// As x_freq increases, you zoom in more and more on the x-axis.\
    /// In general, frequency should always be below 0. 
    pub x_frequency  : f32,

    /// the starting y frequency.\
    /// As y_freq increases, you zoom in more and more on the y-axis.\
    /// In general, frequency should always be below 0. 
    pub y_frequency  : f32,

    /// the starting z frequency.\
    /// As z_freq increases, you zoom in more and more on the z-axis.\
    /// In general, frequency should always be below 0. 
    pub z_frequency  : f32,

    /// The rate of change of the frequency.\
    /// As lacunarity increases, the "variance" decreases.\
    /// This means more hills and valleys, but same overall structure.\
    /// Less than zero lacunarities leads to a decrease in detail.
    pub lacunarity   : f32,

    /// The rate of change of amplitude.\
    /// Persistence increases the "wobble", if you will.\
    /// It's hard to explain, just gotta play with it.\
    /// In general, persistence should stay between 0 and 1. 
    pub persistence  : f32,

    /// The max number this generator can output.
    pub max: f32,

    /// The min number this generator can output.
    pub min: f32,

    /// The permutation the noise algorithm will use to \
    /// inform its number generation. 
    perm: [u8; 512],
    seed: u128,

}

impl Simplex {

    pub fn new(
        octaves: u8, x_frequency: f32, y_frequency: f32, z_frequency: f32,
        lacunarity: f32, persistence: f32, max: f32, min: f32, seed: u128
    ) -> Self {
        Self { octaves, x_frequency, y_frequency, z_frequency,
               lacunarity, persistence, max, min, perm: get_perm(seed), seed }
    }

    /// Change the range field of this noise generator. \
    /// Will cause this gen to produce values in a different range. 
    #[inline]
    pub fn set_range(&mut self, max: f32, min: f32) {
        self.max = max;
        self.min = min;
    }

    pub fn change_seed(&mut self, seed: u128) {
        self.seed = seed;
        self.perm = get_perm(seed);
    }

    /// Generates a single noise value. \
    /// `x` and `y` are the input values, and dictate the algorithm on how to behave. \
    /// This function also applies Fractal Brownian Motion.
    pub fn generate2D (&self, x: f32, y: f32) -> f32 {

        // Create temporary values to hold sums
        let mut output: f32 = 0.0;
        let mut denom : f32 = 0.0;
    
        // temp values to hold starting frequencies.
        let mut xfreq = self.x_frequency;
        let mut yfreq = self.y_frequency;

        // amplitude always set to 1
        let mut amp = 1.0;
    
        // octaves sets how many times we run this part
        for _i in 0..self.octaves {
            // add product of amp and the output of simplex3d to get the noise value for this octave. 
            output += amp * simplex2d(x * xfreq, y * yfreq, &self.perm);
            // add to denom so we can calculate range. 
            denom += amp;

            // multiply lacunarity to frequency.
            xfreq *= self.lacunarity;
            yfreq *= self.lacunarity;

            // multiply amp by persistence. 
            amp *= self.persistence;
        }

        // Calculate range and converted to target range.
        (((output / denom) + 1.0) * (self.max - self.min)) / 2.0 + self.min
    }

    /// Generates a single noise value. \
    /// `x` and `y` are the input values, and dictate the algorithm on how to behave. \
    /// This function also applies Fractal Brownian Motion.
    pub fn generate3D (&self, x: f32, y: f32, z: f32) -> f32 {

        // Create temporary values to hold sums
        let mut output: f32 = 0.0;
        let mut denom : f32 = 0.0;

        // temp values to hold starting frequencies.
        let mut xfreq = self.x_frequency;
        let mut yfreq = self.y_frequency;
        let mut zfreq = self.z_frequency;

        // amplitude always set to 1
        let mut amp = 1.0;

        // octaves sets how many times we run this part
        for _i in 0..self.octaves {
            // add product of amp and the output of simplex3d to get the noise value for this octave. 
            output += amp * simplex3d(x * xfreq, y * yfreq, z * zfreq, &self.perm);
            // add to denom so we can calculate range. 
            denom += amp;

            // multiply lacunarity to frequency.
            xfreq *= self.lacunarity;
            yfreq *= self.lacunarity;
            zfreq *= self.lacunarity;

            // multiply amp by persistence. 
            amp *= self.persistence;
        }

        // Calculate range and converted to target range. 
        (((output / denom) + 1.0) * (self.max - self.min)) / 2.0 + self.min
    }

    /// Same as generate2D, but takes the absolute value.\
    /// To make best use of this, set your min to negative your max.
    #[inline]
    pub fn ridged2D (&self, x: f32, y: f32) -> f32 {
        f32::abs(self.generate2D(x, y))
    }

    /// Same as generate3D, but takes the absolute value.\
    /// To make best use of this, set your min to negative your max.
    #[inline]
    pub fn ridged3D (&self, x: f32, y: f32, z: f32) -> f32 {
        f32::abs(self.generate3D(x, y, z))
    }

    /// Generates a noisemap of values.\
    /// * x_start -> the x offset for the x input values
    /// * y_start -> the y offset for the y input values
    /// 
    /// * map -> A 1-dimensional array with 2-dimensions - x and y. 
    /// * map_width -> the x dimension of the array. 
    /// 
    /// Think of x_start and y_start as the position of the map if it was in coordinate space - make them 0 and 0 if the you just want the values.
    /// 
    /// The input values for the noise function will be every number between x_start and map_width, 
    /// and every number between y_start and map_height, which is calculated using `map.len();`.
    pub fn generate_noisemap2D (&self, x_start: f32, y_start: f32, map: &mut [f32], map_width: usize) {
        for x in 0..map_width {
            for y in 0..(map.len() / map_width) {
                map[x + map_width * y] = self.generate2D(x_start + x as f32, y_start + y as f32);
            }
        }
    }

    /// Generates a noisemap of values.\
    /// * x_start -> the x offset for the x input values
    /// * y_start -> the y offset for the y input values
    /// * z_start -> the z offset for the z input values
    /// 
    /// * map -> A 1-dimensional array with 3-dimensions - x, y, and z. 
    /// * map_width -> the x dimension of the array. 
    /// * map_height -> the y dimension of the array. 
    /// 
    /// Think of x_start, y_start, and z_start as the position of the map if it was in coordinate space - make them 0, 0, 0 if the you just want the values.
    /// 
    /// The input values for the noise function will be every number between x_start and map_width, 
    /// and every number between y_start and map_height, and every number between z_start and map_depth, which is calculated using `map.len();`.
    pub fn generate_noisemap3D (&self, x_start: f32, y_start: f32, z_start: f32, map: &mut [f32], map_width: usize, map_height: usize) {
        for x in 0..map_width {
            for y in 0..map_height {
                for z in 0..(map.len() / (map_width * map_height)) {
                    map[x + map_width * y + map_width * map_height * z] = 
                        self.generate3D(x_start + x as f32, y_start + y as f32, z_start + z as f32);
                }
            }
        }
    }

}

impl Default for Simplex {
    fn default() -> Self {
        Simplex::new(
            3, // octaves
            0.01, // x_freq
            0.01, // y_freq
            0.01, // z_freq
            2.5, // lacunarity
            0.5, // persistence
            255.0, // max
            0.0, // min
            67893402, // Seed
        )
    }
}

impl PartialEq for Simplex {
    fn eq(&self, other: &Self) -> bool {
        self.seed == other.seed
    }
}

unsafe impl Send for Simplex { }
unsafe impl Sync for Simplex { }
