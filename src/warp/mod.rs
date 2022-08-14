
use super::simplex::*;

pub mod gen;
use gen::*;

pub struct DomainWarp {
    simplex1: Simplex,
    simplex2: Simplex,
    simplex3: Simplex,

    warps: [f32; 6],
    weight: f32
}

impl DomainWarp {

    pub fn new (simplex1: Simplex, simplex2: Simplex, simplex3: Simplex, warp_values: [f32; 6], weight: f32) -> Self {
        Self { simplex1, simplex2, simplex3, warps: warp_values, weight }
    }

    pub fn generate2D (&self, x: f32, y: f32) -> f32 {
        domain_warp2d (&self, x, y)
    }

}