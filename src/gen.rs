
extern crate nanorand;
use nanorand::{Pcg64, Rng};

const PERMUTATION: [u8; 512] = [
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69,
    142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219,
    203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
    74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230,
    220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76,
    132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173,
    186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206,
    59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163,
    70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162,
    241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204,
    176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141,
    128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194,
    233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234,
    75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174,
    20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83,
    111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
    63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188,
    159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147,
    118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170,
    213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253,
    19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193,
    238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31,
    181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
    222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
];

pub fn get_perm(seed: u128) -> [u8; 512] {
    let mut rng = Pcg64::new_seed(seed);
    let mut perm = PERMUTATION;
    rng.shuffle(&mut perm);
    perm
}

const F2: f32 = 0.366025403;
const G2: f32 = 0.211324865;

/// ---------------------------------------
/// Generate 2d Noise

pub fn simplex2d (mut x: f32, mut y: f32, perm: &[u8; 512]) -> f32 {
    let mut n0: f32 = 0.0;
    let mut n1: f32 = 0.0;
    let mut n2: f32 = 0.0;

    let s = (x + y) * F2;
    let xs = x + s;
    let ys = y + s;
    let i = fast_floor(xs);
    let j = fast_floor(ys);

    let t: f32 = ((i + j) as f32) * G2;
    let x_0 = i as f32 - t;
    let y_0 = j as f32 - t;
    let x_0 = x - x_0;
    let y_0 = y - y_0;

    let mut i1: i32 = 0;
    let mut j1: i32 = 0;
    if x_0 > y_0 {
        i1 = 1;
        j1 = 0;
    } else {
        i1 = 0;
        j1 = 1;
    }

    let x1 = x_0 - i1 as f32 + G2;
    let y1 = y_0 - j1 as f32 + G2;
    let x2 = x_0 - 1.0 + 2.0 * G2;
    let y2 = y_0 - 1.0 + 2.0 * G2;

    let ii = modulo(i, 256);
    let jj = modulo(j, 256);

    let mut t0 = 0.5 - x_0 * x_0 - y_0 * y_0;
    if t0 < 0.0 {
        n0 = 0.0;
    } else {
        t0 *= t0;
        let temp = perm[jj as usize];
        n0 = t0 * t0 * gradient(perm[(ii + temp as i32) as usize], x_0, y_0);
    }

    let mut t1 = 0.5 - x1 * x1 - y1 * y1;
    if t1 < 0.0 {
        n1 = 0.0;
    } else {
        t1 *= t1;
        let temp = perm[(jj + j1) as usize];
        n1 = t1 * t1 * gradient(perm[(ii + i1 + temp as i32) as usize], x1, y1);
    }

    let mut t2 = 0.5 - x2 * x2 - y2 * y2;
    if t2 < 0.0 {
        n2 = 0.0;
    } else {
        t2 *= t2;
        let temp = perm[(jj + 1) as usize];
        n2 = t2 * t2 * gradient(perm[(ii + 1 + temp as i32) as usize], x2, y2);
    }

    return 40.0 * (n0 + n1 + n2);
}   

/// Quickly finds the floor of a number faster than std can.
fn fast_floor(x: f32) -> i32 {
    if x > 0.0 {
        return x as i32;
    } else {
        return (x as i32) - 1;
    }
}

/// This function is private and is not intended to be used by an end-user.
/// Function for simplex noise algorithm.
/// Calculates Modulo
fn modulo(x: i32, m: i32) -> i32 {
    let a = x % m;
    if 0 > a {
        return a + m;
    } else {
        return a;
    }
}

/// This function is private and is not intended to be used by an end-user.
/// Function for simplex noise algorithm.
/// Calculates gradients.
fn gradient(hash: u8, x: f32, y: f32) -> f32 {
    let h = hash & 7;

    let mut u: f32 = if 4 > h { x } else { y };
    let v: f32 = if 4 > h { y } else { x };

    if h & 1 != 0 {
        u *= -1.0;
    }

    return u + (if h & 2 != 0 { -2.0 * v } else { 2.0 * v });
}

// - 3D NOISE - //

// Simple skewing factors for the 3D case
const F3: f32 = 0.333333333;
const G3: f32 = 0.166666667;

pub fn generate3D (x: f32, y: f32, z: f32, perm: &[u8; 512]) -> f32 {

    let s = (x + y + z) * F3;

    let i = fast_floor(x + s);
    let j = fast_floor(y + s);
    let k = fast_floor(z + s);

    let t = (i + j + k) as f32 * G3;
    let x0 = x - (i as f32 - t);
    let y0 = y - (j as f32 - t);
    let z0 = z - (k as f32 - t);

    let i1: i32;
    let j1: i32;
    let k1: i32;

    let i2: i32;
    let j2: i32;
    let k2: i32;

    if (x0 >= y0)
    {
        if (y0 >= z0)
        { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0; } // X Y Z order
        else if (x0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; } // X Z Y order
        else { i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; } // Z X Y order
    }
    else
    { // x0<y0
        if (y0 < z0) { i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; } // Z Y X order
        else if (x0 < z0) { i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; } // Y Z X order
        else { i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; } // Y X Z order
    }

    let x1 = x0 - i1 as f32 + G3; // Offsets for second corner in (x,y,z) coords
    let y1 = y0 - j1 as f32 + G3;
    let z1 = z0 - k1 as f32 + G3;
    let x2 = x0 - i2 as f32 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
    let y2 = y0 - j2 as f32 + 2.0 * G3;
    let z2 = z0 - k2 as f32 + 2.0 * G3;
    let x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
    let y3 = y0 - 1.0 + 3.0 * G3;
    let z3 = z0 - 1.0 + 3.0 * G3;

    let ii = modulo(i, 256);
    let jj = modulo(j, 256);
    let kk = modulo(k, 256);

    let mut n: f32 = 0.0;

    let mut t = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
    if (t >= 0.0) {
        t *= t;
        n += t * t * gradient_3d(perm[(ii + perm[(jj + perm[kk as usize] as i32) as usize] as i32) as usize] as i32, x0, y0, z0);
    }

    t = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
    if (t >= 0.0)
    {
        t *= t;
        n += t * t * gradient_3d(perm[(ii + i1 + perm[(jj + j1 + perm[(kk + k1) as usize] as i32) as usize] as i32) as usize] as i32, x1, y1, z1);
    }

    t = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
    if (t >= 0.0)
    {
        t *= t;
        n += t * t * gradient_3d(perm[(ii + i2 + perm[(jj + j2 + perm[(kk + k2) as usize] as i32) as usize] as i32) as usize] as i32, x2, y2, z2);
    }

    t = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
    if (t >= 0.0) 
    {
        t *= t;
        n += t * t * gradient_3d(perm[(ii + 1 + perm[(jj + 1 + perm[(kk + 1) as usize] as i32) as usize] as i32) as usize] as i32, x3, y3, z3);
    }

    32.0 * n

}

fn gradient_3d(hash: i32, x: f32, y: f32, z: f32) -> f32 {
    let h = hash & 15;
    let u = if (h < 8) { x } else { y };
    let v = if (h < 4) { y } else { if (h == 12 || h == 14) { x } else { z } };
    (if (h & 1 != 0) { -u } else { u }) + (if (h & 2 != 0) { -v } else { v })
}
