
#![allow(unused_parens)]
#![allow(non_snake_case)]


// just useful stuff
    // p.x = index / width;
    // p.y = index % width;

    // p.i = x + width / y

    // i = x + width*y + width*height*z;

    // x = i % width;
    // y = (i / width)%height;
    // z = i / (width*height);

pub mod simplex;
pub use simplex::*;

pub mod warp;
pub use warp::*;