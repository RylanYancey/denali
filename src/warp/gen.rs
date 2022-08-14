
use super::*;

/*
    // https://github.com/rshipley160/learn-cuda/wiki/Thread-and-Block-Scheduling
    int id = (gridDim.x * blockDim.y * blockDim.x * blockIdx.y) + (blockDim.y * blockDim.x * blockIdx.x) + (threadIdx.y * blockDim.x) + threadIdx.x;

    if (id > size) return;

    int x = id / width;
    int y = id % height;

    float weight = w;

    float qx = generate(x, y, &fractal);
    float qy = generate(y + a, x + b, &fractal);

    float rx = generate(x + weight*qx + c, y + weight*qy + d, &fractal);
    float ry = generate(x + weight*qx + e, y + weight*qy + f, &fractal);

    float value1 = generate(x + weight * rx, y + weight * ry, &fractal1);
    float value2 = generate(x + weight * rx, y + weight * ry, &fractal2);
    float value3 = generate(x + weight * rx, y + weight * ry, &fractal3);

    red[id] = value1;
    green[id] = value2;
    blue[id] = value3;
*/

pub fn domain_warp2d (warp: &DomainWarp, x: f32, y: f32) -> f32 {

    let qx = warp.simplex1.generate2D(x, y);
    let qy = warp.simplex1.generate2D(y + warp.warps[0], x + warp.warps[1]);

    let rx = warp.simplex2.generate2D(x + warp.weight * qx + warp.warps[2], y + warp.weight * qy + warp.warps[3]);
    let ry = warp.simplex2.generate2D(x + warp.weight * qx + warp.warps[4], y + warp.weight * qy + warp.warps[4]);

    warp.simplex3.generate2D(x + warp.weight * rx, y + warp.weight * ry)

}