import { GPU } from "gpu.js";
import { vec3 } from "gl-matrix";
import { Particles } from "./src/particles.js";
import { MACGrid } from "./src/mac-grid.js";
import { createAdvectParticlesKernel } from "./src/simulation/advect-particles.js";
import { createParticleToGridKernel } from "./src/simulation/transfer-particle-to-grid.js";
import { compileKernels } from "./src/simulation/kernels.js";
import { createGridToParticlesKernel } from "./src/simulation/flip.js";

const particles = new Particles(1, {
  min: vec3.fromValues(0.0, 0.0, 0.0),
  max: vec3.fromValues(1.0, 1.0, 1.0),
});
particles.particleBuffer[0] = 0.25;
particles.particleBuffer[1] = 0.25;
particles.particleBuffer[2] = 0.25;
particles.particleBuffer[3] = 0;
particles.particleBuffer[4] = 0;
particles.particleBuffer[5] = 0;

console.log("old particle buffer:");
console.log(particles.particleBuffer);

// create MAC grid with one cell
const grid = new MACGrid(
  {
    min: vec3.fromValues(0.0, 0.0, 0.0),
    max: vec3.fromValues(1.0, 1.0, 1.0),
  },
  1
);

const new_grid = new MACGrid(
  {
    min: vec3.fromValues(0.0, 0.0, 0.0),
    max: vec3.fromValues(1.0, 1.0, 1.0),
  },
  1
);

new_grid.velocity_x[0][0][0] = 3
new_grid.velocity_x[1][0][0] = 4
new_grid.velocity_y[0][0][0] = 1
new_grid.velocity_y[0][1][0] = 2
new_grid.velocity_z[0][0][0] = 5
new_grid.velocity_z[0][0][1] = 6
// console.log(`created grid of size: (${grid.nx}, ${grid.ny}, ${grid.nz})`);


// compile kernels
const gpu = new GPU();
const kernels = compileKernels(gpu, particles, grid);

const flipKernel = createGridToParticlesKernel(gpu, particles.count(), grid.nx, grid.ny, grid.nz, grid.cellSize);
const updatedFLIPParticleBuffer = flipKernel(grid.velocity_x, grid.velocity_y, grid.velocity_z, new_grid.velocity_x, new_grid.velocity_y, new_grid.velocity_z,particles.particleBuffer);

console.log("\nnew particle buffer:");
console.log(updatedFLIPParticleBuffer);
