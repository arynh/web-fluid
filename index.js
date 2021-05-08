import { GPU } from "gpu.js";
import { vec3 } from "gl-matrix";
import { Particles } from "./src/particles.js";
import { MACGrid } from "./src/mac-grid.js";
import { createAdvectParticlesKernel } from "./src/simulation/advect-particles.js";
import { createParticleToGridKernel } from "./src/simulation/transfer-particle-to-grid.js";
import { compileKernels } from "./src/simulation/kernels.js";

const particles = new Particles(1, {
  min: vec3.fromValues(0.0, 0.0, 0.0),
  max: vec3.fromValues(1.0, 1.0, 1.0),
});
particles.particleBuffer[0] = 0.25;
particles.particleBuffer[1] = 0.25;
particles.particleBuffer[2] = 0.25;
particles.particleBuffer[3] = 1;
particles.particleBuffer[4] = 1;
particles.particleBuffer[5] = 1;

console.log("old particle buffer:");
console.log(particles.particleBuffer);

// create MAC grid with one cell
const grid = new MACGrid(
  {
    min: vec3.fromValues(0.0, 0.0, 0.0),
    max: vec3.fromValues(1.0, 1.0, 1.0),
  },
  0.5
);
console.log(`created grid of size: (${grid.nx}, ${grid.ny}, ${grid.nz})`);

// compile kernels
const gpu = new GPU();
const kernels = compileKernels(gpu, particles, grid);
