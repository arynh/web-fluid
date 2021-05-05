import { GPU } from "gpu.js";
import { vec3 } from "gl-matrix";
import { initialize3DArray } from "./src/utils.js";
import { createClearGridKernel } from "./src/simulation/clear-grid-velocities.js";
import { createAddGravityKernel } from "./src/simulation/add-gravity.js";
import {
  createParticleToGridKernel,
  createParticleWeightsKernel,
  createWeightedSumKernel,
  createNewVelocitiesKernel,
} from "./src/simulation/transfer-particle-to-grid.js";
import { Particles } from "./src/particles.js";
import { MACGrid } from "./src/mac-grid.js";

// create one test particle with velocity <1,1,1>
const particles = new Particles(1, {
  min: vec3.fromValues(0.0, 0.0, 0.0),
  max: vec3.fromValues(1.0, 1.0, 1.0),
});
particles.particleBuffer[0] = 0.5;
particles.particleBuffer[1] = 0.5;
particles.particleBuffer[2] = 0.5;
particles.particleBuffer[3] = 2;
particles.particleBuffer[4] = 1;
particles.particleBuffer[5] = 1;
console.log("particle buffer:");
console.log(particles.particleBuffer);
console.log("particle indices:");
console.log(particles.particleIndices);

// create MAC grid with one cell
const grid = new MACGrid(
  {
    min: vec3.fromValues(0.0, 0.0, 0.0),
    max: vec3.fromValues(1.0, 1.0, 1.0),
  },
  1
);
console.log("grid pressure:");
console.log(grid.pressure);
console.log("grid x velocity:");
console.log(grid.velocity_x);

const particleCount = particles.count();
const edgeCount = grid.count[0];

const gpu = new GPU();
const particleWeightsKernel = createParticleWeightsKernel(
  gpu,
  particleCount,
  edgeCount
);
const weightedSumKernel = createWeightedSumKernel(
  gpu,
  particleCount,
  edgeCount
);
const newVelocitiesKernel = createNewVelocitiesKernel(gpu, edgeCount);

const weights = particleWeightsKernel(particles.particleBuffer, 1, 0);
console.log("\nweights:");
console.log(weights);

const sums = weightedSumKernel(weights);
console.log("\nsums:");
console.log(sums);
const newGridVelocities = newVelocitiesKernel(sums);

// const particleToGridKernel = createParticleToGridKernel(gpu, 1, 2);
// const newGridVelocities = particleToGridKernel(particles.particleBuffer, 1, 0);

console.log("\nnew grid velocities:");
console.log(newGridVelocities);
