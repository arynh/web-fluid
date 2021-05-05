import { GPU } from "gpu.js";
import { vec3 } from "gl-matrix";
import { createParticleToGridKernel } from "./src/simulation/transfer-particle-to-grid.js";
import { Particles } from "./src/particles.js";
import { MACGrid } from "./src/mac-grid.js";

// create one test particle with velocity <1,1,1>
const particles = new Particles(10, {
  min: vec3.fromValues(0.0, 0.0, 0.0),
  max: vec3.fromValues(1.0, 1.0, 1.0),
});
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
  0.5
);
grid.addDefaultSolids();
console.log("grid states:");
console.log(grid.voxelStates);

const particleCount = particles.count();

const gpu = new GPU();
const particleToGridKernel = createParticleToGridKernel(
  gpu,
  particleCount,
  grid.count[0],
  grid.count[1] - 1,
  grid.count[2] - 1
);

const newGridVelocitiesX = particleToGridKernel(
  particles.particleBuffer,
  grid.cellSize,
  0
);
const newGridVelocitiesY = particleToGridKernel(
  particles.particleBuffer,
  grid.cellSize,
  1
);
const newGridVelocitiesZ = particleToGridKernel(
  particles.particleBuffer,
  grid.cellSize,
  2
);

console.log("\nnew grid velocities:");
console.log(newGridVelocitiesX);
console.log(newGridVelocitiesY);
console.log(newGridVelocitiesZ);
