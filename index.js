import { GPU } from "gpu.js";
import { vec3 } from "gl-matrix";
import { createClassifyVoxelsKernel } from "./src/simulation/classify-voxels.js";
import { Particles } from "./src/particles.js";
import { MACGrid } from "./src/mac-grid.js";

// create one test particle with velocity <1,1,1>
const particles = new Particles(1, {
  min: vec3.fromValues(0.0, 0.0, 0.0),
  max: vec3.fromValues(1.0, 1.0, 1.0),
});
particles.particleBuffer[0] = 0.4;
particles.particleBuffer[1] = 0.4;
particles.particleBuffer[2] = 0.4;
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
  0.25
);
grid.addDefaultSolids();
console.log("grid states:");
console.log(grid.voxelStates);

const particleCount = particles.count();

const gpu = new GPU();

//create classify voxels kernel
const classifyVoxelsKernel = createClassifyVoxelsKernel(
  gpu,
  particleCount,
  grid.nx,
  grid.ny,
  grid.nz
);

//classify the voxels using the gpu
const classifyVoxels = classifyVoxelsKernel(
  grid.voxelStates,
  particles.particleBuffer,
  grid.cellSize
);

console.log("\nnew voxels:");
console.log(classifyVoxels);
