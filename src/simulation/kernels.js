import { createAddGravityKernel } from "./add-gravity.js";
import { createAdvectParticlesKernel } from "./advect-particles.js";
import { createClassifyVoxelsKernel } from "./classify-voxels.js";
import { createCopyKernel } from "./copy-kernel.js";
import {
  createEnforceBoundaryXKernel,
  createEnforceBoundaryYKernel,
  createEnforceBoundaryZKernel,
} from "./enforce-boundary-conditions.js";
import { createParticleToGridKernel } from "./transfer-particle-to-grid.js";
import { createClearGridKernel } from "./clear-grid-velocities.js";

export const compileKernels = (gpu, particles, grid) => {
  const start = Date.now();

  const pressureSize = [grid.nx, grid.ny, grid.nz];
  const velocityXSize = [grid.nx + 1, grid.ny, grid.nz];
  const velocityYSize = [grid.nx, grid.ny + 1, grid.nz];
  const velocityZSize = [grid.nx, grid.ny, grid.nz + 1];

  // make blank velocity arrays
  const clearVelocityX = createClearGridKernel(gpu, ...velocityXSize);
  const clearVelocityY = createClearGridKernel(gpu, ...velocityYSize);
  const clearVelocityZ = createClearGridKernel(gpu, ...velocityZSize);

  // project particle velocities to the grid
  const particleToXGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityXSize
  );
  const particleToYGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityYSize
  );
  const particleToZGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityZSize
  );

  // copy grid quantities to save
  const copyPressure = createCopyKernel(gpu, ...pressureSize);
  const copyXVelocity = createCopyKernel(gpu, ...velocityXSize);
  const copyYVelocity = createCopyKernel(gpu, ...velocityYSize);
  const copyZVelocity = createCopyKernel(gpu, ...velocityZSize);

  // mark cells as solid, fluid, or air
  const classifyVoxels = createClassifyVoxelsKernel(
    gpu,
    particles.count(),
    ...pressureSize
  );

  // add gravitational influence
  const addGravity = createAddGravityKernel(gpu, ...velocityYSize);

  // enforce boundary conditions
  const enforceXBoundary = createEnforceBoundaryXKernel(gpu, ...velocityXSize);
  const enforceYBoundary = createEnforceBoundaryYKernel(gpu, ...velocityYSize);
  const enforceZBoundary = createEnforceBoundaryZKernel(gpu, ...velocityZSize);

  // do pressure solve
  // TODO: implement these kernels

  // update the velocities of the particles using PIC/FLIP
  // TODO: implement this shader

  // update the positions of the particles
  const advectParticles = createAdvectParticlesKernel(gpu, particles.count());

  const end = Date.now();
  console.log(`Kernels compiled in ${end - start} ms.`);

  return {
    clearVelocityX: clearVelocityX,
    clearVelocityY: clearVelocityY,
    clearVelocityZ: clearVelocityZ,
    particleToXGrid: particleToXGrid,
    particleToYGrid: particleToYGrid,
    particleToZGrid: particleToZGrid,
    copyPressure: copyPressure,
    copyXVelocity: copyXVelocity,
    copyYVelocity: copyYVelocity,
    copyZVelocity: copyZVelocity,
    classifyVoxels: classifyVoxels,
    addGravity: addGravity,
    enforceXBoundary: enforceXBoundary,
    enforceYBoundary: enforceYBoundary,
    enforceZBoundary: enforceZBoundary,
    advectParticles: advectParticles,
  };
};
