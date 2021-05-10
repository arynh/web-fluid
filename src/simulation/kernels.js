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
import { createGridToParticlesKernel } from "./transfer-grid-to-particles.js";

export const compileKernels = (gpu, particles, grid) => {
  const start = Date.now();

  const gridSize = [grid.nx, grid.ny, grid.nz];
  const velocityXSize = [grid.nx + 1, grid.ny, grid.nz];
  const velocityYSize = [grid.nx, grid.ny + 1, grid.nz];
  const velocityZSize = [grid.nx, grid.ny, grid.nz + 1];

  // project particle velocities to the grid
  const particleToXGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityXSize,
    0
  );
  const particleToYGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityYSize,
    1
  );
  const particleToZGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityZSize,
    2
  );

  // copy grid quantities to save
  const copyPressure = createCopyKernel(gpu, ...gridSize);
  const copyXVelocity = createCopyKernel(gpu, ...velocityXSize);
  const copyYVelocity = createCopyKernel(gpu, ...velocityYSize);
  const copyZVelocity = createCopyKernel(gpu, ...velocityZSize);

  // mark cells as solid, fluid, or air
  const classifyVoxels = createClassifyVoxelsKernel(
    gpu,
    particles.count(),
    ...gridSize
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
  const gridToParticles = createGridToParticlesKernel(
    gpu,
    particles.count(),
    ...gridSize,
    grid.cellSize
  );

  // update the positions of the particles
  const advectParticles = createAdvectParticlesKernel(
    gpu,
    particles.count(),
    grid.cellSize
  );

  const end = Date.now();
  console.log(`Kernels compiled in ${end - start} ms.`);

  return {
    particleToXGrid: particleToXGrid.setPipeline(true),
    particleToYGrid: particleToYGrid.setPipeline(true),
    particleToZGrid: particleToZGrid.setPipeline(true),
    copyPressure: copyPressure.setPipeline(true),
    copyXVelocity: copyXVelocity.setPipeline(true),
    copyYVelocity: copyYVelocity.setPipeline(true),
    copyZVelocity: copyZVelocity.setPipeline(true),
    classifyVoxels: classifyVoxels.setPipeline(true),
    addGravity: addGravity.setPipeline(true),
    enforceXBoundary: enforceXBoundary.setPipeline(true),
    enforceYBoundary: enforceYBoundary.setPipeline(true),
    enforceZBoundary: enforceZBoundary.setPipeline(true),
    gridToParticles: gridToParticles,
    advectParticles: advectParticles.setPipeline(true),
  };
};
