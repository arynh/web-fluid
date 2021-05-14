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
import { createNegativeDivergenceKernel } from "./pressure-solve/negative-divergence.js";
import {
  createVelocityXUpdateKernel,
  createVelocityYUpdateKernel,
  createVelocityZUpdateKernel,
} from "./velocity-update.js";
import { FLUID_DENSITY } from "../simulation.js";
import { createJacobiIterationKernel } from "./pressure-solve/jacobi-iteration.js";

export const compileKernels = (gpu, particles, grid) => {
  const start = Date.now();

  const gridSize = [grid.nx, grid.ny, grid.nz];
  const velocityXSize = [grid.nx + 1, grid.ny, grid.nz];
  const velocityYSize = [grid.nx, grid.ny + 1, grid.nz];
  const velocityZSize = [grid.nx, grid.ny, grid.nz + 1];
  const DIMENSION = { X: 0, Y: 1, Z: 2 };

  // project particle velocities to the grid
  const particleToXGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityXSize,
    DIMENSION.X
  );
  const particleToYGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityYSize,
    DIMENSION.Y
  );
  const particleToZGrid = createParticleToGridKernel(
    gpu,
    particles.count(),
    ...velocityZSize,
    DIMENSION.Z
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

  // build negative divergence vector with boundary conditions
  const buildD = createNegativeDivergenceKernel(
    gpu,
    ...gridSize,
    grid.cellSize
  );

  const jacobi = createJacobiIterationKernel(gpu, ...gridSize);

  const pressureSolve = {
    buildD: buildD.setPipeline(true).setImmutable(true),
    jacobi: jacobi.setPipeline(true).setImmutable(true),
  };

  // update grid velocities using the pressure gradient
  const updateVelocityX = createVelocityXUpdateKernel(
    gpu,
    ...velocityXSize,
    FLUID_DENSITY,
    grid.cellSize
  );
  const updateVelocityY = createVelocityYUpdateKernel(
    gpu,
    ...velocityYSize,
    FLUID_DENSITY,
    grid.cellSize
  );
  const updateVelocityZ = createVelocityZUpdateKernel(
    gpu,
    ...velocityZSize,
    FLUID_DENSITY,
    grid.cellSize
  );

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
    grid.cellSize,
    ...gridSize
  );

  const end = Date.now();
  console.log(`Kernels compiled in ${end - start} ms.`);

  return {
    particleToXGrid: particleToXGrid.setPipeline(true).setImmutable(true),
    particleToYGrid: particleToYGrid.setPipeline(true).setImmutable(true),
    particleToZGrid: particleToZGrid.setPipeline(true).setImmutable(true),
    copyPressure: copyPressure.setPipeline(true).setImmutable(true),
    copyXVelocity: copyXVelocity.setPipeline(true).setImmutable(true),
    copyYVelocity: copyYVelocity.setPipeline(true).setImmutable(true),
    copyZVelocity: copyZVelocity.setPipeline(true).setImmutable(true),
    classifyVoxels: classifyVoxels.setPipeline(true).setImmutable(true),
    addGravity: addGravity.setPipeline(true).setImmutable(true),
    enforceXBoundary: enforceXBoundary.setPipeline(true).setImmutable(true),
    enforceYBoundary: enforceYBoundary.setPipeline(true).setImmutable(true),
    enforceZBoundary: enforceZBoundary.setPipeline(true).setImmutable(true),
    gridToParticles: gridToParticles,
    advectParticles: advectParticles.setPipeline(true).setImmutable(true),
    pressureSolve: pressureSolve,
    updateVelocityX: updateVelocityX.setPipeline(true).setImmutable(true),
    updateVelocityY: updateVelocityY.setPipeline(true).setImmutable(true),
    updateVelocityZ: updateVelocityZ.setPipeline(true).setImmutable(true),
  };
};
