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
import {
  createComponentWiseAddKernel,
  createComponentWiseMultiplyKernel,
  createScalarMultiplyKernel,
} from "./pressure-solve/vector-math.js";
import {
  createADiagKernel,
  createAXKernel,
  createAYKernel,
  createAZKernel,
} from "./pressure-solve/build-coefficient-matrix.js";
import { createNegativeDivergenceKernel } from "./pressure-solve/negative-divergence.js";
import {
  createFlattenKernel,
  createUnflattenKernel,
} from "./pressure-solve/flatten.js";
import { createApplyAKernel } from "./pressure-solve/apply-a.js";
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

  // build coefficient matrix
  const buildADiag = createADiagKernel(gpu, ...gridSize, grid.cellSize);
  const buildAX = createAXKernel(gpu, ...gridSize, grid.cellSize);
  const buildAY = createAYKernel(gpu, ...gridSize, grid.cellSize);
  const buildAZ = createAZKernel(gpu, ...gridSize, grid.cellSize);

  // build negative divergence vector with boundary conditions
  const buildD = createNegativeDivergenceKernel(
    gpu,
    ...gridSize,
    grid.cellSize
  );
  const flatten = createFlattenKernel(gpu, ...gridSize);
  const unflatten = createUnflattenKernel(gpu, ...gridSize);

  // compile kernels to do vector operations
  const pcgVectorLength = grid.nx * grid.ny * grid.nz;
  const componentWiseAdd = createComponentWiseAddKernel(gpu, pcgVectorLength);
  const componentWiseMultiply = createComponentWiseMultiplyKernel(
    gpu,
    pcgVectorLength
  );
  // implement dot product's sum portion on the CPU
  const dot = (a, b) =>
    componentWiseMultiply(a, b).reduce((sum, n) => sum + n, 0);
  const scalarMultiply = createScalarMultiplyKernel(gpu, pcgVectorLength);
  const applyA = createApplyAKernel(gpu, pcgVectorLength, ...gridSize);
  const math = {
    componentWiseAdd: componentWiseAdd.setPipeline(true).setImmutable(true),
    dot: dot,
    scalarMultiply: scalarMultiply.setPipeline(true).setImmutable(true),
    applyA: applyA.setPipeline(true).setImmutable(true),
  };

  // PCG methods
  const zeroVector = gpu
    .createKernel(function () {
      return 0;
    })
    .setOutput([pcgVectorLength]);

  const jacobi = createJacobiIterationKernel(gpu, ...gridSize);

  const pressureSolve = {
    buildADiag: buildADiag.setPipeline(true),
    buildAX: buildAX.setPipeline(true),
    buildAY: buildAY.setPipeline(true),
    buildAZ: buildAZ.setPipeline(true),
    buildD: buildD.setPipeline(true),
    flatten: flatten.setPipeline(true).setImmutable(true),
    unflatten: unflatten.setPipeline(true),
    math: math,
    zeroVector: zeroVector.setPipeline(true).setImmutable(true),
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
    particleToXGrid: particleToXGrid.setPipeline(true),
    particleToYGrid: particleToYGrid.setPipeline(true),
    particleToZGrid: particleToZGrid.setPipeline(true),
    copyPressure: copyPressure.setPipeline(true),
    copyXVelocity: copyXVelocity.setPipeline(true),
    copyYVelocity: copyYVelocity.setPipeline(true),
    copyZVelocity: copyZVelocity.setPipeline(true),
    classifyVoxels: classifyVoxels.setPipeline(true).setImmutable(true),
    addGravity: addGravity.setPipeline(true).setImmutable(true),
    enforceXBoundary: enforceXBoundary.setPipeline(true).setImmutable(true),
    enforceYBoundary: enforceYBoundary.setPipeline(true).setImmutable(true),
    enforceZBoundary: enforceZBoundary.setPipeline(true).setImmutable(true),
    gridToParticles: gridToParticles,
    advectParticles: advectParticles.setPipeline(true),
    pressureSolve: pressureSolve,
    updateVelocityX: updateVelocityX.setPipeline(true).setImmutable(true),
    updateVelocityY: updateVelocityY.setPipeline(true).setImmutable(true),
    updateVelocityZ: updateVelocityZ.setPipeline(true).setImmutable(true),
  };
};
