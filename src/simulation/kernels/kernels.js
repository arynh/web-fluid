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
import { createGridToParticlesKernel } from "../transfer-grid-to-particles.js";
import {
  createApplyAKernel,
  createComponentWiseAddKernel,
  createComponentWiseMultiplyKernel,
  createScalarMultiplyKernel,
} from "./pressure-solve/vector-math.js";
import {
  createADiagKernel,
  createAXKernel,
  createAYKernel,
  createAZKernel,
} from "./kernels/pressure-solve/build-coefficient-matrix.js.js";
import { createNegativeDivergenceKernel } from "./kernels/pressure-solve/negative-divergence.js.js";

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
  const applyA = createApplyAKernel(gpu, pcgVectorLength);
  const math = {
    componentWiseAdd: componentWiseAdd,
    dot: dot,
    scalarMultiply: scalarMultiply,
    applyA: applyA,
  };

  // PCG methods
  const zeroVector = gpu
    .createKernel(function () {
      return 0;
    })
    .setOutput([pcgVectorLength]);

  const pressureSolve = {
    buildADiag: buildADiag,
    buildAX: buildAX,
    buildAY: buildAY,
    buildAZ: buildAZ,
    buildD: buildD,
    math: math,
    zeroVector: zeroVector,
  };

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
    gridToParticles: gridToParticles,
    advectParticles: advectParticles,
    pressureSolve: pressureSolve,
  };
};