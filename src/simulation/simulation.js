import { MACGrid } from "./mac-grid.js";
import { Particles } from "./particles.js";
import { solve } from "./pressure-solve.js";
import { compileKernels } from "./kernels/kernels.js";

export const FLUID_DENSITY = 997;
const SOLVER_TOLERANCE = 1e-4;
const SOLVER_ITERATION_LIMIT = 200;

export class Simulation {
  constructor(gpu, config) {
    this.particles = new Particles(
      config.particleDensity,
      config.particleBounds
    );
    this.grid = new MACGrid(
      config.gridBounds,
      2.0 / Math.cbrt(config.particleDensity)
    );
    this.grid.addDefaultSolids();
    this.kernels = compileKernels(gpu, this.particles, this.grid);
  }

  step(dt) {
    let particleBufferCopy = new Float32Array(this.particles.particleBuffer);
    // transfer particle velocities to the grid and interpolate
    this.grid.velocityX = this.kernels.particleToXGrid(
      particleBufferCopy,
      this.grid.cellSize
    );
    this.grid.velocityY = this.kernels.particleToYGrid(
      particleBufferCopy,
      this.grid.cellSize
    );
    this.grid.velocityZ = this.kernels.particleToZGrid(
      particleBufferCopy,
      this.grid.cellSize
    );

    // copy grid values to store the old ones
    this.grid.pressureOld = this.kernels.copyPressure(this.grid.pressure);
    this.grid.velocityXOld = this.kernels.copyXVelocity(this.grid.velocityX);
    this.grid.velocityYOld = this.kernels.copyYVelocity(this.grid.velocityY);
    this.grid.velocityZOld = this.kernels.copyZVelocity(this.grid.velocityZ);

    // mark cells as solid, fluid, or air
    this.grid.voxelStates = this.kernels.classifyVoxels(
      this.grid.voxelStates.toArray(),
      particleBufferCopy,
      this.grid.cellSize
    );

    // perform gravity update
    this.grid.velocityY = this.kernels.addGravity(this.grid.velocityY, dt);

    // enforce boundary conditions
    this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
    this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
    this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

    // do the pressure solve with a zero divergence velocity field
    this.grid.pressure = solve(
      this.kernels.pressureSolve,
      this.grid.voxelStates,
      dt,
      this.grid.velocityX,
      this.grid.velocityY,
      this.grid.velocityZ,
      SOLVER_TOLERANCE,
      SOLVER_ITERATION_LIMIT,
      this.grid.pressure,
      this.grid.pressureOld
    );

    // update the velocity fields with the new pressure gradients
    this.grid.velocityX = this.kernels.updateVelocityX(
      this.grid.velocityX,
      this.grid.pressure,
      this.grid.voxelStates,
      dt
    );
    this.grid.velocityY = this.kernels.updateVelocityY(
      this.grid.velocityY,
      this.grid.pressure,
      this.grid.voxelStates,
      dt
    );
    this.grid.velocityZ = this.kernels.updateVelocityZ(
      this.grid.velocityZ,
      this.grid.pressure,
      this.grid.voxelStates,
      dt
    );

    // enforce boundary conditions
    this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
    this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
    this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

    // update the velocities of the particles
    this.particles.particleBuffer = this.kernels
      .gridToParticles(
        this.grid.velocityXOld,
        this.grid.velocityYOld,
        this.grid.velocityZOld,
        this.grid.velocityX,
        this.grid.velocityY,
        this.grid.velocityZ,
        particleBufferCopy
      )
      .toArray();
    // advect the particles to find their new positions
    this.particles.particleBuffer = this.kernels
      .advectParticles(
        new Float32Array(this.particles.particleBuffer),
        dt,
        this.grid.velocityX,
        this.grid.velocityY,
        this.grid.velocityZ
      )
      .toArray();
  }
}
