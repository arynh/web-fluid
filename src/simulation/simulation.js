import { MACGrid } from "./mac-grid.js";
import { Particles } from "./particles.js";
import { compileKernels } from "./kernels/kernels.js";

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
    this.kernels = compileKernels(gpu, this.particles, this.grid);
  }

  step(dt) {
    let particleBufferCopy = new Float32Array(this.particles.particleBuffer);
    console.log("SIM START");
    // transfer particle velocities to the grid and interpolate
    //console.log("1: " + this.grid.velocityX[0][0][0])
    console.log(this.particles.particleBuffer);
    this.grid.velocityX = this.kernels.particleToXGrid(
      particleBufferCopy,
      this.grid.cellSize
    );
    //console.log("2: " + this.grid.velocityX[0][0][0])
    //console.log(this.particles.particleBuffer)
    this.grid.velocityY = this.kernels.particleToYGrid(
      particleBufferCopy,
      this.grid.cellSize
    );
    //console.log("3: " + this.grid.velocityX[0][0][0])
    this.grid.velocityZ = this.kernels.particleToZGrid(
      particleBufferCopy,
      this.grid.cellSize
    );

    //console.log("4: " + this.grid.velocityX[0][0][1])

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
    //this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
    this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
    //this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

    // do the pressure solve with a zero divergence velocity field
    // TODO: implement this!

    // enforce boundary conditions
    // this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
    // this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
    // this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

    // update the velocities of the particles
    console.log("PRE " + this.particles.particleBuffer);
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
    console.log("POST " + this.particles.particleBuffer);
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
