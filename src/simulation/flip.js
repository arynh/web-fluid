import { ATTRIBUTE_COUNT } from "../particles.js";

// kernel for subracting the new grid velocities from the old grid velocities
export const createGridVelocityDifferenceKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (oldVelocities, newVelocities) {
      return (
        oldVelocities[this.thread.x][this.thread.y][this.thread.z] -
        newVelocities[this.thread.x][this.thread.y][this.thread.z]
      );
    })
    .setOutput([nx, ny, nz]);

// FLIP Kernel
export const createFLIPKernel = (gpu, particleCount, cellSize) =>
  gpu
    .createKernel(function (particles, diffGridVx, diffGridVy, diffGridVz) {
      // mod to figure out which index we are at (0-5) for each particle
      let index_mod = this.thread.x % this.constants.ATTRIBUTE_COUNT;
      // if we are looking at the position just return the position
      if (index_mod === 0 || index_mod === 1 || index_mod === 2) {
        return particles[this.thread.x];
      }
      // get the positions - index changes depending which velocity we are looking at
      let pos_x = particles[this.thread.x - index_mod];
      let pos_y = particles[this.thread.x - index_mod + 1];
      let pos_z = particles[this.thread.x - index_mod + 2];
      // get the lower and upper grid positions
      let grid_lower_x = Math.floor(pos_x / this.constants.CELL_SIZE);
      let grid_upper_x = Math.ceil(pos_x / this.constants.CELL_SIZE);
      let grid_lower_y = Math.floor(pos_y / this.constants.CELL_SIZE);
      let grid_upper_y = Math.ceil(pos_y / this.constants.CELL_SIZE);
      let grid_lower_z = Math.floor(pos_z / this.constants.CELL_SIZE);
      let grid_upper_z = Math.ceil(pos_z / this.constants.CELL_SIZE);
      if (index_mod === 3) {
        // vx
        // get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_x - grid_lower_x) / this.constants.CELL_SIZE;
        return lerp(
          diffGridVx[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVx[grid_upper_x][grid_lower_y][grid_lower_z],
          lerpWeight
        );
      } else if (index_mod === 4) {
        // vy
        // get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_y - grid_lower_y) / this.constants.CELL_SIZE;
        return lerp(
          diffGridy[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVy[grid_lower_x][grid_upper_y][grid_lower_z],
          lerpWeight
        );
      } else if (index_mod === 5) {
        // vz
        // get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_z - grid_lower_z) / this.constants.CELL_SIZE;
        return lerp(
          diffGridVz[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVz[grid_lower_x][grid_lower_y][grid_upper_z],
          lerpWeight
        );
      }
    })
    .addFunction(function lerp(a, b, t) {
      return (1 - t) * a + t * b;
    })
    .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT, CELL_SIZE: cellSize })
    .setOutput([particleCount * ATTRIBUTE_COUNT]);

export const createUpdateVelocitiesKernel = (
  gpu,
  particleCount,
  nx,
  ny,
  nz
) => {
  const velocityXDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx + 1,
    ny,
    nz
  );
  const velocityYDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx,
    ny + 1,
    nz
  );
  const velocityZDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx,
    ny,
    nz + 1
  );
  const flipKernel = createFLIPKernel(gpu, particleCount);
  return gpu.combineKernels(
    velocityXDifference,
    velocityYDifference,
    velocityZDifference,
    flipKernel,
    function (
      oldXVelocity,
      oldYVelocity,
      oldZVelocity,
      newXVelocity,
      newYVelocity,
      newZVelocity,
      particles
    ) {
      return flipKernel(
        particles,
        velocityXDifference(oldXVelocity, newXVelocity),
        velocityYDifference(oldYVelocity, newYVelocity),
        velocityZDifference(oldZVelocity, newZVelocity)
      );
    }
  );
};
