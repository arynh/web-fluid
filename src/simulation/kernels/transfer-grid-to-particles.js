import { ATTRIBUTE_COUNT } from "../particles.js";

// kernel for subracting the new grid velocities from the old grid velocities
const createGridVelocityDifferenceKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (oldVelocities, newVelocities) {
      return (
        newVelocities[this.thread.z][this.thread.y][this.thread.x] -
        oldVelocities[this.thread.z][this.thread.y][this.thread.x]
      );
    })
    .setOutput([nx, ny, nz]);

// FLIP Kernel
const createFLIPKernel = (gpu, particleCount, cellSize) =>
  gpu
    .createKernel(function (
      particles,
      diffGridVx,
      diffGridVy,
      diffGridVz,
      oldVx,
      oldVy,
      oldVz
    ) {
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
        let lerpWeight =
          (pos_x - grid_lower_x * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        return (
          lerp(
            oldVx[grid_lower_z][grid_lower_y][grid_lower_x],
            oldVx[grid_lower_z][grid_lower_y][grid_upper_x],
            lerpWeight
          ) +
          lerp(
            diffGridVx[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVx[grid_lower_z][grid_lower_y][grid_upper_x],
            lerpWeight
          )
        );
      } else if (index_mod === 4) {
        // vy
        // get the lerp weight and return the lerp'd velocity
        let lerpWeight =
          (pos_y - grid_lower_y * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        return (
          lerp(
            oldVy[grid_lower_z][grid_lower_y][grid_lower_x],
            oldVy[grid_lower_z][grid_upper_y][grid_lower_x],
            lerpWeight
          ) +
          lerp(
            diffGridVy[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVy[grid_lower_z][grid_upper_y][grid_lower_x],
            lerpWeight
          )
        );
      } else if (index_mod === 5) {
        // vz
        // get the lerp weight and return the lerp'd velocity
        let lerpWeight =
          (pos_z - grid_lower_z * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        return (
          lerp(
            oldVz[grid_lower_z][grid_lower_y][grid_lower_x],
            oldVz[grid_upper_z][grid_lower_y][grid_lower_x],
            lerpWeight
          ) +
          lerp(
            diffGridVz[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVz[grid_upper_z][grid_lower_y][grid_lower_x],
            lerpWeight
          )
        );
      }
    })
    .addFunction(function lerp(a, b, t) {
      return (1 - t) * a + t * b;
    })
    .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT, CELL_SIZE: cellSize })
    .setOutput([particleCount * ATTRIBUTE_COUNT]);

export const createGridToParticlesKernel = (
  gpu,
  particleCount,
  nx,
  ny,
  nz,
  cellSize
) => {
  const velocityXDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx + 1,
    ny,
    nz
  ).setPipeline(true);

  const velocityYDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx,
    ny + 1,
    nz
  ).setPipeline(true);

  const velocityZDifference = createGridVelocityDifferenceKernel(
    gpu,
    nx,
    ny,
    nz + 1
  ).setPipeline(true);

  const flipKernel = createFLIPKernel(gpu, particleCount, cellSize).setPipeline(
    true
  );

  return (
    oldXVelocity,
    oldYVelocity,
    oldZVelocity,
    newXVelocity,
    newYVelocity,
    newZVelocity,
    particles
  ) =>
    flipKernel(
      particles,
      velocityXDifference(oldXVelocity, newXVelocity),
      velocityYDifference(oldYVelocity, newYVelocity),
      velocityZDifference(oldZVelocity, newZVelocity),
      oldXVelocity,
      oldYVelocity,
      oldZVelocity
    );
};
