export const createAddGravityKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocity_y, dt, voxelStates) {
      // non-physic gravity for the lolz
      return velocity_y[this.thread.z][this.thread.y][this.thread.x] - dt * 6.0;
    })
    .setOutput([nx, ny, nz]);
