export const createAddGravityKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocity_y, dt) {
      return (
        velocity_y[this.thread.z][this.thread.y][this.thread.x] - dt * 1
      );
    })
    .setOutput([nx, ny, nz]);
