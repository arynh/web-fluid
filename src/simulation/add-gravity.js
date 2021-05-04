export const createAddGravityKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocity_y, dt) {
      return (
        velocity_y[this.thread.x][this.thread.y][this.thread.z] - dt * 9.81
      );
    })
    .setOutput([nx, ny, nz]);
