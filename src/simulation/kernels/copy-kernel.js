export const createCopyKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (array) {
      return array[this.thread.x][this.thread.y][this.thread.z];
    })
    .setOutput([nx, ny, nz]);
