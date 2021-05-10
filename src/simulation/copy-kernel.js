export const createCopyKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (array) {
      return array[this.thread.z][this.thread.y][this.thread.x];
    })
    .setOutput([nx, ny, nz]);
