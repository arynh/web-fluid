export const createClearGridKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function () {
      return 0;
    })
    .setOutput([nx, ny, nz]);
