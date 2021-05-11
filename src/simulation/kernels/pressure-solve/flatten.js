/**
 * Map from 3D arrays to 1D vectors.
 */

const createFlattenKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (array) {
      const aux = this.thread.x % (this.constants.NX * this.constants.NY);
      const i = Math.floor(aux / this.constants.NY);
      const j = aux % this.constants.NX;
      const k = Math.floor(
        this.thread.x / (this.constants.NX * this.constants.NY)
      );
      return array[i][j][k];
    })
    .setTactic("precision")
    .setConstants({ NX: nx, NY: ny, NZ: nz })
    .setOutput([nx * ny * nz]);

const createUnflattenKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (vector) {})
    .setTactic("precision")
    .setOutput([nx, ny, nz]);
