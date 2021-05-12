/**
 * Map from 3D arrays to 1D vectors.
 *
 * To unflatten:
 * `(k * ny + i) * nx + j`
 */
export const createFlattenKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (array) {
      const aux = this.thread.x % (this.constants.NX * this.constants.NY);
      const i = Math.floor(aux / this.constants.NY);
      const j = aux % this.constants.NX;
      const k = Math.floor(
        this.thread.x / (this.constants.NX * this.constants.NY)
      );
      return array[k][j][i];
    })
    .setTactic("precision")
    .setConstants({ NX: nx, NY: ny, NZ: nz })
    .setOutput([nx * ny * nz]);

export const createUnflattenKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (flat) {
      return flat[
        (this.thread.z * this.constants.NY + this.thread.x) *
          this.constants.NX +
          this.thread.y
      ];
    })
    .setTactic("precision")
    .setConstants({ NX: nx, NY: ny, NZ: nz })
    .setOutput([nx, ny, nz]);
