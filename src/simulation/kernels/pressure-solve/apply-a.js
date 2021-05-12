import { STATE_ENUM } from "../../mac-grid.js";

/**
 * Produce the matrix-vector product `Ax` from the sparsely stored A and x.
 */
export const createApplyAKernel = (gpu, vectorLength, nx, ny, nz) =>
  gpu
    .createKernel(function (Adiag, Ax, Ay, Az, x, voxelStates) {
      const aux = this.thread.x % (this.constants.NX * this.constants.NY);
      const i = Math.floor(aux / this.constants.NY);
      const j = aux % this.constants.NX;
      const k = Math.floor(
        this.thread.x / (this.constants.NX * this.constants.NY)
      );

      // only consider fluid cells
      if (voxelStates[k][j][i] !== this.constants.FLUID) {
        return 0;
      }

      let vectorIndex = this.thread.x;
      let accumulator = Adiag[k][j][i] * x[vectorIndex];

      // negative x neighbor
      if (voxelStates[k][j][i - 1] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i - 1, j, k);
        accumulator += Ax[k][j][i - 1] * x[vectorIndex];
      }
      // positive x neighbor
      if (voxelStates[k][j][i + 1] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i + 1, j, k);
        accumulator += Ax[k][j][i + 1] * x[vectorIndex];
      }
      // negative y neighbor
      if (voxelStates[k][j - 1][i] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i, j - 1, k);
        accumulator += Ay[k][j - 1][i] * x[vectorIndex];
      }
      // positive y neighbor
      if (voxelStates[k][j + 1][i] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i, j + 1, k);
        accumulator += Ay[k][j + 1][i] * x[vectorIndex];
      }
      // negative z neighbor
      if (voxelStates[k - 1][j][i] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i, j, k - 1);
        accumulator += Az[k - 1][j][i] * x[vectorIndex];
      }
      // positive z neighbor
      if (voxelStates[k + 1][j][i] === this.constants.FLUID) {
        vectorIndex = gridToVectorIndex(i, j, k + 1);
        accumulator += Az[k + 1][j][i] * x[vectorIndex];
      }

      return accumulator;
    })
    .addFunction(function gridToVectorIndex(i, j, k) {
      return (k * this.constants.NY + i) * this.constants.NX + j;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      VECTOR_LENGTH: vectorLength,
      NX: nx,
      NY: ny,
      NZ: nz,
      FLUID: STATE_ENUM.FLUID,
    })
    .setOutput([vectorLength]);
