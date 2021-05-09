import { STATE_ENUM } from "../mac-grid.js";

export const createNegativeDivergenceKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, velocityX, velocityY, velocityZ) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;

      if (voxelStates[i][j][k] !== FLUID) {
        return 0;
      }

      const scale = 1.0 / this.constants.CELL_SIZE;

      let divergence =
        -scale *
        (velocityX[i + 1][j][k] -
          velocityX[i][j][k] +
          velocityY[i][j + 1][k] -
          velocityY[i][j][k] +
          velocityZ[i][j][k + 1] -
          velocityZ[i][j][k]);

      // TODO: add pseudocode 5.4 to adjust RHS for boundary conditions
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID: STATE_ENUM.FLUID,
    })
    .setOutput([nx, ny, nz]);
