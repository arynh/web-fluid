import { STATE_ENUM } from "../../mac-grid.js";

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

      // modifying RHS (divergence) to account for solid velocities
      if (voxelStates[i - 1][j][k] === SOLID) {
        divergence -= scale * velocityX[i][j][k];
      }
      if (voxelStates[i + 1][j][k] === SOLID) {
        divergence += scale * velocityX[i + 1][j][k];
      }

      if (voxelStates[i][j - 1][k] === SOLID) {
        divergence -= scale * velocityY[i][j][k];
      }
      if (voxelStates[i][j + 1][k] === SOLID) {
        divergence += scale * velocityY[i][j + 1][k];
      }

      if (voxelStates[i][j][k - 1] === SOLID) {
        divergence -= scale * velocityZ[i][j][k];
      }

      if (voxelStates[i][j][k + 1] === SOLID) {
        divergence += scale * velocityZ[i][j][k + 1];
      }

      return divergence;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);
