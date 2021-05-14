import { STATE_ENUM } from "../../mac-grid.js";

export const createNegativeDivergenceKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, velocityX, velocityY, velocityZ) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;

      if (voxelStates[k][j][i] !== this.constants.FLUID) {
        return 0;
      }

      const scale = 1.0 / this.constants.CELL_SIZE;

      let divergence =
        -scale *
        (velocityX[k][j][i + 1] -
          velocityX[k][j][i] +
          velocityY[k][j + 1][i] -
          velocityY[k][j][i] +
          velocityZ[k + 1][j][i] -
          velocityZ[k][j][i]);

      /* TODO: experiment with these boundary conditions to see if we can
          lose less energy

      // modifying RHS (divergence) to account for solid velocities
      if (voxelStates[k][j][i - 1] === this.constants.SOLID) {
        divergence -= scale * velocityX[k][j][i];
      }
      if (voxelStates[k][j][i + 1] === this.constants.SOLID) {
        divergence += scale * velocityX[k][j][i + 1];
      }

      if (voxelStates[k][j - 1][i] === this.constants.SOLID) {
        divergence -= scale * velocityY[k][j][i];
      }
      if (voxelStates[k][j + 1][i] === this.constants.SOLID) {
        divergence += scale * velocityY[k][j + 1][i];
      }

      if (voxelStates[k - 1][j][i] === this.constants.SOLID) {
        divergence -= scale * velocityZ[k][j][i];
      }

      if (voxelStates[k + 1][j][i] === this.constants.SOLID) {
        divergence += scale * velocityZ[k + 1][j][i];
      }
      */

      return divergence;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);
