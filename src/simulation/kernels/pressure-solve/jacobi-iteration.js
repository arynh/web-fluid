import { STATE_ENUM } from "../../mac-grid.js";

export const createJacobiIterationKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (negativeDivergence, pressure, voxelStates) {
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;

      if (voxelStates[k][j][i] === this.constants.AIR) {
        return 0;
      } else if (voxelStates[k][j][i] === this.constants.SOLID) {
        return 100;
      }

      const divergenceCenter = negativeDivergence[k][j][i];

      let left = 0.0;
      let right = 0.0;
      let bottom = 0.0;
      let top = 0.0;
      let back = 0.0;
      let front = 0.0;

      // negative x neighbor
      if (voxelStates[k][j][i - 1] === this.constants.FLUID) {
        left = pressure[k][j][i - 1];
      }
      // positive x neighbor
      if (voxelStates[k][j][i + 1] === this.constants.FLUID) {
        right = pressure[k][j][i + 1];
      }
      // negative y neighbor
      if (voxelStates[k][j - 1][i] === this.constants.FLUID) {
        bottom = pressure[k][j - 1][i];
      }
      // positive y neighbor
      if (voxelStates[k][j + 1][i] === this.constants.FLUID) {
        top = pressure[k][j + 1][i];
      }
      // negative z neighbor
      if (voxelStates[k - 1][j][i] === this.constants.FLUID) {
        front = pressure[k - 1][j][i];
      }
      // positive z neighbor
      if (voxelStates[k + 1][j][i] === this.constants.FLUID) {
        back = pressure[k + 1][j][i];
      }

      return (
        (left + right + bottom + top + back + front + divergenceCenter) / 6.0
      );
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      FLUID: STATE_ENUM.FLUID,
      AIR: STATE_ENUM.AIR,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);
