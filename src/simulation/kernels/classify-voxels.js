import { STATE_ENUM } from "../mac-grid.js";

export const createClassifyVoxelsKernel = (gpu, particleCount, nx, ny, nz) =>
  gpu
    .createKernel(function (voxelStates, particles, cellSize) {
      // get spatial location of grid
      let x = cellSize * this.thread.x;
      let y = cellSize * this.thread.y;
      let z = cellSize * this.thread.z;

      let particle_exists = false;
      for (let i = 0; i < this.constants.particleCount; i++) {
        let pos_x = particles[i * 6];
        let pos_y = particles[i * 6 + 1];
        let pos_z = particles[i * 6 + 2];
        // check if there is a particle in that grid
        if (
          pos_x - x <= cellSize &&
          pos_x - x > 0 &&
          pos_y - y <= cellSize &&
          pos_y - y > 0 &&
          pos_z - z <= cellSize &&
          pos_z - z > 0
        ) {
          particle_exists = true;
          break;
        }
      }
      // set any cell with a particle in it that isn't solid to fluid
      // any cell without a particle or solid is air
      if (particle_exists) {
        if (
          voxelStates[this.thread.x][this.thread.y][this.thread.z] !==
          this.constants.SOLID
        ) {
          return this.constants.FLUID;
        } else {
          return this.constants.SOLID;
        }
      } else {
        // if there isn't, and the state is fluid, flip it to air
        if (
          voxelStates[this.thread.x][this.thread.y][this.thread.z] !==
          this.constants.SOLID
        ) {
          return this.constants.AIR;
        } else {
          return this.constants.SOLID;
        }
      }
    })
    .setConstants({
      particleCount: particleCount,
      AIR: STATE_ENUM.AIR,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);
