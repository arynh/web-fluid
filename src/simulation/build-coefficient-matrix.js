import { FLUID_DENSITY } from "./simulation.js";
import { STATE_ENUM } from "../mac-grid.js";

/**
 * Assumption!
 * As long as there is the default solid walls around the fluid, edge cases
 * will be fine. Otherwise, this will need to be modified to support more
 * complex boundary conditions.
 */

export const createADiagKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, dt) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;
      const FLUID = this.constants.FLUID;
      const AIR = this.constants.AIR;

      // only consider fluid cells
      if (voxelStates[i][j][k] !== FLUID) {
        return 0;
      }

      const scale =
        dt /
        (this.constants.FLUID_DENSITY *
          this.constants.CELL_SIZE *
          this.constants.CELL_SIZE);

      let accumulator = 0;

      // negative x neighbor
      if (voxelStates[i - 1][j][k] === FLUID) {
        accumulator += scale;
      }
      // positive x neighbor
      if (
        voxelStates[i + 1][j][k] === FLUID ||
        voxelStates[i + 1][j][k] === AIR
      ) {
        accumulator += scale;
      }

      // negative y neighbor
      if (voxelStates[i][j - 1][k] === FLUID) {
        accumulator += scale;
      }
      // positive y neighbor
      if (
        voxelStates[i][j + 1][k] === FLUID ||
        voxelStates[i][j + 1][k] === AIR
      ) {
        accumulator += scale;
      }

      // negative z neighbor
      if (voxelStates[i][j][k - 1] === FLUID) {
        accumulator += scale;
      }
      // positive z neighbor
      if (
        voxelStates[i][j][k + 1] === FLUID ||
        voxelStates[i][j][k + 1] === AIR
      ) {
        accumulator += scale;
      }

      return accumulator;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID_DENSITY: FLUID_DENSITY,
      AIR: STATE_ENUM.AIR,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);

export const createAXKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, dt) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;
      const FLUID = this.constants.FLUID;

      // only consider fluid cells
      if (voxelStates[i][j][k] !== FLUID) {
        return 0;
      }

      const scale =
        dt /
        (this.constants.FLUID_DENSITY *
          this.constants.CELL_SIZE *
          this.constants.CELL_SIZE);

      let accumulator = 0;
      //positive x neighbor
      if (voxelStates[(i + 1, j, k)] === FLUID) {
        accumulator = -scale;
      }
      return accumulator;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID_DENSITY: FLUID_DENSITY,
      AIR: STATE_ENUM.AIR,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);

export const createAYKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, dt) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;
      const FLUID = this.constants.FLUID;

      // only consider fluid cells
      if (voxelStates[i][j][k] !== FLUID) {
        return 0;
      }

      const scale =
        dt /
        (this.constants.FLUID_DENSITY *
          this.constants.CELL_SIZE *
          this.constants.CELL_SIZE);

      let accumulator = 0;
      //positive y neighbor
      if (voxelStates[(i, j + 1, k)] === FLUID) {
        accumulator = -scale;
      }
      return accumulator;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID_DENSITY: FLUID_DENSITY,
      AIR: STATE_ENUM.AIR,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);

export const createAZKernel = (gpu, nx, ny, nz, cellSize) =>
  gpu
    .createKernel(function (voxelStates, dt) {
      // for brevity
      const i = this.thread.x;
      const j = this.thread.y;
      const k = this.thread.z;
      const FLUID = this.constants.FLUID;

      // only consider fluid cells
      if (voxelStates[i][j][k] !== FLUID) {
        return 0;
      }

      const scale =
        dt /
        (this.constants.FLUID_DENSITY *
          this.constants.CELL_SIZE *
          this.constants.CELL_SIZE);

      let accumulator = 0;
      //positive z neighbor
      if (voxelStates[(i, j, k + 1)] === FLUID) {
        accumulator = -scale;
      }
      return accumulator;
    })
    .setTactic("precision") // vector math should be high precision
    .setConstants({
      CELL_SIZE: cellSize,
      FLUID_DENSITY: FLUID_DENSITY,
      AIR: STATE_ENUM.AIR,
      FLUID: STATE_ENUM.FLUID,
      SOLID: STATE_ENUM.SOLID,
    })
    .setOutput([nx, ny, nz]);
