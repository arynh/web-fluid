import { STATE_ENUM } from "../mac-grid.js";

export const createVelocityXUpdateKernel = (
  gpu,
  nx,
  ny,
  nz,
  fluidDensity,
  cellSize
) =>
  gpu
    .createKernel(function (velocity, pressure, voxelStates, dt) {
      // only consider boundaries which have a fluid cell on at least one side
      if (
        this.thread.x === 0 ||
        this.thread.x === this.constants.NX - 1 ||
        !(
          voxelStates[this.thread.z][this.thread.y][this.thread.x - 1] ===
            this.constants.FLUID ||
          voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
            this.constants.FLUID
        )
      ) {
        return 0;
      }

      const pressureGradient =
        (pressure[this.thread.z][this.thread.y][this.thread.x] -
          pressure[this.thread.z][this.thread.y][this.thread.x - 1]) /
        this.constants.CELL_SIZE;

      return (
        velocity[this.thread.z][this.thread.y][this.thread.x] -
        (dt * pressureGradient) / this.constants.FLUID_DENSITY
      );
    })
    .setConstants({
      FLUID_DENSITY: fluidDensity,
      CELL_SIZE: cellSize,
      NX: nx,
      NY: ny,
      NZ: nz,
      FLUID: STATE_ENUM.FLUID,
    })
    .setOutput([nx, ny, nz]);

export const createVelocityYUpdateKernel = (
  gpu,
  nx,
  ny,
  nz,
  fluidDensity,
  cellSize
) =>
  gpu
    .createKernel(function (velocity, pressure, voxelStates, dt) {
      // only consider boundaries which have a fluid cell on at least one side
      if (
        this.thread.y === 0 ||
        this.thread.y === this.constants.NY - 1 ||
        !(
          voxelStates[this.thread.z][this.thread.y - 1][this.thread.x] ===
            this.constants.FLUID ||
          voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
            this.constants.FLUID
        )
      ) {
        return 0;
      }

      const pressureGradient =
        (pressure[this.thread.z][this.thread.y][this.thread.x] -
          pressure[this.thread.z][this.thread.y - 1][this.thread.x]) /
        this.constants.CELL_SIZE;

      const oldVelocity = velocity[this.thread.z][this.thread.y][this.thread.x];
      const newVelocity =
        oldVelocity - (dt * pressureGradient) / this.constants.FLUID_DENSITY;

      // if (newVelocity > 0) {
      //   debugger;
      // }

      return newVelocity;
    })
    .setConstants({
      FLUID_DENSITY: fluidDensity,
      CELL_SIZE: cellSize,
      NX: nx,
      NY: ny,
      NZ: nz,
      FLUID: STATE_ENUM.FLUID,
    })
    .setOutput([nx, ny, nz]);

export const createVelocityZUpdateKernel = (
  gpu,
  nx,
  ny,
  nz,
  fluidDensity,
  cellSize
) =>
  gpu
    .createKernel(function (velocity, pressure, voxelStates, dt) {
      // only consider boundaries which have a fluid cell on at least one side
      if (
        this.thread.z === 0 ||
        this.thread.z === this.constants.NZ - 1 ||
        !(
          voxelStates[this.thread.z - 1][this.thread.y][this.thread.x] ===
            this.constants.FLUID ||
          voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
            this.constants.FLUID
        )
      ) {
        return 0;
      }

      const pressureGradient =
        (pressure[this.thread.z][this.thread.y][this.thread.x] -
          pressure[this.thread.z - 1][this.thread.y][this.thread.x]) /
        this.constants.CELL_SIZE;

      return (
        velocity[this.thread.z][this.thread.y][this.thread.x] -
        (dt * pressureGradient) / this.constants.FLUID_DENSITY
      );
    })
    .setConstants({
      FLUID_DENSITY: fluidDensity,
      CELL_SIZE: cellSize,
      NX: nx,
      NY: ny,
      NZ: nz,
      FLUID: STATE_ENUM.FLUID,
    })
    .setOutput([nx, ny, nz]);
