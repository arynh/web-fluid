import { vec3 } from "gl-matrix";
import { initialize3DArray } from "./utils.js";

/** Enum for voxel states. */
export const STATE_ENUM = {
  AIR: 0,
  FLUID: 1,
  SOLID: 2,
};

/**
 * Represent a MAC grid and the quantities associated with it. See Bridson,
 * 2015, chapter 2 for details on the structure of the grid. Pressures are
 * stored at the center of each voxel, and normal velocities are stored at the
 * boundaries on the voxels.
 */
export class MACGrid {
  /**
   * Construct a new MAC grid with the given specification. This grid
   * construction is based off of Austin Eng's representation.
   *
   * @param {{min: vec3, max: vec3}} boundaries The bounds of the grid.
   * @param {number} cellSize The width of each voxel.
   */
  constructor(boundaries, cellSize) {
    this.min = boundaries.min;
    this.max = boundaries.max;

    // adjust the max extent to align to an integer number of cells
    vec3.set(
      this.max,
      this.min[0] +
        cellSize * Math.ceil((this.max[0] - this.min[0]) / cellSize),
      this.min[1] +
        cellSize * Math.ceil((this.max[1] - this.min[1]) / cellSize),
      this.min[2] + cellSize * Math.ceil((this.max[2] - this.min[2]) / cellSize)
    );

    this.cellSize = cellSize;
    this.count = vec3.create();
    this.size = vec3.create();
    vec3.sub(this.size, this.max, this.min);
    vec3.scale(this.count, this.size, 1.0 / this.cellSize);
    vec3.add(this.count, this.count, vec3.fromValues(1, 1, 1));
    vec3.floor(this.count, this.count);

    /// initialize
    this.nx = this.count[0] - 1;
    this.ny = this.count[1] - 1;
    this.nz = this.count[2] - 1;
    this.pressure = initialize3DArray(this.nx, this.ny, this.nz);
    this.velocity_x = initialize3DArray(this.nx + 1, this.ny, this.nz);
    this.velocity_y = initialize3DArray(this.nx, this.ny + 1, this.nz);
    this.velocity_z = initialize3DArray(this.nx, this.ny, this.nz + 1);

    // initialize voxel states
    this.voxelStates = initialize3DArray(this.nx, this.ny, this.nz);
  }

  /**
   * Set all voxels on the boundary of the grid to be solids. Leave everything
   * else untouched (by default, air).
   */
  addDefaultSolids() {
    for (let i = 0; i < this.nx; i++) {
      for (let j = 0; j < this.ny; j++) {
        for (let k = 0; k < this.nz; k++) {
          if (
            i === 0 ||
            j === 0 ||
            k === 0 ||
            i === this.nx - 1 ||
            j === this.ny - 1 ||
            k === this.nz - 1
          ) {
            this.voxelStates[i][j][k] = STATE_ENUM.SOLID;
          }
        }
      }
    }
  }
}
