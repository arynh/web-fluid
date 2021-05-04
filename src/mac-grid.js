import { vec3 } from "gl-matrix";
import { initialize3DArray } from "./utils.js";

/**
 * Represent a MAC grid and the quantities associated with it. See Bridson,
 * 2015, chapter 2 for details on the structure of the grid. Pressures are
 * stored at the center of each voxel, and normal velocities are stored at the
 * boundaries on the voxels.
 */
class MACGrid {
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
    vec3.scale(this.count, boundaries._size, 1.0 / this.cellSize);
    vec3.add(this.count, this.count, vec3.fromValues(1, 1, 1));
    vec3.floor(this.count, this.count);

    /// initialize
    let nx = this.count.x;
    let ny = this.count.y;
    let nz = this.count.z;
    this.pressure = initialize3DArray(nx, ny, nz);
    this.velocity_x = initialize3DArray(nx + 1, ny, nz);
    this.velocity_y = initialize3DArray(nx, ny + 1, nz);
    this.velocity_z = initialize3DArray(nx, ny, nz + 1);
  }
}
