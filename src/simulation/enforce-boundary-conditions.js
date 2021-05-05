/**
 * Assumption!
 * This assumes that we're using the default arrangement of having all solid
 * voxels around the edge of our rectangular domain.
 */

export const createEnforceBoundaryXKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.x === 1 || this.thread.x === this.constants.nx - 2) {
        return 0;
      }
      return velocities[this.thread.x][this.thread.y][this.thread.z];
    })
    .setConstants({ nx: nx })
    .setOutput([nx, ny, nz]);

export const createEnforceBoundaryYKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.y === 1 || this.thread.y === this.constants.ny - 2) {
        return 0;
      }
      return velocities[this.thread.x][this.thread.y][this.thread.z];
    })
    .setConstants({ ny: ny })
    .setOutput([nx, ny, nz]);

export const createEnforceBoundaryZKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.z === 1 || this.thread.z === this.constants.nz - 2) {
        return 0;
      }
      return velocities[this.thread.x][this.thread.y][this.thread.z];
    })
    .setConstants({ nz: nz })
    .setOutput([nx, ny, nz]);
