/**
 * Assumption!
 * This assumes that we're using the default arrangement of having all solid
 * voxels around the edge of our rectangular domain.
 */

export const createEnforceBoundaryXKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.x === 0 || this.thread.x === this.constants.nx - 1) {
        return 0;
      }
      return velocities[this.thread.z][this.thread.y][this.thread.x];
    })
    .setConstants({ nx: nx })
    .setOutput([nx, ny, nz]);

export const createEnforceBoundaryYKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.y === 0 || this.thread.y === this.constants.ny - 1) {
        return 0;
      }
      return velocities[this.thread.z][this.thread.y][this.thread.x];
    })
    .setConstants({ ny: ny })
    .setOutput([nx, ny, nz]);

export const createEnforceBoundaryZKernel = (gpu, nx, ny, nz) =>
  gpu
    .createKernel(function (velocities) {
      if (this.thread.z === 0 || this.thread.z === this.constants.nz - 1) {
        return 0;
      }
      return velocities[this.thread.z][this.thread.y][this.thread.x];
    })
    .setConstants({ nz: nz })
    .setOutput([nx, ny, nz]);
