import { ATTRIBUTE_COUNT } from "../particles.js";
//add kernel for subracting the two and sending it to the createFLIPKernel

export const createFLIPKernel = (gpu,particleCount) =>
  gpu
    .addFunction(function lerp(a, b, t) {
      return t * a + (1 - t) * b;
    })
    .createKernel(function (particles,diffGridVx,diffGridVy,diffGridVz) {
      //mod to figure out which index we are at (0-5) for each particle
      let index_mod = this.thread.x % this.constants.ATTRIBUTE_COUNT;
      //if we are looking at the position just return the position
      if (index_mod === 0 || index_mod === 1 || index_mod === 2) {
        return this.thread.x;
      }
      //get the positions - index changes depending which velocity we are looking at
      let pos_x = particles[this.thread.x-index_mod];
      let pos_y = particles[this.thread.x-index_mod+1];
      let pos_z = particles[this.thread.x-index_mod+2];
      //get the lower and upper grid positions
      let grid_lower_x = Math.floor(pos_x/cellSize);
      let grid_upper_x = Math.ceil(pos_x/cellSize);
      let grid_lower_y = Math.floor(pos_y/cellSize);
      let grid_upper_y = Math.ceil(pos_y/cellSize);
      let grid_lower_z = Math.floor(pos_z/cellSize);
      let grid_upper_z = Math.ceil(pos_z/cellSize);
      if (index_mod === 3) {
        //vx
        //get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_x-grid_lower_x)/cellSize;
        return lerp(
          diffGridVx[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVx[grid_upper_x][grid_upper_y][grid_upper_z],
          lerpWeight
        );
      } else if (index_mod === 4) {
        //vy
        //get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_y-grid_lower_y)/cellSize;
        return lerp(
          diffGridy[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVy[grid_upper_x][grid_upper_y][grid_upper_z],
          lerpWeight
        );
      } else if (index_mod === 5) {
        //vz
        //get the lerp weight and return the lerp'd velocity
        let lerpWeight = (pos_z-grid_lower_z)/cellSize;
        return lerp(
          diffGridVz[grid_lower_x][grid_lower_y][grid_lower_z],
          diffGridVz[grid_upper_x][grid_upper_y][grid_upper_z],
          lerpWeight
        );
      }
    })
    .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT,
                  cellSize: cellSize})
    .setOutput([particleCount*ATTRIBUTE_COUNT]);
