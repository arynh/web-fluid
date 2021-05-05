import {STATE_ENUM} from "../mac-grid.js"
//export makes it public so you can import into index.js
export const createClassifyVoxelsKernel = (gpu,nx,ny,nz) =>
gpu
  .createKernel(function (voxelStates,particles,cellSize) {
    //check if there is a particle in that grid
    // get spatial location of grid velocity vector
    let x = this.thread.x * cellSize;
    let y = this.thread.y * cellSize;
    let z = this.thread.z * cellSize;

    var particle_exists = false;
    for (let i = 0; i < particles.length; i++) {
      let pos_x = particles[i*6];
      let pos_y = particles[i*6+1];
      let pos_z = particles[i*6+2];
      if (Math.abs(pos_x-x) < cellSize && Math.abs(pos_y-y) < cellSize && Math.abs(pos_z-z) < cellSize) {
        particle_exists = true;
        break;
      }
    }
    //use cellsize and the xyz to figure out your grid cell dims and use particles pos to compare
    //if there is and the state is air flip it to fluid
    if (particle_exists) {
      if (voxelStates[this.thread.x][this.thread.y][this.thread.z] !== STATE_ENUM.SOLID){
        return STATE_ENUM.FLUID;
      }
    } else {
      //if there isn't and the state is flud flip it to air
      if (voxelStates[this.thread.x][this.thread.y][this.thread.z] !== STATE_ENUM.SOLID){
        return STATE_ENUM.AIR;
      }
    }
  })
  .setOutput([nx, ny, nz]);
