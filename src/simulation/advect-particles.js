import { ATTRIBUTE_COUNT } from "../particles.js";

export const createAdvectParticlesKernel = (gpu, particleCount, cellSize) =>
  gpu
    .createKernel(function (
      particles,
      dt,
      velocityFieldX,
      velocityFieldY,
      velocityFieldZ
    ) {
      // check which position component we're looking at
      if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 0) {
        // get position
        let x = particles[this.thread.x];
        let y = particles[this.thread.x + 1];
        let z = particles[this.thread.x + 2];

        // get x velocity
        let vx = particles[this.thread.x + 3];

        // carry out 2nd order Runge-Kutta solver in one dimension
        let k1 = dt * vx;
        let xIntermediate = x + k1 / 2;

        // interpolate the velocity at the intermediate x value
        let lerpWeight =
          (xIntermediate -
            Math.floor(xIntermediate) * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        let vxIntermediate = lerp(
          velocityFieldX[Math.floor(xIntermediate / this.constants.CELL_SIZE)][
            Math.floor(y / this.constants.CELL_SIZE)
          ][Math.floor(z / this.constants.CELL_SIZE)],
          velocityFieldX[Math.ceil(xIntermediate / this.constants.CELL_SIZE)][
            Math.floor(y / this.constants.CELL_SIZE)
          ][Math.floor(z / this.constants.CELL_SIZE)],
          lerpWeight
        );
        let k2 = dt * vxIntermediate;
        return x + k2;
      } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 1) {
        // get position
        let x = particles[this.thread.x - 1];
        let y = particles[this.thread.x];
        let z = particles[this.thread.x + 1];

        // get y velocity
        let vy = particles[this.thread.x + 3];

        // carry out 2nd order Runge-Kutta solver in one dimension
        let k1 = dt * vy;
        let yIntermediate = y + k1 / 2;

        // interpolate the velocity at the intermediate y value
        let lerpWeight =
          (yIntermediate -
            Math.floor(yIntermediate) * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        let vyIntermediate = lerp(
          velocityFieldY[Math.floor(x / this.constants.CELL_SIZE)][
            Math.floor(yIntermediate / this.constants.CELL_SIZE)
          ][Math.floor(z / this.constants.CELL_SIZE)],
          velocityFieldY[Math.floor(x / this.constants.CELL_SIZE)][
            Math.ceil(yIntermediate / this.constants.CELL_SIZE)
          ][Math.floor(z / this.constants.CELL_SIZE)],
          lerpWeight
        );
        let k2 = dt * vyIntermediate;
        return y + k2;
      } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 2) {
        // get position
        let x = particles[this.thread.x - 2];
        let y = particles[this.thread.x - 1];
        let z = particles[this.thread.x];

        // get z velocity
        let vz = particles[this.thread.x + 3];

        // carry out 2nd order Runge-Kutta solver in one dimension
        let k1 = dt * vz;
        let zIntermediate = z + k1 / 2;

        // interpolate the velocity at the intermediate z value
        let lerpWeight =
          (zIntermediate -
            Math.floor(zIntermediate) * this.constants.CELL_SIZE) /
          this.constants.CELL_SIZE;
        let vzIntermediate = lerp(
          velocityFieldZ[Math.floor(x / this.constants.CELL_SIZE)][
            Math.floor(y / this.constants.CELL_SIZE)
          ][Math.floor(zIntermediate / this.constants.CELL_SIZE)],
          velocityFieldZ[Math.floor(x / this.constants.CELL_SIZE)][
            Math.floor(y / this.constants.CELL_SIZE)
          ][Math.ceil(zIntermediate / this.constants.CELL_SIZE)],
          lerpWeight
        );
        let k2 = dt * vzIntermediate;
        return z + k2;
      } else {
        // don't change the velocities
        return particles[this.thread.x];
      }
    })
    .addFunction(function lerp(a, b, t) {
      return (1 - t) * a + t * b;
    })
    .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT, CELL_SIZE: cellSize })
    .setOutput([ATTRIBUTE_COUNT * particleCount]);
