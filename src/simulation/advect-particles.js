import { ATTRIBUTE_COUNT } from "../particles.js";

export const createAdvectParticlesKernel = (gpu, particleCount) =>
  gpu
    .addFunction(function lerp(a, b, t) {
      return t * a + (1 - t) * b;
    })
    .createKernel(function (
      particles,
      dt,
      cellSize,
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
        let lerpWeight = (xIntermediate - Math.floor(xIntermediate)) / cellSize;
        let vxIntermediate = lerp(
          velocityFieldX[Math.floor(xIntermediate / cellSize)][
            Math.floor(y / cellSize)
          ][Math.floor(z / cellSize)],
          velocityFieldX[Math.ceil(xIntermediate / cellSize)][
            Math.floor(y / cellSize)
          ][Math.floor(z / cellSize)],
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
        let lerpWeight = (yIntermediate - Math.floor(yIntermediate)) / cellSize;
        let vyIntermediate = lerp(
          velocityFieldY[Math.floor(x / cellSize)][
            Math.floor(yIntermediate / cellSize)
          ][Math.floor(z / cellSize)],
          velocityFieldY[Math.floor(x / cellSize)][
            Math.ceil(yIntermediate / cellSize)
          ][Math.floor(z / cellSize)],
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
        let lerpWeight = (zIntermediate - Math.floor(zIntermediate)) / cellSize;
        let vzIntermediate = lerp(
          velocityFieldZ[Math.floor(x / cellSize)][Math.floor(y / cellSize)][
            Math.floor(zIntermediate / cellSize)
          ],
          velocityFieldZ[Math.floor(x / cellSize)][Math.floor(y / cellSize)][
            Math.ceil(zIntermediate / cellSize)
          ],
          lerpWeight
        );
        let k2 = dt * vzIntermediate;
        return z + k2;
      } else {
        // don't change the velocities
        return particles[this.thread.x];
      }
    })
    .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT })
    .setOutput([ATTRIBUTE_COUNT * particleCount]);
