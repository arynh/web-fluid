/**
 * Transfer particle velocities to the grid. See Bridson, section 7.6 on
 * particle methods.
 */

/**
 * How to parallelize?
 *
 *
 * for each particle:
 *  - get position of particle
 *  - find grid edges within dx
 *  for each of these edges:
 *      - add velocity * weight
 *      - where weight is given as triangle function
 */
export const createParticleToGridKernel = (
  gpu,
  particleCount,
  nx,
  ny,
  nz,
  dimension
) =>
  gpu
    .createKernel(function (particles, cellSize) {
      // get spatial location of grid velocity vector
      let x = cellSize * this.thread.x;
      let y = cellSize * this.thread.y;
      let z = cellSize * this.thread.z;

      // declare numerator and denominator of the weighted sum
      let numerator = 0;
      let denominator = 0;
      /* loop through particles to find ones that are close, add their
      velocity contribution to the grid velocity */
      for (
        let particleIndex = 0;
        particleIndex < this.constants.PARTICLE_COUNT;
        particleIndex++
      ) {
        // calculate distance in each dimension
        let distance_x = particles[particleIndex * 6] - x;
        let distance_y = particles[particleIndex * 6 + 1] - y;
        let distance_z = particles[particleIndex * 6 + 2] - z;

        // if it's far, skip it
        if (
          Math.abs(distance_x) > cellSize ||
          Math.abs(distance_y) > cellSize ||
          Math.abs(distance_z) > cellSize
        ) {
          continue;
        }

        // calculate the weight according to the trilinear interpolation
        let weight =
          triangle(distance_x / cellSize) *
          triangle(distance_y / cellSize) *
          triangle(distance_z / cellSize);

        numerator +=
          particles[particleIndex * 6 + 3 + this.constants.DIMENSION] * weight;
        denominator += weight;
      }

      // check for divide by zero
      if (Math.abs(denominator) < 0.0001) {
        return 0;
      }
      return numerator / denominator;
    })
    .addFunction(function triangle(r) {
      let r_magnitude = Math.abs(r);
      if (r_magnitude >= 1) {
        return 0;
      }
      return 1 - r_magnitude;
    })
    .setConstants({ PARTICLE_COUNT: particleCount, DIMENSION: dimension })
    .setOutput([nx, ny, nz]);
