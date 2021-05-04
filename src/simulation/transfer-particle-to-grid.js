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

/** Calulate weights
 * In:  particle buffer
 * Out: weight per particle (dimension: [# particles])
 */
const createParticleWeightsKernel = (gpu, particleCount, edgeCount) =>
  gpu
    .createKernel(function (particleIndices, particles, cellSize, dimension) {
      let position = particles[particleIndices[this.thread.x] + dimension];
      let velocity = particles[particleIndices[this.thread.x] + 3 + dimension];

      // get the indices of the grid boundries that lie around this particle
      let scaledPosition = position / cellSize;
      let lowerEdgeIndex = Math.floor(scaledPosition);
      let upperEdgeIndex = Math.ceil(scaledPosition);

      // TODO: finish this one
      // FIXME: please fix my brain
    })
    .setOutput([particleCount, 2, edgeCount]);

const createWeightedSumKernel = (gpu, particleCount, edgeCount) =>
  gpu
    .createKernel(
      function (weightedSumValues) {
        let sum = 0;
        for (let i = 0; i < particleCount; i++) {
          sum += weightedSumValues[i][this.thread.x][this.thread.y];
        }
        return sum;
      },
      { constants: { particleCount: particleCount } }
    )
    .setOutput([2, edgeCount]);

const createNewVelocitiesKernel = (gpu, edgeCount) =>
  gpu
    .createKernel(function (weightedSumValues) {
      let denominator = weightedSumValues[1][this.thread.x];
      if (Math.abs(denominator) < 0.0001) {
        return 0;
      }
      return weightedSumValues[0][this.thread.x] / denominator;
    })
    .setOutput([edgeCount]);
