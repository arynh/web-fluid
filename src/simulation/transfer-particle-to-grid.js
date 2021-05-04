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
export const createParticleWeightsKernel = (gpu, particleCount, edgeCount) =>
  gpu
    .createKernel(function (particles, cellSize, dimension) {
      let position = particles[this.thread.x * 6 + dimension];
      let velocity = particles[this.thread.x * 6 + 3 + dimension];

      // get the indices of the grid boundries that lie around this particle
      let scaledPosition = position / cellSize;
      let lowerEdgeIndex = Math.floor(scaledPosition);
      let upperEdgeIndex = Math.ceil(scaledPosition);
      // if the current edge isn't relevant, return 0
      if (this.thread.z != lowerEdgeIndex && this.thread.z != upperEdgeIndex) {
        return 0;
      }

      // calculate the weight based on the triangle weighting function
      let weight = 0;
      if (this.thread.z == lowerEdgeIndex) {
        weight = 1 - (position - lowerEdgeIndex * cellSize) / cellSize;
      } else {
        weight = 1 + (position - upperEdgeIndex * cellSize) / cellSize;
      }

      if (this.thread.y == 0) {
        return velocity * weight;
      } else {
        return weight;
      }
    })
    .setOutput([particleCount, 2, edgeCount]);

export const createWeightedSumKernel = (gpu, particleCount, edgeCount) =>
  gpu
    .createKernel(function (weightedSumValues) {
      let sum = 0;
      for (let i = 0; i < this.constants.particleCount; i++) {
        sum += weightedSumValues[i][this.thread.x][this.thread.y];
      }
      return sum;
    })
    .setConstants({ particleCount: particleCount })
    .setOutput([2, edgeCount]);

export const createNewVelocitiesKernel = (gpu, edgeCount) =>
  gpu
    .createKernel(function (weightedSumValues) {
      let denominator = weightedSumValues[1][this.thread.x];
      if (Math.abs(denominator) < 0.0001) {
        return 0;
      }
      return weightedSumValues[0][this.thread.x] / denominator;
    })
    .setOutput([edgeCount]);

export const createParticleToGridKernel = (gpu, particleCount, edgeCount) => {
  const particleWeightsKernel = createParticleWeightsKernel(
    gpu,
    particleCount,
    edgeCount
  );
  const weightedSumKernel = createWeightedSumKernel(
    gpu,
    particleCount,
    edgeCount
  );
  const newVelocitiesKernel = createNewVelocitiesKernel(gpu, edgeCount);

  // pipeline these three kernels
  return gpu.combineKernels(
    particleWeightsKernel,
    weightedSumKernel,
    newVelocitiesKernel,
    function (particles, cellSize, dimension) {
      return newVelocitiesKernel(
        weightedSumKernel(particleWeightsKernel(particles, cellSize, dimension))
      );
    }
  );
};
