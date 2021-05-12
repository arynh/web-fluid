/**
 * Component-wise add two vectors.
 */
export const createComponentWiseAddKernel = (gpu, vectorLength) =>
  gpu
    .createKernel(function (a, b) {
      return a[this.thread.x] + b[this.thread.x];
    })
    .setTactic("precision") // vector math should be high precision
    .setOutput([vectorLength]);

/**
 * Component-wise multiply two vectors.
 */
export const createComponentWiseMultiplyKernel = (gpu, vectorLength) =>
  gpu
    .createKernel(function (a, b) {
      return a[this.thread.x] * b[this.thread.x];
    })
    .setTactic("precision") // vector math should be high precision
    .setOutput([vectorLength]);

/**
 * Muliply a vector `a` by a scalar.
 */
export const createScalarMultiplyKernel = (gpu, vectorLength) =>
  gpu
    .createKernel(function (a, scalar) {
      return scalar * a[this.thread.x];
    })
    .setTactic("precision") // vector math should be high precision
    .setOutput([vectorLength]);
