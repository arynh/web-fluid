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

/**
 * Produce the matrix-vector product `Ax` from the sparsely stored A and x.
 * FIXME: I'm not quite sure how to do this, this is just a guess.
 */
export const createApplyAKernel = (gpu, vectorLength) =>
  gpu
    .createKernel(function (Adiag, Ax, Ay, Az, x) {
      const i = this.thread.x;
      // FIXME: idk how to do this
    })
    .setTactic("precision") // vector math should be high precision
    .setOutput([vectorLength]);
