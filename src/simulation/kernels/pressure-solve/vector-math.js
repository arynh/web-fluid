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
 */
export const createApplyAKernel = (gpu, vectorLength) =>
  gpu
    .createKernel(function (Adiag, Ax, Ay, Az, x) {
      const aux = this.thread.x % (this.constants.NX * this.constants.NY);
      const i = Math.floor(aux / this.constants.NY);
      const j = aux % this.constants.NX;
      const k = Math.floor(
        this.thread.x / (this.constants.NX * this.constants.NY)
      );

      if (this.thread.x === 0) {
      } else if (this.thread.x === 1) {
      } else if (this.thread.x === 2) {
      } else if (this.thread.x === 3) {
      } else if (this.thread.x === 4) {
      } else if (this.thread.x === 5) {
      }
    })
    .addFunction(function () {})
    .setTactic("precision") // vector math should be high precision
    .setConstants({ VECTOR_LENGTH: vectorLength })
    .setOutput([vectorLength]);
