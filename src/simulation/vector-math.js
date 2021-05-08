/**
 * Take the product of a square matrix A and vector x.
 */
const createMatrixVectorProductKernel = (gpu, vectorLength) =>
  gpu
    .setTactic("precision") // vector math should be high precision
    .createKernel(function (A, x) {
      let sum = 0;
      for (let i = 0; i < this.constants.vectorLength; i++) {
        sum += A[this.thread.x][i] * x[i];
      }
      return sum;
    })
    .setConstants({ vectorLength: vectorLength })
    .setOutput([vectorLength]);
