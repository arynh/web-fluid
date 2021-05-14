export const solve = (
  kernels,
  voxelStates,
  velocityX,
  velocityY,
  velocityZ,
  iterationLimit,
  pressureOld
) => {
  let p = pressureOld;

  const d = kernels.buildD(voxelStates, velocityX, velocityY, velocityZ);

  // Jacobi Iteration
  for (let i = 0; i < iterationLimit; i++) {
    p = kernels.jacobi(d, p, voxelStates);
  }

  return p;
};
