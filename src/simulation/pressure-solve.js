export const solve = (
  kernels,
  voxelStates,
  dt,
  velocityX,
  velocityY,
  velocityZ,
  tolerance,
  iterationLimit,
  pressure,
  pressureOld
) => {
  let p = pressureOld;

  const d = kernels.buildD(voxelStates, velocityX, velocityY, velocityZ);

  // JACOBI ITERATION
  for (let i = 0; i < iterationLimit; i++) {
    p = kernels.jacobi(d, p, voxelStates);
  }

  return p;
};
