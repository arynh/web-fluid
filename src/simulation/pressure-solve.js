const solve = (
  kernels,
  voxelStates,
  dt,
  velocityX,
  velocityY,
  velocityZ,
  tolerance,
  iterationLimit
) => {
  // build coefficient matrix
  const Adiag = kernels.buildADiag(voxelStates, dt);
  const Ax = kernels.buildAX(voxelStates, dt);
  const Ay = kernels.buildAY(voxelStates, dt);
  const Az = kernels.buildAZ(voxelStates, dt);
  // build negative divergence vector
  // TODO: d needs to be flattened
  const d = kernels.buildD(voxelStates, velocityX, velocityY, velocityZ);

  // follow PCG algorithm set out in Bridson
  let p = kernels.zeroVector();
  let r = d;
  let z = applyPreconditioner(r); // TODO: implement the preconditioner.
  let s = z;
  let sigma = kernels.math.dot(z, r);

  let iterationCount = 0;
  while (!checkResidual(r, tolerance) && iterationCount < iterationLimit) {
    // z <- As
    z = kernels.math.applyA(Adiag, Ax, Ay, Az, s);

    // alpha <- sigma / (z dot s)
    let alpha = sigma / kernels.math.dot(z, s);

    // p <- p + alpha * s
    p = kernels.math.componentWiseAdd(p, kernels.math.scalarMultiply(s, alpha));

    // r <- r - alpha * z
    r = kernels.math.componentWiseAdd(
      r,
      kernels.math.scalarMultiply(z, -alpha)
    );
    console.log(`error at iteration ${iterationCount}: ${error(r)}`);

    if (checkResidual(r, tolerance)) {
      return p;
    }

    z = applyPreconditioner(r);

    // sigma' <- z dot r
    let _sigma = kernels.math.dot(z, r);

    // beta <- sigma' / sigma
    let beta = _sigma / sigma;

    // s <- z + beta * s
    s = kernels.math.componentWiseAdd(z, kernels.math.scalarMultiply(s, beta));

    // sigma <- sigma'
    sigma = _sigma;

    iterationCount++;
  }

  console.error("Maximum iterations used in PCG solver!");

  return p;
};

const checkResidual = (r, tolerance) => {
  return error(r) <= tolerance;
};

const error = (r) => r.reduce((max, n) => Math.max(max, Math.abs(n)), -1);
