export const solve = (
  kernels,
  voxelStates,
  dt,
  velocityX,
  velocityY,
  velocityZ,
  tolerance,
  iterationLimit
) => {
  const start = Date.now();

  // build coefficient matrix
  const Adiag = kernels.buildADiag(voxelStates, dt);
  const Ax = kernels.buildAX(voxelStates, dt);
  const Ay = kernels.buildAY(voxelStates, dt);
  const Az = kernels.buildAZ(voxelStates, dt);
  // build negative divergence vector
  const d = kernels.flatten(
    kernels.buildD(voxelStates, velocityX, velocityY, velocityZ)
  );
  // console.log("d:");
  // console.log(d.toArray());
  // console.log("error:");
  // console.log(error(d.toArray()));

  // follow PCG algorithm set out in Bridson
  let p = kernels.zeroVector();
  let r = d;

  // check if there is a trivial solution
  if (checkResidual(r.toArray(), tolerance)) {
    const end = Date.now();
    console.log(`Solver took 0 iterations in ${end - start} ms.`);
    const _p = p.toArray();
    free([p, r]);
    return _p;
  }

  let _r = [];
  let z = r; // applyPreconditioner(r); // TODO: implement the preconditioner.
  let s = z;
  let sigma = kernels.math.dot(z, r);

  let iterationCount = 0;
  while (iterationCount++ < iterationLimit) {
    // z <- As
    z = kernels.math.applyA(Adiag, Ax, Ay, Az, s, voxelStates);
    // console.log("A:");
    // console.log(Adiag.toArray());
    // console.log(Ax.toArray());
    // console.log(Ay.toArray());
    // console.log(Az.toArray());
    // console.log("S:");
    // console.log(s.toArray());
    // console.log("z <- As:");
    // console.log(z.toArray());

    // alpha <- sigma / (z dot s)
    let alpha = sigma / kernels.math.dot(z, s);

    // p <- p + alpha * s
    p = kernels.math.componentWiseAdd(p, kernels.math.scalarMultiply(s, alpha));

    // r <- r - alpha * z
    r = kernels.math.componentWiseAdd(
      r,
      kernels.math.scalarMultiply(z, -alpha)
    );

    // transfer r to CPU to do error calculation
    _r = r.toArray();

    // if (iterationCount % 50 === 0) {
    //   console.log(`error at iteration ${iterationCount}: ${error(_r)}`);
    // }

    // if the residual is sufficiently small, return early
    if (checkResidual(_r, tolerance)) {
      const end = Date.now();
      // console.log(
      //   `Solver took ${iterationCount - 1} iterations in ${end - start} ms.`
      // );
      const _p = p.toArray();
      free([p, z, s, r]);
      return _p;
    }
    z = r; // applyPreconditioner(r); // TODO: implement the preconditioner

    // sigma' <- z dot r
    let _sigma = kernels.math.dot(z, r);

    // beta <- sigma' / sigma
    let beta = _sigma / sigma;

    // s <- z + beta * s
    s = kernels.math.componentWiseAdd(z, kernels.math.scalarMultiply(s, beta));

    // sigma <- sigma'
    sigma = _sigma;
  }

  const _p = p.toArray();
  free([p, z, s, r]);

  // console.error("Maximum iterations used in PCG solver!");

  const end = Date.now();
  // console.log(
  //   `Solver took ${iterationCount - 1} iterations in ${end - start} ms.`
  // );

  return _p;
};

const checkResidual = (r, tolerance) => error(r) <= tolerance;
const error = (r) => r.reduce((max, n) => Math.max(max, Math.abs(n)), -1);
const free = (textures) => textures.map((t) => t.delete());
