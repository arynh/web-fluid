(function () {
  'use strict';

  /**
   * Common utilities
   * @module glMatrix
   */
  var ARRAY_TYPE = typeof Float32Array !== 'undefined' ? Float32Array : Array;
  if (!Math.hypot) Math.hypot = function () {
    var y = 0,
        i = arguments.length;

    while (i--) {
      y += arguments[i] * arguments[i];
    }

    return Math.sqrt(y);
  };

  /**
   * 3 Dimensional Vector
   * @module vec3
   */

  /**
   * Creates a new, empty vec3
   *
   * @returns {vec3} a new 3D vector
   */

  function create() {
    var out = new ARRAY_TYPE(3);

    if (ARRAY_TYPE != Float32Array) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
    }

    return out;
  }
  /**
   * Creates a new vec3 initialized with the given values
   *
   * @param {Number} x X component
   * @param {Number} y Y component
   * @param {Number} z Z component
   * @returns {vec3} a new 3D vector
   */

  function fromValues(x, y, z) {
    var out = new ARRAY_TYPE(3);
    out[0] = x;
    out[1] = y;
    out[2] = z;
    return out;
  }
  /**
   * Set the components of a vec3 to the given values
   *
   * @param {vec3} out the receiving vector
   * @param {Number} x X component
   * @param {Number} y Y component
   * @param {Number} z Z component
   * @returns {vec3} out
   */

  function set(out, x, y, z) {
    out[0] = x;
    out[1] = y;
    out[2] = z;
    return out;
  }
  /**
   * Adds two vec3's
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a the first operand
   * @param {ReadonlyVec3} b the second operand
   * @returns {vec3} out
   */

  function add(out, a, b) {
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
    return out;
  }
  /**
   * Subtracts vector b from vector a
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a the first operand
   * @param {ReadonlyVec3} b the second operand
   * @returns {vec3} out
   */

  function subtract(out, a, b) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
    return out;
  }
  /**
   * Math.floor the components of a vec3
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a vector to floor
   * @returns {vec3} out
   */

  function floor(out, a) {
    out[0] = Math.floor(a[0]);
    out[1] = Math.floor(a[1]);
    out[2] = Math.floor(a[2]);
    return out;
  }
  /**
   * Scales a vec3 by a scalar number
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a the vector to scale
   * @param {Number} b amount to scale the vector by
   * @returns {vec3} out
   */

  function scale(out, a, b) {
    out[0] = a[0] * b;
    out[1] = a[1] * b;
    out[2] = a[2] * b;
    return out;
  }
  /**
   * Alias for {@link vec3.subtract}
   * @function
   */

  var sub = subtract;
  /**
   * Perform some operation over an array of vec3s.
   *
   * @param {Array} a the array of vectors to iterate over
   * @param {Number} stride Number of elements between the start of each vec3. If 0 assumes tightly packed
   * @param {Number} offset Number of elements to skip at the beginning of the array
   * @param {Number} count Number of vec3s to iterate over. If 0 iterates over entire array
   * @param {Function} fn Function to call for each vector in the array
   * @param {Object} [arg] additional argument to pass to fn
   * @returns {Array} a
   * @function
   */

  (function () {
    var vec = create();
    return function (a, stride, offset, count, fn, arg) {
      var i, l;

      if (!stride) {
        stride = 3;
      }

      if (!offset) {
        offset = 0;
      }

      if (count) {
        l = Math.min(count * stride + offset, a.length);
      } else {
        l = a.length;
      }

      for (i = offset; i < l; i += stride) {
        vec[0] = a[i];
        vec[1] = a[i + 1];
        vec[2] = a[i + 2];
        fn(vec, vec, arg);
        a[i] = vec[0];
        a[i + 1] = vec[1];
        a[i + 2] = vec[2];
      }

      return a;
    };
  })();

  /**
   * Create a new 3D array of zeros with the given dimensions.
   *
   * @param {number} x x dimension length
   * @param {number} y y dimension length
   * @param {number} z z dimension length
   * @returns The new array, all values set to zero.
   */
  const initialize3DArray = (x, y, z) => {
    let a = [];
    for (let i = 0; i < x; i++) {
      a.push([]);
      for (let j = 0; j < y; j++) {
        a[i].push(new Float32Array(z));
      }
    }
    a.toArray = function() {return this;};
    return a;
  };

  /** Enum for voxel states. */
  const STATE_ENUM = {
    AIR: 0,
    FLUID: 1,
    SOLID: 2,
  };

  /**
   * Represent a MAC grid and the quantities associated with it. See Bridson,
   * 2015, chapter 2 for details on the structure of the grid. Pressures are
   * stored at the center of each voxel, and normal velocities are stored at the
   * boundaries on the voxels.
   */
  class MACGrid {
    /**
     * Construct a new MAC grid with the given specification. This grid
     * construction is based off of Austin Eng's representation.
     *
     * @param {{min: vec3, max: vec3}} boundaries The bounds of the grid.
     * @param {number} cellSize The width of each voxel.
     */
    constructor(boundaries, cellSize) {
      this.min = boundaries.min;
      this.max = boundaries.max;

      // adjust the max extent to align to an integer number of cells
      set(
        this.max,
        this.min[0] +
          cellSize * Math.ceil((this.max[0] - this.min[0]) / cellSize),
        this.min[1] +
          cellSize * Math.ceil((this.max[1] - this.min[1]) / cellSize),
        this.min[2] + cellSize * Math.ceil((this.max[2] - this.min[2]) / cellSize)
      );

      this.cellSize = cellSize;
      this.count = create();
      this.size = create();
      sub(this.size, this.max, this.min);
      scale(this.count, this.size, 1.0 / this.cellSize);
      add(this.count, this.count, fromValues(1, 1, 1));
      floor(this.count, this.count);

      /// initialize
      this.nx = this.count[0] - 1;
      this.ny = this.count[1] - 1;
      this.nz = this.count[2] - 1;
      this.pressure = initialize3DArray(this.nx, this.ny, this.nz);
      this.pressureOld = null;
      this.velocityX = [[["hi"]]]; //initialize3DArray(this.nx + 1, this.ny, this.nz);
      this.velocityXOld = null;
      this.velocityY = null; //initialize3DArray(this.nx, this.ny + 1, this.nz);
      this.velocityYOld = null;
      this.velocityZ = null; //initialize3DArray(this.nx, this.ny, this.nz + 1);
      this.velocityZOld = null;

      // initialize voxel states
      this.voxelStates = initialize3DArray(this.nx, this.ny, this.nz);

      console.log(
        `Created a MAC grid with dimensions (${this.nx}, ${this.ny}, ${this.nz}).`
      );
    }

    /**
     * Set all voxels on the boundary of the grid to be solids. Leave everything
     * else untouched (by default, air).
     */
    addDefaultSolids() {
      for (let i = 0; i < this.nx; i++) {
        for (let j = 0; j < this.ny; j++) {
          for (let k = 0; k < this.nz; k++) {
            if (
              i === 0 ||
              j === 0 ||
              k === 0 ||
              i === this.nx - 1 ||
              j === this.ny - 1 ||
              k === this.nz - 1
            ) {
              this.voxelStates[i][j][k] = STATE_ENUM.SOLID;
            }
          }
        }
      }
    }
  }

  const ATTRIBUTE_COUNT = 6;

  /**
   * Represent the particle cloud.
   */
  class Particles {
    /**
     * Create a set of particles. The density of the particles determines how
     * many are made, and the bounds determine the initial position of the
     * particles. They start out evenly distributed throughout this box.
     *
     * @param {number} density
     * @param {{min: vec3, max: vec3}} bounds The initial minimum and maximum
     * extent of the box of particles.
     */
    constructor(density, bounds) {
      this.particleBuffer = [];
      this.particleIndices = [];
      let particle_counter = 0;
      let gap_between = 1 / Math.cbrt(density);
      for (let x = bounds.min[0]; x < bounds.max[0]; x += gap_between) {
        for (let y = bounds.min[1]; y < bounds.max[1]; y += gap_between) {
          for (let z = bounds.min[2]; z < bounds.max[2]; z += gap_between) {
            // push initial particle quantities
            this.particleBuffer.push(x); // initial position
            this.particleBuffer.push(y);
            this.particleBuffer.push(z);
            this.particleBuffer.push(0); // initial velocity
            this.particleBuffer.push(0);
            this.particleBuffer.push(0);

            // push particles indices
            this.particleIndices.push(particle_counter++ * ATTRIBUTE_COUNT);
          }
        }
      }
      console.log(`Created ${this.count()} particles.`);
    }

    /**
     * @returns The count of particles in the cloud.
     */
    count() {
      return this.particleBuffer.length / ATTRIBUTE_COUNT;
    }

    /**
     * Get the ith particle from the buffer.
     *
     * @param {number} i The index of the particle to retrieve.
     * @returns The particle as an object.
     */
    get(i) {
      if (i < 0 || i >= this.count()) {
        console.error("Index out of bounds in particle buffer!");
        return null;
      }

      return {
        x_position: this.particleBuffer[ATTRIBUTE_COUNT * i],
        y_position: this.particleBuffer[ATTRIBUTE_COUNT * i + 1],
        z_position: this.particleBuffer[ATTRIBUTE_COUNT * i + 2],
        x_velocity: this.particleBuffer[ATTRIBUTE_COUNT * i + 3],
        y_velocity: this.particleBuffer[ATTRIBUTE_COUNT * i + 4],
        z_velocity: this.particleBuffer[ATTRIBUTE_COUNT * i + 5],
      };
    }
  }

  const createAddGravityKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocity_y, dt) {
        return (
          velocity_y[this.thread.z][this.thread.y][this.thread.x] - dt * 9.81
        );
      })
      .setOutput([nx, ny, nz]);

  const createAdvectParticlesKernel = (
    gpu,
    particleCount,
    cellSize,
    nx,
    ny,
    nz
  ) =>
    gpu
      .createKernel(function (
        particles,
        dt,
        velocityFieldX,
        velocityFieldY,
        velocityFieldZ
      ) {
        // check which position component we're looking at
        if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 0) {
          // get position
          let x = particles[this.thread.x];
          let y = particles[this.thread.x + 1];
          let z = particles[this.thread.x + 2];

          // get x velocity
          let vx = particles[this.thread.x + 3];

          // carry out 2nd order Runge-Kutta solver in one dimension
          let k1 = dt * vx;
          let xIntermediate = x + k1 / 2;

          // interpolate the velocity at the intermediate x value
          let lerpWeight =
            (xIntermediate -
              Math.floor(xIntermediate) * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          let vxIntermediate = lerp(
            velocityFieldX[Math.floor(z / this.constants.CELL_SIZE)][
              Math.floor(y / this.constants.CELL_SIZE)
            ][Math.floor(xIntermediate / this.constants.CELL_SIZE)],
            velocityFieldX[Math.floor(z / this.constants.CELL_SIZE)][
              Math.floor(y / this.constants.CELL_SIZE)
            ][Math.ceil(xIntermediate / this.constants.CELL_SIZE)],
            lerpWeight
          );
          let k2 = dt * vxIntermediate;
          let projectedPosition = x + k2;
          if (projectedPosition < this.constants.CELL_SIZE) {
            projectedPosition = this.constants.CELL_SIZE;
          } else if (
            projectedPosition >
            (this.constants.NX - 1) * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              (this.constants.NX - 1) * this.constants.CELL_SIZE;
          }
          return projectedPosition;
        } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 1) {
          // get position
          let x = particles[this.thread.x - 1];
          let y = particles[this.thread.x];
          let z = particles[this.thread.x + 1];

          // get y velocity
          let vy = particles[this.thread.x + 3];

          // carry out 2nd order Runge-Kutta solver in one dimension
          let k1 = dt * vy;
          let yIntermediate = y + k1 / 2;

          // interpolate the velocity at the intermediate y value
          let lerpWeight =
            (yIntermediate -
              Math.floor(yIntermediate) * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          let vyIntermediate = lerp(
            velocityFieldY[Math.floor(z / this.constants.CELL_SIZE)][
              Math.floor(yIntermediate / this.constants.CELL_SIZE)
            ][Math.floor(x / this.constants.CELL_SIZE)],
            velocityFieldY[Math.floor(z / this.constants.CELL_SIZE)][
              Math.ceil(yIntermediate / this.constants.CELL_SIZE)
            ][Math.floor(x / this.constants.CELL_SIZE)],
            lerpWeight
          );
          let k2 = dt * vyIntermediate;
          let projectedPosition = y + k2;
          if (projectedPosition < this.constants.CELL_SIZE) {
            projectedPosition = this.constants.CELL_SIZE;
          } else if (
            projectedPosition >
            (this.constants.NY - 1) * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              (this.constants.NY - 1) * this.constants.CELL_SIZE;
          }
          return projectedPosition;
        } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 2) {
          // get position
          let x = particles[this.thread.x - 2];
          let y = particles[this.thread.x - 1];
          let z = particles[this.thread.x];

          // get z velocity
          let vz = particles[this.thread.x + 3];

          // carry out 2nd order Runge-Kutta solver in one dimension
          let k1 = dt * vz;
          let zIntermediate = z + k1 / 2;

          // interpolate the velocity at the intermediate z value
          let lerpWeight =
            (zIntermediate -
              Math.floor(zIntermediate) * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          let vzIntermediate = lerp(
            velocityFieldZ[Math.floor(zIntermediate / this.constants.CELL_SIZE)][
              Math.floor(y / this.constants.CELL_SIZE)
            ][Math.floor(x / this.constants.CELL_SIZE)],
            velocityFieldZ[Math.ceil(zIntermediate / this.constants.CELL_SIZE)][
              Math.floor(y / this.constants.CELL_SIZE)
            ][Math.floor(x / this.constants.CELL_SIZE)],
            lerpWeight
          );
          let k2 = dt * vzIntermediate;
          let projectedPosition = z + k2;
          if (projectedPosition < this.constants.CELL_SIZE) {
            projectedPosition = this.constants.CELL_SIZE;
          } else if (
            projectedPosition >
            (this.constants.NZ - 1) * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              (this.constants.NZ - 1) * this.constants.CELL_SIZE;
          }
          return projectedPosition;
        } else {
          // don't change the velocities
          return particles[this.thread.x];
        }
      })
      .addFunction(function lerp(a, b, t) {
        return (1 - t) * a + t * b;
      })
      .setConstants({
        ATTRIBUTE_COUNT: ATTRIBUTE_COUNT,
        CELL_SIZE: cellSize,
        NX: nx,
        NY: ny,
        NZ: nz,
      })
      .setOutput([ATTRIBUTE_COUNT * particleCount]);

  const createClassifyVoxelsKernel = (gpu, particleCount, nx, ny, nz) =>
    gpu
      .createKernel(function (voxelStates, particles, cellSize) {
        // get spatial location of grid
        let x = cellSize * this.thread.x;
        let y = cellSize * this.thread.y;
        let z = cellSize * this.thread.z;

        let particle_exists = false;
        for (let i = 0; i < this.constants.particleCount; i++) {
          let pos_x = particles[i * 6];
          let pos_y = particles[i * 6 + 1];
          let pos_z = particles[i * 6 + 2];
          // check if there is a particle in that grid
          if (
            pos_x - x <= cellSize &&
            pos_x - x > 0 &&
            pos_y - y <= cellSize &&
            pos_y - y > 0 &&
            pos_z - z <= cellSize &&
            pos_z - z > 0
          ) {
            particle_exists = true;
            break;
          }
        }
        // set any cell with a particle in it that isn't solid to fluid
        // any cell without a particle or solid is air
        if (particle_exists) {
          if (
            voxelStates[this.thread.z][this.thread.y][this.thread.x] !==
            this.constants.SOLID
          ) {
            return this.constants.FLUID;
          } else {
            return this.constants.SOLID;
          }
        } else {
          // if there isn't, and the state is fluid, flip it to air
          if (
            voxelStates[this.thread.z][this.thread.y][this.thread.x] !==
            this.constants.SOLID
          ) {
            return this.constants.AIR;
          } else {
            return this.constants.SOLID;
          }
        }
      })
      .setConstants({
        particleCount: particleCount,
        AIR: STATE_ENUM.AIR,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createCopyKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (array) {
        return array[this.thread.z][this.thread.y][this.thread.x];
      })
      .setOutput([nx, ny, nz]);

  /**
   * Assumption!
   * This assumes that we're using the default arrangement of having all solid
   * voxels around the edge of our rectangular domain.
   */

  const createEnforceBoundaryXKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocities) {
        if (this.thread.x === 1 || this.thread.x === this.constants.nx - 2) {
          return 0;
        }
        return velocities[this.thread.z][this.thread.y][this.thread.x];
      })
      .setConstants({ nx: nx })
      .setOutput([nx, ny, nz]);

  const createEnforceBoundaryYKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocities) {
        if (this.thread.y === 1 || this.thread.y === this.constants.ny - 2) {
          return 0;
        }
        return velocities[this.thread.z][this.thread.y][this.thread.x];
      })
      .setConstants({ ny: ny })
      .setOutput([nx, ny, nz]);

  const createEnforceBoundaryZKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocities) {
        if (this.thread.z === 1 || this.thread.z === this.constants.nz - 2) {
          return 0;
        }
        return velocities[this.thread.z][this.thread.y][this.thread.x];
      })
      .setConstants({ nz: nz })
      .setOutput([nx, ny, nz]);

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
  const createParticleToGridKernel = (
    gpu,
    particleCount,
    nx,
    ny,
    nz,
    dimension
  ) =>
    gpu
      .createKernel(function (particles, cellSize) {
        // get spatial location of grid velocity vector
        let x = cellSize * this.thread.x;
        let y = cellSize * this.thread.y;
        let z = cellSize * this.thread.z;

        // declare numerator and denominator of the weighted sum
        let numerator = 0;
        let denominator = 0;
        /* loop through particles to find ones that are close, add their
           velocity contribution to the grid velocity */
        for (
          let particleIndex = 0;
          particleIndex < this.constants.PARTICLE_COUNT;
          particleIndex++
        ) {
          // calculate distance in each dimension
          let distance_x = particles[particleIndex * 6] - x;
          let distance_y = particles[particleIndex * 6 + 1] - y;
          let distance_z = particles[particleIndex * 6 + 2] - z;

          // if it's far, skip it
          if (
            Math.abs(distance_x) > cellSize ||
            Math.abs(distance_y) > cellSize ||
            Math.abs(distance_z) > cellSize
          ) {
            continue;
          }

          // calculate the weight according to the trilinear interpolation
          let weight =
            triangle(distance_x / cellSize) *
            triangle(distance_y / cellSize) *
            triangle(distance_z / cellSize);

          numerator +=
            particles[particleIndex * 6 + 3 + this.constants.DIMENSION] * weight;
          denominator += weight;
        }

        // check for divide by zero
        if (Math.abs(denominator) < 0.0001) {
          return 0;
        }
        return numerator / denominator;
      })
      .addFunction(function triangle(r) {
        let r_magnitude = Math.abs(r);
        if (r_magnitude >= 1) {
          return 0;
        }
        return 1 - r_magnitude;
      })
      .setConstants({ PARTICLE_COUNT: particleCount, DIMENSION: dimension })
      .setOutput([nx, ny, nz]);

  // kernel for subracting the new grid velocities from the old grid velocities
  const createGridVelocityDifferenceKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (oldVelocities, newVelocities) {
        return (
          newVelocities[this.thread.z][this.thread.y][this.thread.x] -
          oldVelocities[this.thread.z][this.thread.y][this.thread.x]
        );
      })
      .setOutput([nx, ny, nz]);

  // FLIP Kernel
  const createFLIPKernel = (gpu, particleCount, cellSize) =>
    gpu
      .createKernel(function (particles, diffGridVx, diffGridVy, diffGridVz) {
        // mod to figure out which index we are at (0-5) for each particle
        let index_mod = this.thread.x % this.constants.ATTRIBUTE_COUNT;
        // if we are looking at the position just return the position
        if (index_mod === 0 || index_mod === 1 || index_mod === 2) {
          return particles[this.thread.x];
        }
        // get the positions - index changes depending which velocity we are looking at
        let pos_x = particles[this.thread.x - index_mod];
        let pos_y = particles[this.thread.x - index_mod + 1];
        let pos_z = particles[this.thread.x - index_mod + 2];
        // get the lower and upper grid positions
        let grid_lower_x = Math.floor(pos_x / this.constants.CELL_SIZE);
        let grid_upper_x = Math.ceil(pos_x / this.constants.CELL_SIZE);
        let grid_lower_y = Math.floor(pos_y / this.constants.CELL_SIZE);
        let grid_upper_y = Math.ceil(pos_y / this.constants.CELL_SIZE);
        let grid_lower_z = Math.floor(pos_z / this.constants.CELL_SIZE);
        let grid_upper_z = Math.ceil(pos_z / this.constants.CELL_SIZE);
        if (index_mod === 3) {
          // vx
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_x - grid_lower_x * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
            diffGridVx[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVx[grid_lower_z][grid_lower_y][grid_upper_x],
            lerpWeight
          );
        } else if (index_mod === 4) {
          // vy
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_y - grid_lower_y * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
            diffGridVy[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVy[grid_lower_z][grid_upper_y][grid_lower_x],
            lerpWeight
          );
        } else if (index_mod === 5) {
          // vz
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_z - grid_lower_z * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
            diffGridVz[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVz[grid_upper_z][grid_lower_y][grid_lower_x],
            lerpWeight
          );
        }
      })
      .addFunction(function lerp(a, b, t) {
        return (1 - t) * a + t * b;
      })
      .setConstants({ ATTRIBUTE_COUNT: ATTRIBUTE_COUNT, CELL_SIZE: cellSize })
      .setOutput([particleCount * ATTRIBUTE_COUNT]);

  const createGridToParticlesKernel = (
    gpu,
    particleCount,
    nx,
    ny,
    nz,
    cellSize
  ) => {
    const velocityXDifference = createGridVelocityDifferenceKernel(
      gpu,
      nx + 1,
      ny,
      nz
    ).setPipeline(true);

    const velocityYDifference = createGridVelocityDifferenceKernel(
      gpu,
      nx,
      ny + 1,
      nz
    ).setPipeline(true);

    const velocityZDifference = createGridVelocityDifferenceKernel(
      gpu,
      nx,
      ny,
      nz + 1
    ).setPipeline(true);

    const flipKernel = createFLIPKernel(gpu, particleCount, cellSize).setPipeline(
      true
    );

    return (
      oldXVelocity,
      oldYVelocity,
      oldZVelocity,
      newXVelocity,
      newYVelocity,
      newZVelocity,
      particles
    ) =>
      flipKernel(
        particles,
        velocityXDifference(oldXVelocity, newXVelocity),
        velocityYDifference(oldYVelocity, newYVelocity),
        velocityZDifference(oldZVelocity, newZVelocity)
      );
  };

  /**
   * Component-wise add two vectors.
   */
  const createComponentWiseAddKernel = (gpu, vectorLength) =>
    gpu
      .createKernel(function (a, b) {
        return a[this.thread.x] + b[this.thread.x];
      })
      .setTactic("precision") // vector math should be high precision
      .setOutput([vectorLength]);

  /**
   * Component-wise multiply two vectors.
   */
  const createComponentWiseMultiplyKernel = (gpu, vectorLength) =>
    gpu
      .createKernel(function (a, b) {
        return a[this.thread.x] * b[this.thread.x];
      })
      .setTactic("precision") // vector math should be high precision
      .setOutput([vectorLength]);

  /**
   * Muliply a vector `a` by a scalar.
   */
  const createScalarMultiplyKernel = (gpu, vectorLength) =>
    gpu
      .createKernel(function (a, scalar) {
        return scalar * a[this.thread.x];
      })
      .setTactic("precision") // vector math should be high precision
      .setOutput([vectorLength]);

  /**
   * Assumption!
   * As long as there is the default solid walls around the fluid, edge cases
   * will be fine. Otherwise, this will need to be modified to support more
   * complex boundary conditions.
   */

  const createADiagKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, dt) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;
        const FLUID = this.constants.FLUID;
        const AIR = this.constants.AIR;

        // only consider fluid cells
        if (voxelStates[k][j][i] !== FLUID) {
          return 0;
        }

        const scale =
          dt /
          (this.constants.FLUID_DENSITY *
            this.constants.CELL_SIZE *
            this.constants.CELL_SIZE);

        let accumulator = 0;

        // negative x neighbor
        if (voxelStates[k][j][i - 1] === FLUID) {
          accumulator += scale;
        }
        // positive x neighbor
        if (
          voxelStates[k][j][i + 1] === FLUID ||
          voxelStates[k][j][i + 1] === AIR
        ) {
          accumulator += scale;
        }

        // negative y neighbor
        if (voxelStates[k][j - 1][i] === FLUID) {
          accumulator += scale;
        }
        // positive y neighbor
        if (
          voxelStates[k][j + 1][i] === FLUID ||
          voxelStates[k][j + 1][i] === AIR
        ) {
          accumulator += scale;
        }

        // negative z neighbor
        if (voxelStates[k - 1][j][i] === FLUID) {
          accumulator += scale;
        }
        // positive z neighbor
        if (
          voxelStates[k + 1][j][i] === FLUID ||
          voxelStates[k + 1][j][i] === AIR
        ) {
          accumulator += scale;
        }

        return accumulator;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID_DENSITY: FLUID_DENSITY,
        AIR: STATE_ENUM.AIR,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createAXKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, dt) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;
        const FLUID = this.constants.FLUID;

        // only consider fluid cells
        if (voxelStates[k][j][i] !== FLUID) {
          return 0;
        }

        const scale =
          dt /
          (this.constants.FLUID_DENSITY *
            this.constants.CELL_SIZE *
            this.constants.CELL_SIZE);

        let accumulator = 0;
        //positive x neighbor
        if (voxelStates[k][j][i + 1] === FLUID) {
          accumulator = -scale;
        }
        return accumulator;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID_DENSITY: FLUID_DENSITY,
        AIR: STATE_ENUM.AIR,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createAYKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, dt) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;
        const FLUID = this.constants.FLUID;

        // only consider fluid cells
        if (voxelStates[k][j][i] !== FLUID) {
          return 0;
        }

        const scale =
          dt /
          (this.constants.FLUID_DENSITY *
            this.constants.CELL_SIZE *
            this.constants.CELL_SIZE);

        let accumulator = 0;
        //positive y neighbor
        if (voxelStates[k][j + 1][i] === FLUID) {
          accumulator = -scale;
        }
        return accumulator;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID_DENSITY: FLUID_DENSITY,
        AIR: STATE_ENUM.AIR,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createAZKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, dt) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;
        const FLUID = this.constants.FLUID;

        // only consider fluid cells
        if (voxelStates[k][j][i] !== FLUID) {
          return 0;
        }

        const scale =
          dt /
          (this.constants.FLUID_DENSITY *
            this.constants.CELL_SIZE *
            this.constants.CELL_SIZE);

        let accumulator = 0;
        //positive z neighbor
        if (voxelStates[k + 1][j][i] === FLUID) {
          accumulator = -scale;
        }
        return accumulator;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID_DENSITY: FLUID_DENSITY,
        AIR: STATE_ENUM.AIR,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createNegativeDivergenceKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, velocityX, velocityY, velocityZ) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;

        if (voxelStates[k][j][i] !== FLUID) {
          return 0;
        }

        const scale = 1.0 / this.constants.CELL_SIZE;

        let divergence =
          -scale *
          (velocityX[k][j][i + 1] -
            velocityX[k][j][i] +
            velocityY[k][j + 1][i] -
            velocityY[k][j][i] +
            velocityZ[k + 1][j][i] -
            velocityZ[k][j][i]);

        // modifying RHS (divergence) to account for solid velocities
        if (voxelStates[k][j][i - 1] === SOLID) {
          divergence -= scale * velocityX[k][j][i];
        }
        if (voxelStates[k][j][i + 1] === SOLID) {
          divergence += scale * velocityX[k][j][i + 1];
        }

        if (voxelStates[k][j - 1][i] === SOLID) {
          divergence -= scale * velocityY[k][j][i];
        }
        if (voxelStates[k][j + 1][i] === SOLID) {
          divergence += scale * velocityY[k][j + 1][i];
        }

        if (voxelStates[k - 1][j][i] === SOLID) {
          divergence -= scale * velocityZ[k][j][i];
        }

        if (voxelStates[k + 1][j][i] === SOLID) {
          divergence += scale * velocityZ[k + 1][j][i];
        }

        return divergence;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  /**
   * Map from 3D arrays to 1D vectors.
   *
   * To unflatten:
   * `(k * ny + i) * nx + j`
   */
  const createFlattenKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (array) {
        const aux = this.thread.x % (this.constants.NX * this.constants.NY);
        const i = Math.floor(aux / this.constants.NY);
        const j = aux % this.constants.NX;
        const k = Math.floor(
          this.thread.x / (this.constants.NX * this.constants.NY)
        );
        return array[k][j][i];
      })
      .setTactic("precision")
      .setConstants({ NX: nx, NY: ny, NZ: nz })
      .setOutput([nx * ny * nz]);

  /**
   * Produce the matrix-vector product `Ax` from the sparsely stored A and x.
   */
  const createApplyAKernel = (gpu, vectorLength, nx, ny, nz) =>
    gpu
      .createKernel(function (Adiag, Ax, Ay, Az, x, voxelStates) {
        const aux = this.thread.x % (this.constants.NX * this.constants.NY);
        const i = Math.floor(aux / this.constants.NY);
        const j = aux % this.constants.NX;
        const k = Math.floor(
          this.thread.x / (this.constants.NX * this.constants.NY)
        );

        // only consider fluid cells
        if (voxelStates[k][j][i] !== this.constants.FLUID) {
          return 0;
        }

        let vectorIndex = this.thread.x;
        let accumulator = Adiag[k][j][i] * x[vectorIndex];

        // negative x neighbor
        if (voxelStates[k][j][i - 1] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i - 1, j, k);
          accumulator += Ax[k][j][i - 1] * x[vectorIndex];
        }
        // positive x neighbor
        if (voxelStates[k][j][i + 1] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i + 1, j, k);
          accumulator += Ax[k][j][i + 1] * x[vectorIndex];
        }
        // negative y neighbor
        if (voxelStates[k][j - 1][i] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i, j - 1, k);
          accumulator += Ay[k][j - 1][i] * x[vectorIndex];
        }
        // positive y neighbor
        if (voxelStates[k][j + 1][i] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i, j + 1, k);
          accumulator += Ay[k][j + 1][i] * x[vectorIndex];
        }
        // negative z neighbor
        if (voxelStates[k - 1][j][i] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i, j, k - 1);
          accumulator += Az[k - 1][j][i] * x[vectorIndex];
        }
        // positive z neighbor
        if (voxelStates[k + 1][j][i] === this.constants.FLUID) {
          vectorIndex = gridToVectorIndex(i, j, k + 1);
          accumulator += Az[k + 1][j][i] * x[vectorIndex];
        }

        return accumulator;
      })
      .addFunction(function gridToVectorIndex(i, j, k) {
        return (k * this.constants.NY + i) * this.constants.NX + j;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        VECTOR_LENGTH: vectorLength,
        NX: nx,
        NY: ny,
        NZ: nz,
        FLUID: STATE_ENUM.FLUID,
      })
      .setOutput([vectorLength]);

  const compileKernels = (gpu, particles, grid) => {
    const start = Date.now();

    const gridSize = [grid.nx, grid.ny, grid.nz];
    const velocityXSize = [grid.nx + 1, grid.ny, grid.nz];
    const velocityYSize = [grid.nx, grid.ny + 1, grid.nz];
    const velocityZSize = [grid.nx, grid.ny, grid.nz + 1];
    const DIMENSION = { X: 0, Y: 1, Z: 2 };

    // project particle velocities to the grid
    const particleToXGrid = createParticleToGridKernel(
      gpu,
      particles.count(),
      ...velocityXSize,
      DIMENSION.X
    );
    const particleToYGrid = createParticleToGridKernel(
      gpu,
      particles.count(),
      ...velocityYSize,
      DIMENSION.Y
    );
    const particleToZGrid = createParticleToGridKernel(
      gpu,
      particles.count(),
      ...velocityZSize,
      DIMENSION.Z
    );

    // copy grid quantities to save
    const copyPressure = createCopyKernel(gpu, ...gridSize);
    const copyXVelocity = createCopyKernel(gpu, ...velocityXSize);
    const copyYVelocity = createCopyKernel(gpu, ...velocityYSize);
    const copyZVelocity = createCopyKernel(gpu, ...velocityZSize);

    // mark cells as solid, fluid, or air
    const classifyVoxels = createClassifyVoxelsKernel(
      gpu,
      particles.count(),
      ...gridSize
    );

    // add gravitational influence
    const addGravity = createAddGravityKernel(gpu, ...velocityYSize);

    // enforce boundary conditions
    const enforceXBoundary = createEnforceBoundaryXKernel(gpu, ...velocityXSize);
    const enforceYBoundary = createEnforceBoundaryYKernel(gpu, ...velocityYSize);
    const enforceZBoundary = createEnforceBoundaryZKernel(gpu, ...velocityZSize);

    // do pressure solve

    // build coefficient matrix
    const buildADiag = createADiagKernel(gpu, ...gridSize, grid.cellSize);
    const buildAX = createAXKernel(gpu, ...gridSize, grid.cellSize);
    const buildAY = createAYKernel(gpu, ...gridSize, grid.cellSize);
    const buildAZ = createAZKernel(gpu, ...gridSize, grid.cellSize);

    // build negative divergence vector with boundary conditions
    const buildD = createNegativeDivergenceKernel(
      gpu,
      ...gridSize,
      grid.cellSize
    );
    const flatten = createFlattenKernel(gpu, ...gridSize);

    // compile kernels to do vector operations
    const pcgVectorLength = grid.nx * grid.ny * grid.nz;
    const componentWiseAdd = createComponentWiseAddKernel(gpu, pcgVectorLength);
    const componentWiseMultiply = createComponentWiseMultiplyKernel(
      gpu,
      pcgVectorLength
    );
    // implement dot product's sum portion on the CPU
    const dot = (a, b) =>
      componentWiseMultiply(a, b).reduce((sum, n) => sum + n, 0);
    const scalarMultiply = createScalarMultiplyKernel(gpu, pcgVectorLength);
    const applyA = createApplyAKernel(gpu, pcgVectorLength, ...gridSize);
    const math = {
      componentWiseAdd: componentWiseAdd.setPipeline(true),
      dot: dot,
      scalarMultiply: scalarMultiply.setPipeline(true),
      applyA: applyA.setPipeline(true),
    };

    // PCG methods
    const zeroVector = gpu
      .createKernel(function () {
        return 0;
      })
      .setOutput([pcgVectorLength]);

    const pressureSolve = {
      buildADiag: buildADiag.setPipeline(true),
      buildAX: buildAX.setPipeline(true),
      buildAY: buildAY.setPipeline(true),
      buildAZ: buildAZ.setPipeline(true),
      buildD: buildD.setPipeline(true),
      flatten: flatten.setPipeline(true),
      math: math,
      zeroVector: zeroVector.setPipeline(true),
    };

    // update the velocities of the particles using PIC/FLIP
    const gridToParticles = createGridToParticlesKernel(
      gpu,
      particles.count(),
      ...gridSize,
      grid.cellSize
    );

    // update the positions of the particles
    const advectParticles = createAdvectParticlesKernel(
      gpu,
      particles.count(),
      grid.cellSize,
      ...gridSize
    );

    const end = Date.now();
    console.log(`Kernels compiled in ${end - start} ms.`);

    return {
      particleToXGrid: particleToXGrid.setPipeline(true),
      particleToYGrid: particleToYGrid.setPipeline(true),
      particleToZGrid: particleToZGrid.setPipeline(true),
      copyPressure: copyPressure.setPipeline(true),
      copyXVelocity: copyXVelocity.setPipeline(true),
      copyYVelocity: copyYVelocity.setPipeline(true),
      copyZVelocity: copyZVelocity.setPipeline(true),
      classifyVoxels: classifyVoxels.setPipeline(true),
      addGravity: addGravity.setPipeline(true),
      enforceXBoundary: enforceXBoundary.setPipeline(true),
      enforceYBoundary: enforceYBoundary.setPipeline(true),
      enforceZBoundary: enforceZBoundary.setPipeline(true),
      gridToParticles: gridToParticles,
      advectParticles: advectParticles.setPipeline(true),
      pressureSolve: pressureSolve,
    };
  };

  const FLUID_DENSITY = 997;

  class Simulation {
    constructor(gpu, config) {
      this.particles = new Particles(
        config.particleDensity,
        config.particleBounds
      );
      this.grid = new MACGrid(
        config.gridBounds,
        2.0 / Math.cbrt(config.particleDensity)
      );
      this.grid.addDefaultSolids();
      this.kernels = compileKernels(gpu, this.particles, this.grid);
    }

    step(dt) {
      let particleBufferCopy = new Float32Array(this.particles.particleBuffer);
      // transfer particle velocities to the grid and interpolate
      this.grid.velocityX = this.kernels.particleToXGrid(
        particleBufferCopy,
        this.grid.cellSize
      );
      this.grid.velocityY = this.kernels.particleToYGrid(
        particleBufferCopy,
        this.grid.cellSize
      );
      this.grid.velocityZ = this.kernels.particleToZGrid(
        particleBufferCopy,
        this.grid.cellSize
      );

      // copy grid values to store the old ones
      this.grid.pressureOld = this.kernels.copyPressure(this.grid.pressure);
      this.grid.velocityXOld = this.kernels.copyXVelocity(this.grid.velocityX);
      this.grid.velocityYOld = this.kernels.copyYVelocity(this.grid.velocityY);
      this.grid.velocityZOld = this.kernels.copyZVelocity(this.grid.velocityZ);

      // mark cells as solid, fluid, or air
      this.grid.voxelStates = this.kernels.classifyVoxels(
        this.grid.voxelStates.toArray(),
        particleBufferCopy,
        this.grid.cellSize
      );

      // perform gravity update
      this.grid.velocityY = this.kernels.addGravity(this.grid.velocityY, dt);

      // enforce boundary conditions
      //this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
      this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
      //this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

      // do the pressure solve with a zero divergence velocity field
      // TODO: implement this!

      // enforce boundary conditions
      // this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
      // this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
      // this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

      // update the velocities of the particles
      this.particles.particleBuffer = this.kernels
        .gridToParticles(
          this.grid.velocityXOld,
          this.grid.velocityYOld,
          this.grid.velocityZOld,
          this.grid.velocityX,
          this.grid.velocityY,
          this.grid.velocityZ,
          particleBufferCopy
        )
        .toArray();
      // advect the particles to find their new positions
      this.particles.particleBuffer = this.kernels
        .advectParticles(
          new Float32Array(this.particles.particleBuffer),
          dt,
          this.grid.velocityX,
          this.grid.velocityY,
          this.grid.velocityZ
        )
        .toArray();
    }
  }

  /*
   * Copyright 2010, Google Inc.
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are
   * met:
   *
   *     * Redistributions of source code must retain the above copyright
   * notice, this list of conditions and the following disclaimer.
   *     * Redistributions in binary form must reproduce the above
   * copyright notice, this list of conditions and the following disclaimer
   * in the documentation and/or other materials provided with the
   * distribution.
   *     * Neither the name of Google Inc. nor the names of its
   * contributors may be used to endorse or promote products derived from
   * this software without specific prior written permission.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   */

  function RayMarchingEffect(resolution) {
    var ext = gl.getExtension("OES_texture_float");
    if (!ext) {
      alert("this machine or browser does not support OES_texture_float");
      return;
    }

    var arrays = tdl.primitives.createCube(1.0);
    var program = tdl.programs.loadProgramFromScriptTags("ray_vs", "ray_fs");
    var textures = [new tdl.textures.ExternalTexture2D()];

    var model = new tdl.models.Model(program, arrays, textures);

    var size = resolution;

    var size3 = size * size * size;
    var max_tex_dim = 16384;
    if (size3 > max_tex_dim * 4) {
      alert("Resolution too high! Something's wrong.");
    }

    var field = new Float32Array(max_tex_dim * 4);

    var tex = textures[0].texture;
    var tex_level = 0;
    var tex_width = max_tex_dim;
    var tex_height = 1;

    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    var startTime = Date.now() / 1000;
    var lastTime = startTime;

    const gpu = new GPU();
    const sim = new Simulation(gpu, {
      particleDensity: 10000,
      particleBounds: {
        min: fromValues(0.3, 0.3, 0.3),
        max: fromValues(0.7, 0.7, 0.7),
      },
      gridBounds: {
        min: fromValues(0.1, 0.1, 0.1),
        max: fromValues(0.9, 0.9, 0.9),
      },
    });

    const fillField = gpu
      .createKernel(function (balls, n, size, radius) {
        let z = Math.floor(this.thread.x / (size * size));
        let y = Math.floor(this.thread.x / size) % size;
        let x = this.thread.x % size;
        var z_w = z / size;
        var y_w = y / size;
        var x_w = x / size;
        let closest = 0;
        let best = 100000;
        // if n too big will need to change loopmaxiterations
        for (let i = 0; i < n; ++i) {
          let cur =
            (x_w - balls[i][0]) * (x_w - balls[i][0]) +
            (y_w - balls[i][1]) * (y_w - balls[i][1]) +
            (z_w - balls[i][2]) * (z_w - balls[i][2]);
          if (cur < best) {
            best = cur;
            closest = i;
          }
        }

        return (
          Math.sqrt(
            (x_w - balls[closest][0]) * (x_w - balls[closest][0]) +
            (y_w - balls[closest][1]) * (y_w - balls[closest][1]) +
            (z_w - balls[closest][2]) * (z_w - balls[closest][2])
          ) - radius
        );
      })
      .setLoopMaxIterations(10000)
      .setPipeline(true)
      .setOutput([max_tex_dim * 4]);

    const smooth = gpu
      .createKernel(function (field, size) {
        let z_c = Math.floor(this.thread.x / (size * size));
        let y_c = Math.floor(this.thread.x / size) % size;
        let x_c = this.thread.x % size;

        let sum = 0;
        let count = 0.001; // weight must be nonzero
        let r = 2;
        for (let z_o = -r; z_o <= r; ++z_o) {
          for (let y_o = -r; y_o <= r; ++y_o) {
            for (let x_o = -r; x_o <= r; ++x_o) {
              let x = x_c + x_o;
              let y = y_c + y_o;
              let z = z_c + z_o;

              // Weighted by e^(-c * r^2)
              let w = Math.pow(
                2.71,
                -1.5 * Math.sqrt(x_o * x_o + y_o * y_o + z_o * z_o)
              );

              if (
                x < 0 ||
                x > size - 1 ||
                y < 0 ||
                y > size - 1 ||
                z < 0 ||
                z > size - 1
              ) {
                // cheaper than continue, will try to read invalid data
                w = 0.0;
              }

              sum += field[z * size * size + y * size + x] * w;
              count += w;
            }
          }
        }

        return sum / count;
      })
      .setPipeline(true)
      .setOutput([max_tex_dim * 4]);

    this.render = function () {
      gl.clearColor(0.0, 0.0, 0.0, 1);
      gl.clearDepth(1.0);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      gl.enable(gl.DEPTH_TEST);
      gl.enable(gl.CULL_FACE);

      let curTime = Date.now() / 1000;
      let deltaTime = curTime - lastTime;
      lastTime = curTime;

      let localTime = Date.now() / 1000 - startTime;

      if (deltaTime > 1 / 20) {
        deltaTime = 1 / 20;
      }

      // step the simulation forwards
      sim.step(deltaTime);

      var uniformsConst = {
        u_field: textures[0],
        time: localTime,
      };

      model.drawPrep(uniformsConst);

      {
        // Set firstDraw = false to only draw 1 frame
        
        // Sine wave water
        /*
        let balls = [];
        let n = 30;
        let radius = 0.04;
        
        for (let x = 0; x < n; ++x) {
          for (let z = 0; z < n; ++z) {
            let xp = (x+0.5) / n;
            let zp = (z+0.5) / n;
            let r = Math.sqrt((xp-0.5) * (xp-0.5) + (zp-0.5) * (zp-0.5));
            //let y = 0.1*((Math.sin(40 * r - 1.5*time) + 1) / 2) / Math.abs(10*(Math.max(r, 0.013))) + 0.05;
            let y = 0.3 * Math.pow(Math.cos(10 * r - 1 * time), 2) / Math.max(10*r, 0.5) + 0.05;
            balls.push([xp, y, zp]);
          }
        }
        */

        let balls = [];
        let radius = 0.04;
        for (let i = 0; i < sim.particles.particleBuffer.length; i += 6) {
          balls.push([
            sim.particles.particleBuffer[i],
            sim.particles.particleBuffer[i + 1],
            sim.particles.particleBuffer[i + 2],
          ]);
        }

        // Swap comment to see with / without smoothing
        //field = fillField(balls, balls.length, size, radius).toArray();
        field = smooth(
          fillField(balls, balls.length, size, radius),
          size
        ).toArray();
      }

      // Send the field to GPU, issue draw
      imm.begin(gl.TRIANGLES, program);

      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texImage2D(
        gl.TEXTURE_2D,
        tex_level,
        gl.RGBA,
        tex_width,
        tex_height,
        0,
        gl.RGBA,
        gl.FLOAT,
        field
      );

      gl.activeTexture(gl.TEXTURE0);

      imm.quad2d(-1, -1, 2, 2, 1);
      imm.end();
    };
  }

  window.RayMarchingEffect = RayMarchingEffect;

}());
//# sourceMappingURL=bundle.js.map
