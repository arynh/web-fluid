(function () {
  'use strict';

  /**
   * Common utilities
   * @module glMatrix
   */
  // Configuration Constants
  var EPSILON = 0.000001;
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
   * 3x3 Matrix
   * @module mat3
   */

  /**
   * Creates a new identity mat3
   *
   * @returns {mat3} a new 3x3 matrix
   */

  function create$3() {
    var out = new ARRAY_TYPE(9);

    if (ARRAY_TYPE != Float32Array) {
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
      out[5] = 0;
      out[6] = 0;
      out[7] = 0;
    }

    out[0] = 1;
    out[4] = 1;
    out[8] = 1;
    return out;
  }

  /**
   * 3 Dimensional Vector
   * @module vec3
   */

  /**
   * Creates a new, empty vec3
   *
   * @returns {vec3} a new 3D vector
   */

  function create$2() {
    var out = new ARRAY_TYPE(3);

    if (ARRAY_TYPE != Float32Array) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
    }

    return out;
  }
  /**
   * Calculates the length of a vec3
   *
   * @param {ReadonlyVec3} a vector to calculate length of
   * @returns {Number} length of a
   */

  function length(a) {
    var x = a[0];
    var y = a[1];
    var z = a[2];
    return Math.hypot(x, y, z);
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
   * Normalize a vec3
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a vector to normalize
   * @returns {vec3} out
   */

  function normalize$2(out, a) {
    var x = a[0];
    var y = a[1];
    var z = a[2];
    var len = x * x + y * y + z * z;

    if (len > 0) {
      //TODO: evaluate use of glm_invsqrt here?
      len = 1 / Math.sqrt(len);
    }

    out[0] = a[0] * len;
    out[1] = a[1] * len;
    out[2] = a[2] * len;
    return out;
  }
  /**
   * Calculates the dot product of two vec3's
   *
   * @param {ReadonlyVec3} a the first operand
   * @param {ReadonlyVec3} b the second operand
   * @returns {Number} dot product of a and b
   */

  function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }
  /**
   * Computes the cross product of two vec3's
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a the first operand
   * @param {ReadonlyVec3} b the second operand
   * @returns {vec3} out
   */

  function cross(out, a, b) {
    var ax = a[0],
        ay = a[1],
        az = a[2];
    var bx = b[0],
        by = b[1],
        bz = b[2];
    out[0] = ay * bz - az * by;
    out[1] = az * bx - ax * bz;
    out[2] = ax * by - ay * bx;
    return out;
  }
  /**
   * Transforms the vec3 with a quat
   * Can also be used for dual quaternions. (Multiply it with the real part)
   *
   * @param {vec3} out the receiving vector
   * @param {ReadonlyVec3} a the vector to transform
   * @param {ReadonlyQuat} q quaternion to transform with
   * @returns {vec3} out
   */

  function transformQuat(out, a, q) {
    // benchmarks: https://jsperf.com/quaternion-transform-vec3-implementations-fixed
    var qx = q[0],
        qy = q[1],
        qz = q[2],
        qw = q[3];
    var x = a[0],
        y = a[1],
        z = a[2]; // var qvec = [qx, qy, qz];
    // var uv = vec3.cross([], qvec, a);

    var uvx = qy * z - qz * y,
        uvy = qz * x - qx * z,
        uvz = qx * y - qy * x; // var uuv = vec3.cross([], qvec, uv);

    var uuvx = qy * uvz - qz * uvy,
        uuvy = qz * uvx - qx * uvz,
        uuvz = qx * uvy - qy * uvx; // vec3.scale(uv, uv, 2 * w);

    var w2 = qw * 2;
    uvx *= w2;
    uvy *= w2;
    uvz *= w2; // vec3.scale(uuv, uuv, 2);

    uuvx *= 2;
    uuvy *= 2;
    uuvz *= 2; // return vec3.add(out, a, vec3.add(out, uv, uuv));

    out[0] = x + uvx + uuvx;
    out[1] = y + uvy + uuvy;
    out[2] = z + uvz + uuvz;
    return out;
  }
  /**
   * Alias for {@link vec3.subtract}
   * @function
   */

  var sub = subtract;
  /**
   * Alias for {@link vec3.length}
   * @function
   */

  var len = length;
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
    var vec = create$2();
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
   * 4 Dimensional Vector
   * @module vec4
   */

  /**
   * Creates a new, empty vec4
   *
   * @returns {vec4} a new 4D vector
   */

  function create$1() {
    var out = new ARRAY_TYPE(4);

    if (ARRAY_TYPE != Float32Array) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
    }

    return out;
  }
  /**
   * Normalize a vec4
   *
   * @param {vec4} out the receiving vector
   * @param {ReadonlyVec4} a vector to normalize
   * @returns {vec4} out
   */

  function normalize$1(out, a) {
    var x = a[0];
    var y = a[1];
    var z = a[2];
    var w = a[3];
    var len = x * x + y * y + z * z + w * w;

    if (len > 0) {
      len = 1 / Math.sqrt(len);
    }

    out[0] = x * len;
    out[1] = y * len;
    out[2] = z * len;
    out[3] = w * len;
    return out;
  }
  /**
   * Perform some operation over an array of vec4s.
   *
   * @param {Array} a the array of vectors to iterate over
   * @param {Number} stride Number of elements between the start of each vec4. If 0 assumes tightly packed
   * @param {Number} offset Number of elements to skip at the beginning of the array
   * @param {Number} count Number of vec4s to iterate over. If 0 iterates over entire array
   * @param {Function} fn Function to call for each vector in the array
   * @param {Object} [arg] additional argument to pass to fn
   * @returns {Array} a
   * @function
   */

  (function () {
    var vec = create$1();
    return function (a, stride, offset, count, fn, arg) {
      var i, l;

      if (!stride) {
        stride = 4;
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
        vec[3] = a[i + 3];
        fn(vec, vec, arg);
        a[i] = vec[0];
        a[i + 1] = vec[1];
        a[i + 2] = vec[2];
        a[i + 3] = vec[3];
      }

      return a;
    };
  })();

  /**
   * Quaternion
   * @module quat
   */

  /**
   * Creates a new identity quat
   *
   * @returns {quat} a new quaternion
   */

  function create() {
    var out = new ARRAY_TYPE(4);

    if (ARRAY_TYPE != Float32Array) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
    }

    out[3] = 1;
    return out;
  }
  /**
   * Set a quat to the identity quaternion
   *
   * @param {quat} out the receiving quaternion
   * @returns {quat} out
   */

  function identity(out) {
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
    out[3] = 1;
    return out;
  }
  /**
   * Sets a quat from the given angle and rotation axis,
   * then returns it.
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyVec3} axis the axis around which to rotate
   * @param {Number} rad the angle in radians
   * @returns {quat} out
   **/

  function setAxisAngle(out, axis, rad) {
    rad = rad * 0.5;
    var s = Math.sin(rad);
    out[0] = s * axis[0];
    out[1] = s * axis[1];
    out[2] = s * axis[2];
    out[3] = Math.cos(rad);
    return out;
  }
  /**
   * Multiplies two quat's
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyQuat} a the first operand
   * @param {ReadonlyQuat} b the second operand
   * @returns {quat} out
   */

  function multiply(out, a, b) {
    var ax = a[0],
        ay = a[1],
        az = a[2],
        aw = a[3];
    var bx = b[0],
        by = b[1],
        bz = b[2],
        bw = b[3];
    out[0] = ax * bw + aw * bx + ay * bz - az * by;
    out[1] = ay * bw + aw * by + az * bx - ax * bz;
    out[2] = az * bw + aw * bz + ax * by - ay * bx;
    out[3] = aw * bw - ax * bx - ay * by - az * bz;
    return out;
  }
  /**
   * Performs a spherical linear interpolation between two quat
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyQuat} a the first operand
   * @param {ReadonlyQuat} b the second operand
   * @param {Number} t interpolation amount, in the range [0-1], between the two inputs
   * @returns {quat} out
   */

  function slerp(out, a, b, t) {
    // benchmarks:
    //    http://jsperf.com/quaternion-slerp-implementations
    var ax = a[0],
        ay = a[1],
        az = a[2],
        aw = a[3];
    var bx = b[0],
        by = b[1],
        bz = b[2],
        bw = b[3];
    var omega, cosom, sinom, scale0, scale1; // calc cosine

    cosom = ax * bx + ay * by + az * bz + aw * bw; // adjust signs (if necessary)

    if (cosom < 0.0) {
      cosom = -cosom;
      bx = -bx;
      by = -by;
      bz = -bz;
      bw = -bw;
    } // calculate coefficients


    if (1.0 - cosom > EPSILON) {
      // standard case (slerp)
      omega = Math.acos(cosom);
      sinom = Math.sin(omega);
      scale0 = Math.sin((1.0 - t) * omega) / sinom;
      scale1 = Math.sin(t * omega) / sinom;
    } else {
      // "from" and "to" quaternions are very close
      //  ... so we can do a linear interpolation
      scale0 = 1.0 - t;
      scale1 = t;
    } // calculate final values


    out[0] = scale0 * ax + scale1 * bx;
    out[1] = scale0 * ay + scale1 * by;
    out[2] = scale0 * az + scale1 * bz;
    out[3] = scale0 * aw + scale1 * bw;
    return out;
  }
  /**
   * Creates a quaternion from the given 3x3 rotation matrix.
   *
   * NOTE: The resultant quaternion is not normalized, so you should be sure
   * to renormalize the quaternion yourself where necessary.
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyMat3} m rotation matrix
   * @returns {quat} out
   * @function
   */

  function fromMat3(out, m) {
    // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
    // article "Quaternion Calculus and Fast Animation".
    var fTrace = m[0] + m[4] + m[8];
    var fRoot;

    if (fTrace > 0.0) {
      // |w| > 1/2, may as well choose w > 1/2
      fRoot = Math.sqrt(fTrace + 1.0); // 2w

      out[3] = 0.5 * fRoot;
      fRoot = 0.5 / fRoot; // 1/(4w)

      out[0] = (m[5] - m[7]) * fRoot;
      out[1] = (m[6] - m[2]) * fRoot;
      out[2] = (m[1] - m[3]) * fRoot;
    } else {
      // |w| <= 1/2
      var i = 0;
      if (m[4] > m[0]) i = 1;
      if (m[8] > m[i * 3 + i]) i = 2;
      var j = (i + 1) % 3;
      var k = (i + 2) % 3;
      fRoot = Math.sqrt(m[i * 3 + i] - m[j * 3 + j] - m[k * 3 + k] + 1.0);
      out[i] = 0.5 * fRoot;
      fRoot = 0.5 / fRoot;
      out[3] = (m[j * 3 + k] - m[k * 3 + j]) * fRoot;
      out[j] = (m[j * 3 + i] + m[i * 3 + j]) * fRoot;
      out[k] = (m[k * 3 + i] + m[i * 3 + k]) * fRoot;
    }

    return out;
  }
  /**
   * Alias for {@link quat.multiply}
   * @function
   */

  var mul = multiply;
  /**
   * Normalize a quat
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyQuat} a quaternion to normalize
   * @returns {quat} out
   * @function
   */

  var normalize = normalize$1;
  /**
   * Sets a quaternion to represent the shortest rotation from one
   * vector to another.
   *
   * Both vectors are assumed to be unit length.
   *
   * @param {quat} out the receiving quaternion.
   * @param {ReadonlyVec3} a the initial vector
   * @param {ReadonlyVec3} b the destination vector
   * @returns {quat} out
   */

  (function () {
    var tmpvec3 = create$2();
    var xUnitVec3 = fromValues(1, 0, 0);
    var yUnitVec3 = fromValues(0, 1, 0);
    return function (out, a, b) {
      var dot$1 = dot(a, b);

      if (dot$1 < -0.999999) {
        cross(tmpvec3, xUnitVec3, a);
        if (len(tmpvec3) < 0.000001) cross(tmpvec3, yUnitVec3, a);
        normalize$2(tmpvec3, tmpvec3);
        setAxisAngle(out, tmpvec3, Math.PI);
        return out;
      } else if (dot$1 > 0.999999) {
        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1;
        return out;
      } else {
        cross(tmpvec3, a, b);
        out[0] = tmpvec3[0];
        out[1] = tmpvec3[1];
        out[2] = tmpvec3[2];
        out[3] = 1 + dot$1;
        return normalize(out, out);
      }
    };
  })();
  /**
   * Performs a spherical linear interpolation with two control points
   *
   * @param {quat} out the receiving quaternion
   * @param {ReadonlyQuat} a the first operand
   * @param {ReadonlyQuat} b the second operand
   * @param {ReadonlyQuat} c the third operand
   * @param {ReadonlyQuat} d the fourth operand
   * @param {Number} t interpolation amount, in the range [0-1], between the two inputs
   * @returns {quat} out
   */

  (function () {
    var temp1 = create();
    var temp2 = create();
    return function (out, a, b, c, d, t) {
      slerp(temp1, a, d, t);
      slerp(temp2, b, c, t);
      slerp(out, temp1, temp2, 2 * t * (1 - t));
      return out;
    };
  })();
  /**
   * Sets the specified quaternion with values corresponding to the given
   * axes. Each axis is a vec3 and is expected to be unit length and
   * perpendicular to all other specified axes.
   *
   * @param {ReadonlyVec3} view  the vector representing the viewing direction
   * @param {ReadonlyVec3} right the vector representing the local "right" direction
   * @param {ReadonlyVec3} up    the vector representing the local "up" direction
   * @returns {quat} out
   */

  (function () {
    var matr = create$3();
    return function (out, view, right, up) {
      matr[0] = right[0];
      matr[3] = right[1];
      matr[6] = right[2];
      matr[1] = up[0];
      matr[4] = up[1];
      matr[7] = up[2];
      matr[2] = -view[0];
      matr[5] = -view[1];
      matr[8] = -view[2];
      return normalize(out, fromMat3(out, matr));
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
      this.count = create$2();
      this.size = create$2();
      sub(this.size, this.max, this.min);
      scale(this.count, this.size, 1.0 / this.cellSize);
      add(this.count, this.count, fromValues(1, 1, 1));
      floor(this.count, this.count);

      /// initialize
      this.nx = this.count[0] - 1;
      this.ny = this.count[1] - 1;
      this.nz = this.count[2] - 1;
      this.pressure = initialize3DArray(this.nx, this.ny, this.nz);
      this.pressureOld = initialize3DArray(this.nx, this.ny, this.nz);
      this.velocityX = null;
      this.velocityXOld = null;
      this.velocityY = null;
      this.velocityYOld = null;
      this.velocityZ = null;
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

  const solve = (
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

  const createAddGravityKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocity_y, dt, voxelStates) {
        // non-physic gravity for the lolz
        return velocity_y[this.thread.z][this.thread.y][this.thread.x] - dt * 6.0;
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
          particles[this.thread.x + 1];
          particles[this.thread.x + 2];

          // get x velocity
          let vx = particles[this.thread.x + 3];

          return x + dt * vx;

          /* FIXME: Runge-Kutta solver does not work yet

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
          if (projectedPosition <= 0) {
            projectedPosition = 0.01;
          } else if (
            projectedPosition >=
            this.constants.NX * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              this.constants.NX * this.constants.CELL_SIZE - 0.01;
          }
          return projectedPosition;

          */
        } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 1) {
          // get position
          particles[this.thread.x - 1];
          let y = particles[this.thread.x];
          particles[this.thread.x + 1];

          // get y velocity
          let vy = particles[this.thread.x + 3];

          return y + dt * vy;

          /* FIXME: Runge-Kutta solver does not work yet

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
          if (projectedPosition <= 0) {
            projectedPosition = 0.01;
          } else if (
            projectedPosition >=
            this.constants.NY * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              this.constants.NY * this.constants.CELL_SIZE - 0.01;
          }
          return projectedPosition;

          */
        } else if (this.thread.x % this.constants.ATTRIBUTE_COUNT === 2) {
          // get position
          particles[this.thread.x - 2];
          particles[this.thread.x - 1];
          let z = particles[this.thread.x];

          // get z velocity
          let vz = particles[this.thread.x + 3];

          return z + dt * vz;

          /* FIXME: Runge-Kutta solver does not work yet

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
          if (projectedPosition <= 0) {
            projectedPosition = 0.01;
          } else if (
            projectedPosition >=
            this.constants.NZ * this.constants.CELL_SIZE
          ) {
            projectedPosition =
              this.constants.NZ * this.constants.CELL_SIZE - 0.01;
          }
          return projectedPosition;
          
          */
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
        if (this.thread.x === 0 || this.thread.x === this.constants.nx - 1) {
          return 0;
        }
        return velocities[this.thread.z][this.thread.y][this.thread.x];
      })
      .setConstants({ nx: nx })
      .setOutput([nx, ny, nz]);

  const createEnforceBoundaryYKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocities) {
        if (this.thread.y === 0 || this.thread.y === this.constants.ny - 1) {
          return 0;
        }
        return velocities[this.thread.z][this.thread.y][this.thread.x];
      })
      .setConstants({ ny: ny })
      .setOutput([nx, ny, nz]);

  const createEnforceBoundaryZKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (velocities) {
        if (this.thread.z === 0 || this.thread.z === this.constants.nz - 1) {
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
      .createKernel(function (
        particles,
        diffGridVx,
        diffGridVy,
        diffGridVz,
        oldVx,
        oldVy,
        oldVz
      ) {
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
        let alpha = 0.95;
        if (index_mod === 3) {
          // vx
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_x - grid_lower_x * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
          lerp(
            diffGridVx[grid_lower_z][grid_lower_y][grid_lower_x]+oldVx[grid_lower_z][grid_lower_y][grid_lower_x],
            diffGridVx[grid_lower_z][grid_lower_y][grid_upper_x]+oldVx[grid_lower_z][grid_lower_y][grid_upper_x],
            lerpWeight
          ),(
            lerp(
              oldVx[grid_lower_z][grid_lower_y][grid_lower_x],
              oldVx[grid_lower_z][grid_lower_y][grid_upper_x],
              lerpWeight
            ) +
            lerp(
              diffGridVx[grid_lower_z][grid_lower_y][grid_lower_x],
              diffGridVx[grid_lower_z][grid_lower_y][grid_upper_x],
              lerpWeight
            )
          ),alpha);
        } else if (index_mod === 4) {
          // vy
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_y - grid_lower_y * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
            lerp(
              diffGridVy[grid_lower_z][grid_lower_y][grid_lower_x]+oldVy[grid_lower_z][grid_lower_y][grid_lower_x],
              diffGridVy[grid_lower_z][grid_upper_y][grid_lower_x]+oldVy[grid_lower_z][grid_upper_y][grid_lower_x],
              lerpWeight)
            ,
            (lerp(
              oldVy[grid_lower_z][grid_lower_y][grid_lower_x],
              oldVy[grid_lower_z][grid_upper_y][grid_lower_x],
              lerpWeight
            ) +
            lerp(
              diffGridVy[grid_lower_z][grid_lower_y][grid_lower_x],
              diffGridVy[grid_lower_z][grid_upper_y][grid_lower_x],
              lerpWeight
            )
          ),alpha);
        } else if (index_mod === 5) {
          // vz
          // get the lerp weight and return the lerp'd velocity
          let lerpWeight =
            (pos_z - grid_lower_z * this.constants.CELL_SIZE) /
            this.constants.CELL_SIZE;
          return lerp(
            lerp(
              diffGridVz[grid_lower_z][grid_lower_y][grid_lower_x]+oldVz[grid_lower_z][grid_lower_y][grid_lower_x],
              diffGridVz[grid_upper_z][grid_lower_y][grid_lower_x]+oldVz[grid_upper_z][grid_lower_y][grid_lower_x],
              lerpWeight)
            ,(
            lerp(
              oldVz[grid_lower_z][grid_lower_y][grid_lower_x],
              oldVz[grid_upper_z][grid_lower_y][grid_lower_x],
              lerpWeight
            ) +
            lerp(
              diffGridVz[grid_lower_z][grid_lower_y][grid_lower_x],
              diffGridVz[grid_upper_z][grid_lower_y][grid_lower_x],
              lerpWeight
            )
          ),alpha);
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
        velocityZDifference(oldZVelocity, newZVelocity),
        oldXVelocity,
        oldYVelocity,
        oldZVelocity
      );
  };

  const createNegativeDivergenceKernel = (gpu, nx, ny, nz, cellSize) =>
    gpu
      .createKernel(function (voxelStates, velocityX, velocityY, velocityZ) {
        // for brevity
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;

        if (voxelStates[k][j][i] !== this.constants.FLUID) {
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

        /* TODO: experiment with these boundary conditions to see if we can
            lose less energy

        // modifying RHS (divergence) to account for solid velocities
        if (voxelStates[k][j][i - 1] === this.constants.SOLID) {
          divergence -= scale * velocityX[k][j][i];
        }
        if (voxelStates[k][j][i + 1] === this.constants.SOLID) {
          divergence += scale * velocityX[k][j][i + 1];
        }

        if (voxelStates[k][j - 1][i] === this.constants.SOLID) {
          divergence -= scale * velocityY[k][j][i];
        }
        if (voxelStates[k][j + 1][i] === this.constants.SOLID) {
          divergence += scale * velocityY[k][j + 1][i];
        }

        if (voxelStates[k - 1][j][i] === this.constants.SOLID) {
          divergence -= scale * velocityZ[k][j][i];
        }

        if (voxelStates[k + 1][j][i] === this.constants.SOLID) {
          divergence += scale * velocityZ[k + 1][j][i];
        }
        */

        return divergence;
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        CELL_SIZE: cellSize,
        FLUID: STATE_ENUM.FLUID,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

  const createVelocityXUpdateKernel = (
    gpu,
    nx,
    ny,
    nz,
    fluidDensity,
    cellSize
  ) =>
    gpu
      .createKernel(function (velocity, pressure, voxelStates, dt) {
        // only consider boundaries which have a fluid cell on at least one side
        if (
          this.thread.x === 0 ||
          this.thread.x === this.constants.NX - 1 ||
          !(
            voxelStates[this.thread.z][this.thread.y][this.thread.x - 1] ===
              this.constants.FLUID ||
            voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
              this.constants.FLUID
          )
        ) {
          return 0;
        }

        const pressureGradient =
          (pressure[this.thread.z][this.thread.y][this.thread.x] -
            pressure[this.thread.z][this.thread.y][this.thread.x - 1]) /
          this.constants.CELL_SIZE;

        const oldVelocity = velocity[this.thread.z][this.thread.y][this.thread.x];
        const newVelocity =
          oldVelocity - (dt * pressureGradient) / this.constants.FLUID_DENSITY;

        return newVelocity;
      })
      .setConstants({
        FLUID_DENSITY: fluidDensity,
        CELL_SIZE: cellSize,
        NX: nx,
        NY: ny,
        NZ: nz,
        FLUID: STATE_ENUM.FLUID,
      })
      .setOutput([nx, ny, nz]);

  const createVelocityYUpdateKernel = (
    gpu,
    nx,
    ny,
    nz,
    fluidDensity,
    cellSize
  ) =>
    gpu
      .createKernel(function (velocity, pressure, voxelStates, dt) {
        // only consider boundaries which have a fluid cell on at least one side
        if (
          this.thread.y === 0 ||
          this.thread.y === this.constants.NY - 1 ||
          !(
            voxelStates[this.thread.z][this.thread.y - 1][this.thread.x] ===
              this.constants.FLUID ||
            voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
              this.constants.FLUID
          )
        ) {
          return 0;
        }

        const pressureGradient =
          (pressure[this.thread.z][this.thread.y][this.thread.x] -
            pressure[this.thread.z][this.thread.y - 1][this.thread.x]) /
          this.constants.CELL_SIZE;

        const oldVelocity = velocity[this.thread.z][this.thread.y][this.thread.x];
        const newVelocity =
          oldVelocity - (dt * pressureGradient) / this.constants.FLUID_DENSITY;

        return newVelocity;
      })
      .setConstants({
        FLUID_DENSITY: fluidDensity,
        CELL_SIZE: cellSize,
        NX: nx,
        NY: ny,
        NZ: nz,
        FLUID: STATE_ENUM.FLUID,
      })
      .setOutput([nx, ny, nz]);

  const createVelocityZUpdateKernel = (
    gpu,
    nx,
    ny,
    nz,
    fluidDensity,
    cellSize
  ) =>
    gpu
      .createKernel(function (velocity, pressure, voxelStates, dt) {
        // only consider boundaries which have a fluid cell on at least one side
        if (
          this.thread.z === 0 ||
          this.thread.z === this.constants.NZ - 1 ||
          !(
            voxelStates[this.thread.z - 1][this.thread.y][this.thread.x] ===
              this.constants.FLUID ||
            voxelStates[this.thread.z][this.thread.y][this.thread.x] ===
              this.constants.FLUID
          )
        ) {
          return 0;
        }

        const pressureGradient =
          (pressure[this.thread.z][this.thread.y][this.thread.x] -
            pressure[this.thread.z - 1][this.thread.y][this.thread.x]) /
          this.constants.CELL_SIZE;

        return (
          velocity[this.thread.z][this.thread.y][this.thread.x] -
          (dt * pressureGradient) / this.constants.FLUID_DENSITY
        );
      })
      .setConstants({
        FLUID_DENSITY: fluidDensity,
        CELL_SIZE: cellSize,
        NX: nx,
        NY: ny,
        NZ: nz,
        FLUID: STATE_ENUM.FLUID,
      })
      .setOutput([nx, ny, nz]);

  const createJacobiIterationKernel = (gpu, nx, ny, nz) =>
    gpu
      .createKernel(function (negativeDivergence, pressure, voxelStates) {
        const i = this.thread.x;
        const j = this.thread.y;
        const k = this.thread.z;

        if (voxelStates[k][j][i] === this.constants.AIR) {
          return 0;
        }

        const divergenceCenter = negativeDivergence[k][j][i];

        let left = 0.0;
        let right = 0.0;
        let bottom = 0.0;
        let top = 0.0;
        let back = 0.0;
        let front = 0.0;

        // negative x neighbor
        if (voxelStates[k][j][i - 1] === this.constants.FLUID) {
          left = pressure[k][j][i - 1];
        }
        // positive x neighbor
        if (voxelStates[k][j][i + 1] === this.constants.FLUID) {
          right = pressure[k][j][i + 1];
        }
        // negative y neighbor
        if (voxelStates[k][j - 1][i] === this.constants.FLUID) {
          bottom = pressure[k][j - 1][i];
        }
        // positive y neighbor
        if (voxelStates[k][j + 1][i] === this.constants.FLUID) {
          top = pressure[k][j + 1][i];
        }
        // negative z neighbor
        if (voxelStates[k - 1][j][i] === this.constants.FLUID) {
          front = pressure[k - 1][j][i];
        }
        // positive z neighbor
        if (voxelStates[k + 1][j][i] === this.constants.FLUID) {
          back = pressure[k + 1][j][i];
        }

        return (
          (left + right + bottom + top + back + front + divergenceCenter) / 6.0
        );
      })
      .setTactic("precision") // vector math should be high precision
      .setConstants({
        FLUID: STATE_ENUM.FLUID,
        AIR: STATE_ENUM.AIR,
        SOLID: STATE_ENUM.SOLID,
      })
      .setOutput([nx, ny, nz]);

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

    // build negative divergence vector with boundary conditions
    const buildD = createNegativeDivergenceKernel(
      gpu,
      ...gridSize,
      grid.cellSize
    );

    const jacobi = createJacobiIterationKernel(gpu, ...gridSize);

    const pressureSolve = {
      buildD: buildD.setPipeline(true).setImmutable(true),
      jacobi: jacobi.setPipeline(true).setImmutable(true),
    };

    // update grid velocities using the pressure gradient
    const updateVelocityX = createVelocityXUpdateKernel(
      gpu,
      ...velocityXSize,
      FLUID_DENSITY,
      grid.cellSize
    );
    const updateVelocityY = createVelocityYUpdateKernel(
      gpu,
      ...velocityYSize,
      FLUID_DENSITY,
      grid.cellSize
    );
    const updateVelocityZ = createVelocityZUpdateKernel(
      gpu,
      ...velocityZSize,
      FLUID_DENSITY,
      grid.cellSize
    );

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
      particleToXGrid: particleToXGrid.setPipeline(true).setImmutable(true),
      particleToYGrid: particleToYGrid.setPipeline(true).setImmutable(true),
      particleToZGrid: particleToZGrid.setPipeline(true).setImmutable(true),
      copyPressure: copyPressure.setPipeline(true).setImmutable(true),
      copyXVelocity: copyXVelocity.setPipeline(true).setImmutable(true),
      copyYVelocity: copyYVelocity.setPipeline(true).setImmutable(true),
      copyZVelocity: copyZVelocity.setPipeline(true).setImmutable(true),
      classifyVoxels: classifyVoxels.setPipeline(true).setImmutable(true),
      addGravity: addGravity.setPipeline(true).setImmutable(true),
      enforceXBoundary: enforceXBoundary.setPipeline(true).setImmutable(true),
      enforceYBoundary: enforceYBoundary.setPipeline(true).setImmutable(true),
      enforceZBoundary: enforceZBoundary.setPipeline(true).setImmutable(true),
      gridToParticles: gridToParticles,
      advectParticles: advectParticles.setPipeline(true).setImmutable(true),
      pressureSolve: pressureSolve,
      updateVelocityX: updateVelocityX.setPipeline(true).setImmutable(true),
      updateVelocityY: updateVelocityY.setPipeline(true).setImmutable(true),
      updateVelocityZ: updateVelocityZ.setPipeline(true).setImmutable(true),
    };
  };

  const FLUID_DENSITY = 3.75;
  const SOLVER_ITERATION_LIMIT = 50;

  class Simulation {
    constructor(gpu, config) {
      this.particles = new Particles(
        config.particleDensity,
        config.particleBounds
      );
      this.grid = new MACGrid(config.gridBounds, 0.025);
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
      this.grid.velocityY = this.kernels.addGravity(
        this.grid.velocityY,
        dt,
        this.grid.voxelStates
      );

      // enforce boundary conditions
      this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
      this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
      this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

      // do the pressure solve with a zero divergence velocity field
      this.grid.pressure = solve(
        this.kernels.pressureSolve,
        this.grid.voxelStates,
        this.grid.velocityX,
        this.grid.velocityY,
        this.grid.velocityZ,
        SOLVER_ITERATION_LIMIT,
        this.grid.pressureOld
      );

      // update the velocity fields with the new pressure gradients
      this.grid.velocityX = this.kernels.updateVelocityX(
        this.grid.velocityX,
        this.grid.pressure,
        this.grid.voxelStates,
        dt
      );
      this.grid.velocityY = this.kernels.updateVelocityY(
        this.grid.velocityY,
        this.grid.pressure,
        this.grid.voxelStates,
        dt
      );
      this.grid.velocityZ = this.kernels.updateVelocityZ(
        this.grid.velocityZ,
        this.grid.pressure,
        this.grid.voxelStates,
        dt
      );

      // enforce boundary conditions
      this.grid.velocityX = this.kernels.enforceXBoundary(this.grid.velocityX);
      this.grid.velocityY = this.kernels.enforceYBoundary(this.grid.velocityY);
      this.grid.velocityZ = this.kernels.enforceZBoundary(this.grid.velocityZ);

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

  function RayMarchingEffect(resolution, density, fakeWater) {
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

    var currentlyPressedKeys = {};
    var eyePt = fromValues(0.5, 1.2, 2);
    var viewDir = fromValues(0.0, -0.9, -1.5);
    /** @global Rotation around up */
    var upRot = create();
    /** @global Rotation around right */
    var rightRot = create();
    /** @global Final change in rotation */
    var finalRot = create();
    /** @global Direction of up in world coordinates */
    var upDir = fromValues(0.0, 1.0, 0.0);
    /** @global Direction of right in world coordinates */
    var rightDir = fromValues(1.0, 0.0, 0.0);

    const gpu = new GPU();
    const sim = new Simulation(gpu, {
      particleDensity: density,
      particleBounds: {
        min: fromValues(0.15, 0.15, 0.15),
        max: fromValues(0.8, 0.75, 0.35),
      },
      gridBounds: {
        min: fromValues(0.1, 0.0, 0.1),
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
      .createKernel(function (field, size, coefficient) {
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

              // Weighted by e^(-r^2 / c)
              let w = Math.pow(
                2.71,
                (-1 * Math.sqrt(x_o * x_o + y_o * y_o + z_o * z_o)) / coefficient
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

    // Key down watcher
    this.handleKeyDown = function (event) {
      var k = event.key.toLowerCase();
      currentlyPressedKeys[k] = true;
      if (k == "arrowup" || k == "arrowleft" || k == "arrowright" || k == "arrowdown") {
        event.preventDefault();
      }
    };

    // Key up watcher
    this.handleKeyUp = function (event) {
      currentlyPressedKeys[event.key.toLowerCase()] = false;
    };

    // Free camera movement
    this.freecamMovement = function (deltaTime) {
      // Rotation speeds
      var rotDegPerSec = 120;
      var rotRadPerSec = rotDegPerSec * 3.14159 / 180;
      var viewDirScaled = fromValues(0.0, 0.0, 0.0);
      var rightDirScaled = fromValues(0.0, 0.0, 0.0);
      scale(viewDirScaled, viewDir, 2 * deltaTime);
      scale(rightDirScaled, rightDir, 2 * deltaTime);

      // Move eye relative to forward
      if (currentlyPressedKeys["w"]) {
        add(eyePt, eyePt, viewDirScaled);
      }
      if (currentlyPressedKeys["s"]) {
        sub(eyePt, eyePt, viewDirScaled);
      }
      if (currentlyPressedKeys["a"]) {
        sub(eyePt, eyePt, rightDirScaled);
      }
      if (currentlyPressedKeys["d"]) {
        add(eyePt, eyePt, rightDirScaled);
      }

      // Rotate eye freely
      var upAmt = 0;
      var rightAmt = 0;
      if (currentlyPressedKeys["arrowleft"]) {
        upAmt += rotRadPerSec * deltaTime;
      }
      if (currentlyPressedKeys["arrowright"]) {
        upAmt -= rotRadPerSec * deltaTime;
      }
      if (currentlyPressedKeys["arrowup"]) {
        rightAmt += rotRadPerSec * deltaTime;
      }
      if (currentlyPressedKeys["arrowdown"]) {
        rightAmt -= rotRadPerSec * deltaTime;
      }

      // Quaternion rotation for eye
      setAxisAngle(upRot, [0, 1, 0], upAmt);
      setAxisAngle(rightRot, rightDir, rightAmt);
      identity(finalRot);
      mul(finalRot, upRot, finalRot);
      mul(finalRot, rightRot, finalRot);
      transformQuat(viewDir, viewDir, finalRot);
      transformQuat(upDir, upDir, finalRot);
      transformQuat(rightDir, rightDir, finalRot);
    };

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

      // step the simulation forwards
      deltaTime = Math.min(deltaTime, 1 / 60);
      if (!fakeWater) {
        sim.step(deltaTime / 2);
      }

      var uniformsConst = {
        u_field: textures[0],
        time: localTime,
        eye: eyePt,
        forward: viewDir,
        automatic: !window.useCamera
      };

      this.freecamMovement(deltaTime);

      model.drawPrep(uniformsConst);

      {
        // Set firstDraw = false to only draw 1 frame
        let balls = [];
        let radius = window.radiusSlider.value / 100;

        if (fakeWater) {
          // Sine wave water
          let n = 30;

          for (let x = 0; x < n; ++x) {
            for (let z = 0; z < n; ++z) {
              let xp = (x + 0.5) / n;
              let zp = (z + 0.5) / n;
              let r = Math.sqrt((xp - 0.5) * (xp - 0.5) + (zp - 0.5) * (zp - 0.5));
              //let y = 0.1*((Math.sin(40 * r - 1.5*time) + 1) / 2) / Math.abs(10*(Math.max(r, 0.013))) + 0.05;
              let y = 0.3 * Math.pow(Math.cos(10 * r - 1 * curTime), 2) / Math.max(10 * r, 0.5) + 0.05;
              balls.push([xp, y, zp]);
            }
          }
        } else {
          for (let i = 0; i < sim.particles.particleBuffer.length; i += 6) {
            balls.push([
              sim.particles.particleBuffer[i],
              sim.particles.particleBuffer[i + 1],
              sim.particles.particleBuffer[i + 2],
            ]);
          }
        }

        // Swap comment to see with / without smoothing
        //field = fillField(balls, balls.length, size, radius).toArray();
        field = smooth(
          fillField(balls, balls.length, size, radius),
          size,
          Math.max(0.001, window.smoothSlider.value)
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
