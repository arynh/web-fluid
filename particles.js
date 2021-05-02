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
   * @param {{min: vec3, max: vec3}} bounds The minimum and maximum extent of
   *    the box of particles.
   */
  constructor(density, bounds) {
    this.particleBuffer = [];
    let gap_between = 1 / Math.cbrt(density);
    for (let x = bounds.min.x; x < bounds.max.x; x += gap_between) {
      for (let y = bounds.min.y; y < bounds.max.y; y += gap_between) {
        for (let z = bounds.min.z; z < bounds.max.z; z += gap_between) {
          this.particleBuffer.push(x); // initial position
          this.particleBuffer.push(y);
          this.particleBuffer.push(z);
          this.particleBuffer.push(0); // initial velocity
          this.particleBuffer.push(0);
          this.particleBuffer.push(0);
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
