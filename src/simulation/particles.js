export const ATTRIBUTE_COUNT = 6;

/**
 * Represent the particle cloud.
 */
export class Particles {
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
