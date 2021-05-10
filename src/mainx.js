import { vec3 } from "gl-matrix";
import { Simulation } from "./simulation/simulation.js";

// create GPU
const gpu = new GPU();
const sim = new Simulation(gpu, {
  particleDensity: 1,
  particleBounds: {
    min: vec3.fromValues(0.5, 0.5, 0.5),
    max: vec3.fromValues(1, 1, 1),
  },
  gridBounds: {
    min: vec3.fromValues(0, 0, 0),
    max: vec3.fromValues(1, 1, 1),
  },
});
console.log(sim.particles.get(0));
sim.step(1 / 10);
console.log(sim.particles.get(0));
