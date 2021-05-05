# web-fluid

A web-based, 3D, real-time fluid simulation.

## Outline of Simulation

For every simulation loop:

- [x] clear the grid velocities
- [x] project particle velocities to the grid
- [x] copy the grid to store the version from the previous time step
- [x] mark cells as solid, fluid, or air
- [x] perform gravity update
- [ ] enforce boundary conditions
- [ ] do the pressure solve
- [ ] extrapolate velocity
- [ ] enforce boundary conditions
- [ ] update the velocities of the particles
- [ ] update the positions of the particles

## Installation

To install this project, clone the repository and run `yarn` in the root directory. Then you can run `yarn start` to start the local development server.
