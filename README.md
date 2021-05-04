# web-fluid

A web-based, 3D, real-time fluid simulation.

## Outline of Simulation

For every simulation loop:

- clear the grid velocities
- project particle velocities to the grid
- copy the grid to store the version from the previous time step
- mark cells, not entirely clear what this step does
- perform gravity update
- enforce boundary conditions
- do the pressure solve
- extrapolate velocity
- enforce boundary conditions
- update the velocities of the particles
- update the positions of the particles

## Installation

To install this project, clone the repository and run `yarn` in the root directory. Then you can run `yarn start` to start the local development server.
