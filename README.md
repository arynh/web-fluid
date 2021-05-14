# web-fluid

A web-based, 3D, real-time fluid simulation.

Final Project for CS 419, Production Computer Graphics, at UIUC.

Authors:

- Idan Raiter
- Curie Hong
- Rishi Pandey
- Aryn Harmon

## Outline of Simulation

For every simulation loop:

- [x] clear the grid velocities
- [x] project particle velocities to the grid
- [x] copy the grid to store the version from the previous time step
- [x] mark cells as solid, fluid, or air
- [x] perform gravity update
- [x] enforce boundary conditions
- [x] do the pressure solve and generate new grid velocities
- [x] enforce boundary conditions
- [x] update the velocities of the particles
- [x] update the positions of the particles

## Installation

To install this project, clone the repository and run `yarn` in the root directory. Then you can run `yarn dev` to start the local development server.
