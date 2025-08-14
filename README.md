# LAMMPS Extension: `fix_wall_reflect_cylinder`

## Overview
This extension adds a new wall reflection fix for **cylindrical geometries** in [LAMMPS](https://www.lammps.org/).  
It enables particles to **elastically reflect** from an infinite cylindrical boundary, preserving kinetic energy and momentum perpendicular to the wall.  

Potential use cases include:
- Confined flows in cylindrical channels
- Nanotube simulations

---

## Motivation
While LAMMPS provides existing wall fixes for planes and spheres, a direct implementation for cylindrical reflections was missing.  
This feature fills that gap by implementing a customizable cylindrical wall model with adjustable radius and center position.

---

## Files Added
- `fix_wall_reflect_cylinder.h` — class declaration, parameter parsing, and LAMMPS integration hooks.
- `fix_wall_reflect_cylinder.cpp` — implementation of the reflection algorithm and force application logic.

---

## Usage
Add the fix in your LAMMPS input script:

```
fix <name> <group> wall/reflectcylinder <radius> <center_y> <center_z>
```

## Example
```
fix wall all wall/reflectcylinder 15.0 7.5 7.5
```

This creates a cylindrical wall of radius 15.0 units, centered at (y=7.5, z=7.5), affecting all atoms in the simulation.

## Installation

1. Clone this fork or copy the two new files into your local src/ directory of LAMMPS:
```
git clone https://github.com/<your-username>/lammps.git
cd lammps/src
```
2. Rebuild LAMMPS as usual:
```
mkdir build
cd build
cmake ../cmake -D BUILD_MPI=ON
```
3. Build LAMMPS:
cmake --build . -j 8 

## Example Script
A minimal example will be provided in:
```
examples/wall_reflect_cylinder/reflect_dpd.lmps
```
This script demonstrates particle reflection inside a cylindrical boundary using Dissipative Particle Dynamics (DPD).

## Notes
- Developed and tested with LAMMPS version X.Y.Z.
- Not yet merged into upstream LAMMPS (standalone feature branch).

## ToDo
- [ ] Add example script and visualization
- [ ] Add test to verify reflection behavior
- [ ] Benchmark performance impact for large system