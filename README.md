# Molecular Dynamics Simulation

A Fortran-based molecular dynamics simulation program for simulating Lennard-Jones particles, with Python-based analysis tools.

## Project Structure

```
.
├── bin/           # Compiled executables
├── obj/           # Object files and module files
├── output/        # Simulation output files
├── src/           # Source files
│   ├── parameters.f90    # Simulation parameters and constants
│   ├── particles.f90     # Particle system and RDF calculations
│   ├── forces.f90        # Force calculations
│   ├── integrator.f90    # Time integration
│   ├── io.f90           # Input/output operations
│   └── main.f90         # Main program
├── analyze_results.py    # Analysis and visualization script
├── Makefile             # Build system
└── README.md            # This file
```

## Features

- Lennard-Jones potential for particle interactions
- Velocity Verlet integration
- Optional Berendsen thermostat for temperature control
- Neighbor list optimization for force calculations
- Periodic boundary conditions
- Radial Distribution Function (RDF) calculation
- Real-time energy, temperature, and pressure monitoring
- Trajectory output in XYZ format
- Interactive visualization and analysis tools

## Requirements

### Compilation
- Fortran compiler (gfortran recommended)
- Make build system

### Analysis
- Python 3.x
- NumPy
- Matplotlib

## Building the Project

To compile the simulation program:

```bash
make
```

This will create the executable `bin/md_sim.exe` and all necessary object files in the `obj/` directory.

## Running the Simulation

To run the simulation:

```bash
make run
```

This will execute the simulation and store all output files in the `output/` directory:
- `simulation_info.txt`: Contains all simulation parameters
- `energies.dat`: Time evolution of energies, temperature, and pressure
- `trajectory.xyz`: Particle positions in XYZ format
- `rdf_*.dat`: Radial distribution function at different time steps

## Analysis

To analyze the simulation results:

```bash
python analyze_results.py
```

This will:
1. Display statistical analysis of the simulation
2. Show interactive plots of:
   - Energy evolution (potential, kinetic, total)
   - Temperature evolution
   - Pressure evolution
   - Energy conservation
   - Radial Distribution Function (initial vs final)

## Output Files

### simulation_info.txt
Contains all simulation parameters including:
- Number of particles
- Box size
- Temperature
- Time step
- Number of steps
- Lennard-Jones parameters
- Thermostat settings
- Neighbor list parameters
- RDF calculation parameters

### energies.dat
Tab-separated file with columns:
1. Time step
2. Potential energy
3. Kinetic energy
4. Total energy
5. Temperature
6. Pressure

### trajectory.xyz
Standard XYZ format containing:
- Number of particles
- Comment line with step number, total energy, temperature, and pressure
- Particle positions (x, y, z)

### rdf_*.dat
Radial distribution function data with columns:
1. Distance (r)
2. g(r)

## Cleaning Up

To clean all compiled files and simulation outputs:

```bash
make clean
```

This will remove:
- All object files (`obj/*.o`)
- All module files (`obj/*.mod`)
- The executable (`bin/md_sim.exe`)
- All output files (`output/*`)

## Notes

- The simulation uses reduced units where σ = ε = m = 1
- Periodic boundary conditions are applied in all directions
- The neighbor list is updated based on particle displacements
- The RDF is calculated and output at specified intervals
- All output files are automatically saved in the `output/` directory 