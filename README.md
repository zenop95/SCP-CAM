# Sequential Convex Program for Collision Avoidance Manoeuvres (SCP CAM)

A complete software framework for computing collision avoidance manoeuvres (CAMs) using Sequential Convex Programming (SCP). Handles single and multiple secondary objects in both short- and long-term conjunction scenarios across LEO and GEO orbits.

## Overview

SCP CAM uses successive convexification to optimize trajectories that avoid collisions while minimizing fuel consumption. The system supports:

- **Single and multiple secondary objects**
- **LEO and GEO orbits**
- **Short-term and long-term conjunctions**
- **Station-keeping and operational constraints**
- **Advanced environment models** (atmosphere, SRP, third-body perturbations)

## System Requirements

- **MATLAB**: R2022a or newer (tested up to R2025b)
- **MOSEK**: 10.0 (optimization solver)
- **WSL (Ubuntu)** with the following:
  - DACE (Differential Algebra Computational Engine)
  - CSPICE (Navigation Toolkit)
  - Astrotools (Astrodynamics library)
  - nlohmann-json3-dev
  - libeigen3-dev, libjsoncpp-dev, libdlib-dev

## Installation

### 1. Install Dependencies

Run the dependency checker script on Windows:

```bash
check_deps.bat
```

This script automatically:
- Detects MOSEK installation 
- Detects MATLAB installation 
- Detects WSL installation
- Installs/verifies DACE in WSL
- Installs/verifies CSPICE
- Installs/verifies Astrotools and pre-requisites
- Installs/verifies JSON

**Note:** Ensure you have a valid MOSEK license before proceeding.

### 2. Build the Project

Run the build script:

```bash
build_scp.bat
```

This compiles the C++ backend using CMake within WSL and generates the necessary binaries.

## Quick Start

1. **Configure the scenario** in `userParams.m`:
   ```matlab
   up.orbit      = 'leo';              % or 'geo'
   up.scenario   = 'single-short';     % scenario type
   up.uMax       = 1e-8;               % max thrust [km/s²]
   up.iterMaxMaj = 10;                 % major iterations
   ```

2. **Run the main script**:
   ```matlab
   main
   ```

3. **View results** via generated figures in postprocessing.

## Configuration Files

### `userParams.m`
User-facing parameters for scenario definition:
- Orbit selection (LEO/GEO)
- Scenario type (single/multiple, short/long-term)
- Discretization and time windows
- Propulsion limits
- Optimizer stopping criteria
- Environmental model flags

### `defaultParams.m`
Internal defaults for:
- Orbital mechanics (GM constant)
- Objective weights
- NLI (nonlinear inequality) limits
- SMD gradient constraints
- Operational toggles

### `scenarioParams.m`
Builds the unscaled scenario by merging user and default parameters, then:
- Generates initial orbits
- Sets up discretization and time grids
- Configures collision avoidance windows

### `scalingParams.m`
Applies nondimensionalization and derives scaled limits.

### `simProperties.m`
Master configuration function that orchestrates all parameter assembly and consistency checks.

## High-Level Algorithm Flow

1. **Initialize**: Set up MATLAB path and load toolboxes
2. **Configure**: Assemble full configuration via `simProperties()`
3. **Optimize**: Run `runSCP()` for major-iteration loop with successive convexification
4. **Validate**: Check solution feasibility and constraint satisfaction
5. **Postprocess**: Generate visualizations and reports via `mainPostProcess()`

## Output

- **pp.majorIter**: Per-iteration data (trajectories, costs, convergence metrics)
- **Figures**: Trajectory plots, risk profiles, state evolution (via `selectFigs()`)
- **Logs**: Installation and build logs for troubleshooting

## Troubleshooting

- **MOSEK not found**: Install MOSEK or verify license file location
- **WSL user detection failed**: Ensure Ubuntu WSL distribution is installed and accessible
- **DACE/CSPICE build errors**: Check `install_log.txt` for detailed error messages
- **Astrotools missing**: Verify source repository is available and re-run `check_deps.bat`

## Author

**Zeno Pavanello**  
E-mail: zeno.pavanello@polimi.it  
Date: 2022–2026

## Compatibility

- MATLAB R2022a–R2025b
- MOSEK 10.0
- Windows 10/11 with WSL2 (Ubuntu)
- GCC 64-bit (Linux backend)

## Credits
If you use this work, please cite the following papers:

-Pavanello, Zeno, Laura Pirovano, and Roberto Armellin. ‘Long-Term Fuel-Optimal Collision Avoidance Maneuvers with Station-Keeping Constraints’. Journal of Guidance, Control, and Dynamics 47, no. 9 (2024): 1855–71. https://doi.org/10.2514/1.G007839.

-Pavanello, Zeno, Laura Pirovano, Roberto Armellin, Andrea De Vittori, and Pierluigi Di Lizia. ‘Collision Avoidance Maneuvers Optimization in the Presence of Multiple Encounters’. Journal of Guidance, Control, and Dynamics 49, no. 1 (2026): 130–48. https://doi.org/10.2514/1.G008332.

