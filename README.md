# DFT Quantum Dot Dielectric Constant Estimator and FEA runner for Bai Lab Experiment (Top Capacitor Plate -> HBN -> 2D QD Pattern -> HBN -> Graphene -> HBN -> Bottom Plate)

A comprehensive simulation framework for quantum dot arrays in capacitor stacks to estimate dielectric constants via Density Functional Theory (DFT) and finite element analysis (FEA).

## Overview

This project simulates the electrostatic behavior of quantum dot arrays embedded in heterostructured capacitor stacks. The workflow combines:
- **DFT calculations** using ASE+GPAW for accurate quantum dot polarizability
- **Finite Element Analysis** using scikit-fem for 3D Poisson equation solving
- **Visualization tools** for potential mapping and field analysis

## Stack Architecture

The simulated heterostructure consists of (bottom to top):
1. **hBN bottom layer** - Dielectric substrate
2. **2D quantum dot array** - CdSe dots with radius-dependent εᵣ profile
3. **hBN spacer layer** - Isolation between dots and readout
4. **Graphene readout plane** - 2D conductor for potential sensing
5. **hBN top layer** - Top dielectric
6. **Capacitor plates** - Boundary conditions at z=0 and z=Lz

## Core Components

### Scripts
- **`scripts/DFT.py`** - DFT calculations using ASE+GPAW for quantum dot polarizability
- **`scripts/FEA.py`** - 3D Poisson solver using scikit-fem with adaptive mesh
- **`scripts/run.py`** - Main workflow orchestrator with parameter management

### Demos
- **`demos/quick_demo.py`** - Fast demonstration without DFT (uses predefined εᵣ profile)
- **`demos/demoDFT.py`** - Full DFT workflow demonstration with visualization

### Data
- **`data/`** - Sample quantum dot structures (CdSe, etc.) in XYZ format

## Quick Start

### Fast Demo (No DFT Required)
```bash
python demos/quick_demo.py
```
**Outputs:**
- `graphene_potential_2D.png` - 2D potential map with quantum dot overlay
- `graphene_potential_3D.png` - 3D surface plot of potential
- `graphene_slice.csv` - Raw potential data at graphene plane
- `summary.txt` - Field statistics and simulation parameters

### Full DFT Workflow
```bash
python demos/demoDFT.py
```
**Requirements:** ASE+GPAW installed
**Uses:** Real quantum dot structures from `data/` directory
**Additional Outputs:**
- `demo_radial_profile.png` - εᵣ vs radius plot
- `demo_2D_profile.png` - 2D dielectric cross-section

## Advanced Usage

### Manual Parameter Control
```bash
python scripts/run.py \
  --epsr-profile "14,12,10,8,6,4" \
  --pattern hexagonal \
  --qd-rows 3 --qd-cols 4 \
  --qd-radius-nm 2.5 \
  --spacing-factor 3.0 \
  --Lx-nm 30 --Ly-nm 30 --Lz-nm 20 \
  --nx 64 --ny 64 --nz 48 \
  --E-Vnm 0.15 \
  --workdir custom_output
```

### DFT-Based εᵣ Profile Generation
```bash
python scripts/run.py \
  --xyz data/Cd68_OPT.xyz \
  --dft-shells 6 \
  --dft-E-Vnm 0.01 \
  --pattern hexagonal \
  --qd-rows 4 --qd-cols 5
```

## Key Parameters

### Geometry
- `--Lx-nm`, `--Ly-nm`, `--Lz-nm` - Simulation box dimensions
- `--qd-radius-nm` - Quantum dot radius
- `--spacing-factor` - Dot spacing as multiple of diameter
- `--pattern` - Array pattern: `single`, `grid`, `hexagonal`

### Physics
- `--epsr-profile` - Comma-separated radial εᵣ values (center to edge)
- `--E-Vnm` - Applied electric field strength
- `--eps-hbn` - hBN dielectric constant (default: 6.93)
- `--lateral-bc` - Boundary conditions: `natural` or `cutoff`

### Numerical
- `--nx`, `--ny`, `--nz` - Mesh resolution
- `--dft-shells` - Number of radial shells for DFT analysis

## Installation

### Prerequisites
- Python 3.8+
- Git

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Optional: DFT Support
For full DFT functionality, ensure GPAW is properly configured:
```bash
# Additional GPAW setup may be required depending on your system
pip install gpaw
```

## Output Files

### Standard Outputs
- **`graphene_slice.csv`** - Potential values at graphene plane (x_nm, y_nm, phi_V)
- **`summary.txt`** - Simulation summary with min/max potentials and parameters
- **`graphene_potential_2D.png`** - 2D heatmap with quantum dot positions
- **`graphene_potential_3D.png`** - 3D surface visualization

### DFT Outputs
- **`epsr_results.jsonl`** - Detailed DFT results including polarizability
- **`demo_radial_profile.png`** - Radial dielectric profile
- **`demo_2D_profile.png`** - 2D cross-sectional view

## System Requirements

### Dependencies
```
numpy>=1.21.0      # Numerical computations
matplotlib>=3.5.0  # Visualization
scikit-fem>=8.0.0  # Finite element solver
ase>=3.22.0        # Atomic simulation environment
gpaw>=22.8.0       # DFT calculations (optional)
```

### Compatibility
- **Operating Systems:** Windows, Linux, macOS
- **Python:** 3.8+
- **Memory:** 4GB+ recommended for larger simulations
- **Storage:** ~100MB for installation + output space

## Contributing

This project is actively maintained. For questions, issues, or contributions, please refer to the repository's issue tracker.

## License

See LICENSE file for details.
