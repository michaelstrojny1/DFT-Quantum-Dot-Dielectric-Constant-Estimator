# Changelog

## v2.1 - Enhanced Workflow Release (2025-09-12)

### New Features
- **Enhanced visualization pipeline** with 2D/3D potential plotting
- **Improved parameter management** with comprehensive CLI options
- **Hexagonal array patterns** with configurable geometry
- **Radial dielectric profiling** for quantum dots
- **Advanced boundary condition handling** (natural/cutoff)
- **Automated output management** with structured file organization

### Core Improvements
- **Streamlined DFT workflow** with optimized shell-based εᵣ extraction
- **Enhanced FEA solver** with better mesh handling and convergence
- **Comprehensive demo suite** including quick_demo.py and demoDFT.py
- **Robust error handling** and validation throughout pipeline
- **Performance optimizations** for larger simulation domains

### Technical Updates
- Updated scikit-fem integration with modern API
- Improved ASE+GPAW compatibility for latest versions
- Enhanced matplotlib visualization with publication-quality plots
- Better memory management for large-scale simulations
- Comprehensive parameter validation and error reporting

### Documentation
- **Complete README overhaul** with detailed usage examples
- **Comprehensive parameter documentation** with physics explanations
- **Installation guide improvements** with dependency management
- **Output file documentation** with format specifications

### Dependencies
- numpy>=1.21.0 (numerical computations)
- matplotlib>=3.5.0 (visualization)
- scikit-fem>=8.0.0 (finite element solver)
- ase>=3.22.0 (atomic simulation environment)
- gpaw>=22.8.0 (DFT calculations, optional)

### Compatibility
- Python 3.8+
- Windows/Linux/macOS
- Memory: 4GB+ recommended
- Storage: ~100MB + output space

---

## v2.0 - Production Release

### Core Features
- DFT calculations using ASE+GPAW
- Finite element Poisson solver
- Quantum dot array simulation
- Real CdSe structures included

### Dependencies
- numpy, matplotlib, scikit-fem, ase, gpaw

### Compatibility
- Python 3.8+
- Windows/Linux/macOS
