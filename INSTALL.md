# Installation Guide

## Quick Installation

### 1. Clone Repository
```bash
git clone https://github.com/michaelstrojny1/DFT-Quantum-Dot-Dielectric-Constant-Estimator.git
cd DFT-Quantum-Dot-Dielectric-Constant-Estimator
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Test Installation
```bash
python demos/quick_demo.py
```

## Detailed Setup

### Prerequisites
- **Python 3.8+** (recommended: 3.9 or 3.10)
- **pip** package manager
- **Git** for repository cloning

### Core Dependencies
```bash
numpy>=1.21.0      # Numerical computations
matplotlib>=3.5.0  # Visualization and plotting
scikit-fem>=8.0.0  # Finite element method solver
```

### Optional Dependencies (for DFT)
```bash
ase>=3.22.0        # Atomic Simulation Environment
gpaw>=22.8.0       # Grid-based Projector Augmented Wave DFT
```

### Platform-Specific Notes

#### Windows
- Ensure Python is added to PATH
- Consider using Anaconda/Miniconda for easier dependency management
- GPAW may require additional setup for MPI support

#### Linux/macOS
- Standard pip installation should work
- For GPAW, you may need development headers:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install python3-dev libopenmpi-dev
  
  # macOS (with Homebrew)
  brew install open-mpi
  ```

### Virtual Environment (Recommended)
```bash
python -m venv dft_env
source dft_env/bin/activate  # Linux/macOS
# or
dft_env\Scripts\activate     # Windows

pip install -r requirements.txt
```

## Verification

### Test Core Functionality
```bash
# Quick demo (no DFT required)
python demos/quick_demo.py

# Check output files
ls quick_demo_output/
```

### Test DFT Functionality (if installed)
```bash
# Full DFT demo
python demos/demoDFT.py

# Check DFT outputs
ls demo_*.png epsr_results.json
```

## Project Structure
```
DFT-Quantum-Dot-Dielectric-Constant-Estimator/
├── scripts/           # Core computational modules
│   ├── DFT.py        # DFT calculations with ASE+GPAW
│   ├── FEA.py        # Finite element Poisson solver
│   └── run.py        # Main workflow orchestrator
├── demos/            # Example workflows and demonstrations
│   ├── quick_demo.py # Fast demo without DFT
│   └── demoDFT.py    # Full DFT workflow demo
├── data/             # Sample quantum dot structures
│   ├── Cd68_OPT.xyz  # Optimized CdSe quantum dot
│   └── dot_demo.xyz  # Demo structure
├── requirements.txt  # Python dependencies
├── README.md         # Comprehensive documentation
├── CHANGELOG.md      # Version history
└── INSTALL.md        # This installation guide
```

## Troubleshooting

### Common Issues

#### Import Errors
```bash
# If scikit-fem import fails
pip install --upgrade scikit-fem

# If matplotlib backend issues
export MPLBACKEND=Agg  # Linux/macOS
set MPLBACKEND=Agg     # Windows
```

#### GPAW Installation Issues
```bash
# Try conda installation instead
conda install -c conda-forge gpaw

# Or build from source (advanced)
pip install gpaw --no-binary gpaw
```

#### Memory Issues
- Reduce mesh resolution (`--nx`, `--ny`, `--nz` parameters)
- Use smaller simulation domains
- Close other applications to free memory

### Getting Help
- Check the repository's issue tracker
- Ensure all dependencies meet minimum version requirements
- Verify Python version compatibility (3.8+)
