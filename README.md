# DFT Quantum Dot Dielectric Constant Estimator

Simulates quantum dot arrays in capacitor stacks to estimate dielectric constants via DFT and finite element analysis.

## Stack Structure
1. hBN stack (bottom)
2. 2D quantum-dot array with radius-dependent εᵣ
3. hBN spacer
4. Graphene readout plane
5. hBN stack (top)
6. Capacitor plates at z=0 and z=Lz

---

## Core Scripts
* `scripts/DFT.py` - DFT calculations using ASE+GPAW
* `scripts/FEA.py` - Poisson solver using scikit-fem
* `scripts/run.py` - Main workflow
* `demos/quick_demo.py` - Fast demo without DFT
* `demos/demoDFT.py` - DFT demo

---

## Quick Start

```bash
python demos/quick_demo.py
```

Outputs: CSV data, summary stats, 2D/3D potential plots

---

## DFT Workflow
```bash
python demos/demoDFT.py
```
Requires ASE+GPAW. Uses real quantum dot structures from `data/` directory.

---

## Manual Usage
```bash
python scripts/run.py --epsr-profile "14,12,10,8,6" --pattern hexagonal --qd-rows 3 --qd-cols 4
```

## Dependencies
```
numpy
matplotlib
ase
gpaw
scikit-fem
```

## Installation
```bash
pip install -r requirements.txt
```
