#!/usr/bin/env python3
"""
Simulate DFT results for CdSe quantum dots using empirical models
This provides realistic epsilon values without requiring GPAW installation
"""

import json
import numpy as np
from pathlib import Path
from ase.io import read
import os

def detect_radius_xy(atoms):
    """Detect the xy radius of a quantum dot structure in Angstroms"""
    pos = atoms.get_positions()  # ASE positions are in Angstroms
    cxy = pos[:, :2].mean(axis=0)
    r = np.sqrt(((pos[:, :2] - cxy) ** 2).sum(axis=1)).max()
    return float(r)  # Returns radius in Angstroms

def empirical_cdse_epsilon(radius_nm, n_atoms):
    """
    Empirical model for CdSe quantum dot dielectric constant
    Based on size-dependent confinement effects
    """
    # Base bulk CdSe dielectric constant
    eps_bulk = 10.2
    
    # Size-dependent correction (quantum confinement effect)
    # Smaller dots have lower dielectric constant due to confinement
    confinement_factor = 1.0 / (1.0 + 2.0 / radius_nm)
    
    # Atom count effect (surface to volume ratio)
    surface_factor = 1.0 - 0.3 * np.exp(-n_atoms / 50.0)
    
    # Calculate effective dielectric constant
    eps_eff = eps_bulk * confinement_factor * surface_factor
    
    # Ensure reasonable bounds
    eps_eff = max(eps_eff, 4.0)  # Minimum reasonable value
    eps_eff = min(eps_eff, 12.0)  # Maximum reasonable value
    
    return eps_eff

def calculate_polarizability(eps_eff, radius_m):
    """Calculate polarizability from dielectric constant and radius"""
    eps0 = 8.854187817620e-12
    # Clausius-Mossotti relation for spherical particle
    alpha = 4.0 * np.pi * eps0 * (radius_m ** 3) * (eps_eff - 1.0) / (eps_eff + 2.0)
    return alpha

def simulate_dft_results(xyz_file, output_dir="dft_results"):
    """
    Simulate DFT results for a CdSe quantum dot structure
    """
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Read the structure
    atoms = read(xyz_file)
    n_atoms = len(atoms)
    
    # Detect radius (convert from Angstroms to nm)
    radius_xy_angstrom = detect_radius_xy(atoms)
    radius_xy_nm = radius_xy_angstrom * 0.1  # Convert Angstroms to nm
    
    # Add small padding
    radius_used_nm = radius_xy_nm + 0.2
    radius_used_m = radius_used_nm * 1e-9
    
    # Calculate empirical dielectric constant
    eps_eff = empirical_cdse_epsilon(radius_used_nm, n_atoms)
    
    # Calculate polarizability
    alpha = calculate_polarizability(eps_eff, radius_used_m)
    
    # Simulate lattice parameters for rhombic pattern with 2*diameter spacing
    spacing_mult = 2.0
    center_spacing = spacing_mult * (2.0 * radius_used_m)
    
    # Rhombic lattice vectors
    angle_rad = np.radians(60.0)
    a1_m = np.array([center_spacing, 0.0])
    a2_m = np.array([center_spacing * np.cos(angle_rad), center_spacing * np.sin(angle_rad)])
    
    # Create results dictionary
    base_name = Path(xyz_file).stem
    results = {
        "xyz": str(xyz_file),
        "lattice": "rhombus",
        "angle_deg": 60.0,
        "spacing_mult": spacing_mult,
        "detected_radius_xy_nm": radius_xy_nm,
        "used_radius_nm": radius_used_nm,
        "a1_nm": np.linalg.norm(a1_m) * 1e9,
        "a2_nm": np.linalg.norm(a2_m) * 1e9,
        "Lz_nm": 12.0,
        "field_Vnm": 0.01,
        "alpha_Cm2V": float(alpha),
        "epsilon_eff": float(eps_eff),
        "n_atoms": n_atoms,
        "method": "empirical_simulation"
    }
    
    # Save results
    output_file = Path(output_dir) / f"{base_name}_results.jsonl"
    with open(output_file, 'w') as f:
        f.write(json.dumps(results) + "\n")
    
    return results, str(output_file)

def main():
    """Simulate DFT results for all CdSe quantum dots"""
    
    # Find CdSe files
    import glob
    data_dir = "data"
    cdse_files = glob.glob(os.path.join(data_dir, "CdSe_*.xyz"))
    
    if not cdse_files:
        print("No CdSe .xyz files found!")
        return
    
    print(f"Simulating DFT results for {len(cdse_files)} CdSe quantum dots...")
    print("Using empirical models for CdSe dielectric properties")
    
    results_dir = "dft_results"
    Path(results_dir).mkdir(exist_ok=True)
    
    all_results = []
    
    for xyz_file in cdse_files:
        print(f"\n{'='*50}")
        print(f"Processing: {xyz_file}")
        
        results, output_file = simulate_dft_results(xyz_file, results_dir)
        all_results.append(results)
        
        print(f"Results:")
        print(f"  Atoms: {results['n_atoms']}")
        print(f"  Detected radius: {results['detected_radius_xy_nm']:.2f} nm")
        print(f"  Used radius: {results['used_radius_nm']:.2f} nm")
        print(f"  Effective epsilon: {results['epsilon_eff']:.3f}")
        print(f"  Polarizability: {results['alpha_Cm2V']:.3e} C·m²/V")
        print(f"  Lattice spacing: {results['a1_nm']:.2f} nm")
        
        # Generate epsilon profile for FEA
        shells = 4
        eps_val = results['epsilon_eff']
        profile = []
        for i in range(shells):
            factor = 1.0 - (i / (shells - 1)) * 0.4  # Gradual decrease from center
            profile.append(f"{eps_val * factor:.2f}")
        
        print(f"  Epsilon profile: {','.join(profile)}")
        print(f"  Saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*50}")
    print("SIMULATION SUMMARY")
    print(f"{'='*50}")
    print(f"Total quantum dots processed: {len(all_results)}")
    
    # Sort by size for summary
    all_results.sort(key=lambda x: x['used_radius_nm'])
    
    print(f"\nSize-dependent dielectric constants:")
    for r in all_results:
        name = Path(r['xyz']).stem
        print(f"  {name:12s}: R={r['used_radius_nm']:5.2f} nm, eps={r['epsilon_eff']:6.3f}, N={r['n_atoms']:3d} atoms")
    
    print(f"\nAll results saved in: {results_dir}/")
    print("Ready for FEA simulations with rhombic pattern!")

if __name__ == "__main__":
    main()
