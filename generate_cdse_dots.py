#!/usr/bin/env python3
"""
Generate CdSe quantum dot structures for DFT calculations
Creates different sized CdSe quantum dots with wurtzite structure
"""

import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import bulk
from ase.cluster import wulff_construction
import os

def create_cdse_dot(diameter_nm, structure='wurtzite'):
    """
    Create a CdSe quantum dot with specified diameter
    
    Parameters:
    diameter_nm: Diameter in nanometers
    structure: 'wurtzite' or 'zincblende'
    """
    
    # CdSe lattice parameters (wurtzite)
    if structure == 'wurtzite':
        # Wurtzite CdSe: a=4.30 Å, c=7.01 Å
        a = 4.30
        c = 7.01
        cdse_bulk = bulk('CdSe', 'wurtzite', a=a, c=c)
    else:
        # Zincblende CdSe: a=6.08 Å
        a = 6.08
        cdse_bulk = bulk('CdSe', 'zincblende', a=a)
    
    # Create a spherical cluster
    radius_angstrom = diameter_nm * 10 / 2  # Convert nm to Angstrom and get radius
    
    # Build a simple spherical quantum dot by replicating and cutting
    # First, create a larger supercell
    n_cells = int(radius_angstrom / a) + 2
    supercell = cdse_bulk.repeat((n_cells, n_cells, n_cells))
    
    # Get center of mass
    center = supercell.get_center_of_mass()
    
    # Select atoms within radius
    positions = supercell.get_positions()
    distances = np.linalg.norm(positions - center, axis=1)
    mask = distances <= radius_angstrom
    
    # Create the quantum dot
    dot_positions = positions[mask]
    dot_symbols = [supercell.get_chemical_symbols()[i] for i in range(len(supercell)) if mask[i]]
    
    quantum_dot = Atoms(symbols=dot_symbols, positions=dot_positions)
    
    # Center the quantum dot
    quantum_dot.center()
    
    return quantum_dot

def create_small_cdse_cluster(n_atoms=68):
    """
    Create a small CdSe cluster with approximately n_atoms
    Tries to maintain Cd:Se ratio close to 1:1
    """
    # Start with a small wurtzite unit cell
    a = 4.30
    c = 7.01
    cdse_bulk = bulk('CdSe', 'wurtzite', a=a, c=c)
    
    # Estimate supercell size needed
    atoms_per_cell = len(cdse_bulk)
    n_cells_estimate = int(np.cbrt(n_atoms / atoms_per_cell)) + 1
    
    # Create supercell
    supercell = cdse_bulk.repeat((n_cells_estimate, n_cells_estimate, n_cells_estimate))
    
    # Create spherical cut
    center = supercell.get_center_of_mass()
    positions = supercell.get_positions()
    distances = np.linalg.norm(positions - center, axis=1)
    
    # Sort by distance and take closest n_atoms
    sorted_indices = np.argsort(distances)
    
    # Try to maintain stoichiometry
    selected_indices = []
    n_cd = 0
    n_se = 0
    target_each = n_atoms // 2
    
    for idx in sorted_indices:
        symbol = supercell.get_chemical_symbols()[idx]
        if symbol == 'Cd' and n_cd < target_each:
            selected_indices.append(idx)
            n_cd += 1
        elif symbol == 'Se' and n_se < target_each:
            selected_indices.append(idx)
            n_se += 1
        
        if len(selected_indices) >= n_atoms:
            break
    
    # Create the cluster
    cluster_positions = positions[selected_indices]
    cluster_symbols = [supercell.get_chemical_symbols()[i] for i in selected_indices]
    
    cluster = Atoms(symbols=cluster_symbols, positions=cluster_positions)
    cluster.center()
    
    return cluster

def main():
    """Generate multiple CdSe quantum dots of different sizes"""
    
    # Create data directory if it doesn't exist
    data_dir = "data"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    # Generate quantum dots of different sizes
    dot_specs = [
        ("CdSe_small", 1.5),   # 1.5 nm diameter
        ("CdSe_medium", 2.5),   # 2.5 nm diameter
        ("CdSe_large", 3.5),    # 3.5 nm diameter
    ]
    
    print("Generating CdSe quantum dots...")
    
    for name, diameter in dot_specs:
        dot = create_cdse_dot(diameter)
        filename = os.path.join(data_dir, f"{name}.xyz")
        
        # Write with energy comment (placeholder)
        write(filename, dot, format='xyz', comment=f"CdSe quantum dot, diameter={diameter} nm")
        
        n_atoms = len(dot)
        n_cd = sum(1 for s in dot.get_chemical_symbols() if s == 'Cd')
        n_se = sum(1 for s in dot.get_chemical_symbols() if s == 'Se')
        
        print(f"Created {name}: {n_atoms} atoms (Cd: {n_cd}, Se: {n_se}), diameter: {diameter} nm")
        print(f"  Saved to: {filename}")
    
    # Also create a smaller cluster similar to the existing Cd68 file
    print("\nGenerating small CdSe clusters...")
    
    cluster_sizes = [34, 68, 102]  # Different cluster sizes
    
    for size in cluster_sizes:
        cluster = create_small_cdse_cluster(size)
        filename = os.path.join(data_dir, f"CdSe_{size}.xyz")
        
        write(filename, cluster, format='xyz', comment=f"CdSe cluster with ~{size} atoms")
        
        n_atoms = len(cluster)
        n_cd = sum(1 for s in cluster.get_chemical_symbols() if s == 'Cd')
        n_se = sum(1 for s in cluster.get_chemical_symbols() if s == 'Se')
        
        print(f"Created CdSe_{size}: {n_atoms} atoms (Cd: {n_cd}, Se: {n_se})")
        print(f"  Saved to: {filename}")
    
    print("\nAll CdSe quantum dot structures have been generated!")
    print("Files saved in the 'data' directory")

if __name__ == "__main__":
    main()
