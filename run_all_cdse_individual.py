#!/usr/bin/env python3
"""
Run FEA simulations for all CdSe quantum dots individually with monitoring
"""

import subprocess
import sys
import json
import glob
import os
from pathlib import Path
import time

def run_individual_fea(dot_name):
    """Run FEA for a single quantum dot with monitoring"""
    
    print(f"\n{'='*50}")
    print(f"Processing {dot_name}")
    print(f"{'='*50}")
    
    # Read DFT results
    dft_file = f"dft_results/{dot_name}_results.jsonl"
    if not os.path.exists(dft_file):
        print(f"ERROR: DFT results file not found: {dft_file}")
        return False
    
    # Check if file has content
    if os.path.getsize(dft_file) == 0:
        print(f"ERROR: DFT results file is empty (calculation likely failed): {dft_file}")
        return False
    
    with open(dft_file, 'r') as f:
        line = f.readline().strip()
        if not line:
            print(f"ERROR: DFT results file is empty: {dft_file}")
            return False
        dft_data = json.loads(line)
    
    # Extract parameters
    eps_eff = dft_data['epsilon_eff']
    radius_nm = dft_data['used_radius_nm']
    n_atoms = dft_data['n_atoms']
    
    # Generate epsilon profile
    shells = 4
    profile = []
    for i in range(shells):
        factor = 1.0 - (i / (shells - 1)) * 0.4
        profile.append(f"{eps_eff * factor:.2f}")
    eps_profile = ",".join(profile)
    
    # Output directory
    output_dir = Path("cdse_fea_output") / dot_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Parameters:")
    print(f"  Atoms: {n_atoms}")
    print(f"  Radius: {radius_nm:.2f} nm")
    print(f"  Epsilon: {eps_eff:.3f}")
    print(f"  Profile: {eps_profile}")
    
    # FEA command with rhombic pattern (2*diameter spacing)
    fea_cmd = [
        sys.executable, "scripts/run.py",
        "--epsr-profile", eps_profile,
        "--pattern", "hexagonal",
        "--qd-rows", "2",  # 2x2 array for faster execution
        "--qd-cols", "2", 
        "--qd-radius-nm", str(radius_nm),
        "--spacing-factor", "2.0",  # 2*diameter spacing as requested
        "--hbn-bot-nm", "2.0",
        "--hbn-mid-nm", "2.5", 
        "--hbn-top-nm", "2.0",
        "--Lx-nm", "18",  # Compact simulation box
        "--Ly-nm", "18",
        "--Lz-nm", "12",
        "--nx", "36",  # Moderate resolution for speed
        "--ny", "36",
        "--nz", "24",
        "--E-Vnm", "0.15",
        "--lateral-bc", "natural",
        "--workdir", str(output_dir)
    ]
    
    print(f"Running FEA simulation...")
    start_time = time.time()
    
    try:
        result = subprocess.run(fea_cmd, capture_output=True, text=True, check=True, timeout=180)
        
        elapsed = time.time() - start_time
        print(f"SUCCESS: Completed in {elapsed:.1f} seconds")
        
        # Check outputs
        expected_files = ["graphene_potential_2D.png", "graphene_potential_3D.png", 
                         "graphene_slice.csv", "summary.txt"]
        
        all_present = True
        for filename in expected_files:
            filepath = output_dir / filename
            if filepath.exists():
                size_kb = filepath.stat().st_size / 1024
                print(f"  Generated: {filename} ({size_kb:.1f} KB)")
            else:
                print(f"  MISSING: {filename}")
                all_present = False
        
        return all_present
        
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT: Simulation exceeded 3 minutes")
        return False
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Simulation failed (exit code {e.returncode})")
        if e.stderr:
            print(f"Error details: {e.stderr}")
        return False

def main():
    """Run FEA simulations for all CdSe quantum dots individually"""
    
    # Find all DFT result files
    dft_dir = Path("dft_results")
    if not dft_dir.exists():
        print("ERROR: No DFT results directory found!")
        print("Please run simulate_cdse_dft.py first.")
        return
    
    # Get list of quantum dots
    quantum_dots = []
    for jsonl_file in dft_dir.glob("CdSe_*_results.jsonl"):
        dot_name = jsonl_file.stem.replace('_results', '')
        quantum_dots.append(dot_name)
    
    if not quantum_dots:
        print("ERROR: No CdSe quantum dot results found!")
        return
    
    print(f"Found {len(quantum_dots)} CdSe quantum dots to process:")
    for dot in sorted(quantum_dots):
        print(f"  - {dot}")
    
    # Process each quantum dot individually
    successful = []
    failed = []
    
    for dot_name in sorted(quantum_dots):
        success = run_individual_fea(dot_name)
        
        if success:
            successful.append(dot_name)
        else:
            failed.append(dot_name)
        
        # Small delay between simulations
        time.sleep(1)
    
    # Final summary
    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Total quantum dots: {len(quantum_dots)}")
    print(f"Successful simulations: {len(successful)}")
    print(f"Failed simulations: {len(failed)}")
    
    if successful:
        print(f"\nSuccessful:")
        for dot in successful:
            print(f"  SUCCESS: {dot}")
    
    if failed:
        print(f"\nFailed:")
        for dot in failed:
            print(f"  FAILED: {dot}")
    
    print(f"\nConfiguration Summary:")
    print(f"  - Pattern: Hexagonal (rhombic-like arrangement)")
    print(f"  - Spacing: 2 × diameter between quantum dot centers")
    print(f"  - Array: 2×2 quantum dots per simulation")
    print(f"  - Each dot uses its own DFT-calculated dielectric profile")
    print(f"\nResults saved in: cdse_fea_output/")

if __name__ == "__main__":
    main()
