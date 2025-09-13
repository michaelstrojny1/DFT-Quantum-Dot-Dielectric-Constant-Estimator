#!/usr/bin/env python3
"""
Run real DFT calculations using ORCA for all CdSe quantum dots
Uses PBEh-3c method for fast but accurate DFT calculations
"""

import subprocess
import sys
import os
import json
from pathlib import Path
import time

def run_orca_dft(xyz_file, output_dir="dft_results"):
    """Run ORCA DFT calculation for a single quantum dot"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Get base filename for output
    base_name = Path(xyz_file).stem
    output_json = Path(output_dir) / f"{base_name}_results.jsonl"
    
    print(f"\n{'='*60}")
    print(f"Running ORCA DFT calculation for: {xyz_file}")
    print(f"Method: PBEh-3c (fast composite DFT)")
    print(f"Output: {output_json}")
    print(f"{'='*60}")
    
    # ORCA command
    orca_cmd = [
        sys.executable, "scripts/DFT_orca.py",
        xyz_file,
        "--lattice", "rhombus",
        "--spacing-mult", "2.0",
        "--angle-deg", "60.0", 
        "--Lz-nm", "12.0",
        "--radius-pad-nm", "0.2",
        "--nprocs", str(min(4, os.cpu_count() or 4)),  # Use up to 4 cores
        "--maxcore-mb", "2000",  # 2GB per core
        "--out-json", str(output_json)
    ]
    
    start_time = time.time()
    
    try:
        result = subprocess.run(orca_cmd, capture_output=True, text=True, check=True)
        
        elapsed = time.time() - start_time
        print(f"SUCCESS: DFT calculation completed in {elapsed:.1f} seconds")
        
        # Print output
        if result.stdout:
            print(result.stdout)
        
        # Parse and display results
        if output_json.exists():
            with open(output_json, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    print(f"\nResults:")
                    print(f"  Detected radius: {data['detected_radius_xy_nm']:.2f} nm")
                    print(f"  Used radius: {data['used_radius_nm']:.2f} nm")
                    print(f"  Polarizability: {data['alpha_A3']:.2f} Å³")
                    print(f"  Effective epsilon: {data['epsilon_eff']:.3f}")
        
        return True, str(output_json)
        
    except subprocess.CalledProcessError as e:
        print(f"FAILED: DFT calculation failed!")
        print(f"Error code: {e.returncode}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return False, None

def main():
    """Run ORCA DFT calculations for all CdSe quantum dots"""
    
    # Check if ORCA is available
    orca_exe = os.environ.get("ORCA_EXE", "orca")
    
    # Test if ORCA is accessible
    try:
        subprocess.run([orca_exe, "--version"], capture_output=True, check=False)
        print(f"ORCA found at: {orca_exe}")
    except FileNotFoundError:
        print("ERROR: ORCA not found!")
        print("\nInstallation instructions:")
        print("1. Download ORCA from: https://orcaforum.kofo.mpg.de")
        print("2. Extract to a folder (e.g., C:\\ORCA)")  
        print("3. Set environment variable:")
        print('   $env:ORCA_EXE = "C:\\ORCA\\orca.exe"')
        print("4. Run this script again")
        return
    
    # Find all CdSe .xyz files
    import glob
    data_dir = "data"
    cdse_files = sorted(glob.glob(os.path.join(data_dir, "CdSe_*.xyz")))
    
    if not cdse_files:
        print("No CdSe .xyz files found in the data directory!")
        return
    
    print(f"\nFound {len(cdse_files)} CdSe quantum dot files:")
    for f in cdse_files:
        size = Path(f).stat().st_size
        with open(f, 'r') as xf:
            n_atoms = int(xf.readline().strip())
        print(f"  - {Path(f).name}: {n_atoms} atoms")
    
    # Sort by file size to run smaller ones first
    cdse_files.sort(key=lambda x: Path(x).stat().st_size)
    
    # Run DFT for each file
    successful_runs = []
    failed_runs = []
    
    print(f"\nStarting ORCA DFT calculations (PBEh-3c method)...")
    print("Note: Larger dots will take longer. Starting with smallest first.")
    
    for xyz_file in cdse_files:
        success, output_file = run_orca_dft(xyz_file)
        
        if success:
            successful_runs.append((xyz_file, output_file))
        else:
            failed_runs.append(xyz_file)
    
    # Summary
    print(f"\n{'='*60}")
    print("ORCA DFT CALCULATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total files processed: {len(cdse_files)}")
    print(f"Successful calculations: {len(successful_runs)}")
    print(f"Failed calculations: {len(failed_runs)}")
    
    if successful_runs:
        print(f"\nSuccessful runs:")
        for xyz_file, output_file in successful_runs:
            print(f"  SUCCESS: {Path(xyz_file).name} -> {output_file}")
    
    if failed_runs:
        print(f"\nFailed runs:")
        for xyz_file in failed_runs:
            print(f"  FAILED: {Path(xyz_file).name}")
    
    if successful_runs:
        print(f"\nAll DFT results saved in: dft_results/")
        print("Ready to run FEA simulations with real DFT data!")
        print("\nNext step: python run_all_cdse_individual.py")

if __name__ == "__main__":
    main()
