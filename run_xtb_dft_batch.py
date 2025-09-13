#!/usr/bin/env python3
"""
Run real quantum mechanical calculations using xTB for all CdSe quantum dots
Uses GFN2-xTB method for fast but accurate polarizability calculations
"""

import subprocess
import sys
import os
import json
from pathlib import Path
import time
import glob

def run_xtb_dft(xyz_file, output_dir="dft_results"):
    """Run xTB calculation for a single quantum dot"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Get base filename for output
    base_name = Path(xyz_file).stem
    output_json = Path(output_dir) / f"{base_name}_results.jsonl"
    
    print(f"\n{'='*60}")
    print(f"Running xTB calculation for: {xyz_file}")
    print(f"Method: GFN2-xTB (real quantum mechanics)")
    print(f"Output: {output_json}")
    print(f"{'='*60}")
    
    # xTB command
    xtb_cmd = [
        sys.executable, "scripts/DFT_xtb.py",
        xyz_file,
        "--lattice", "rhombus",
        "--spacing-mult", "2.0",
        "--angle-deg", "60.0", 
        "--Lz-nm", "12.0",
        "--radius-pad-nm", "0.2",
        "--method", "2",  # GFN2-xTB
        "--out-json", str(output_json)
    ]
    
    start_time = time.time()
    
    try:
        result = subprocess.run(xtb_cmd, capture_output=True, text=True, check=True)
        
        elapsed = time.time() - start_time
        print(f"SUCCESS: xTB calculation completed in {elapsed:.1f} seconds")
        
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
        print(f"FAILED: xTB calculation failed!")
        print(f"Error code: {e.returncode}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return False, None

def main():
    """Run xTB calculations for all CdSe quantum dots"""
    
    # Check if xTB is available
    xtb_exe = r"C:\ORCA\xtb\xtb-6.6.1\bin\xtb.exe"
    
    if not Path(xtb_exe).exists():
        print("ERROR: xTB not found!")
        print(f"Expected at: {xtb_exe}")
        return
    
    print(f"xTB found at: {xtb_exe}")
    
    # Find all CdSe .xyz files
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
    
    # Run xTB for each file
    successful_runs = []
    failed_runs = []
    
    print(f"\nStarting xTB calculations (GFN2-xTB method)...")
    print("Note: Larger dots will take longer. Starting with smallest first.")
    
    for xyz_file in cdse_files:
        success, output_file = run_xtb_dft(xyz_file)
        
        if success:
            successful_runs.append((xyz_file, output_file))
        else:
            failed_runs.append(xyz_file)
    
    # Summary
    print(f"\n{'='*60}")
    print("xTB CALCULATION SUMMARY")
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
        print("Ready to run FEA simulations with real quantum mechanical data!")
        print("\nNext step: python run_all_cdse_individual.py")

if __name__ == "__main__":
    main()
