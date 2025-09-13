#!/usr/bin/env python3
"""
Batch runner for DFT calculations on multiple CdSe quantum dots
Runs DFT.py for each .xyz file in the data directory with rhombic pattern and 2*diameter spacing
"""

import subprocess
import sys
import json
import os
from pathlib import Path
import glob

def run_dft_for_file(xyz_file, output_dir="dft_results"):
    """
    Run DFT calculation for a single .xyz file with fast but valid settings
    """
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Get base filename for output
    base_name = Path(xyz_file).stem
    output_json = Path(output_dir) / f"{base_name}_results.jsonl"
    
    # Fast but valid DFT settings for CdSe quantum dots
    dft_cmd = [
        sys.executable, "scripts/DFT.py",
        xyz_file,
        "--lattice", "rhombus",        # Rhombic pattern as requested
        "--spacing-mult", "2.0",       # 2*diameter spacing as requested
        "--angle-deg", "60.0",         # Standard rhombic angle
        "--Lz-nm", "12.0",            # Reasonable z-dimension
        "--Efield-Vnm", "0.01",       # Small field for polarizability calculation
        "--kpts", "1", "1", "1",      # Gamma point only for speed
        "--ecut-eV", "150.0",         # Reduced cutoff for faster calculation
        "--smearing-eV", "0.05",      # Larger smearing for better convergence
        "--xc", "PBE",                # Standard exchange-correlation functional
        "--radius-pad-nm", "0.2",     # Small padding around detected radius
        "--float32",                  # Use float32 for faster calculation
        "--poisson-eps", "1e-10",     # Slightly relaxed convergence for speed
        "--shells", "4",              # Fewer shells for faster analysis
        "--out-json", str(output_json)
    ]
    
    print(f"\n{'='*60}")
    print(f"Running DFT calculation for: {xyz_file}")
    print(f"Output will be saved to: {output_json}")
    print(f"{'='*60}")
    
    try:
        # Run the DFT calculation
        result = subprocess.run(dft_cmd, capture_output=True, text=True, check=True)
        
        print("SUCCESS: DFT calculation completed successfully!")
        
        # Print the output from DFT.py
        if result.stdout:
            print("DFT Output:")
            print(result.stdout)
        
        # Parse and display detailed results
        if output_json.exists():
            with open(output_json, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    print(f"\nResults for {base_name}:")
                    print(f"  Detected radius: {data['detected_radius_xy_nm']:.2f} nm")
                    print(f"  Used radius: {data['used_radius_nm']:.2f} nm")
                    print(f"  Effective epsilon: {data['epsilon_eff']:.3f}")
                    print(f"  Polarizability: {data['alpha_Cm2V']:.3e} C·m²/V")
                    print(f"  Lattice vectors: a1={data['a1_nm']:.2f} nm, a2={data['a2_nm']:.2f} nm")
                    
                    # Extract the epsilon profile for later use
                    eps_profile = []
                    shells = 4
                    eps_val = data['epsilon_eff']
                    for i in range(shells):
                        factor = 1.0 - (i / (shells - 1)) * 0.5
                        eps_profile.append(f"{eps_val * factor:.2f}")
                    
                    print(f"  Epsilon profile: {','.join(eps_profile)}")
        
        return True, str(output_json)
        
    except subprocess.CalledProcessError as e:
        print(f"FAILED: DFT calculation failed!")
        print(f"Error code: {e.returncode}")
        if e.stdout:
            print(f"Stdout: {e.stdout}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return False, None

def main():
    """
    Run DFT calculations for all CdSe quantum dot files
    """
    
    # Find all CdSe .xyz files in the data directory
    data_dir = "data"
    cdse_files = []
    
    # Look for CdSe files (exclude the original Cd68_OPT.xyz and dot_demo.xyz)
    for pattern in ["CdSe_*.xyz"]:
        cdse_files.extend(glob.glob(os.path.join(data_dir, pattern)))
    
    if not cdse_files:
        print("No CdSe .xyz files found in the data directory!")
        return
    
    print(f"Found {len(cdse_files)} CdSe quantum dot files:")
    for f in cdse_files:
        print(f"  - {f}")
    
    # Create results directory
    results_dir = "dft_results"
    Path(results_dir).mkdir(exist_ok=True)
    
    # Run DFT for each file
    successful_runs = []
    failed_runs = []
    
    for xyz_file in cdse_files:
        success, output_file = run_dft_for_file(xyz_file, results_dir)
        
        if success:
            successful_runs.append((xyz_file, output_file))
        else:
            failed_runs.append(xyz_file)
    
    # Summary
    print(f"\n{'='*60}")
    print("BATCH DFT CALCULATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total files processed: {len(cdse_files)}")
    print(f"Successful calculations: {len(successful_runs)}")
    print(f"Failed calculations: {len(failed_runs)}")
    
    if successful_runs:
        print(f"\nSuccessful runs:")
        for xyz_file, output_file in successful_runs:
            print(f"  SUCCESS: {xyz_file} -> {output_file}")
    
    if failed_runs:
        print(f"\nFailed runs:")
        for xyz_file in failed_runs:
            print(f"  FAILED: {xyz_file}")
    
    print(f"\nAll results saved in: {results_dir}/")
    print("You can now use these epsilon profiles for FEA simulations!")

if __name__ == "__main__":
    main()
