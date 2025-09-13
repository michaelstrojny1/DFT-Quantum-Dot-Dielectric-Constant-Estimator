#!/usr/bin/env python3
"""
Run FEA simulations for CdSe quantum dots with rhombic pattern and 2*diameter spacing
Uses the DFT-calculated epsilon profiles from our simulation
"""

import subprocess
import sys
import json
from pathlib import Path
import glob

def load_dft_results():
    """Load all DFT results from the dft_results directory"""
    results = {}
    dft_dir = Path("dft_results")
    
    if not dft_dir.exists():
        print("No DFT results found! Run simulate_cdse_dft.py first.")
        return {}
    
    for jsonl_file in dft_dir.glob("*_results.jsonl"):
        with open(jsonl_file, 'r') as f:
            for line in f:
                data = json.loads(line)
                # Extract base name (e.g., CdSe_medium from CdSe_medium_results.jsonl)
                base_name = jsonl_file.stem.replace('_results', '')
                results[base_name] = data
    
    return results

def generate_epsilon_profile(eps_eff, shells=4):
    """Generate radial epsilon profile from effective epsilon"""
    profile = []
    for i in range(shells):
        # Gradual decrease from center to edge
        factor = 1.0 - (i / (shells - 1)) * 0.4
        profile.append(f"{eps_eff * factor:.2f}")
    return ",".join(profile)

def run_fea_simulation(dot_name, dft_data, output_dir="cdse_fea_output"):
    """Run FEA simulation for a specific quantum dot configuration"""
    
    # Create output directory
    fea_output = Path(output_dir) / dot_name
    fea_output.mkdir(parents=True, exist_ok=True)
    
    # Extract parameters from DFT data
    eps_eff = dft_data['epsilon_eff']
    radius_nm = dft_data['used_radius_nm']
    
    # Generate epsilon profile
    eps_profile = generate_epsilon_profile(eps_eff)
    
    # Calculate spacing with 2*diameter factor as requested
    spacing_factor = 2.0  # 2*diameter spacing
    
    # FEA simulation parameters optimized for rhombic pattern
    fea_cmd = [
        sys.executable, "scripts/run.py",
        "--epsr-profile", eps_profile,
        "--pattern", "hexagonal",  # Hexagonal gives better rhombic-like arrangement
        "--qd-rows", "3",
        "--qd-cols", "3", 
        "--qd-radius-nm", str(radius_nm),
        "--spacing-factor", str(spacing_factor),
        "--hbn-bot-nm", "2.0",
        "--hbn-mid-nm", "2.5", 
        "--hbn-top-nm", "2.0",
        "--Lx-nm", str(max(25, radius_nm * spacing_factor * 4)),  # Dynamic box size
        "--Ly-nm", str(max(25, radius_nm * spacing_factor * 4)),
        "--Lz-nm", "15",
        "--nx", "48",  # Reduced for faster execution
        "--ny", "48",
        "--nz", "36",
        "--E-Vnm", "0.15",  # Reasonable field strength
        "--lateral-bc", "natural",
        "--workdir", str(fea_output)
    ]
    
    print(f"\n{'='*60}")
    print(f"Running FEA simulation for: {dot_name}")
    print(f"  Radius: {radius_nm:.2f} nm")
    print(f"  Effective epsilon: {eps_eff:.3f}")
    print(f"  Epsilon profile: {eps_profile}")
    print(f"  Spacing factor: {spacing_factor} (2*diameter)")
    print(f"  Output directory: {fea_output}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(fea_cmd, capture_output=True, text=True, check=True)
        
        print("SUCCESS: FEA simulation completed!")
        
        # Check for output files
        expected_files = [
            "graphene_potential_2D.png",
            "graphene_potential_3D.png", 
            "graphene_slice.csv",
            "summary.txt"
        ]
        
        found_files = []
        for filename in expected_files:
            filepath = fea_output / filename
            if filepath.exists():
                found_files.append(filename)
                size_kb = filepath.stat().st_size / 1024
                print(f"  Generated: {filename} ({size_kb:.1f} KB)")
        
        if len(found_files) == len(expected_files):
            print(f"  All expected output files generated successfully!")
        else:
            print(f"  Warning: Only {len(found_files)}/{len(expected_files)} files generated")
        
        return True, str(fea_output)
        
    except subprocess.CalledProcessError as e:
        print(f"FAILED: FEA simulation failed!")
        print(f"Error code: {e.returncode}")
        if e.stdout:
            print(f"Stdout: {e.stdout}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return False, None

def main():
    """Run FEA simulations for all CdSe quantum dots"""
    
    # Load DFT results
    dft_results = load_dft_results()
    
    if not dft_results:
        print("No DFT results available. Please run simulate_cdse_dft.py first.")
        return
    
    print(f"Found DFT results for {len(dft_results)} quantum dots:")
    for name, data in dft_results.items():
        print(f"  {name}: R={data['used_radius_nm']:.2f} nm, eps={data['epsilon_eff']:.3f}")
    
    # Run FEA simulations
    successful_runs = []
    failed_runs = []
    
    for dot_name, dft_data in dft_results.items():
        success, output_dir = run_fea_simulation(dot_name, dft_data)
        
        if success:
            successful_runs.append((dot_name, output_dir))
        else:
            failed_runs.append(dot_name)
    
    # Summary
    print(f"\n{'='*60}")
    print("FEA SIMULATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total quantum dots processed: {len(dft_results)}")
    print(f"Successful FEA simulations: {len(successful_runs)}")
    print(f"Failed FEA simulations: {len(failed_runs)}")
    
    if successful_runs:
        print(f"\nSuccessful simulations:")
        for dot_name, output_dir in successful_runs:
            print(f"  SUCCESS: {dot_name} -> {output_dir}")
            print(f"    - 2D potential map: {output_dir}/graphene_potential_2D.png")
            print(f"    - 3D potential plot: {output_dir}/graphene_potential_3D.png")
            print(f"    - Raw data: {output_dir}/graphene_slice.csv")
    
    if failed_runs:
        print(f"\nFailed simulations:")
        for dot_name in failed_runs:
            print(f"  FAILED: {dot_name}")
    
    print(f"\nAll FEA results saved in: cdse_fea_output/")
    print("Each quantum dot has its own subdirectory with visualization files!")
    
    # Show rhombic pattern info
    print(f"\nRhombic Pattern Configuration:")
    print(f"  - Pattern: Hexagonal arrangement (rhombic-like)")
    print(f"  - Spacing: 2 × diameter between dot centers")
    print(f"  - Array size: 3×3 quantum dots")
    print(f"  - Each dot uses size-dependent dielectric profile")

if __name__ == "__main__":
    main()
