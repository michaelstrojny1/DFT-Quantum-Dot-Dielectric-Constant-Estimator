#!/usr/bin/env python3
"""
Run a single CdSe FEA simulation with monitoring and error handling
"""

import subprocess
import sys
import json
from pathlib import Path

def run_single_fea(dot_name="CdSe_large"):
    """Run FEA for a single quantum dot with full monitoring"""
    
    # Load DFT results
    dft_file = Path("dft_results") / f"{dot_name}_results.jsonl"
    
    if not dft_file.exists():
        print(f"DFT results not found: {dft_file}")
        return False
    
    with open(dft_file, 'r') as f:
        dft_data = json.loads(f.readline())
    
    # Extract parameters
    eps_eff = dft_data['epsilon_eff']
    radius_nm = dft_data['used_radius_nm']
    
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
    
    # FEA command with rhombic pattern (2*diameter spacing)
    fea_cmd = [
        sys.executable, "scripts/run.py",
        "--epsr-profile", eps_profile,
        "--pattern", "hexagonal",
        "--qd-rows", "2",  # Smaller array for testing
        "--qd-cols", "2", 
        "--qd-radius-nm", str(radius_nm),
        "--spacing-factor", "2.0",  # 2*diameter spacing as requested
        "--hbn-bot-nm", "2.0",
        "--hbn-mid-nm", "2.5", 
        "--hbn-top-nm", "2.0",
        "--Lx-nm", "20",  # Fixed size for consistency
        "--Ly-nm", "20",
        "--Lz-nm", "15",
        "--nx", "40",  # Moderate resolution
        "--ny", "40",
        "--nz", "30",
        "--E-Vnm", "0.15",
        "--lateral-bc", "natural",
        "--workdir", str(output_dir)
    ]
    
    print(f"Running FEA simulation for {dot_name}")
    print(f"  Radius: {radius_nm:.2f} nm")
    print(f"  Epsilon: {eps_eff:.3f}")
    print(f"  Profile: {eps_profile}")
    print(f"  Output: {output_dir}")
    print(f"  Command: {' '.join(fea_cmd)}")
    
    try:
        result = subprocess.run(fea_cmd, capture_output=True, text=True, check=True, timeout=300)
        
        print("SUCCESS: FEA simulation completed!")
        
        # Check outputs
        expected_files = ["graphene_potential_2D.png", "graphene_potential_3D.png", 
                         "graphene_slice.csv", "summary.txt"]
        
        for filename in expected_files:
            filepath = output_dir / filename
            if filepath.exists():
                size_kb = filepath.stat().st_size / 1024
                print(f"  SUCCESS: {filename} ({size_kb:.1f} KB)")
            else:
                print(f"  MISSING: {filename}")
        
        return True
        
    except subprocess.TimeoutExpired:
        print("TIMEOUT: FEA simulation took too long (>5 minutes)")
        return False
    except subprocess.CalledProcessError as e:
        print(f"ERROR: FEA simulation failed (exit code {e.returncode})")
        if e.stdout:
            print(f"Stdout: {e.stdout}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return False

if __name__ == "__main__":
    run_single_fea()
