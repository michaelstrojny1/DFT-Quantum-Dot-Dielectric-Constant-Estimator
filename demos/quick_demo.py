#!/usr/bin/env python3
"""
Quick demo using user chosen epsilon profile for quantum dot (no DFT required)
Tests the FEA solver without DFT
"""

import subprocess, sys, os
from pathlib import Path

def main():
    demo_cmd = [
        sys.executable, "scripts/run.py",
        "--epsr-profile", "14,12,10,8,6,4",  # Manual shell profile
        "--pattern", "hexagonal",
        "--qd-rows", "3", "--qd-cols", "4", 
        "--qd-radius-nm", "2.5",
        "--spacing-factor", "3.0",
        "--hbn-bot-nm", "2.0", "--hbn-mid-nm", "2.5", "--hbn-top-nm", "2.0",
        "--Lx-nm", "30", "--Ly-nm", "30", "--Lz-nm", "20",
        "--nx", "64", "--ny", "64", "--nz", "48",  # Reduced grid for speed
        "--E-Vnm", "0.15",
        "--lateral-bc", "natural",
        "--workdir", "quick_demo_output"
    ]
    
    try:
        result = subprocess.run(demo_cmd, check=True, capture_output=True, text=True)
        
        output_dir = Path("quick_demo_output")
        if output_dir.exists():
            for f in output_dir.iterdir():
                if f.is_file():
                    print(f"  - {f.name} ({f.stat().st_size} bytes)")
        
        print("\nVisualization files:")
        print("  - graphene_potential_2D.png (2D potential with QD overlay)")
        print("  - graphene_potential_3D.png (3D surface plot)")
        print("  - graphene_slice.csv (raw potential data)")
        print("  - summary.txt (field statistics)")
        
    except subprocess.CalledProcessError as e:
        print("ERROR: Demo failed!")
        if result.stderr:
            print(result.stderr)
        return 1
    except FileNotFoundError:
        print("ERROR: Could not find run.py script")
        print("Make sure you're in the correct directory")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
