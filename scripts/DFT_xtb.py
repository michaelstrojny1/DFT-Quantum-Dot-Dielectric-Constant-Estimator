#!/usr/bin/env python3
"""
Real quantum mechanical dielectric constant estimator using xTB
Uses GFN2-xTB method for fast but accurate polarizability calculations
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from ase.io import read

# Physical constants
EPS0 = 8.8541878128e-12  # F/m
E_CHARGE = 1.602176634e-19  # C

def detect_radius_xy(atoms):
    """Detect a robust 3D radius of the quantum dot.

    Historically this used the max XY extent, which can under-estimate the
    true size when the dot is elongated along z or rotated. We switch to the
    95th percentile of 3D distances from the center to be orientation- and
    outlier-robust.
    """
    pos = atoms.get_positions()
    center = np.mean(pos, axis=0)
    distances = np.linalg.norm(pos - center, axis=1)
    r95 = np.percentile(distances, 95)
    return float(r95) * 1e-10  # Convert Å to m

def epsilon_from_alpha(alpha_SI, R_m):
    """Convert polarizability to effective dielectric constant using Clausius-Mossotti"""
    V_sphere = (4.0 / 3.0) * math.pi * R_m**3
    alpha_norm = alpha_SI / (4.0 * math.pi * EPS0 * V_sphere)
    eps_eff = (1.0 + 2.0 * alpha_norm) / (1.0 - alpha_norm)
    return eps_eff

def build_rhombus_vectors(center_spacing_m: float, angle_deg: float = 60.0):
    """Build rhombic lattice vectors"""
    angle_rad = math.radians(angle_deg)
    a1 = np.array([center_spacing_m, 0.0])
    a2 = np.array([center_spacing_m * math.cos(angle_rad), center_spacing_m * math.sin(angle_rad)])
    return a1, a2

def build_rect_vectors(ax_m: float, ay_m: float):
    """Build rectangular lattice vectors"""
    return np.array([ax_m, 0.0]), np.array([0.0, ay_m])

def center_xy_in_cell(atoms, a1, a2, Lz):
    """Center atoms in unit cell"""
    atoms_copy = atoms.copy()
    pos = atoms_copy.get_positions()
    
    # Center in xy
    center_xy = np.mean(pos[:, :2], axis=0)
    cell_center_xy = 0.5 * (a1[:2] + a2[:2])
    shift_xy = cell_center_xy - center_xy
    pos[:, :2] += shift_xy
    
    # Center in z
    center_z = np.mean(pos[:, 2])
    pos[:, 2] += (Lz * 1e10 / 2.0) - center_z  # Convert m to Å
    
    atoms_copy.set_positions(pos)
    
    # Set cell
    cell = np.zeros((3, 3))
    cell[0, :2] = a1[:2] * 1e10  # Convert m to Å
    cell[1, :2] = a2[:2] * 1e10
    cell[2, 2] = Lz * 1e10
    atoms_copy.set_cell(cell)
    atoms_copy.set_pbc([True, True, True])
    
    return atoms_copy

# xTB helpers
XTB_EXE = r"C:\ORCA\xtb\xtb-6.6.1\bin\xtb.exe"

def run_xtb_polarizability(atoms, method="2", charge: int = 0, uhf: int = 0):
    """Run xTB calculation to get polarizability with robust fallbacks"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Write structure file
        xyz_file = tmpdir_path / "structure.xyz"
        atoms.write(str(xyz_file))

        # Define a series of attempts with varying settings for robustness
        attempts = [
            {"etemp": 500,  "iterations": 500,  "timeout": 1200, "use_property": False},
            {"etemp": 500,  "iterations": 800,  "timeout": 1500, "use_property": False},
            {"etemp": 800,  "iterations": 1200, "timeout": 1800, "use_property": True},
            {"etemp": 1000, "iterations": 1500, "timeout": 2400, "use_property": True},
        ]

        last_err = None
        for idx, opts in enumerate(attempts, start=1):
            # Build xTB command
            cmd = [
                XTB_EXE,
                str(xyz_file),
                f"--gfn{method}",
                "--chrg", str(charge),
                "--uhf", str(uhf),
                "--etemp", str(opts["etemp"]),
                "--iterations", str(opts["iterations"]),
                "--verbose",
            ]
            if opts.get("use_property"):
                cmd.insert(-1, "--property")

            try:
                result = subprocess.run(
                    cmd,
                    cwd=str(tmpdir_path),
                    capture_output=True,
                    text=True,
                    encoding='utf-8',
                    errors='ignore',
                    check=True,
                    timeout=opts["timeout"],
                )

                # Parse polarizability from output
                output = result.stdout

                alpha_patterns = [
                    r"Mol\.\s+α\(0\)\s+/au\s*:\s*([0-9]+\.?[0-9]*)",
                    r"isotropic\s+polarizability\s+α\s*=\s*([0-9]+\.?[0-9]*)\s*Bohr\^3",
                    r"polarizability\s+α\s*=\s*([0-9]+\.?[0-9]*)\s*Bohr\^3",
                    r"α\(iso\)\s*=\s*([0-9]+\.?[0-9]*)\s*Bohr\^3",
                ]

                alpha_bohr3 = None
                for pattern in alpha_patterns:
                    match = re.search(pattern, output, re.IGNORECASE)
                    if match:
                        alpha_bohr3 = float(match.group(1))
                        break

                if alpha_bohr3 is None:
                    # Try to find it in any .out file produced by xTB
                    prop_files = list(tmpdir_path.glob("*.out"))
                    for prop_file in prop_files:
                        with open(prop_file, 'r', encoding='utf-8', errors='ignore') as f:
                            prop_content = f.read()
                            for pattern in alpha_patterns:
                                match = re.search(pattern, prop_content, re.IGNORECASE)
                                if match:
                                    alpha_bohr3 = float(match.group(1))
                                    break
                        if alpha_bohr3 is not None:
                            break

                if alpha_bohr3 is None:
                    last_err = RuntimeError("Could not parse polarizability from xTB output")
                    continue

                # Convert Bohr^3 to Å^3
                bohr_to_angstrom = 0.529177210903
                alpha_A3 = alpha_bohr3 * (bohr_to_angstrom**3)

                return alpha_A3, output

            except subprocess.TimeoutExpired:
                last_err = RuntimeError("xTB calculation timed out")
                continue
            except subprocess.CalledProcessError as e:
                # Keep last error detail for reporting
                msg = e.stderr or e.stdout or ""
                last_err = RuntimeError(f"xTB calculation failed: {msg.strip()}")
                continue

        # If we reach here, all attempts failed
        if last_err is None:
            last_err = RuntimeError("xTB calculation failed for unknown reasons")
        raise last_err

def main():
    """Main function for xTB-based DFT calculations"""
    
    ap = argparse.ArgumentParser(description="Real DFT dielectric constant estimation using xTB")
    ap.add_argument("xyz", nargs="+", help="XYZ files to process")
    ap.add_argument("--lattice", choices=["rhombus", "rect"], default="rhombus")
    ap.add_argument("--spacing-mult", type=float, default=2.0)
    ap.add_argument("--angle-deg", type=float, default=60.0)
    ap.add_argument("--Lz-nm", type=float, default=12.0)
    ap.add_argument("--radius-pad-nm", type=float, default=0.2)
    ap.add_argument("--dot-radius-nm", type=float, default=None)
    ap.add_argument("--method", choices=["1", "2"], default="2", help="GFN method (1 or 2)")
    ap.add_argument("--out-json", default="xtb_results.jsonl")
    
    args = ap.parse_args()
    
    # Check if xTB is available
    if not Path(XTB_EXE).exists():
        print(f"ERROR: xTB not found at {XTB_EXE}")
        return 1
    
    Lz = args.Lz_nm * 1e-9
    results = []
    
    for xyz_path in args.xyz:
        p = Path(xyz_path)
        if not p.exists():
            print(f"[WARN] Missing file: {p}", file=sys.stderr)
            continue
        
        print(f"\nProcessing: {p.name}")
        print(f"Method: GFN{args.method}-xTB")
        
        # Read structure
        atoms_raw = read(str(p))
        n_atoms = len(atoms_raw)
        print(f"Atoms: {n_atoms}")
        
        # Detect radius
        R_xy_det_m = detect_radius_xy(atoms_raw)
        R_used_m = (args.dot_radius_nm * 1e-9) if args.dot_radius_nm else (
            R_xy_det_m + args.radius_pad_nm * 1e-9
        )
        
        # Build lattice
        center_spacing = args.spacing_mult * (2.0 * R_used_m)
        if args.lattice == "rhombus":
            a1, a2 = build_rhombus_vectors(center_spacing, angle_deg=args.angle_deg)
        else:
            a1, a2 = build_rect_vectors(center_spacing, center_spacing)
        
        # Center in cell
        atoms = center_xy_in_cell(atoms_raw, a1, a2, Lz)
        
        try:
            # Run xTB calculation
            print("Running xTB calculation...")
            alpha_A3, xtb_output = run_xtb_polarizability(atoms, method=args.method)
            
            # Convert to SI polarizability
            alpha_SI = 4.0 * math.pi * EPS0 * (alpha_A3 * 1e-30)
            
            # Calculate effective dielectric constant
            eps_eff = epsilon_from_alpha(alpha_SI, R_used_m)
            
            print(f"Polarizability: {alpha_A3:.2f} Å³")
            print(f"Effective epsilon: {eps_eff:.3f}")
            
            # Store results
            rec = {
                "xyz": str(p),
                "engine": "xTB",
                "method": f"GFN{args.method}-xTB",
                "lattice": args.lattice,
                "angle_deg": args.angle_deg if args.lattice == "rhombus" else None,
                "spacing_mult": args.spacing_mult,
                "detected_radius_xy_nm": R_xy_det_m * 1e9,
                "used_radius_nm": R_used_m * 1e9,
                "a1_nm": float(np.linalg.norm(a1) * 1e9),
                "a2_nm": float(np.linalg.norm(a2) * 1e9),
                "Lz_nm": args.Lz_nm,
                "alpha_A3": float(alpha_A3),
                "alpha_SI": float(alpha_SI),
                "epsilon_eff": float(eps_eff),
                "n_atoms": n_atoms
            }
            results.append(rec)
            
        except Exception as e:
            print(f"ERROR processing {p.name}: {e}")
            continue
    
    # Write results
    with open(args.out_json, "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    
    # Print summary
    print(f"\n{'='*50}")
    print("xTB DFT RESULTS SUMMARY")
    print(f"{'='*50}")
    for r in results:
        print(f"{Path(r['xyz']).name:15} R={r['used_radius_nm']:.2f}nm  eps={r['epsilon_eff']:.3f}")
    
    print(f"Results saved to: {args.out_json}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
