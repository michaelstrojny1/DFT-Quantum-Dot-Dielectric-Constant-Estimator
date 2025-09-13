#!/usr/bin/env python3
"""
DFT_orca.py

Real DFT dielectric constant estimator using ORCA on Windows/Linux/macOS.
- Computes static polarizability via ORCA's %elprop Polar True (linear response)
- Converts to an effective dielectric constant using Clausius–Mossotti
- Supports multiple XYZ files; rhombic lattice spacing like scripts/DFT.py

Prereqs:
- Install ORCA and ensure the orca executable is on PATH or set ORCA_EXE env var
- ORCA 5.x recommended
- Method/basis: B3LYP D3BJ with def2-SVP (with def2/J and RIJCOSX for speed)
- For heavy atoms (Cd, Se), def2-SVP will automatically use def2-ECP in ORCA

Example:
  python scripts/DFT_orca.py data/CdSe_34.xyz data/CdSe_68.xyz \
      --lattice rhombus --spacing-mult 2.0 --Lz-nm 12 --out-json epsr_results_orca.jsonl

Outputs:
- JSON lines with keys: xyz, epsilon_eff, alpha_A3, used_radius_nm, etc.
- Final printed summary per dot
"""

import argparse, json, os, re, shutil, subprocess, sys, math, tempfile
from pathlib import Path
import numpy as np
from ase.io import read
from ase import Atoms

EPS0 = 8.854187817620e-12

# ---------- Geometry helpers (no GPAW required) ----------

def center_xy_in_cell(atoms: Atoms, a1, a2, Lz):
    a = atoms.copy()
    cell = np.array([[a1[0], a2[0], 0.0],
                     [a1[1], a2[1], 0.0],
                     [0.0,    0.0,   Lz ]], dtype=float)
    a.set_cell(cell, scale_atoms=False)
    a.pbc = (True, True, False)
    pos = a.get_positions()
    xy_center = pos[:, :2].mean(axis=0)
    z_center = pos[:, 2].mean()
    cell_xy_center = 0.5*(cell[0, :2] + cell[1, :2])
    shift_xy = cell_xy_center - xy_center
    shift_z = (Lz/2.0) - z_center
    shift = np.array([shift_xy[0], shift_xy[1], shift_z])
    a.set_positions(pos + shift)
    return a


def detect_radius_xy(atoms: Atoms) -> float:
    pos = atoms.get_positions()  # in Angstrom
    cxy = pos[:, :2].mean(axis=0)
    r = np.sqrt(((pos[:, :2] - cxy) ** 2).sum(axis=1)).max()  # Angstrom
    return float(r) * 1e-10  # meters


def build_rhombus_vectors(side_m: float, angle_deg: float = 60.0):
    ang = math.radians(angle_deg)
    a1 = np.array([side_m, 0.0])
    a2 = np.array([side_m * math.cos(ang), side_m * math.sin(ang)])
    return a1, a2


def build_rect_vectors(ax_m: float, ay_m: float):
    return np.array([ax_m, 0.0]), np.array([0.0, ay_m])


# ---------- ORCA helpers ----------

ORCA_EXE = os.environ.get("ORCA_EXE", "orca")

ORCA_HEADER = """! PBEh-3c SP NormalSCF SlowConv
%pal
  nprocs %d
end
%maxcore %d
%scf
  MaxIter 200
end
%elprop
  Polar True
end
"""

# Notes:
# - PBEh-3c: Fast composite DFT method (still real DFT, ~10-100x faster than B3LYP/def2)
#   Includes: PBEh functional, modified def2-mSVP basis, D3BJ dispersion, gCP correction
#   Excellent for large systems like quantum dots with hundreds of atoms
# - NormalSCF for speed; SlowConv to aid convergence on metallic-like clusters
# - Polar True computes static (zero-frequency) polarizability via response

POLAR_REGEXS = [
    re.compile(r"Isotropic\s+polarizability\s*:\s*([0-9]+\.?[0-9]*)\s*A\^3", re.IGNORECASE),
    re.compile(r"Average\s+polarizability\s*:\s*([0-9]+\.?[0-9]*)\s*A\^3", re.IGNORECASE),
    re.compile(r"Isotropic\s+polarizability\s*:\s*([0-9]+\.?[0-9]*)\s*au", re.IGNORECASE),
]

AU_POL_TO_A3 = (0.529177210903)**3  # 1 a.u. pol = a0^3 in Ang^3


def write_orca_input(inp_path: Path, atoms: Atoms, nprocs: int, maxcore_mb: int):
    lines = []
    lines.append(ORCA_HEADER % (nprocs, maxcore_mb))
    lines.append("* xyz 0 1")
    for sym, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"{sym:2s}  {x: .8f}  {y: .8f}  {z: .8f}")
    lines.append("*")
    inp_path.write_text("\n".join(lines))


def run_orca(inp_path: Path, out_path: Path) -> None:
    try:
        with open(out_path, "w", encoding="utf-8", errors="ignore") as fout:
            subprocess.run([ORCA_EXE, str(inp_path)], stdout=fout, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"ORCA failed for {inp_path}: exit {e.returncode}")


def parse_polarizability_A3(out_text: str) -> float | None:
    # Try regexes
    for rx in POLAR_REGEXS:
        m = rx.search(out_text)
        if m:
            val = float(m.group(1))
            # If the match was in a.u., convert to A^3
            if "au" in rx.pattern:
                return val * AU_POL_TO_A3
            return val
    # Fallback: try to find a polarizability table and average diagonal (A^3)
    # Heuristic for robustness
    try:
        lines = out_text.splitlines()
        idx = None
        for i, L in enumerate(lines):
            if "POLARIZABILITY TENSOR" in L.upper() and ("A**3" in L.upper() or "ANGSTROM^3" in L.upper()):
                idx = i
                break
        if idx is not None:
            block = "\n".join(lines[idx: idx+20])
            nums = re.findall(r"([-+]?[0-9]*\.?[0-9]+)", block)
            vals = [float(x) for x in nums[:9]]  # 3x3
            if len(vals) >= 9:
                # average of diagonal
                a_iso = (vals[0] + vals[4] + vals[8]) / 3.0
                return a_iso
    except Exception:
        pass
    return None


def epsilon_from_alpha(alpha_SI: float, R_m: float) -> float:
    f = alpha_SI / (4.0 * math.pi * EPS0 * (R_m ** 3))
    if f >= 0.999:
        return float("inf")
    if f <= -0.999:
        return -1.0
    return (1.0 + 2.0 * f) / (1.0 - f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz", nargs="+", help="One or more XYZ files")
    ap.add_argument("--lattice", choices=["rhombus", "rect"], default="rhombus")
    ap.add_argument("--spacing-mult", type=float, default=2.0)
    ap.add_argument("--angle-deg", type=float, default=60.0)
    ap.add_argument("--Lz-nm", type=float, default=12.0)
    ap.add_argument("--radius-pad-nm", type=float, default=0.2)
    ap.add_argument("--dot-radius-nm", type=float, default=None)
    ap.add_argument("--nprocs", type=int, default=max(1, os.cpu_count() or 1))
    ap.add_argument("--maxcore-mb", type=int, default=2000)
    ap.add_argument("--out-json", default="epsr_results_orca.jsonl")
    args = ap.parse_args()

    # Quick check for ORCA availability
    if shutil.which(ORCA_EXE) is None:
        print(f"ERROR: Could not find ORCA executable '{ORCA_EXE}'. Set ORCA_EXE env var or add to PATH.")
        sys.exit(1)

    Lz = args.Lz_nm * 1e-9
    results = []

    with open(args.out_json, "w", encoding="utf-8") as jf:
        for xyz_path in args.xyz:
            p = Path(xyz_path)
            if not p.exists():
                print(f"[WARN] Missing file: {p}", file=sys.stderr)
                continue

            atoms_raw = read(str(p))
            R_xy_det_m = detect_radius_xy(atoms_raw)  # meters
            R_used_m = (args.dot_radius_nm * 1e-9) if args.dot_radius_nm else (
                R_xy_det_m + args.radius_pad_nm * 1e-9
            )

            center_spacing = args.spacing_mult * (2.0 * R_used_m)
            if args.lattice == "rhombus":
                a1, a2 = build_rhombus_vectors(center_spacing, angle_deg=args.angle_deg)
            else:
                a1, a2 = build_rect_vectors(center_spacing, center_spacing)
            atoms = center_xy_in_cell(atoms_raw, a1, a2, Lz)

            # Prepare ORCA run in temp dir
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                inp = td / (p.stem + ".inp")
                out = td / (p.stem + ".out")
                write_orca_input(inp, atoms, args.nprocs, args.maxcore_mb)
                run_orca(inp, out)
                out_text = out.read_text(errors="ignore")

            alpha_A3 = parse_polarizability_A3(out_text)
            if alpha_A3 is None:
                print(f"ERROR: Could not parse polarizability from ORCA output for {p.name}", file=sys.stderr)
                continue

            # Convert alpha volume (A^3) to SI polarizability: alpha_SI = 4*pi*eps0 * alpha_vol
            alpha_SI = 4.0 * math.pi * EPS0 * (alpha_A3 * 1e-30)
            eps_eff = epsilon_from_alpha(alpha_SI, R_used_m)

            rec = {
                "xyz": str(p),
                "engine": "ORCA",
                "method": "PBEh-3c (fast composite DFT)",
                "nprocs": args.nprocs,
                "maxcore_mb": args.maxcore_mb,
                "lattice": args.lattice,
                "angle_deg": args.angle_deg if args.lattice == "rhombus" else None,
                "spacing_mult": args.spacing_mult,
                "detected_radius_xy_nm": R_xy_det_m * 1e9,
                "used_radius_nm": R_used_m * 1e9,
                "a1_nm": float(np.linalg.norm(a1) * 1e9),
                "a2_nm": float(np.linalg.norm(a2) * 1e9),
                "Lz_nm": args.Lz_nm,
                "alpha_A3": float(alpha_A3),
                "epsilon_eff": float(eps_eff),
            }
            jf.write(json.dumps(rec) + "\n")
            results.append(rec)
            print(f"{p.name}\tR={rec['used_radius_nm']:.3f} nm\talpha={alpha_A3:.3f} A^3\teps_eff={eps_eff:.6f}")

    if not results:
        sys.exit(2)


if __name__ == "__main__":
    main()
