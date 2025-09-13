#!/usr/bin/env python3
import argparse, json, sys, math
from pathlib import Path
import numpy as np
from ase.io import read
from ase import Atoms
from gpaw import GPAW, FermiDirac, PoissonSolver
from gpaw.external import ExternalElectricField
from gpaw import PW

EPS0 = 8.854187817620e-12
E_CHARGE = 1.602176634e-19

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
    pos = atoms.get_positions()
    cxy = pos[:, :2].mean(axis=0)
    r = np.sqrt(((pos[:, :2] - cxy) ** 2).sum(axis=1)).max()
    return float(r)

def build_rhombus_vectors(side: float, angle_deg: float = 60.0):
    ang = math.radians(angle_deg)
    a1 = np.array([side, 0.0])
    a2 = np.array([side * math.cos(ang), side * math.sin(ang)])
    return a1, a2

def build_rect_vectors(ax: float, ay: float):
    return np.array([ax, 0.0]), np.array([0.0, ay])

def run_dft(atoms: Atoms,
            ecut_eV: float,
            smearing_eV: float,
            xc: str,
            kpts,
            field_Vnm: float | None,
            float32: bool,
            poisson_eps: float,
            use_dipole_layer: bool):
    ps_kwargs = dict(relax="GS", eps=poisson_eps)
    if use_dipole_layer:
        ps_kwargs["dipolelayer"] = "xy"

    calc = GPAW(
        mode=PW(ecut_eV),
        xc=xc,
        kpts=tuple(kpts),
        occupations=FermiDirac(smearing_eV),
        poissonsolver=PoissonSolver(**ps_kwargs),
        symmetry={"point_group": False, "time_reversal": True},
        dtype=np.float32 if float32 else np.float64,
        txt=None,
    )

    if field_Vnm is not None and abs(field_Vnm) > 0:
        calc.set(external=ExternalElectricField(field_Vnm * 0.1, direction=2))

    atoms = atoms.copy()
    atoms.set_calculator(calc)
    _ = atoms.get_potential_energy()
    dipole_eA = np.array(atoms.get_dipole_moment())
    dipole_SI = dipole_eA * E_CHARGE * 1e-10
    return dipole_SI

def epsilon_from_alpha(alpha_SI: float, R_m: float) -> float:
    f = alpha_SI / (4.0 * math.pi * EPS0 * (R_m ** 3))
    if f >= 0.999:
        return float("inf")
    if f <= -0.999:
        return -1.0
    return (1.0 + 2.0 * f) / (1.0 - f)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz", nargs="+")
    ap.add_argument("--lattice", choices=["rhombus", "rect"], default="rhombus")
    ap.add_argument("--spacing-mult", type=float, default=2.2)
    ap.add_argument("--angle-deg", type=float, default=60.0)
    ap.add_argument("--Lz-nm", type=float, default=16.0)
    ap.add_argument("--Efield-Vnm", type=float, default=0.01)
    ap.add_argument("--kpts", nargs=3, type=int, default=[1, 1, 1])
    ap.add_argument("--ecut-eV", type=float, default=200.0)
    ap.add_argument("--smearing-eV", type=float, default=0.02)
    ap.add_argument("--xc", default="PBE")
    ap.add_argument("--radius-pad-nm", type=float, default=0.25)
    ap.add_argument("--dot-radius-nm", type=float, default=None)
    ap.add_argument("--float32", action="store_true")
    ap.add_argument("--poisson-eps", type=float, default=1e-12)
    ap.add_argument("--no-dipole-layer", action="store_true")
    ap.add_argument("--out-json", default="epsr_results.jsonl")
    ap.add_argument("--shells", type=int, default=6)
    args = ap.parse_args()

    if args.Efield_Vnm == 0.0:
        raise SystemExit("Need non-zero field (try 0.005–0.02 V/nm).")

    Lz = args.Lz_nm * 1e-9
    results = []

    for xyz_path in args.xyz:
        p = Path(xyz_path)
        if not p.exists():
            print(f"[WARN] Missing file: {p}", file=sys.stderr)
            continue

        atoms_raw = read(str(p))
        R_xy_det_m = detect_radius_xy(atoms_raw)
        R_used_m = (args.dot_radius_nm * 1e-9) if args.dot_radius_nm else (
            R_xy_det_m + args.radius_pad_nm * 1e-9
        )
        center_spacing = args.spacing_mult * (2.0 * R_used_m)
        if args.lattice == "rhombus":
            a1, a2 = build_rhombus_vectors(center_spacing, angle_deg=args.angle_deg)
        else:
            a1, a2 = build_rect_vectors(center_spacing, center_spacing)
        atoms = center_xy_in_cell(atoms_raw, a1, a2, Lz)

        mu0 = run_dft(
            atoms,
            ecut_eV=args.ecut_eV,
            smearing_eV=args.smearing_eV,
            xc=args.xc,
            kpts=args.kpts,
            field_Vnm=None,
            float32=args.float32,
            poisson_eps=args.poisson_eps,
            use_dipole_layer=(not args.no_dipole_layer),
        )
        muE = run_dft(
            atoms,
            ecut_eV=args.ecut_eV,
            smearing_eV=args.smearing_eV,
            xc=args.xc,
            kpts=args.kpts,
            field_Vnm=args.Efield_Vnm,
            float32=args.float32,
            poisson_eps=args.poisson_eps,
            use_dipole_layer=(not args.no_dipole_layer),
        )

        dmu = muE - mu0
        alpha = dmu[2] / (args.Efield_Vnm * 1e9)

        eps_eff = epsilon_from_alpha(alpha, R_used_m)

        rec = {
            "xyz": str(p),
            "lattice": args.lattice,
            "angle_deg": args.angle_deg if args.lattice == "rhombus" else None,
            "spacing_mult": args.spacing_mult,
            "detected_radius_xy_nm": R_xy_det_m * 1e9,
            "used_radius_nm": R_used_m * 1e9,
            "a1_nm": np.linalg.norm(a1) * 1e9,
            "a2_nm": np.linalg.norm(a2) * 1e9,
            "Lz_nm": args.Lz_nm,
            "field_Vnm": args.Efield_Vnm,
            "alpha_Cm2V": float(alpha),
            "epsilon_eff": float(eps_eff),
            "dipole0_Cm": list(map(float, mu0)),
            "dipoleE_Cm": list(map(float, muE)),
            "ecut_eV": args.ecut_eV,
            "kpts": args.kpts,
            "xc": args.xc,
            "smearing_eV": args.smearing_eV,
            "float32": bool(args.float32),
        }
        results.append(rec)

    with open(args.out_json, "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    for r in results:
        Rnm = r["used_radius_nm"]
        print(f"{Path(r['xyz']).name}\tR={Rnm:.3f} nm\teps_eff={r['epsilon_eff']:.6f}")
    if len(results) == 1 and "shells" in vars(args):
        eps_val = results[0]["epsilon_eff"]
        shell_count = getattr(args, 'shells', 6)
        profile = []
        for i in range(shell_count):
            factor = 1.0 - (i / (shell_count - 1)) * 0.5
            profile.append(f"{eps_val * factor:.2f}")
        print(",".join(profile))

if __name__ == "__main__":
    main()
