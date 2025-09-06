import argparse, json
import numpy as np
from ase.io import read
from ase import Atoms
from gpaw import GPAW, FermiDirac, PoissonSolver
from gpaw.external import ExternalElectricField

EPS0 = 8.854187817620e-12  # F/m
BOHR = 5.29177210903e-11
E_CHARGE = 1.602176634e-19

def center_in_cell(atoms: Atoms, Lx, Ly, Lz):
    a = atoms.copy()
    a.set_cell([[Lx,0,0],[0,Ly,0],[0,0,Lz]], scale_atoms=False)
    a.pbc = (True, True, True)
    pos = a.get_positions()
    shift = np.array([Lx/2, Ly/2, Lz/2]) - pos.mean(axis=0)
    a.set_positions(pos + shift)
    return a

def detect_radius(atoms: Atoms):
    pos = atoms.get_positions()
    c = pos.mean(axis=0)
    return float(np.sqrt(((pos - c)**2).sum(axis=1)).max())

def run_dft(atoms: Atoms, ecut_eV, kpts, smearing_eV, xc, field_Vnm=None):
    calc = GPAW(
        mode="fd",
        ecut=ecut_eV,
        xc=xc,
        kpts=tuple(kpts),
        occupations=FermiDirac(smearing_eV),
        poissonsolver=PoissonSolver(relax="GS", eps=1e-12),
        symmetry={"point_group": False, "time_reversal": True},
        txt=None,
    )
    if field_Vnm is not None and abs(field_Vnm) > 0:
        calc.set(external=ExternalElectricField(field_Vnm*0.1, direction=2))
    atoms = atoms.copy()
    atoms.set_calculator(calc)
    _ = atoms.get_potential_energy()
    dipole = atoms.get_dipole_moment()
    dipole_SI = np.array(dipole) * E_CHARGE * 1e-10
    return dipole_SI

def epsilon_from_alpha(alpha, R):
    f = alpha / (4*np.pi*EPS0*R**3)
    return (1 + 2*f) / (1 - f)

def main():
    ap = argparse.ArgumentParser(description="Finite-field DFT → εr of spherical dot")
    ap.add_argument("xyz", help=".xyz of spherical quantum dot")
    ap.add_argument("--Lx-nm", type=float, default=16.0)
    ap.add_argument("--Ly-nm", type=float, default=16.0)
    ap.add_argument("--Lz-nm", type=float, default=24.0)
    ap.add_argument("--Efield-Vnm", type=float, default=0.01)
    ap.add_argument("--kpts", nargs=3, type=int, default=[3,3,1])
    ap.add_argument("--ecut-eV", type=float, default=500.0)
    ap.add_argument("--smearing-eV", type=float, default=0.01)
    ap.add_argument("--xc", default="PBE")
    ap.add_argument("--dot-radius-nm", type=float, default=None)
    ap.add_argument("--radius-pad-nm", type=float, default=0.25)
    ap.add_argument("--out-json", default="epsr_result.json")
    args = ap.parse_args()

    if args.Efield_Vnm == 0.0:
        raise SystemExit("Need non-zero field (try 0.005–0.02 V/nm).")

    Lx, Ly, Lz = args.Lx_nm*1e-9, args.Ly_nm*1e-9, args.Lz_nm*1e-9
    atoms = center_in_cell(read(args.xyz), Lx, Ly, Lz)

    R_detect = detect_radius(atoms)
    R = args.dot_radius_nm*1e-9 if args.dot_radius_nm else (R_detect + args.radius_pad_nm*1e-9)
    if 2*R > min(Lx,Ly,Lz): raise SystemExit("Dot too large for cell.")

    mu0 = run_dft(atoms, args.ecut_eV, args.kpts, args.smearing_eV, args.xc, field_Vnm=None)
    muE = run_dft(atoms, args.ecut_eV, args.kpts, args.smearing_eV, args.xc, field_Vnm=args.Efield_Vnm)

    dmu = muE - mu0
    alpha = dmu[2] / (args.Efield_Vnm*1e9)
    eps_eff = epsilon_from_alpha(alpha, R)

    result = {
        "xyz": args.xyz,
        "detected_radius_nm": R_detect*1e9,
        "used_radius_nm": R*1e9,
        "field_Vnm": args.Efield_Vnm,
        "alpha_Cm2V": float(alpha),
        "epsilon_eff": float(eps_eff)
    }
    with open(args.out_json,"w") as f: json.dump(result,f,indent=2)
    print(f"{eps_eff:.6f}")

if __name__=="__main__":
    main()
