import argparse, subprocess, sys, csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def run_dft_get_profile(dft_py, xyz_path, shells, efield_vnm):
    cmd = [sys.executable, dft_py, xyz_path, "--shells", str(shells),
           "--Efield-Vnm", str(efield_vnm)]
    print("[runner] DFT:", " ".join(cmd))
    out = subprocess.check_output(cmd, text=True).strip()
    profile = out.splitlines()[-1].strip()
    if not any(ch.isdigit() for ch in profile):
        raise RuntimeError("DFT output did not end with a numeric comma list.")
    print(f"[runner] eps_r shells from DFT: {profile}")
    return profile

def compute_spacing(spacing_nm, spacing_factor, qd_radius_nm):
    if spacing_nm is not None:
        return float(spacing_nm)
    return float(spacing_factor) * (2.0 * float(qd_radius_nm))

def centers(pattern, Lx, Ly, zc, spacing, rows, cols, hex_angle_deg=60.0):
    P = []
    if pattern == "single":
        P = [(0.5 * Lx, 0.5 * Ly, zc)]
    elif pattern == "grid":
        sx = spacing; sy = spacing
        x0 = 0.5 * (Lx - sx * (cols - 1))
        y0 = 0.5 * (Ly - sy * (rows - 1))
        for r in range(rows):
            for c in range(cols):
                P.append((x0 + c * sx, y0 + r * sy, zc))
    elif pattern == "hexagonal":
        a = spacing
        ay = a * np.sin(np.deg2rad(hex_angle_deg))
        cols = max(cols, 2); rows = max(rows, 2)
        w = a * (cols - 1) + 0.5 * a
        h = ay * (rows - 1)
        x0 = 0.5 * (Lx - w)
        y0 = 0.5 * (Ly - h)
        for r in range(rows):
            for c in range(cols):
                xs = x0 + c * a + (0.5 * a if (r % 2) else 0.0)
                ys = y0 + r * ay
                P.append((xs, ys, zc))
    else:
        P = [(0.5 * Lx, 0.5 * Ly, zc)]
    return np.asarray(P)

def load_graphene_csv(csv_path):
    xs, ys, ph = [], [], []
    with open(csv_path, "r") as f:
        r = csv.DictReader(f)
        for row in r:
            xs.append(float(row["x_nm"]))
            ys.append(float(row["y_nm"]))
            ph.append(float(row["phi_V"]))
    xs = np.asarray(xs); ys = np.asarray(ys); ph = np.asarray(ph)
    nx = len(np.unique(xs))
    ny = len(np.unique(ys))
    return xs.reshape(ny, nx), ys.reshape(ny, nx), ph.reshape(ny, nx)

def plot_outputs(workdir, Lx_nm, Ly_nm, qd_centers_nm, qd_radius_nm):
    csv_path = Path(workdir) / "graphene_slice.csv"
    if not csv_path.exists():
        print("[runner] warning: graphene_slice.csv not found, skipping plots")
        return
    X, Y, Phi = load_graphene_csv(csv_path)

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    im = ax.imshow(Phi, origin="lower",
                   extent=[0, Lx_nm, 0, Ly_nm], cmap="RdBu_r")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
    ax.set_title("Graphene potential (2D) with QD overlay")
    cbar = fig.colorbar(im, ax=ax); cbar.set_label("V")
    for (cx, cy, _) in qd_centers_nm:
        circ = plt.Circle((cx, cy), qd_radius_nm, fill=False, color="k", lw=0.8)
        ax.add_patch(circ)
        ax.plot([cx], [cy], "k.", ms=2)
    fig.tight_layout()
    fig.savefig(Path(workdir) / "graphene_potential_2D.png", dpi=160)

    fig3 = plt.figure(figsize=(7.8, 6.0))
    ax3 = fig3.add_subplot(111, projection="3d")
    ax3.plot_surface(X, Y, Phi, linewidth=0, antialiased=True, cmap="viridis")
    ax3.set_xlabel("x (nm)"); ax3.set_ylabel("y (nm)"); ax3.set_zlabel("V")
    ax3.set_title("Graphene potential (3D)")
    fig3.tight_layout()
    fig3.savefig(Path(workdir) / "graphene_potential_3D.png", dpi=160)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epsr-profile", type=str)
    ap.add_argument("--xyz", type=str)
    ap.add_argument("--dft", default="scripts/DFT.py")
    ap.add_argument("--dft-shells", type=int, default=6)
    ap.add_argument("--dft-E-Vnm", type=float, default=0.01)
    ap.add_argument("--Lx-nm", type=float, default=40.0)
    ap.add_argument("--Ly-nm", type=float, default=40.0)
    ap.add_argument("--Lz-nm", type=float, default=30.0)
    ap.add_argument("--hbn-bot-nm", type=float, default=3.0)
    ap.add_argument("--hbn-mid-nm", type=float, default=3.0)
    ap.add_argument("--hbn-top-nm", type=float, default=3.0)
    ap.add_argument("--graphene-z-nm", type=float, default=6.0)
    ap.add_argument("--pattern", choices=["single", "grid", "hexagonal"], default="hexagonal")
    ap.add_argument("--qd-rows", type=int, default=4)
    ap.add_argument("--qd-cols", type=int, default=5)
    ap.add_argument("--qd-radius-nm", type=float, default=3.0)
    ap.add_argument("--spacing-nm", type=float, default=None)
    ap.add_argument("--spacing-factor", type=float, default=4.0)
    ap.add_argument("--hex-angle", type=float, default=60.0)
    ap.add_argument("--E-Vnm", type=float, default=0.2)
    ap.add_argument("--lateral-bc", choices=["natural", "cutoff"], default="natural")
    ap.add_argument("--eps-hbn", type=float, default=6.93)
    ap.add_argument("--eps-vac", type=float, default=1.0)
    ap.add_argument("--nx", type=int, default=96)
    ap.add_argument("--ny", type=int, default=96)
    ap.add_argument("--nz", type=int, default=72)
    ap.add_argument("--fea", default="scripts/FEA.py")
    ap.add_argument("--workdir", default="fea_out")

    args = ap.parse_args()

    if args.epsr_profile:
        profile = args.epsr_profile.strip()
        print(f"[runner] Using provided eps_r shells: {profile}")
    elif args.xyz:
        profile = run_dft_get_profile(args.dft, args.xyz, args.dft_shells, args.dft_E_Vnm)
    else:
        sys.exit("Error: provide either --epsr-profile or --xyz (to run DFT.py).")

    spacing_nm = compute_spacing(args.spacing_nm, args.spacing_factor, args.qd_radius_nm)

    topV = args.E_Vnm * args.Lz_nm
    fea_cmd = [
        sys.executable, args.fea,
        "--nx", str(args.nx), "--ny", str(args.ny), "--nz", str(args.nz),
        "--Lx-nm", str(args.Lx_nm), "--Ly-nm", str(args.Ly_nm), "--Lz-nm", str(args.Lz_nm),
        "--hbn-bot-nm", str(args.hbn_bot_nm),
        "--hbn-mid-nm", str(args.hbn_mid_nm),
        "--hbn-top-nm", str(args.hbn_top_nm),
        "--eps-hbn", str(args.eps_hbn),
        "--eps-vac", str(args.eps_vac),
        "--graphene-z-nm", str(args.graphene_z_nm),
        "--pattern", args.pattern,
        "--qd-rows", str(args.qd_rows),
        "--qd-cols", str(args.qd_cols),
        "--spacing-nm", str(spacing_nm),
        "--qd-radius-nm", str(args.qd_radius_nm),
        "--epsr-profile", profile,
        "--bottom_V", "0.0",
        "--top_V", str(topV),
        "--lateral-bc", args.lateral_bc,
        "--workdir", args.workdir,
    ]
    print("[runner] FEA:", " ".join(fea_cmd))
    subprocess.run(fea_cmd, check=True)

    Lx_m = args.Lx_nm * 1e-9; Ly_m = args.Ly_nm * 1e-9
    zc_m = args.hbn_bot_nm * 1e-9 + args.hbn_mid_nm * 1e-9 + args.qd_radius_nm * 1e-9
    C = centers(args.pattern, Lx_m, Ly_m, zc_m,
                spacing_nm * 1e-9, args.qd_rows, args.qd_cols, args.hex_angle)
    C_nm = C.copy(); C_nm[:, 0:2] *= 1e9
    plot_outputs(args.workdir, args.Lx_nm, args.Ly_nm, C_nm, args.qd_radius_nm)

    print("[runner] Done. See:", args.workdir,
          "(graphene_slice.csv, summary.txt, graphene_potential_2D.png, graphene_potential_3D.png)")

if __name__ == "__main__":
    main()
