import argparse, numpy as np
from pathlib import Path
from skfem import MeshHex, ElementHex1, Basis, asm, solve
from skfem.helpers import grad, dot

EPS0 = 8.854187817620e-12

def linspace_nodes(L, n):
    return np.linspace(0.0, L, n + 1)

def generate_centers(pattern, Lx, Ly, zc, spacing, rows, cols, hex_angle_deg=60.0):
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
    return np.array(P)

def build_epsfun(Lx, Ly, Lz, hbn_bot, hbn_mid, hbn_top, eps_hbn,
                 centers, R, epsr_shells, eps_vac):
    shells = np.array(epsr_shells, dtype=float)
    nshell = len(shells)
    def epsfun(x):
        X, Y, Z = x
        epsr = np.full_like(X, float(eps_vac))
        mask_bot = (Z <= hbn_bot)
        mask_mid = (Z >= hbn_bot) & (Z <= (hbn_bot + hbn_mid))
        mask_top = (Z >= (Lz - hbn_top))
        epsr[mask_bot | mask_mid | mask_top] = eps_hbn
        if centers.size > 0:
            for cx, cy, cz in centers:
                rr = np.sqrt((X - cx)**2 + (Y - cy)**2 + (Z - cz)**2)
                inside = rr <= R
                if not np.any(inside): 
                    continue
                si = np.minimum((rr[inside] / R * nshell).astype(int), nshell - 1)
                epsr[inside] = shells[si]
        return EPS0 * epsr
    return epsfun

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=96)
    ap.add_argument("--ny", type=int, default=96)
    ap.add_argument("--nz", type=int, default=72)
    ap.add_argument("--Lx-nm", type=float, default=40.0)
    ap.add_argument("--Ly-nm", type=float, default=40.0)
    ap.add_argument("--Lz-nm", type=float, default=30.0)
    ap.add_argument("--hbn-bot-nm", type=float, default=3.0)
    ap.add_argument("--hbn-mid-nm", type=float, default=3.0)
    ap.add_argument("--hbn-top-nm", type=float, default=3.0)
    ap.add_argument("--eps-hbn", type=float, default=6.93)
    ap.add_argument("--eps-vac", type=float, default=1.0)
    ap.add_argument("--bottom_V", type=float, default=0.0)
    ap.add_argument("--top_V", type=float, default=None)
    ap.add_argument("--E-Vnm", type=float, default=0.2)
    ap.add_argument("--graphene-z-nm", type=float, default=6.0)
    ap.add_argument("--pattern", choices=["single", "grid", "hexagonal"], default="hexagonal")
    ap.add_argument("--qd-rows", type=int, default=4)
    ap.add_argument("--qd-cols", type=int, default=5)
    ap.add_argument("--spacing-nm", type=float, default=12.0)
    ap.add_argument("--qd-radius-nm", type=float, default=3.0)
    ap.add_argument("--epsr-profile", type=str, required=True)
    ap.add_argument("--lateral-bc", choices=["natural", "cutoff"], default="natural")
    ap.add_argument("--workdir", default="fea_out")
    args = ap.parse_args()

    Lx = args.Lx_nm * 1e-9; Ly = args.Ly_nm * 1e-9; Lz = args.Lz_nm * 1e-9
    hbn_bot = args.hbn_bot_nm * 1e-9
    hbn_mid = args.hbn_mid_nm * 1e-9
    hbn_top = args.hbn_top_nm * 1e-9
    zg = args.graphene_z_nm * 1e-9
    spacing = args.spacing_nm * 1e-9
    R = args.qd_radius_nm * 1e-9
    shells = [float(s) for s in args.epsr_profile.split(",")]
    eps_hbn = args.eps_hbn
    eps_vac = args.eps_vac

    if args.top_V is None:
        args.top_V = args.bottom_V + args.E_Vnm * args.Lz_nm

    zc = hbn_bot + hbn_mid + R
    centers = generate_centers(args.pattern, Lx, Ly, zc, spacing, args.qd_rows, args.qd_cols)

    xs = linspace_nodes(Lx, args.nx)
    ys = linspace_nodes(Ly, args.ny)
    zs = linspace_nodes(Lz, args.nz)
    m = MeshHex.init_tensor(xs, ys, zs)

    e = ElementHex1()
    basis = Basis(m, e)

    epsfun = build_epsfun(Lx, Ly, Lz, hbn_bot, hbn_mid, hbn_top, eps_hbn, centers, R, shells, eps_vac)
    eps_node = epsfun(m.p)
    eps_q = basis.interpolate(eps_node)

    from skfem import BilinearForm, LinearForm
    @BilinearForm
    def a(u, v, w):
        return w.eps * dot(grad(u), grad(v))
    @LinearForm
    def l(v, w):
        return 0.0 * v

    A = asm(a, basis, eps=eps_q)
    b = asm(l, basis)

    tol = 1e-12
    z_coords = m.p[2]
    bottom_mask = np.abs(z_coords - 0.0) < tol
    top_mask = np.abs(z_coords - Lz) < tol
    bottom_nodes = np.where(bottom_mask)[0]
    top_nodes = np.where(top_mask)[0]
    
    dofs = {}
    dofs.update({int(i): float(args.bottom_V) for i in bottom_nodes})
    dofs.update({int(i): float(args.top_V) for i in top_nodes})

    if args.lateral_bc == "cutoff":
        lateral_nodes = []
        lateral_nodes.extend(np.where(np.isclose(m.p[0], 0.0, atol=tol))[0])
        lateral_nodes.extend(np.where(np.isclose(m.p[0], Lx, atol=tol))[0])
        lateral_nodes.extend(np.where(np.isclose(m.p[1], 0.0, atol=tol))[0])
        lateral_nodes.extend(np.where(np.isclose(m.p[1], Ly, atol=tol))[0])
        lateral_nodes = [i for i in lateral_nodes if i < basis.N]
        dofs.update({int(i): 0.0 for i in lateral_nodes})

    if dofs:
        # Use interior DOF method that works with scikit-fem 11.0.0
        all_dofs = set(range(basis.N))
        boundary_set = set(dofs.keys())
        interior_dofs = np.array(list(all_dofs - boundary_set), dtype=np.intp)
        
        # Build full solution vector with boundary values
        x_full = np.zeros(basis.N)
        for idx, val in dofs.items():
            x_full[idx] = val
            
        from skfem import solve, enforce
        A_enforced, b_enforced = enforce(A, b, I=interior_dofs, x=x_full)
        phi = solve(A_enforced, b_enforced)
    else:
        phi = solve(A, b)

    node_z = m.p[2]
    uniq_z = np.unique(np.round(node_z, 15))
    k = int(np.argmin(np.abs(uniq_z - zg)))
    zsel = uniq_z[k]
    sel = np.where(np.isclose(node_z, zsel, atol=1e-12))[0]

    nxn, nyn = args.nx + 1, args.ny + 1
    layer_offset = k * (nxn * nyn)
    layer_nodes = np.arange(layer_offset, layer_offset + nxn * nyn)
    px = m.p[0, layer_nodes].reshape((nyn, nxn))
    py = m.p[1, layer_nodes].reshape((nyn, nxn))
    pphi = phi[layer_nodes].reshape((nyn, nxn))

    outdir = Path(args.workdir); outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "graphene_slice.csv", "w") as f:
        f.write("x_nm,y_nm,phi_V\n")
        for j in range(nyn):
            for i in range(nxn):
                f.write(f"{px[j,i]*1e9:.6f},{py[j,i]*1e9:.6f},{pphi[j,i]:.8e}\n")

    with open(outdir / "summary.txt", "w") as f:
        f.write(f"phi_min,phi_max,{phi.min():.6e},{phi.max():.6e}\n")
        f.write(f"Lx_nm,Ly_nm,Lz_nm,{Lx*1e9:.3f},{Ly*1e9:.3f},{Lz*1e9:.3f}\n")
        f.write(f"hbn_bot/mid/top_nm,{hbn_bot*1e9:.3f},{hbn_mid*1e9:.3f},{hbn_top*1e9:.3f}\n")
        f.write(f"qd_R_nm,{R*1e9:.3f}\n")
        f.write(f"epsr_shells,{','.join(map(str, shells))}\n")
        f.write(f"centers_N,{len(centers)}\n")

if __name__ == "__main__":
    main()
