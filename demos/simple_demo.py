#!/usr/bin/env python3
# simple_demo.py - Simplified quantum dot simulation for debugging

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def simple_potential_demo():
    """Create a simple 2D potential map to demonstrate the workflow."""
    
    # Grid setup (in nm)
    Lx_nm = 40.0
    Ly_nm = 40.0
    nx = 96
    ny = 96
    
    # Create coordinate grids
    x_nm = np.linspace(0, Lx_nm, nx)
    y_nm = np.linspace(0, Ly_nm, ny)
    X, Y = np.meshgrid(x_nm, y_nm)
    
    # QD parameters
    qd_centers = [(20.0, 15.0), (25.0, 25.0), (15.0, 30.0)]  # (x, y) in nm
    qd_radius = 3.0  # nm
    
    # Create potential field
    phi = np.zeros_like(X)
    
    # Background electric field (V/nm * nm = V)
    E_field = 0.2  # V/nm
    Lz_nm = 30.0
    background_potential = E_field * 20  # Approximate field effect
    phi += background_potential * (Y / Ly_nm)  # Linear field
    
    # Add QD effects (simplified dielectric response)
    for cx, cy in qd_centers:
        # Distance from each QD center
        r = np.sqrt((X - cx)**2 + (Y - cy)**2)
        
        # Simple Gaussian potential perturbation
        qd_mask = r <= qd_radius
        phi[qd_mask] += 0.5 * np.exp(-(r[qd_mask] / qd_radius)**2)
    
    # Add some noise for realism
    phi += 0.1 * np.random.normal(0, 0.1, phi.shape)
    
    return X, Y, phi, qd_centers

def plot_and_save_demo():
    """Generate demo plots and save them."""
    X, Y, phi, qd_centers = simple_potential_demo()
    
    # Create output directory
    outdir = Path("out_demo")
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Save CSV data
    with open(outdir / "graphene_slice.csv", "w") as f:
        f.write("x_nm,y_nm,phi_V\\n")
        for j in range(phi.shape[0]):
            for i in range(phi.shape[1]):
                f.write(f"{X[j,i]:.6f},{Y[j,i]:.6f},{phi[j,i]:.8e}\\n")
    
    # Save summary
    with open(outdir / "summary.txt", "w") as f:
        f.write(f"phi_min,phi_max,{phi.min():.6e},{phi.max():.6e}\\n")
        f.write(f"Lx_nm,Ly_nm,Lz_nm,40.000,40.000,30.000\\n")
        f.write(f"hbn_bot/mid/top_nm,3.000,3.000,3.000\\n")
        f.write(f"qd_R_nm,3.000\\n")
        f.write(f"epsr_shells,14,12,10,8,6\\n")
        f.write(f"centers_N,{len(qd_centers)}\\n")
    
    # 2D potential plot with QD overlay
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    im = ax.imshow(phi, origin="lower", extent=[0, 40, 0, 40], cmap="RdBu_r")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("Graphene potential (2D) with QD overlay")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("V")
    
    # Add QD circles
    for cx, cy in qd_centers:
        circ = plt.Circle((cx, cy), 3.0, fill=False, color="k", lw=0.8)
        ax.add_patch(circ)
        ax.plot([cx], [cy], "k.", ms=2)
    
    fig.tight_layout()
    fig.savefig(outdir / "graphene_potential_2D.png", dpi=160)
    
    # 3D surface plot
    fig3 = plt.figure(figsize=(7.8, 6.0))
    ax3 = fig3.add_subplot(111, projection="3d")
    ax3.plot_surface(X, Y, phi, linewidth=0, antialiased=True, cmap="viridis")
    ax3.set_xlabel("x (nm)")
    ax3.set_ylabel("y (nm)")
    ax3.set_zlabel("V")
    ax3.set_title("Graphene potential (3D)")
    fig3.tight_layout()
    fig3.savefig(outdir / "graphene_potential_3D.png", dpi=160)
    
    print("Demo completed successfully!")
    print(f"Output files saved to: {outdir}")
    print("Files created:")
    print("  - graphene_slice.csv (potential data)")
    print("  - summary.txt (simulation parameters)")
    print("  - graphene_potential_2D.png (2D plot with QD overlay)")
    print("  - graphene_potential_3D.png (3D surface plot)")

if __name__ == "__main__":
    plot_and_save_demo()
