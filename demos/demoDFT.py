#!/usr/bin/env python3
"""
Demo: run DFT → extract radial εr profile, then visualize it.

- Uses fast DFT settings for speed/demo.
- Plots:
   * εr vs radius (center → edge)
   * 2D cross-section of the dielectric sphere
"""

import subprocess, json
import numpy as np
import matplotlib.pyplot as plt
import os

# Path to your QD .xyz file
xyz = r"C:\Users\micha\BaiLabDFT\Cd68_OPT.xyz"

# Fast/demo DFT parameters
dft_cmd = [
    "python", "DFT.py", xyz,
    "--Lx-nm", "16", "--Ly-nm", "16", "--Lz-nm", "24",
    "--Efield-Vnm", "0.005",         # small field
    "--shells", "10",
    "--ecut-eV", "200",              # lower cutoff
    "--kpts", "2", "2", "1",
    "--sigma-vox", "0.5",
    "--out-json", "demo_epsr.json",
    "--out-txt", "demo_epsr.txt"
]
print("[demo] Running DFT:", " ".join(dft_cmd))
sub = subprocess.run(dft_cmd, capture_output=True, text=True)
if sub.returncode != 0:
    print("[ERROR] DFT failed:\n", sub.stderr)
    exit(1)

# Load profile
with open("demo_epsr.json") as f:
    info = json.load(f)
profile = np.array(info["epsr_profile_center_to_edge"])
edges = np.array(info["edges_nm"])

# Plot εr vs radius
r_centers = 0.5*(edges[:-1] + edges[1:])
plt.figure(figsize=(6,4))
plt.plot(r_centers, profile, marker='o')
plt.xlabel("Radius (nm)")
plt.ylabel("εᵣ")
plt.title("Radial Dielectric Profile")
plt.grid(True)
plt.tight_layout()
plt.savefig("demo_radial_profile.png", dpi=150)
print("[demo] Saved εᵣ vs radius → demo_radial_profile.png")

# Generate and plot 2D slice through center
nx, ny, nz = len(profile)*4, len(profile)*4, nz  # fine grid
edges3d = np.linspace(0, info["used_radius_nm"], nx)
X, Y = np.meshgrid(edges3d, edges3d)
R_xy = np.sqrt(X**2 + Y**2)

# Assign εᵣ by shell
eps2d = np.zeros_like(R_xy)
for i, r0 in enumerate(aux := edges[1:]):
    mask = R_xy <= r0
    eps2d[mask] = profile[i]

plt.figure(figsize=(5,4))
plt.imshow(eps2d, extent=[0, edges3d.max(), 0, edges3d.max()],
           origin='lower', cmap='viridis')
plt.colorbar(label="εᵣ")
plt.title("2D Cross-section (top view)")
plt.xlabel("x (nm)"); plt.ylabel("y (nm)")
plt.tight_layout()
plt.savefig("demo_eps2d_slice.png", dpi=150)
print("[demo] Saved 2D εᵣ slice → demo_eps2d_slice.png")

print("[demo] Done.")
