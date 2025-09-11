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
from pathlib import Path

# Path to your QD .xyz file
xyz = "data/dot_demo.xyz"
if not Path(xyz).exists():
    exit(1)

dft_cmd = [
    "python", "scripts/DFT.py", xyz,
    "--Lz-nm", "16", "--Efield-Vnm", "0.005", "--shells", "6",
    "--ecut-eV", "200", "--kpts", "1", "1", "1", "--float32",
    "--spacing-mult", "2.2", "--out-json", "demo_epsr.json"
]
sub = subprocess.run(dft_cmd, capture_output=True, text=True)
if sub.returncode != 0:
    exit(1)

with open("demo_epsr.json") as f:
    lines = f.read().strip().split('\n')
    info = json.loads(lines[0])

epsilon_eff = info["epsilon_eff"]
radius_nm = info["used_radius_nm"]

if len(lines) > 1:
    shell_line = None
    for line in lines:
        if ',' in line and all(c.replace('.','').replace(',','').isdigit() or c==',' for c in line.strip()):
            shell_line = line.strip()
            break
    profile = np.array([float(x) for x in shell_line.split(',')]) if shell_line else np.array([epsilon_eff] * 6)
else:
    profile = np.array([epsilon_eff] * 6)

edges = np.linspace(0, radius_nm, len(profile) + 1)
r_centers = 0.5*(edges[:-1] + edges[1:])
plt.figure(figsize=(6,4))
plt.plot(r_centers, profile, marker='o')
plt.xlabel("Radius (nm)")
plt.ylabel("εᵣ")
plt.title("Radial Dielectric Profile")
plt.grid(True)
plt.tight_layout()
plt.savefig("demo_radial_profile.png", dpi=150)

nx, ny = len(profile)*4, len(profile)*4
edges3d = np.linspace(0, radius_nm, nx)
X, Y = np.meshgrid(edges3d, edges3d)
R_xy = np.sqrt(X**2 + Y**2)
eps2d = np.ones_like(R_xy)
for i, r0 in enumerate(edges[1:]):
    mask = R_xy <= r0
    eps2d[mask] = profile[i]

plt.figure(figsize=(6,4))
im = plt.imshow(eps2d, origin="lower", extent=[0, radius_nm, 0, radius_nm], cmap="viridis")
plt.colorbar(im, label="εᵣ")
plt.xlabel("x (nm)"); plt.ylabel("y (nm)")
plt.title("2D Dielectric Cross-Section")
plt.tight_layout()
plt.savefig("demo_2D_profile.png", dpi=150)

print("[demo] Saved 2D εᵣ slice → demo_2D_profile.png")

print("[demo] Done.")
