# Stack we simulate

From bottom (z=0) to top (z=Lz) in capacitor with gap L

1. **hBN stakc (bottom)**
2. **2D quantum-dot array** (each dot is a dielectric sphere with radius-dependent εr you supply or compute via DFT.py)
3. **hBN stack (middle spacer)**
4. **Graphene sheet** (readout plane; we just sample φ at this z, no self-consistent graphene in this repo)
5. **hBN stack (top)**
6. **Capacitor plates** at z=0 (bottom\_V) and z=Lz (top\_V = E·Lz)

We simulate the periodic system via unit cell. We can apply these boundary conditions: **periodic** (default) or **cutoff** (Dirichlet at x/y sides).

---

## Files

* `FEA.py` — Poisson solver
  * Input: quantum dot radius, pattern of dots, spacing between dot centers, thickness of HBN stacks, incident electric field, list of uniform dielectric constant intervals spanning radius
  * Output `graphene_slice.csv` (potential on graphene), `summary.txt` (max and avg displacemnt field in both HBN stacks).
* `DFT.py` — **Optional** estimates a radial εr profile from a dot `.xyz` file using Density Functional Theory

  * Uses ASE+GPAW (you can swap for your own DFT; only the **comma list** matters here).
* `run.py` — Run workflow

---

## Quick start (demo / mock)

Run with a **manual εr shell profile** (center→edge):

```bash
# Hexagonal array, spacing = 4×diameter (default), 3 nm radius dots
python run.py \
  --epsr-profile 14,12,10,8,6 \
  --pattern hexagonal --qd-rows 4 --qd-cols 5 \
  --qd-radius-nm 3 --spacing-factor 4 \
  --hbn-bot-nm 3 --hbn-mid-nm 3 --hbn-top-nm 3 \
  --Lx-nm 40 --Ly-nm 40 --Lz-nm 30 \
  --E-Vnm 0.2 \
  --lateral-bc natural \
  --workdir out_demo
```

Outputs in `out_demo/`:

* `graphene_slice.csv` (x\_nm, y\_nm, phi\_V)
* `summary.txt`
* `graphene_potential_2D.png` (with QD overlay)
* `graphene_potential_3D.png`

---

## Full workflow (with DFT)

Requires ASE+GPAW installed.

```bash
# Use .xyz to infer the radial εr shells via DFT helper, then run FEA
python run.py \
  --xyz dot.xyz \
  --dft-shells 6 --dft-E-Vnm 0.01 \
  --pattern grid --qd-rows 3 --qd-cols 3 \
  --qd-radius-nm 3 --spacing-factor 3.5 \
  --hbn-bot-nm 3 --hbn-mid-nm 3 --hbn-top-nm 3 \
  --Lx-nm 40 --Ly-nm 40 --Lz-nm 30 \
  --E-Vnm 0.2 \
  --lateral-bc natural \
  --workdir out_dft
```

Tip: if you **already** know a good shell list from your own DFT, just pass `--epsr-profile` and skip `--xyz`.

---

## Inputs (most detailed)

* **Geometry / grid**

  * `--Lx-nm, --Ly-nm, --Lz-nm` (box size), `--nx, --ny, --nz` (grid)
* **Stack**

  * `--hbn-bot-nm, --hbn-mid-nm, --hbn-top-nm` (thicknesses)
  * `--graphene-z-nm` (z where φ is sampled; typically bottom+mid+hints)
  * `--eps-hbn` (default 6.93), `--eps-vac` (default 1.0)
* **Field / BC**

  * `--E-Vnm` → sets `top_V = E * Lz`
  * `--lateral-bc {natural|cutoff}`
* **Dots / pattern**

  * `--qd-radius-nm`
  * `--spacing-nm` **or** `--spacing-factor` (spacing = factor × diameter if `--spacing-nm` omitted)
  * `--pattern {single|grid|hexagonal}`, `--qd-rows`, `--qd-cols`, `--hex-angle`
* **Dielectric profile**

  * **Manual**: `--epsr-profile "14,12,10,8,6"`
  * **From DFT**: `--xyz dot.xyz`, `--dft-shells`, `--dft-E-Vnm`

---

## Outputs

* `graphene_slice.csv` — sampled potential at graphene plane (z = `graphene-z-nm`)
* `graphene_potential_2D.png` — 2D φ with **QD circles** overlaid
* `graphene_potential_3D.png` — 3D surface of φ(x,y)
* `summary.txt` — quick min/max and small stats

---

## Example configs

### 1) Minimal demo (hex pattern, manual εr shells)

```bash
python run.py \
  --epsr-profile 16,14,11,8,6 \
  --pattern hexagonal --qd-rows 5 --qd-cols 6 \
  --qd-radius-nm 3 --spacing-factor 4 \
  --hbn-bot-nm 2.5 --hbn-mid-nm 3 --hbn-top-nm 2.5 \
  --Lx-nm 48 --Ly-nm 48 --Lz-nm 30 \
  --E-Vnm 0.15 --lateral-bc natural \
  --workdir out_hex_demo
```
