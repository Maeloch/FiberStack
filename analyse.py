"""
analyze_pile.py — exploitation de pile.npz : porosité locale par tranches,
distribution des orientations, occupation surfacique.

Usage :
    python analyze_pile.py                # -> pile.npz par défaut
    python analyze_pile.py pile3.npz      # fichier explicite

Conventions : unités (mm, mg, ms). |cos(theta)| proche de 0 = fibre
COUCHÉE (axe horizontal) ; proche de 1 = fibre VERTICALE. L'élévation
par rapport à l'horizontale = 90° - arccos(|cos theta|).
"""

import sys
import numpy as np

npz_path = sys.argv[1] if len(sys.argv) > 1 else "pile.npz"
data = np.load(npz_path)
centers    = data["centers"]      # (N,3)
directions = data["directions"]  # (N,3) unitaires
lengths    = data["lengths"]      # (N,)
radii      = data["radii"]        # (N,)
core       = data["core"]          # (N,) bool

MEASURE_SIDE = 7.0     # mm, même constante que settle_si4
half_m = MEASURE_SIDE / 2

# ---------- 1) volumes et porosité globale (zone de mesure) ----------
volumes = np.pi * radii**2 * lengths
core_mask = core & (np.abs(centers[:,0]) <= half_m) & (np.abs(centers[:,1]) <= half_m)
z_max = centers[core_mask][:,2].max() + lengths[core_mask].max()/2
domain_volume = MEASURE_SIDE**2 * z_max
porosity = 1 - volumes[core_mask].sum() / domain_volume
print(f"fichier : {npz_path}")
print(f"fibres core : {core_mask.sum()}")
print(f"hauteur du tas (core) : {z_max:.3f} mm")
print(f"porosité globale : {porosity*100:.2f} %")

# ---------- 2) porosité par tranche verticale ----------
# Approximation : le volume de chaque fibre est réparti sur les tranches
# couvertes par son axe (proxy, surestime légèrement le solide aux bords).
n_slices = 20
z_edges = np.linspace(0, z_max, n_slices + 1)
slice_vol = (MEASURE_SIDE**2) * (z_max / n_slices)
c = centers[core_mask]; l = lengths[core_mask]; v = volumes[core_mask]
e1 = c - 0.5*l[:,None]*directions[core_mask]
e2 = c + 0.5*l[:,None]*directions[core_mask]
zmin_f = np.minimum(e1[:,2], e2[:,2]); zmax_f = np.maximum(e1[:,2], e2[:,2])
solid_by_slice = np.zeros(n_slices)
for k in range(n_slices):
    lo, hi = z_edges[k], z_edges[k+1]
    overlap = np.minimum(zmax_f, hi) - np.maximum(zmin_f, lo)
    frac = np.clip(overlap / np.maximum(zmax_f - zmin_f, 1e-12), 0, 1)
    solid_by_slice[k] = (v * frac).sum()
poro_by_slice = 1 - solid_by_slice / slice_vol
print("\nprofil de porosité (tranche, z, porosité) :")
for k in range(n_slices):
    print(f"  {k:2d}  z={0.5*(z_edges[k]+z_edges[k+1]):6.3f} mm  "
          f"porosité={poro_by_slice[k]*100:7.2f} %")

# ---------- 3) orientations ----------
# |cos(theta)| : 0 = fibre à plat, 1 = fibre verticale.
u = directions[core_mask]
cos_theta = np.abs(u[:,2])
elevation = 90.0 - np.degrees(np.arccos(np.clip(cos_theta, 0, 1)))
print(f"\nélévation par rapport à l'horizontale : "
      f"médiane={np.median(elevation):.2f}°, "
      f"90%={np.quantile(elevation,0.9):.2f}°")
print("un C569 réel attendu : quelques degrés de médiane, une queue "
      "jusqu'à 15-20° (chevauchements en toit)")

# ---------- 4) occupation surfacique (footprint) ----------
# fraction du plan couverte par les projections : proxy du "voile"
cell = 0.05  # mm
nx = int(MEASURE_SIDE / cell)
occ = np.zeros((nx, nx), dtype=bool)
e1c = centers[core_mask] - 0.5*lengths[core_mask,None]*directions[core_mask]
e2c = centers[core_mask] + 0.5*lengths[core_mask,None]*directions[core_mask]
for a, b in zip(e1c, e2c):
    ts = np.linspace(0, 1, 50)
    pts = a[None,:] + ts[:,None]*(b - a)[None,:]
    ixs = ((pts[:,0] + half_m) / cell).astype(int)
    iys = ((pts[:,1] + half_m) / cell).astype(int)
    ok = (ixs>=0)&(ixs<nx)&(iys>=0)&(iys<nx)
    occ[ixs[ok], iys[ok]] = True
print(f"\noccupation surfacique (maille {cell*1000:.0f} µm) : "
      f"{occ.mean()*100:.1f} %")
print("complément du 'voile' vu du dessus — monte vers 60-80 % à "
      "densité C569 (~19 % de fraction solide)")
