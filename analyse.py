"""
analyse.py — exploitation des fichiers pile .npz : porosité globale,
profil de porosité par tranches, orientations, occupation surfacique.

v2.1 — correctif occupation : échantillonnage ADAPTATIF le long de
  l'axe (au moins 2 points par maille, plancher 50). Les 50 points
  fixes sous-échantillonnaient les fibres longues : à L = 10 mm un
  point tous les 200 µm, à L = 10√2 tous les 283 µm — la majorité des
  mailles traversées n'était pas marquée. Les occupations crossply
  v1 (39,1/63,8/86,6 %) étaient des artefacts de ce biais : la courbe
  occupation-φ doit être REMESURÉE. Les tas settle (L <= 1,3 mm)
  restent à 50 points : références inchangées.

v2 — revue de code :
  - MEASURE_SIDE importé de structure.py (fini le littéral 7.0 : si le
    domaine change, l'analyse suit au lieu de mesurer une zone fantôme)
  - métadonnées du npz affichées au chargement (route, seed, n_total...)
  - hauteur EXACTE (extension verticale réelle de chaque capsule) :
    re-baselinage fait — 300 fibres : 0,150 mm ; 1000 : 1,002 mm
  - fibres ponctuelles en z affectées à la tranche de leur axe
    (v1 les comptait zéro : toute structure parfaitement horizontale
    était invisible du profil — le crossply l'a révélé)
  - avertissement explicite quand tranche < 3 diamètres : la métrique
    par tranches est HORS DOMAINE pour les structures denses ; la
    métrique voxel la remplacera (voxelize.py, à venir)

Usage :
    python analyse.py                    # -> pile.npz par défaut
    python analyse.py pile_1000_0.npz
"""

import sys
import numpy as np

from structure import load_pile, MEASURE_SIDE

npz_path = sys.argv[1] if len(sys.argv) > 1 else "pile.npz"
data, meta = load_pile(npz_path)
centers    = data["centers"]      # (N,3)
directions = data["directions"]   # (N,3) unitaires
lengths    = data["lengths"]      # (N,)
radii      = data["radii"]        # (N,)
core       = data["core"]         # (N,) bool
half_m = MEASURE_SIDE / 2

# ---------- 1) volumes et porosité globale (zone de mesure) ----------
# NB : les fibres qui coupent la frontière de la zone versent leur
# volume intégral (surcomptage ~1.4-1.7 selon la structure — mesuré
# sur crossply). Pour les structures denses, la porosité de confiance
# est celle du module voxel ; cette valeur reste un proxy comparable
# entre structures de même géométrie.
volumes = np.pi * radii**2 * lengths
core_mask = core & (np.abs(centers[:,0]) <= half_m) \
                  & (np.abs(centers[:,1]) <= half_m)
z_ext = 0.5 * lengths * np.abs(directions[:, 2])
z_max = (centers[:, 2] + z_ext + radii)[core_mask].max()
domain_volume = MEASURE_SIDE**2 * z_max
porosity = 1 - volumes[core_mask].sum() / domain_volume
print(f"fibres core : {core_mask.sum()}")
print(f"hauteur du tas (core) : {z_max:.3f} mm")
print(f"porosité globale : {porosity*100:.2f} %")

# ---------- 2) porosité par tranche verticale ----------
n_slices = 20
slice_h = z_max / n_slices
d_mean = 2 * radii[core_mask].mean()
if slice_h < 3 * d_mean:
    print(f"\nAVERTISSEMENT : tranches de {slice_h*1e3:.0f} µm < 3 diamètres "
          f"({d_mean*1e3:.0f} µm) — le profil par tranches est HORS DOMAINE "
          f"DE VALIDITÉ, les valeurs solides sont surcomptées. "
          f"Structure dense : utiliser la métrique voxel.")
z_edges = np.linspace(0, z_max, n_slices + 1)
slice_vol = (MEASURE_SIDE**2) * slice_h
c = centers[core_mask]; l = lengths[core_mask]; v = volumes[core_mask]
e1 = c - 0.5*l[:,None]*directions[core_mask]
e2 = c + 0.5*l[:,None]*directions[core_mask]
zmin_f = np.minimum(e1[:,2], e2[:,2]); zmax_f = np.maximum(e1[:,2], e2[:,2])
span = zmax_f - zmin_f
degenerate = span < 1e-9             # fibres ponctuelles en z (horizontales)
mid_z = 0.5 * (zmin_f + zmax_f)
solid_by_slice = np.zeros(n_slices)
for k in range(n_slices):
    lo, hi = z_edges[k], z_edges[k+1]
    frac = np.zeros(len(c))
    frac[~degenerate] = np.clip(
        (np.minimum(zmax_f[~degenerate], hi)
         - np.maximum(zmin_f[~degenerate], lo))
        / np.maximum(span[~degenerate], 1e-12), 0, 1)
    frac[degenerate] = ((mid_z[degenerate] >= lo)
                       & (mid_z[degenerate] < hi))
    solid_by_slice[k] = (v * frac).sum()
poro_by_slice = 1 - solid_by_slice / slice_vol
print("\nprofil de porosité (tranche, z, porosité) :")
for k in range(n_slices):
    print(f"  {k:2d}  z={0.5*(z_edges[k]+z_edges[k+1]):6.3f} mm  "
          f"porosité={poro_by_slice[k]*100:7.2f} %")

# ---------- 3) orientations ----------
u = directions[core_mask]
cos_theta = np.abs(u[:,2])
elevation = 90.0 - np.degrees(np.arccos(np.clip(cos_theta, 0, 1)))
print(f"\nélévation par rapport à l'horizontale : "
      f"médiane={np.median(elevation):.2f}°, "
      f"90%={np.quantile(elevation,0.9):.2f}°")
print("un C569 réel attendu : quelques degrés de médiane, une queue "
      "jusqu'à 15-20° (chevauchements en toit)")

# ---------- 4) occupation surfacique (footprint) ----------
cell = 0.05  # mm
nx = int(MEASURE_SIDE / cell)
occ = np.zeros((nx, nx), dtype=bool)
e1c = centers[core_mask] - 0.5*lengths[core_mask,None]*directions[core_mask]
e2c = centers[core_mask] + 0.5*lengths[core_mask,None]*directions[core_mask]
for a, b in zip(e1c, e2c):
    seg = b - a
    L_f = float(np.linalg.norm(seg))
    # v2.1 : échantillonnage adaptatif — au moins 2 points par maille,
    # plancher 50 (les tas settle, L <= 1,3 mm, restent à 50 points,
    # leurs références d'occupation ne changent pas)
    n_s = max(50, int(np.ceil(2.0 * L_f / cell)) + 1)
    ts = np.linspace(0, 1, n_s)
    pts = a[None,:] + ts[:,None]*seg[None,:]
    ixs = ((pts[:,0] + half_m) / cell).astype(int)
    iys = ((pts[:,1] + half_m) / cell).astype(int)
    ok = (ixs>=0)&(ixs<nx)&(iys>=0)&(iys<nx)
    occ[ixs[ok], iys[ok]] = True
print(f"\noccupation surfacique (maille {cell*1000:.0f} µm) : "
      f"{occ.mean()*100:.1f} %")
print("proxy du 'voile' vu de dessus — marque les cellules traversées")
print("par l'AXE (rayon ignoré) : comparable entre structures simulées,")
print("mais PAS directement à une occupation µCT (fraction solide par")
print("pixel) — ne pas opposer ce chiffre à une image sans conversion.")
