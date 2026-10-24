"""
analyse.py — exploitation des fichiers pile .npz : porosité globale,
profil de porosité par tranches, orientations, occupation surfacique.

v2 — revue de code :
  - MEASURE_SIDE importé de structure.py (fini le littéral 7.0 : si le
    domaine change, l'analyse suit au lieu de mesurer une zone fantôme)
  - métadonnées du npz affichées au chargement (route, seed, n_total...)
  - hauteur EXACTE (extension verticale réelle de chaque capsule) :
    les hauteurs v2 ne sont PAS comparables aux logs v5, qui gonflaient
    de L/2 par le pire-cas vertical — re-baselinage à faire au Memo
  - fibres ponctuelles en z affectées à la tranche de leur axe
    (v1 les comptait zéro : toute structure parfaitement horizontale
    était invisible du profil — le crossply l'a révélé)
  - avertissement explicite quand tranche < 3 diamètres : la métrique
    par tranches est HORS DOMAINE pour les structures denses ; la
    métrique voxel la remplacera (voxelize.py, à venir)

v2.2 — tuiles pressMT : la zone de mesure est lue dans les
métadonnées du npz (tile_side, burn_in) au lieu du MEASURE_SIDE du
monolithe. Avant ce patch, les porosités des tuiles étaient calculées
sur une zone 7x7 fantôme (99,8 % affichés) alors que le cœur utile
fait 2,5x2,5 ; toute porosité « globale » de tuile sortie avant v2.2
est à rejeter.

Usage :
    python analyse.py                    # -> pile.npz par défaut
    python analyse.py pile_1000_0.npz
    python analyse.py pressedMT_1000_2_t00.npz
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

# zone de mesure : monolithe -> MEASURE_SIDE ; tuile pressMT -> zone
# utile (tile_side - 2*burn_in), mêmes bornes que pressmt.export_tile
if int(meta.get("pressmt", 0)) == 1 and "tile_side" in meta:
    zone_side = float(meta["tile_side"]) - 2.0 * float(meta.get("burn_in", 0.0))
    print(f"zone de mesure pressMT : {zone_side:.2f} mm de côté "
          f"(tuile {float(meta['tile_side']):.2f} - 2x"
          f"{float(meta.get('burn_in', 0.0)):.2f} de rebut)")
elif "core_half" in meta:
    zone_side = 2.0 * float(meta["core_half"])
    print(f"zone de mesure (core_half, route D) : {zone_side:.2f} mm de côté")
else:
    zone_side = MEASURE_SIDE
half_m = zone_side / 2

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
domain_volume = zone_side**2 * z_max
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
slice_vol = (zone_side**2) * slice_h
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
nx = int(zone_side / cell)
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
print("proxy du 'voile' vu de dessus — marque les cellules traversées")
print("par l'AXE (rayon ignoré) : comparable entre structures simulées,")
print("mais PAS directement à une occupation µCT (fraction solide par")
print("pixel) — ne pas opposer ce chiffre à une image sans conversion.")
