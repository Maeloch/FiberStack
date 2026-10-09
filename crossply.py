"""
crossply.py — stratifié croisé : couches de fibres parallèles, petite
rotation et petite translation de l'axe à chaque couche (casser la
symétrie du réseau). Matériau de RÉFÉRENCE étalonné : pas de physique,
la porosité est contrôlée analytiquement par l'espacement s.

φ analytique = π·E[r²] / (s · p), p = FILTER_THICKNESS / n_layers.

v2 — revue de code :
  - L = DOMAIN_SIDE·√2 : couverture des coins à TOUTE rotation
    (v1 laissait les coins de la marge découverts aux couches ~45°)
  - reach = half·√2 + s : les lignes de centres couvrent la projection
    diagonale du domaine à toute rotation
  - assertion p >= FIBER_DIAMETER_HIGH, i.e. n_layers <= 12 : au-delà
    les couches s'interpénètrent PAR CONSTRUCTION (v1 le laissait
    faire en silence — les runs 14c étaient hors domaine de validité)
  - deux passes : géométrie des couches d'abord, diamètres
    précalculés ensuite (fini la réallocation à rallonge de v1)
  - export via le contrat unifié (route="crossply", spacing,
    n_layers, seed, n_total)

NB : « fraction solide comptée » suit la convention d'analyse.py
(volume intégral des fibres core, fibres coupant la frontière
comptées entières) — surcomptage attendu ; en v2 (L = 10√2) le ratio
comptée/analytique grimpe vers ~2,5 (v1 : ~1,66). La métrique de
référence pour les structures denses reste la métrique voxel.

Usage :
    python crossply.py <espacement_mm> <n_layers> [seed]
    python crossply.py 0.15 12
"""

import math
import sys

import numpy as np

from structure import (DOMAIN_SIDE, MEASURE_SIDE, FILTER_THICKNESS,
                       FIBER_DIAMETER_LOW, FIBER_DIAMETER_HIGH,
                       top_exact, save_pile)

CODE_VERSION = 2

GOLDEN = 0.61803398875      # décalage irrationnel de l'axe entre couches
LAYER_TURN_DEG = 5.0        # rotation de l'axe à chaque couche

# Lois de tirage : copie CONFORME de settle.py. Si l'une change,
# changer les DEUX — ou hoister dans structure.py.
_d_mu = 0.5 * (math.log(FIBER_DIAMETER_LOW) + math.log(FIBER_DIAMETER_HIGH))
_d_sigma = (math.log(FIBER_DIAMETER_HIGH) - math.log(FIBER_DIAMETER_LOW)) / 4.0


def sample_diameters(rng, n):
    out = np.empty(n)
    filled = 0
    while filled < n:
        draw = rng.lognormal(mean=_d_mu, sigma=_d_sigma, size=(n - filled) * 2)
        valid = draw[(draw >= FIBER_DIAMETER_LOW) & (draw <= FIBER_DIAMETER_HIGH)]
        take = min(len(valid), n - filled)
        out[filled:filled + take] = valid[:take]
        filled += take
    return out


def build_crossply(spacing, n_layers, seed=0, layer_turn_deg=LAYER_TURN_DEG):
    """Retourne (fibers, pitch, n_per_layer). Chaque fibre : dict
    center/direction/length/radius/core (mêmes clés que settle)."""
    pitch = FILTER_THICKNESS / n_layers
    # tolérance 1e-9 : 0.3/12 = 0.024999999999999998 en IEEE-754, alors
    # que le littéral 25e-3 vaut 0.025000000000000001 — sans epsilon,
    # l'égalité exacte pas = d_max (le cas 12 couches, valide par
    # construction) échouait sur des poussières de bits
    if pitch < FIBER_DIAMETER_HIGH - 1e-9:
        raise ValueError(
            f"n_layers={n_layers} -> pas vertical {pitch*1e3:.1f} µm < "
            f"d_max {FIBER_DIAMETER_HIGH*1e3:.0f} µm : les couches "
            f"s'interpénètrent par construction (maximum 12 couches)")

    half = DOMAIN_SIDE / 2.0
    L = DOMAIN_SIDE * math.sqrt(2)          # v2 : couverture des coins
    reach = half * math.sqrt(2) + spacing    # v2 : projection diagonale
    n_per_layer = int(math.ceil(2.0 * reach / spacing))

    # -------- passe 1 : géométrie des couches (aucun tirage) --------
    layers = []
    for k in range(n_layers):
        theta = math.radians(k * layer_turn_deg)
        u_dir = np.array([math.cos(theta), math.sin(theta), 0.0])
        n_dir = np.array([-math.sin(theta), math.cos(theta), 0.0])
        # petite translation de l'axe de rotation, couche par couche
        shift = (((k * GOLDEN) % 1.0) - 0.5) * spacing
        z_k = (k + 0.5) * pitch
        layers.append((k, u_dir, n_dir, shift, z_k))

    n_total = n_layers * n_per_layer

    # -------- passe 2 : diamètres précalculés en bloc --------
    rng = np.random.default_rng(seed)
    diameters = sample_diameters(rng, n_total)

    half_m = MEASURE_SIDE / 2.0
    fibers = []
    idx = 0
    for (k, u_dir, n_dir, shift, z_k) in layers:
        for j in range(n_per_layer):
            r = float(diameters[idx]) / 2.0
            idx += 1
            # décalage perpendiculaire (le réseau) + jitter axial
            # (casser la collinéarité parfaite, sans toucher à φ)
            u = (j - (n_per_layer - 1) / 2.0) * spacing + shift
            a = rng.uniform(-0.25, 0.25) * spacing
            center = a * u_dir + u * n_dir
            center = np.array([center[0], center[1], z_k])
            core = (abs(center[0]) <= half_m) and (abs(center[1]) <= half_m)
            fibers.append(dict(center=center, direction=u_dir.copy(),
                               length=L, radius=r, core=core))
    return fibers, pitch, n_per_layer


def main():
    if len(sys.argv) < 3:
        print("usage : python crossply.py <espacement_mm> <n_layers> [seed]")
        sys.exit(1)
    spacing = float(sys.argv[1])
    n_layers = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    fibers, pitch, n_per_layer = build_crossply(spacing, n_layers, seed)

    centers = np.array([f["center"] for f in fibers])
    directions = np.array([f["direction"] for f in fibers])
    lengths = np.array([f["length"] for f in fibers])
    radii = np.array([f["radius"] for f in fibers])
    core = np.array([f["core"] for f in fibers])

    z_max = max(top_exact(c, d, L, r)
                for c, d, L, r, ck in zip(centers, directions, lengths,
                                           radii, core) if ck)
    vols = np.pi * radii ** 2 * lengths
    phi_counted = vols[core].sum() / (MEASURE_SIDE ** 2 * z_max)
    # analytique : fraction solide réelle du réseau (volume in-domaine)
    phi_exact = np.pi * (radii ** 2).mean() / (spacing * pitch)

    fname = f"crossply_{int(spacing*1000)}um_{n_layers}c.npz"
    save_pile(fname, centers, directions, lengths, radii, core,
              route="crossply", code_version=CODE_VERSION,
              spacing=spacing, n_layers=n_layers, seed=seed,
              n_total=len(fibers))

    print(f"=== {len(fibers)} fibres ({n_per_layer}/couche x {n_layers}), "
          f"{int(core.sum())} core, hauteur {z_max:.3f} mm, "
          f"fraction solide comptée {phi_counted*100:.1f} % "
          f"(convention analyse, surcomptée), "
          f"analytique {phi_exact*100:.1f} % -> {fname} ===")


if __name__ == "__main__":
    main()
