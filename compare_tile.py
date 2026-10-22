"""
compare_tile.py — certification de la route tuiles contre le monolithe.

Question : le cœur d'une tuile pressMT (zone utile 2,5×2,5 par
défaut) produit-il la même matière que la zone centrale de même côté
d'un monolithe press, à même seed et même box-z ?

DISCIPLINE OPÉRATEUR (le npz ne stocke pas box-z) : pour que la
comparaison isolse l'effet des murs, le monolithe doit être lancé avec
--box-z 1.5 ET le même seed que les tuiles —
    python press.py --random 1000 --seed 2 --target 0.3 --box-z 1.5
À box égal, brouillard volumique (fibres/mm³) et nombre d'incréments
(49) sont identiques ; il ne reste comme différence que la distance
des murs (±2,5 tuile vs ±5,6 monolithe) — c'est l'objet du test. Si
l'écart de fraction solide persiste, le balayage tile_side 5 -> 7,5
mesure la longueur de sanglage.

Métriques par zone :
  - fraction solide φ = V_fibres / (side² × z_max). Biais de
    surcomptage aux frontières (fibres coupant la zone versent tout
    leur volume) identique des deux côtés : comparable.
  - hauteur z_max, élévation des axes (médiane, q90)
  - profil de porosité en 20 tranches verticales, abscisse en
    profondeur RELATIVE z/z_max — la forme se compare même si les
    hauteurs diffèrent

Usage :
    python compare_tile.py --monolith pressed_random_1000_2.npz \
        --tiles "pressedMT_1000_2_t*.npz"
    # zone centrale du monolithe = côté de la zone utile des tuiles
    # (ou --side 2.5 pour l'imposer)

Système d'unités : (mm, mg, ms). Nécessite : numpy, structure.py.
"""

import argparse
import glob

import numpy as np

from structure import load_pile


def zone_stats(path, half):
    """Stats de la zone carrée |x|,|y| <= half d'un npz, + profil de
    porosité en profondeur relative (logique analyse.py v2.2)."""
    data, meta = load_pile(path)
    centers = data["centers"]
    directions = data["directions"]
    lengths = data["lengths"]
    radii = data["radii"]
    m = ((np.abs(centers[:, 0]) <= half)
         & (np.abs(centers[:, 1]) <= half))
    c, d = centers[m], directions[m]
    l, r = lengths[m], radii[m]
    v = np.pi * r ** 2 * l
    z_ext = 0.5 * l * np.abs(d[:, 2])
    z_max = float((c[:, 2] + z_ext + r).max())
    side = 2.0 * half
    phi = float(v.sum() / (side ** 2 * z_max))
    elev = np.degrees(np.arcsin(np.clip(np.abs(d[:, 2]), 0, 1)))

    n_slices = 20
    e1 = c - 0.5 * l[:, None] * d
    e2 = c + 0.5 * l[:, None] * d
    zmin_f = np.minimum(e1[:, 2], e2[:, 2])
    zmax_f = np.maximum(e1[:, 2], e2[:, 2])
    span = zmax_f - zmin_f
    degen = span < 1e-9              # fibres ponctuelles en z
    mid_z = 0.5 * (zmin_f + zmax_f)
    edges = np.linspace(0.0, z_max, n_slices + 1)
    solid = np.zeros(n_slices)
    for k in range(n_slices):
        lo, hi = edges[k], edges[k + 1]
        frac = np.zeros(len(c))
        frac[~degen] = np.clip(
            (np.minimum(zmax_f[~degen], hi)
             - np.maximum(zmin_f[~degen], lo))
            / np.maximum(span[~degen], 1e-12), 0, 1)
        frac[degen] = (mid_z[degen] >= lo) & (mid_z[degen] < hi)
        solid[k] = (v * frac).sum()
    slice_vol = side ** 2 * (z_max / n_slices)
    poro = 1.0 - solid / slice_vol
    return dict(path=path, n=int(m.sum()), half=half, z_max=z_max,
                phi=phi, med=float(np.median(elev)),
                q90=float(np.quantile(elev, 0.9)), poro=poro,
                seed=meta.get("seed"), meta=meta)


def main():
    ap = argparse.ArgumentParser(
        description="Certification tuile pressMT vs monolithe press.")
    ap.add_argument("--monolith", required=True,
                    help="npz du monolithe (route press, pas pressmt)")
    ap.add_argument("--tiles", required=True,
                    help="motif glob des npz de tuiles, ex "
                         "'pressedMT_1000_2_t*.npz'")
    ap.add_argument("--side", type=float, default=0.0,
                    help="côté de la zone du monolithe (défaut : zone "
                         "utile des tuiles, lue dans leurs métadonnées)")
    args = ap.parse_args()

    tile_paths = sorted(glob.glob(args.tiles))
    if not tile_paths:
        raise SystemExit(f"aucune tuile ne matche {args.tiles!r}")
    tiles = [zone_stats(p, None if args.side else _tile_half(p))
             for p in tile_paths]
    side = args.side or 2.0 * tiles[0]["half"]
    mono = zone_stats(args.monolith, side / 2.0)

    # avertissements de protocole
    seeds = {t["seed"] for t in tiles}
    if mono["seed"] not in seeds or len(seeds) > 1:
        print(f"ATTENTION : seeds différents (monolithe {mono['seed']}, "
              f"tuiles {sorted(seeds)}) — la comparaison n'isole rien")
    print("RAPPEL protocole : le monolithe doit avoir été lancé avec "
          "le MÊME box-z que les tuiles (1,5) — le npz ne le stocke pas")
    print()

    print(f"zone de comparaison : {side:.2f} mm de côté")
    hdr = (f"{'structure':<28}{'n':>5}{'φ (%)':>9}{'z_max (mm)':>12}"
           f"{'élev méd':>10}{'q90':>7}")
    print(hdr)
    print("-" * len(hdr))
    print(f"{'monolithe (zone centrale)':<28}{mono['n']:>5}"
          f"{mono['phi']*100:>9.2f}{mono['z_max']:>12.3f}"
          f"{mono['med']:>10.2f}{mono['q90']:>7.2f}")
    for t in tiles:
        name = t["path"].rsplit("/", 1)[-1]
        print(f"{name:<28}{t['n']:>5}{t['phi']*100:>9.2f}"
              f"{t['z_max']:>12.3f}{t['med']:>10.2f}{t['q90']:>7.2f}")
    phis = np.array([t["phi"] for t in tiles])
    print("-" * len(hdr))
    print(f"tuiles : φ = {phis.mean()*100:.2f} ± {phis.std()*100:.2f} % ; "
          f"monolithe : φ = {mono['phi']*100:.2f} % ; "
          f"rapport tuile/monolithe = {phis.mean()/mono['phi']:.2f}")
    print("rapport ~1 = route tuiles certifiée ; sinon : effet murs "
          "(balayage tile_side 5 -> 7,5 pour mesurer sa longueur)")

    # profils côte à côte, abscisse relative z/z_max
    print(f"\nprofil de porosité (profondeur relative, 20 tranches) :")
    cols = [f"{'z/z_max':>8}{'mono':>8}"] + \
           [f"{t['path'].rsplit('/', 1)[-1][10:-4]:>8}" for t in tiles]
    print("".join(cols))
    for k in range(20):
        zr = (k + 0.5) / 20.0
        row = f"{zr:>8.2f}{mono['poro'][k]*100:>8.2f}"
        row += "".join(f"{t['poro'][k]*100:>8.2f}" for t in tiles)
        print(row)
    print("NB : tranches plus fines que ~3 diamètres = surcomptage "
          "solide (avertissement analyse.py) ; comparer les formes, "
          "pas les valeurs absolues des tranches.")


def _tile_half(path):
    """Demi-côté de la zone utile d'un npz pressmt."""
    data = np.load(path, allow_pickle=False)
    meta = {k: data[k].item() for k in data.files if data[k].ndim == 0}
    if int(meta.get("pressmt", 0)) != 1 or "tile_side" not in meta:
        raise SystemExit(f"{path} : pas une tuile pressMT et pas de "
                          f"--side donné")
    return (float(meta["tile_side"])
            - 2.0 * float(meta.get("burn_in", 0.0))) / 2.0


if __name__ == "__main__":
    main()