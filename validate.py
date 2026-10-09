"""
validate.py — contrôle qualité d'un .npz FiberStack. Verdict final
VERT (exploitable) / ROUGE (ne pas utiliser, investiguer).

Le protocole manuel a attrapé 5 incidents ; ce script en attrape le
maximum automatiquement, pour qu'un run de nuit se lise en 10 secondes
le matin.

Contrôles :
  - intégrité : pas de NaN/inf, directions unitaires, longueurs/rayons
    positifs et dans les bornes du matériau
  - invariant n_total : jamais plus de fibres que demandé
  - cohérence du flag core (recalculé depuis les centres)
  - hauteur / profondeur : rien sous le plancher au-delà du biais
    documenté de pénétration Bullet (~1 rayon)
  - pénétration max entre fibres (biais documenté ~20-25 µm : le seuil
    est 4 rayons max, au-delà c'est une interpénétration de placement,
    pas le biais du solveur)
  - stats core : élévation médiane/q90, z min/max (informatif)
  - --tiles  : histogramme 4x4 des centres ; échec si bloc vide ou
               déséquilibre max/min >= 2 (symptôme du bug historique
               du tir xy des tuiles)
  - --crossply : ratio comptée/analytique ; fourchette [2.1, 2.6]
               resserrée après mesure : 2.33/2.33/2.35 sur les runs v2
               s=0.30/0.15/0.08. Les fichiers v1 (L=10, ~1.66)
               échouent désormais — attendu, ils sont dépréciés.

Usage :
    python validate.py pile_1000_0.npz
    python validate.py pile_tiles_4000_0.npz --tiles
    python validate.py crossply_150um_12c.npz --crossply
"""

import argparse

import numpy as np

from structure import (DOMAIN_SIDE, MEASURE_SIDE, FILTER_THICKNESS,
                       FIBER_DIAMETER_HIGH,
                       top_exact, bottom_exact, seg_dist, load_pile)


def run_checks(path, check_tiles=False, check_crossply=False):
    data, meta = load_pile(path)
    centers = data["centers"]
    directions = data["directions"]
    lengths = data["lengths"]
    radii = data["radii"]
    core = data["core"]
    n = len(centers)

    results = []   # (nom, ok, détail)

    def check(name, ok, detail=""):
        results.append((name, bool(ok), detail))
        print(("PASS  " if ok else "ÉCHEC") + f"  {name}" +
              (f"  [{detail}]" if detail else ""))

    # ---------- 1) intégrité des tableaux ----------
    arrays_ok = all(np.all(np.isfinite(a)) for a in
                    (centers, directions, lengths, radii))
    check("tableaux finis (pas de NaN/inf)", arrays_ok)

    norms = np.linalg.norm(directions, axis=1)
    check("directions unitaires",
          np.all(np.abs(norms - 1.0) < 1e-6),
          f"max |‖u‖-1| = {np.abs(norms-1.0).max():.2e}")

    r_max = float(radii.max())
    check("rayons dans les bornes du matériau",
          (radii > 0).all()
          and r_max <= 0.5 * FIBER_DIAMETER_HIGH * 1.01,
          f"r_max = {r_max*1e3:.2f} µm "
          f"(borne {0.5*FIBER_DIAMETER_HIGH*1e3:.1f} µm)")
    check("longueurs positives", (lengths > 0).all(),
          f"L min/max = {lengths.min():.3f}/{lengths.max():.3f} mm")

    # ---------- 2) invariant n_total ----------
    if "n_total" in meta:
        n_total = int(meta["n_total"])
        check("invariant n_total (jamais plus que demandé)",
              n <= n_total,
              f"{n} fibres / {n_total} demandées"
              + ("" if n == n_total else
                 f" — {n_total - n} non déposées (perdues/abandonnées)"))
    else:
        print("PASS   (pas de n_total dans les métadonnées — ancien format)")

    # ---------- 3) cohérence du flag core ----------
    half_m = MEASURE_SIDE / 2.0
    core_re = ((np.abs(centers[:, 0]) <= half_m)
               & (np.abs(centers[:, 1]) <= half_m))
    mism = int((core_re != core).sum())
    check("flag core cohérent avec les centres", mism == 0,
          f"{mism} divergences" if mism else f"{int(core.sum())} core")

    # ---------- 4) hauteur / plancher ----------
    z_min = min(bottom_exact(centers[i], directions[i], lengths[i], radii[i])
                for i in range(n))
    check("rien sous le plancher (tolérance biais Bullet)",
          z_min >= -0.05,
          f"z_min = {z_min:.4f} mm")

    # ---------- 5) pénétration entre fibres ----------
    # hash spatial ; seuil = 4 rayons max (biais documenté ~1 rayon)
    cell = 0.5
    cells = {}
    for idx in range(n):
        c = centers[idx]
        reach = 0.5 * lengths[idx] + radii[idx]
        ix0 = int((c[0] - reach) // cell); ix1 = int((c[0] + reach) // cell)
        iy0 = int((c[1] - reach) // cell); iy1 = int((c[1] + reach) // cell)
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                cells.setdefault((ix, iy), []).append(idx)
    max_pen, n_pairs = 0.0, 0
    seen = set()
    for members in cells.values():
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                i, j = members[a], members[b]
                if i >= j or (i, j) in seen:
                    continue
                seen.add((i, j))
                if abs(centers[i][2] - centers[j][2]) > (
                        0.5 * (lengths[i] + lengths[j]) + radii[i] + radii[j]):
                    continue
                p1 = centers[i] - 0.5 * lengths[i] * directions[i]
                q1 = centers[i] + 0.5 * lengths[i] * directions[i]
                p2 = centers[j] - 0.5 * lengths[j] * directions[j]
                q2 = centers[j] + 0.5 * lengths[j] * directions[j]
                d = seg_dist(p1, q1, p2, q2)
                max_pen = max(max_pen, radii[i] + radii[j] - d)
                n_pairs += 1
    check("pénétration <= 4 rayons max",
          max_pen <= 4.0 * r_max,
          f"max {max_pen*1e3:.1f} µm sur {n_pairs} paires "
          f"(biais documenté ~20-25 µm)")

    # ---------- 6) stats core (informatif) ----------
    cm = core if core.dtype == bool else core.astype(bool)
    if cm.sum() > 0:
        cc = centers[cm]; cd = directions[cm]
        elev = np.degrees(np.arcsin(np.clip(np.abs(cd[:, 2]), 0, 1)))
        z_max_core = max(top_exact(cc[i], cd[i], lengths[k], radii[k])
                        for i, k in enumerate(np.where(cm)[0]))
        print(f"info  élévation core : médiane={np.median(elev):.2f}°, "
              f"q90={np.quantile(elev, 0.9):.2f}° ; "
              f"z_max core={z_max_core:.3f} mm")

    # ---------- 7) --tiles : équilibre des tuiles ----------
    if check_tiles:
        half = DOMAIN_SIDE / 2.0
        hist, _, _ = np.histogram2d(centers[:, 0], centers[:, 1],
                                   bins=4, range=[[-half, half],
                                                  [-half, half]])
        print("info  histogramme 4x4 :")
        for row in hist:
            print("      " + " ".join(f"{v:5.0f}" for v in row))
        check("tuiles : aucun bloc vide", (hist > 0).all(),
              f"min={hist.min():.0f}")
        check("tuiles : équilibre max/min < 2", hist.max() < 2 * hist.min(),
              f"max/min = {hist.max()/max(hist.min(),1):.2f}")

    # ---------- 8) --crossply : comptée vs analytique ----------
    if check_crossply:
        if "spacing" not in meta or "n_layers" not in meta:
            check("crossply : métadonnées spacing/n_layers présentes",
                  False, "absentes — npz d'une autre route ?")
        else:
            spacing = float(meta["spacing"])
            n_layers = int(meta["n_layers"])
            pitch = FILTER_THICKNESS / n_layers
            z_max = max(top_exact(centers[i], directions[i], lengths[i],
                                  radii[i]) for i in range(n))
            phi_c = (np.pi * radii[cm] ** 2 * lengths[cm]).sum() \
                    / (MEASURE_SIDE ** 2 * z_max)
            phi_a = np.pi * (radii ** 2).mean() / (spacing * pitch)
            ratio = phi_c / phi_a
            check("crossply : ratio comptée/analytique dans [2.1, 2.6]",
                  2.1 <= ratio <= 2.6,
                  f"comptée {phi_c*100:.1f} % / analytique {phi_a*100:.1f} % "
                  f"= {ratio:.2f} (v2 mesuré 2.33/2.33/2.35 ; v1 L=10 "
                  f"~1.66 -> rouge attendu, fichiers dépréciés)")

    # ---------- verdict ----------
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print()
    if n_fail == 0:
        print("VERDICT : VERT — structure exploitable")
    else:
        print(f"VERDICT : ROUGE — {n_fail} échec(s), "
              f"ne pas utiliser sans investiguer")
    return n_fail


def main():
    ap = argparse.ArgumentParser(description="Contrôle qualité d'un npz FiberStack")
    ap.add_argument("npz", help="fichier .npz à valider")
    ap.add_argument("--tiles", action="store_true",
                    help="contrôler l'équilibre des 4x4 tuiles")
    ap.add_argument("--crossply", action="store_true",
                    help="contrôler le ratio comptée/analytique")
    args = ap.parse_args()
    raise SystemExit(1 if run_checks(args.npz, args.tiles, args.crossply) else 0)


if __name__ == "__main__":
    main()
