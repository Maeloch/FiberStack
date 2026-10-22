"""
routed.py — route D : germes booléens. Générateur géométrique
INSTANTANÉ d'un lit de fibres : centres uniformes, directions
contrôlées, AUCUN contrôle de chevauchement, fraction solide visée
par arithmétique. Le « simulateur sans physique » : ce que le modèle
booléen de la littérature filtre (germes poissonniens de grains
aléatoires) fait depuis toujours — pas une presse, une référence.

À QUOI ÇA SERT (et à quoi ça ne sert pas) :
  + tester la chaîne voxel/écoulement à densité C569 en secondes,
    sans attendre une nuit de presse ;
  + contrôle négatif des routes physiques : la différence
    routed vs press mesure CE QUE LA PHYSIQUE APPORTE (arches,
    orientation aplatie par la plaque, cordons de contacts) — si
    les stats pore/porosité sont identiques, la presse ne sert à
    rien, et il faut le savoir ;
  + banc d'étalonnage des outils (validate, analyse, voxelize) sur
    des structures de géométrie connue analytiquement ;
  - MAIS : pas de contacts = pas de gorges de pores réalistes à
    l'échelle des contacts, pas d'arches, élévation artificielle
    (paramètre, pas résultat). Interdit comme substitut du C569 :
    c'est un modèle, pas une simulation.

LE PIÈGE DU DOUBLE-COMPTAGE (ce que le collègue paie sans le savoir)
: avec chevauchements, Σ V_fibres / V_boîte n'est PAS la fraction
solide — les intersections sont comptées deux fois. Pour des germes
uniformes indépendants, P(point couvert) = 1 - (1 - V̄_fibre/V)^n
(processus binomial, exact par fibre i.i.d.) : à ΣV/V = 19 %, la
vraie fraction solide n'est que 17,3 % — 1,7 point de porosité
d'erreur. Deux modes :
  --mode union (DÉFAUT) : n est inversé pour que la vraie fraction
     solide (union des volumes) atteigne φ : n = ln(1-φ)/ln(1-V̄/V) ;
  --mode sum : n = φ·V/V̄, la « compacité » du collègue (ΣV/V),
     avec l'union réelle affichée pour rappel.

Matériau : MÊMES échantillonneurs que press (lognormal d, truncnorm
L) — les structures D sont comparables aux structures A/B, et
validate contrôle les bornes matérielles normalement.

Export : contrat unifié (route=\"routed\", métadonnées core_half,
phi_sum, phi_union). validate ≥ v2.2.1 : pénétration EXEMPTÉE pour
cette route (chevauchements par construction), plancher et
longueurs contrôlés, flag core relu depuis core_half. analyse ≥
v2.2.1 : zone de mesure = core_half.

Usage :
    python routed.py --phi 0.19 --seed 0
    python routed.py --phi 0.19 --seed 1 --side 5 --margin 0.65
    python routed.py --phi 0.10 --mode sum --elev-max 8

Système d'unités : (mm, mg, ms) — hérité, même si la route D ne
simule aucune force. Nécessite : numpy, structure.py, press.py.
"""

import argparse
import math
import os

import numpy as np

from structure import (DOMAIN_SIDE, MEASURE_SIDE, FILTER_THICKNESS,
                       FIBER_LENGTH_MEAN, FIBER_LENGTH_STD, save_pile)
from press import sample_diameters, sample_lengths

CODE_VERSION = 1
N_ESTIMATE = 20000          # tirages pour estimer V̄_fibre
MAX_DIR_TRIES = 500         # redirections d'axe par fibre


def sample_direction(rng, elev_max_deg):
    """Azimut uniforme, élévation uniforme dans ±elev_max degrés.
    Isotropie plane parfaite + dispersion verticale PARAMÉTRÉE (5°
    par défaut ≈ stats des mats pressés : médiane ~2,5°, q90 ~4,5°).
    À ne pas confondre avec un résultat : ici c'est un ingrédient."""
    phi = rng.uniform(0.0, 2.0 * np.pi)
    el = math.radians(rng.uniform(-elev_max_deg, elev_max_deg))
    ce = math.cos(el)
    return np.array([ce * math.cos(phi), ce * math.sin(phi), math.sin(el)])


def main():
    ap = argparse.ArgumentParser(
        description="Route D : lit de fibres par germes booléens "
                    "(chevauchements admis, compacité visée).")
    ap.add_argument("--phi", type=float, default=0.19,
                    help="fraction solide visée (défaut : C569 = 0,19)")
    ap.add_argument("--mode", choices=("union", "sum"), default="union",
                    help="union (défaut) : corrige le double-comptage ; "
                         "sum : φ = ΣV/V_boîte, sans correction")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--side", type=float, default=DOMAIN_SIDE,
                    help="côté du domaine (défaut : 10 mm)")
    ap.add_argument("--thickness", type=float, default=FILTER_THICKNESS,
                    help="épaisseur de la tranche (défaut : 0,3 mm)")
    ap.add_argument("--margin", type=float, default=None,
                    help="rebut par côté pour le flag core (défaut : "
                         "1,5 mm à side 10, sinon L_max/2 = 0,65 mm)")
    ap.add_argument("--elev-max", type=float, default=5.0,
                    help="dispersion d'élévation ±degrés (défaut : 5°)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not (0.0 < args.phi < 1.0):
        raise SystemExit("--phi doit être une fraction solide ]0 ; 1[")

    half = args.side / 2.0
    l_max = FIBER_LENGTH_MEAN + 2.0 * FIBER_LENGTH_STD
    if args.side <= l_max:
        raise SystemExit(f"domaine {args.side} mm <= L_max {l_max} mm : "
                         f"une fibre ne tient pas à plat — élargir --side")
    if args.margin is None:
        margin = max((args.side - MEASURE_SIDE) / 2.0, l_max / 2.0)
    else:
        margin = args.margin
    if margin < l_max / 2.0:
        print(f"AVERTISSEMENT : marge {margin:.2f} < L_max/2 {l_max/2:.2f} "
              f"— la bande de raréfaction près des murs déborde dans le core")

    # volume moyen d'une fibre : estimateur sur flux dédié (le flux
    # principal reste déterministe, indépendant de l'estimation)
    est = np.random.default_rng((args.seed, 0x0D))
    d_est = sample_diameters(est, N_ESTIMATE)
    l_est = sample_lengths(est, N_ESTIMATE)
    v_mean = float(np.mean(np.pi * (d_est / 2.0) ** 2 * l_est))

    v_box = args.side ** 2 * args.thickness
    if args.mode == "union":
        # P(couvert) = 1-(1-V̄/V)^n = phi  ->  n exact (binomial)
        n = int(math.ceil(math.log1p(-args.phi)
                          / math.log1p(-v_mean / v_box)))
    else:
        n = int(math.ceil(args.phi * v_box / v_mean))

    rng = np.random.default_rng(args.seed)
    lengths = sample_lengths(rng, n)
    diameters = sample_diameters(rng, n)

    centers, directions, radii = [], [], []
    n_steer = 0
    for k in range(n):
        L = float(lengths[k])
        r = float(diameters[k]) / 2.0
        u = None
        for _ in range(MAX_DIR_TRIES):
            cand = sample_direction(rng, args.elev_max)
            ez = 0.5 * L * abs(cand[2]) + r
            if ez <= 0.5 * args.thickness:
                u = cand
                break
            n_steer += 1
        if u is None:
            raise SystemExit(f"élévation max {args.elev_max}° impossible "
                             f"à loger dans {args.thickness} mm "
                             f"(L={L:.2f}) — baisser --elev-max")
        cz = rng.uniform(ez, args.thickness - ez)
        ex = 0.5 * L * abs(u[0])
        ey = 0.5 * L * abs(u[1])
        # tir UNIFORME dans le rectangle amputé des demi-envergures :
        # la capsule entière est dans le domaine, pas de clip vers le
        # mur (le clip concentrerait la densité en bande), et la
        # raréfaction près des murs est mangée par la marge core
        cx = rng.uniform(-(half - ex), half - ex)
        cy = rng.uniform(-(half - ey), half - ey)
        centers.append([cx, cy, cz])
        directions.append(u)
        radii.append(r)

    lengths = np.asarray(lengths, dtype=float)
    radii = np.asarray(radii, dtype=float)
    centers = np.asarray(centers, dtype=float)
    directions = np.asarray(directions, dtype=float)

    vols = np.pi * radii ** 2 * lengths
    phi_sum = float(vols.sum() / v_box)
    v_pop = float(vols.mean())
    phi_union = float(1.0 - (1.0 - v_pop / v_box) ** n)
    core_half = half - margin
    core = ((np.abs(centers[:, 0]) <= core_half)
            & (np.abs(centers[:, 1]) <= core_half))
    frac_core_exp = ((2.0 * core_half) / args.side) ** 2

    print(f"route D : {n} fibres, domaine {args.side:.1f}x{args.side:.1f}"
          f"x{args.thickness} mm, seed {args.seed}, mode {args.mode}")
    print(f"  φ ΣV/V (brut, style collègue) : {phi_sum*100:.2f} %")
    print(f"  φ union (vraie fraction solide, binomial) : "
          f"{phi_union*100:.2f} %  (double-comptage = "
          f"{(phi_sum - phi_union)*100:.2f} points)")
    print(f"  directions redessinées : {n_steer} (élévation ne logeant "
          f"pas dans la tranche)")
    print(f"  core : {int(core.sum())} fibres dans {2*core_half:.2f} mm "
          f"de côté ({frac_core_exp*100:.1f} % attendus)")

    out = os.path.join(
        args.outdir,
        f"routed_phi{int(round(args.phi * 100))}_{args.seed}.npz")
    save_pile(out, centers, directions, lengths, radii, core,
              route="routed", code_version=CODE_VERSION,
              seed=args.seed, n_total=n, n_actual=n,
              phi_sum=phi_sum, phi_union=phi_union, mode=args.mode,
              side=args.side, thickness=args.thickness,
              margin=margin, core_half=core_half,
              elev_max=args.elev_max)
    print(f"=== {out} ===")


if __name__ == "__main__":
    main()