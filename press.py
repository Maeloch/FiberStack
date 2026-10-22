"""
press.py — compression d'une structure de fibres vers FILTER_THICKNESS.

Deux routes :
  A) python press.py pile_tiles_4000_0.npz --target 0.3
     tas issu de settle/tiles, chargé VIVANT (corps dynamiques),
     pré-relaxation sous gravité, puis descente de plaque.
  B) python press.py --random 4000 --seed 1 --target 0.3
     suspension homogène générée en place (l'idéalisation de la
     phase liquide de fabrication : fibres réparties dans un volume,
     puis pressées — pas de dépôt gravitaire du tout).

v2.3 — CAGE À LA TAILLE DU TAS :
     le lanceur des runs seed 1 (v2.1/v2.2) n'était PAS la plaque :
     deux cinématiques radicalement différentes (25 µm/incrément vs
     0,125 µm/pas, splitImpulse) donnaient le même nuage de débris à
     1 mm près. Les fibres (L ~1 mm, centres dans +/-5 mm) débordent
     du domaine par une extrémité pour ~13 % d'entre elles — la cage
     à faces fixes +/-5,0 les ENTAILLAIT au chargement, et le solveur
     les éjectait pendant la pré-relaxation (invisible : aucun
     relevé après pré-relaxation). Désormais la cage se construit
     APRÈS chargement, juste au-delà de l'étendue réelle du tas
     (rayon d'action + 0,1 mm), la plaque passe à 1,5x domaine pour
     couvrir la cage, et un relevé d'état suit la pré-relaxation —
     une éjection y sera visible AVANT le premier incrément.
     NB : c'est probablement la tare dominante des runs v1 aussi
     (1829/4000 hors domaine). NB damping : pybullet applique
     pow(1-d, dt_brut) — LIQUID_DAMPING=0.05 ~ 5 %/ms, PAS inerte.

v2.2 — plaque cinématique par pas (~0,125 µm) + useSplitImpulse :
     retire l'injection de vitesse de Baumgarte propre à la plaque
     téléportée (ERP x pénétration / pas ~ 2,5 mm/ms). Correctif
     nécessaire mais non suffisant — voir v2.3.

v2.1 — route B EN APESANTEUR pendant la presse : une suspension qui
     tombe n'est plus homogène ; gravité restaurée uniquement pour
     le test de consolidation final.

v2 — P1 rampe de frottement sur les fibres ; P2 gardes tous les 20
     incréments ; P3 --seed ; P4 échecs de placement comptés ;
     P5 hauteur exacte à l'export.

Tout verdict inter-routes exige ces runs v2.3, 3 seeds minimum (S6).
Les runs v1 restent interprétables comme « presse sèche raide » avec
entaillage de cage + lanceur de plaque en prime.

Système d'unités : (mm, mg, ms).
Nécessite : pip install pybullet scipy numpy
"""

import argparse
import math
import time

import numpy as np
from scipy.stats import truncnorm
import pybullet as p

from structure import (DOMAIN_SIDE, MEASURE_SIDE, BUFFER_MARGIN,
                       CELLULOSE_DENSITY, FIBER_DIAMETER_LOW,
                       FIBER_DIAMETER_HIGH, FIBER_LENGTH_MEAN,
                       FIBER_LENGTH_STD,
                       top_exact, bottom_exact, quat_align_z,
                       quat_to_direction, seg_dist, save_pile, load_pile)

# ------------------------------------------------------------- paramètres

CODE_VERSION = 2

# presse
INCR = 0.025               # mm par incrément de descente
STEPS_PER_INCR = 200
PRELAX_STEPS = 1500        # route A : sous gravité ; route B : sans gravité
PAUSE_STEPS = 400          # P2 : pause plaque-fixe en cas de dérive
GUARD_EVERY = 20           # P2 : contrôle tous les N incréments
VMAX_LIMIT = 0.25          # mm/ms = 4x la vitesse de descente de la plaque
MAX_PAUSES = 3             # P2 : abort si la dérive persiste
FINAL_SETTLE_MIN = 1000    # minimum de stabilisation finale
FINAL_SETTLE_MAX = 10000    # plafond si le silence cinématique ne vient pas
VQUIET = 0.01              # mm/ms : seuil de repos pour l'export

# plaque et parois
PLATE_HALF_T = 0.05        # demi-épaisseur de la plaque
WALL_H = 8.0               # demi-hauteur des parois (parois 0..16 mm)
WALL_GAP = 0.1             # cage : étendue du tas + ce garde (mm)
WALL_FRICTION = 0.1        # parois fixes, glissantes — la rampe est
                           # pour les FIBRES, pas pour la cage
PLATE_FRICTION = 0.05      # v3 : plaque glissante CONSTANTE — la rampe
                           # ne s'applique plus à la plaque. NB Bullet :
                           # friction effective d'un contact = PRODUIT des
                           # coefficients des deux corps. Donc fibre-plaque
                           # = lat x 0.05, fibre-fibre = lat² (0,0225 -> 1,0
                           # aux défauts), parois/fond = lat x 0,1. Toute
                           # valeur nominale se lit au carré.
                           
# rampe de frottement fibre-fibre (P1)
F0_DEFAULT = 0.15          # latéral en début de presse
F1_DEFAULT = 1.0           # latéral en fin de presse (= settle « à fond »)
ROLL_START, ROLL_END = 0.02, 0.30

# route B : suspension homogène
BOX_Z_DEFAULT = 8.0        # hauteur du volume de suspension (mm)
CLEAR = 0.01               # garde anti-chevauchement au placement (mm)
MAX_PLACE_TRIES = 20       # redéploiements en z avant échec (P4)
LIQUID_DAMPING = 0.05      # ~5 %/ms (pybullet : pow(1-d, dt_brut)) —
                           # viscosité de phase liquide

# solveur (identique à settle.py, + splitImpulse propre à la presse)
DT = 2e-3
NUM_SUBSTEPS = 4
DAMPING = 1.0 - 0.1 ** (1.0 / 1000.0)
GRAVITY = -9.81e-3         # mm/ms²

# Lois de tirage : copie CONFORME de settle.py (mêmes lois, même schéma
# RNG). Si l'une change, changer les DEUX — ou hoister dans structure.py.
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


def sample_lengths(rng, n):
    return truncnorm.rvs(-2.0, 2.0, loc=FIBER_LENGTH_MEAN,
                         scale=FIBER_LENGTH_STD, size=n, random_state=rng)


# ------------------------------------------------------------ monde

def set_fiber_friction(client, body, lat, roll, spin, damping=DAMPING):
    p.changeDynamics(body, -1, lateralFriction=lat, rollingFriction=roll,
                     spinningFriction=spin, frictionAnchor=False,
                     restitution=0.0,
                     linearDamping=damping, angularDamping=damping,
                     physicsClientId=client)


def make_fiber_body(client, center, direction, length, radius, mass,
                    lat, roll, spin, damping=DAMPING):
    quat = quat_align_z(direction)
    col = p.createCollisionShape(p.GEOM_CAPSULE, radius=radius,
                                 height=max(length - 2 * radius, radius * 0.1),
                                 physicsClientId=client)
    body = p.createMultiBody(baseMass=mass, baseCollisionShapeIndex=col,
                             basePosition=list(center), baseOrientation=quat,
                             physicsClientId=client)
    set_fiber_friction(client, body, lat, roll, spin, damping)
    return body


def pile_extent(fibers):
    """Étendue réelle du tas, extrémités de capsules comprises. Les
    fibres (L ~1 mm) débordent quasi toutes du domaine des centres :
    c'est cette étendue-là que la cage doit encadrer, pas ±5,0."""
    ext_x = ext_y = 0.0
    for f in fibers:
        c, u, L, r = f["center"], f["direction"], f["length"], f["radius"]
        ext_x = max(ext_x, abs(c[0]) + 0.5 * L * abs(u[0]) + r)
        ext_y = max(ext_y, abs(c[1]) + 0.5 * L * abs(u[1]) + r)
    return ext_x, ext_y


def build_plate(client, z_bottom):
    """Plaque statique (masse 0), déplacée à CHAQUE PAS de
    INCR/STEPS_PER_INCR (~0,125 µm) — plus jamais de 25 µm d'un coup
    (biais de vitesse de Baumgarte). 1,5x domaine : couvre la cage
    élargie (v2.3)."""
    half = DOMAIN_SIDE * 0.75
    col = p.createCollisionShape(p.GEOM_BOX,
                                 halfExtents=[half, half, PLATE_HALF_T],
                                 physicsClientId=client)
    plate = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                               basePosition=[0.0, 0.0, z_bottom + PLATE_HALF_T],
                               physicsClientId=client)
    return plate


def build_wells(client, ext_x, ext_y):
    """Cage de 4 parois, faces internes à l'étendue du tas + WALL_GAP
    (v2.3 : construite APRÈS chargement — une cage à +/-5,0 fixes
    entaillait les ~13 % de fibres qui débordent du domaine, et le
    solveur les éjectait ; c'était le lanceur des runs seed 1).
    Friction FIXE 0,1 : la cage guide, elle ne tord pas."""
    face_x = ext_x + WALL_GAP
    face_y = ext_y + WALL_GAP
    wall_len = max(face_x, face_y) + 0.3
    bodies = []
    specs = [([wall_len, 0.05, WALL_H], [0.0,  face_y + 0.05, WALL_H]),
             ([wall_len, 0.05, WALL_H], [0.0, -(face_y + 0.05), WALL_H]),
             ([0.05, wall_len, WALL_H], [ face_x + 0.05, 0.0, WALL_H]),
             ([0.05, wall_len, WALL_H], [-(face_x + 0.05), 0.0, WALL_H])]
    for half_ext, pos in specs:
        col = p.createCollisionShape(p.GEOM_BOX, halfExtents=half_ext,
                                      physicsClientId=client)
        wall = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                                 basePosition=pos, physicsClientId=client)
        p.changeDynamics(wall, -1, lateralFriction=WALL_FRICTION,
                         rollingFriction=0.0, spinningFriction=0.0,
                         restitution=0.0, physicsClientId=client)
        bodies.append(wall)
    print(f"cage : faces internes x=±{face_x:.2f}, y=±{face_y:.2f} mm "
          f"(étendue du tas + {WALL_GAP} mm)")
    return bodies


# ------------------------------------------------------------ route A : chargement vivant

def load_npz_live(client, npz_path, lat, roll, spin):
    """Reprend un npz settle/tiles en corps DYNAMIQUES (masse réelle).
    Les fichiers route="press" sont refusés : re-presser un tas déjà
    pressé n'a pas de sens physique ici."""
    data, meta = load_pile(npz_path, routes=("settle", "tiles"))
    fibers = []
    for c, u, L, r in zip(data["centers"], data["directions"],
                          data["lengths"], data["radii"]):
        c, u, L, r = np.array(c), np.array(u), float(L), float(r)
        mass = CELLULOSE_DENSITY * math.pi * r ** 2 * L
        body = make_fiber_body(client, c, u, L, r, mass, lat, roll, spin)
        fibers.append(dict(body=body, center=c, direction=u,
                           length=L, radius=r))
    return fibers, meta


# ------------------------------------------------------------ route B : suspension homogène

def _seg_ends(center, direction, length):
    return (np.asarray(center) - 0.5 * length * np.asarray(direction),
            np.asarray(center) + 0.5 * length * np.asarray(direction))


def _hash_add(hashcells, idx, center, length, radius, cell=1.0):
    reach = 0.5 * length + radius
    ix0 = int(math.floor((center[0] - reach) / cell))
    ix1 = int(math.floor((center[0] + reach) / cell))
    iy0 = int(math.floor((center[1] - reach) / cell))
    iy1 = int(math.floor((center[1] + reach) / cell))
    iz0 = int(math.floor((center[2] - reach) / cell))
    iz1 = int(math.floor((center[2] + reach) / cell))
    for ix in range(ix0, ix1 + 1):
        for iy in range(iy0, iy1 + 1):
            for iz in range(iz0, iz1 + 1):
                hashcells.setdefault((ix, iy, iz), []).append(idx)


def _touches_any(placed, hashcells, center, direction, length, radius,
                 cell=1.0):
    """z-check (P4) : le candidat ne doit interpénétrer AUCUNE fibre
    déjà placée — une interpénétration initiale = explosion Bullet.
    Portée = demi-longueur + rayon du candidat + ~0,7 mm adverses."""
    reach = 0.5 * length + radius + 0.7
    p1, q1 = _seg_ends(center, direction, length)
    ix0 = int(math.floor((center[0] - reach) / cell))
    ix1 = int(math.floor((center[0] + reach) / cell))
    iy0 = int(math.floor((center[1] - reach) / cell))
    iy1 = int(math.floor((center[1] + reach) / cell))
    iz0 = int(math.floor((center[2] - reach) / cell))
    iz1 = int(math.floor((center[2] + reach) / cell))
    for ix in range(ix0, ix1 + 1):
        for iy in range(iy0, iy1 + 1):
            for iz in range(iz0, iz1 + 1):
                for j in hashcells.get((ix, iy, iz), ()):
                    f = placed[j]
                    p2, q2 = _seg_ends(f["center"], f["direction"], f["length"])
                    if seg_dist(p1, q1, p2, q2) < radius + f["radius"] + CLEAR:
                        return True
    return False


def random_suspension(client, n, seed, box_z, lat, roll, spin, damping,
                      verbose=True):
    """Route B : n fibres réparties dans un volume 10x10xbox_z.
    Grille 2D brassée (une fibre par cellule, jitter dans la cellule)
    pour la répartition xy ; z tiré au hasard, avec z-check
    anti-chevauchement. Retourne (fibers, n_failed) — l'appelant
    décide de l'abort (P4 : > 1 % d'échecs = bug ou densité absurde).
    NB : les fibres débordent du domaine par leurs extrémités — la
    cage est construite après coup à leur étendue (v2.3)."""
    rng = np.random.default_rng(seed)
    diameters = sample_diameters(rng, n)
    lengths = sample_lengths(rng, n)
    half = DOMAIN_SIDE / 2.0
    g = int(math.ceil(math.sqrt(n)))
    cellsize = DOMAIN_SIDE / g
    order = np.arange(g * g)
    rng.shuffle(order)
    grid = [(k // g, k % g) for k in order]

    placed = []
    hashcells = {}
    n_failed = 0
    for k in range(n):
        r = float(diameters[k]) / 2.0
        L = float(lengths[k])
        gi, gj = grid[k % len(grid)]
        cx = -half + (gi + rng.uniform(0.15, 0.85)) * cellsize
        cy = -half + (gj + rng.uniform(0.15, 0.85)) * cellsize
        phi = rng.uniform(0.0, 2.0 * np.pi)
        # horizontales, azimut uniforme (l'hypothèse de fabrication :
        # pas de rotation à gérer, la presse s'en chargera)
        direction = np.array([np.cos(phi), np.sin(phi), 0.0])
        center = None
        for _ in range(MAX_PLACE_TRIES):
            z = rng.uniform(r, box_z - r)
            cand = np.array([cx, cy, z])
            if not _touches_any(placed, hashcells, cand, direction, L, r):
                center = cand
                break
        if center is None:
            n_failed += 1
            if verbose:
                print(f"  <-- placement échoué (fibre {k}) : "
                      f"{MAX_PLACE_TRIES} tirages z en chevauchement")
            continue
        mass = CELLULOSE_DENSITY * math.pi * r ** 2 * L
        body = make_fiber_body(client, center, direction, L, r, mass,
                               lat, roll, spin, damping)
        f = dict(body=body, center=center, direction=direction,
                 length=L, radius=r)
        placed.append(f)
        _hash_add(hashcells, len(placed) - 1, center, L, r)
    return placed, n_failed


# ------------------------------------------------------------ presse

def _survey(client, fibers):
    """Relit l'état courant des corps (position/orientation/vitesse)
    dans chaque dict fibre — les orientations BOUGENT pendant la
    presse, les direction initiales sont périmées."""
    vmax = 0.0
    z_top = 0.0
    for f in fibers:
        pos, quat = p.getBasePositionAndOrientation(f["body"],
                                                     physicsClientId=client)
        lin, _ = p.getBaseVelocity(f["body"], physicsClientId=client)
        d = quat_to_direction(quat)
        f["pos"], f["dir"] = pos, d
        vmax = max(vmax, float(np.linalg.norm(lin)))
        z_top = max(z_top, top_exact(pos, d, f["length"], f["radius"]))
    return vmax, z_top


def run_press(client, fibers, target, frict0, frict1, verbose=True,
              zero_g=False, damping=DAMPING):
    """Descente de plaque : INCR total par incrément, mais déplacée à
    CHAQUE PAS (v2.2) de INCR/STEPS_PER_INCR. Rampe P1 sur toutes les
    fibres à chaque incrément, gardes P2 tous les GUARD_EVERY
    incréments.

    zero_g=True (route B) : presse en apesanteur (phase liquide) ;
    gravité restaurée juste avant le settle final — le test de
    consolidation : le mat doit tenir la cible sous son propre poids.

    v2.4 : export REPOSÉ — on attend le silence cinématique
    (vmax < VQUIET) au lieu d'un compte de pas fixe. Le run 1000
    seed 1 a été exporté 2 ms après l'arrêt de la plaque : les fibres
    suivaient encore à ~0,07 mm/ms et ont continué sur leur élan
    (sommet 0,302 -> 0,224 mm). Sans ça, les hauteurs finales ne
    sont comparables ni entre routes ni entre runs."""
    z_start = max(top_exact(f["center"], f["direction"], f["length"], f["radius"])
                  for f in fibers) + 0.02
    n_incr = int(math.ceil((z_start - target) / INCR))
    total_steps = n_incr * STEPS_PER_INCR
    dz = (z_start - target) / total_steps     # ~0,125 µm/pas, cible exacte
    plate = build_plate(client, z_start)
    print(f"presse : {n_incr} incréments de {INCR*1e3:.0f} µm x "
          f"{STEPS_PER_INCR} pas (plaque {dz*1e3:.3f} µm/pas), "
          f"rampe {frict0}->{frict1} "
          f"(rolling {ROLL_START}->{ROLL_END})"
          + ("  [APESANTEUR — phase liquide]" if zero_g else ""))

    prev_above = 0
    n_pauses = 0
    t0 = time.perf_counter()
    t_step = 0
    for k in range(1, n_incr + 1):
        frac = k / n_incr
        
        # ---- P1 : la rampe s'applique aux FIBRES uniquement (v3)
        lat = frict0 + (frict1 - frict0) * frac
        roll = ROLL_START + (ROLL_END - ROLL_START) * frac
        for f in fibers:
            set_fiber_friction(client, f["body"], lat, roll, roll, damping)
        # v3 : plaque en friction constante basse — elle presse,
        # elle n'agrippe pas. (v2 : elle recevait la rampe et finissait
        # à frict1, collant les fibres en fin de course.)
        p.changeDynamics(plate, -1, lateralFriction=PLATE_FRICTION,
                         rollingFriction=0.0, spinningFriction=0.0,
                         restitution=0.0, physicsClientId=client)
                         
        for _ in range(STEPS_PER_INCR):
            z_bot = z_start - dz * (t_step + 1)
            p.resetBasePositionAndOrientation(plate,
                                              [0.0, 0.0, z_bot + PLATE_HALF_T],
                                              [0.0, 0.0, 0.0, 1.0],
                                              physicsClientId=client)
            p.stepSimulation(physicsClientId=client)
            t_step += 1

        # ---- P2 : gardes
        if k % GUARD_EVERY == 0 or k == n_incr:
            vmax, z_top = _survey(client, fibers)
            n_above = sum(1 for f in fibers
                          if bottom_exact(f["pos"], f["dir"], f["length"],
                                          f["radius"]) > z_bot + 0.02)
            if verbose:
                print(f"  incr {k}/{n_incr}  plaque z={z_bot:.3f} mm  "
                      f"sommet tas z={z_top:.3f} mm  vmax={vmax:.3f} mm/ms  "
                      f"au-dessus={n_above}  ({time.perf_counter()-t0:.0f} s)")
            # dérive = explosion (vmax) ou escalade progressive de la plaque
            drift = vmax > VMAX_LIMIT or n_above > prev_above
            if drift:
                n_pauses += 1
                print(f"  <-- DÉRIVE (vmax={vmax:.3f}, au-dessus "
                      f"{prev_above}->{n_above}) : pause {PAUSE_STEPS} pas, "
                      f"plaque fixe ({n_pauses}/{MAX_PAUSES})")
                if n_pauses >= MAX_PAUSES:
                    raise RuntimeError(
                        "presse divergente : la dérive persiste après "
                        f"{MAX_PAUSES} pauses — matériau ou paramètres à revoir")
                for _ in range(PAUSE_STEPS):
                    p.stepSimulation(physicsClientId=client)
            else:
                n_pauses = 0
            prev_above = n_above

    # ---- test de consolidation (route B) : gravité restaurée
    if zero_g:
        p.setGravity(0, 0, GRAVITY, physicsClientId=client)
        print(f"gravité restaurée : test de consolidation — "
              f"le mat doit tenir la cible ; export à l'état reposé")
    # ---- v2.4 : stabilisation jusqu'au repos (min/max encadrés)
    waited = 0
    while waited < FINAL_SETTLE_MIN:
        p.stepSimulation(physicsClientId=client)
        waited += 1
    while waited < FINAL_SETTLE_MAX:
        if waited % 500 == 0:
            v, _ = _survey(client, fibers)
            if v < VQUIET:
                print(f"  repos atteint à {waited} pas "
                      f"(vmax={v:.4f} mm/ms)")
                break
        p.stepSimulation(physicsClientId=client)
        waited += 1
    else:
        print(f"  AVERTISSEMENT : pas de silence cinématique après "
              f"{FINAL_SETTLE_MAX} pas — export d'un état agité")
    return plate


# ------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Compression d'une structure de fibres vers FILTER_THICKNESS")
    ap.add_argument("input", nargs="?", default=None,
                    help="route A : npz d'un tas settle/tiles")
    ap.add_argument("--target", type=float, required=True,
                    help="épaisseur finale visée (mm), ex. 0.3")
    ap.add_argument("--random", type=int, default=None, metavar="N",
                    help="route B : générer N fibres en suspension homogène")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed de la suspension (route B)")
    ap.add_argument("--box-z", type=float, default=BOX_Z_DEFAULT,
                    help="hauteur du volume de suspension (route B, mm)")
    ap.add_argument("--frict0", type=float, default=F0_DEFAULT,
                    help="frottement fibre-fibre en début de presse")
    ap.add_argument("--frict1", type=float, default=F1_DEFAULT,
                    help="frottement fibre-fibre en fin de presse")
    ap.add_argument("--no-walls", action="store_true",
                    help="désactiver la cage (debug uniquement)")
    ap.add_argument("--out", default=None, help="fichier de sortie")
    args = ap.parse_args()

    if args.random is None and args.input is None:
        ap.error("donner un npz (route A) ou --random N (route B)")
    if args.random is not None and args.input is not None:
        ap.error("route A et route B sont exclusives")

    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, GRAVITY, physicsClientId=client)  # mm/ms²
    # splitImpulse (v2.2) : la correction de position passe par une
    # pseudo-vitesse séparée — la pénétration résiduelle n'injecte plus
    # de vitesse réelle dans les fibres. NB pybullet : le paramètre
    # s'appelle useSplitImpulse (splitImpulse est rejeté par l'API).
    # Propre à la presse : settle.py ne touche pas à ce paramètre
    # (régression bit-exacte à préserver).
    p.setPhysicsEngineParameter(fixedTimeStep=DT, numSubSteps=NUM_SUBSTEPS,
                                numSolverIterations=50,
                                useSplitImpulse=1, physicsClientId=client)

    half = DOMAIN_SIDE / 2.0
    half_floor = 2.0 * half            # plancher 2x le domaine
    floor_col = p.createCollisionShape(p.GEOM_BOX,
                                       halfExtents=[half_floor, half_floor, 1e-2],
                                       physicsClientId=client)
    floor = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=floor_col,
                              basePosition=[0, 0, -1e-2], physicsClientId=client)
    p.changeDynamics(floor, -1, lateralFriction=0.1, rollingFriction=0.0,
                     spinningFriction=0.0, restitution=0.0,
                     physicsClientId=client)

    zero_g = False
    source = None
    if args.random is not None:
        # ---------------- route B : suspension homogène ----------------
        print(f"route B : suspension homogène, {args.random} fibres, "
              f"seed {args.seed}, volume {DOMAIN_SIDE}x{DOMAIN_SIDE}"
              f"x{args.box_z} mm — presse en apesanteur (phase liquide)")
        fibers, n_failed = random_suspension(client, args.random, args.seed,
                                             args.box_z, args.frict0,
                                             ROLL_START, ROLL_START,
                                             LIQUID_DAMPING)
        n_total = args.random
        n_actual = len(fibers)
        if n_failed > max(1, int(0.01 * n_total)):
            raise RuntimeError(
                f"P4 : {n_failed} échecs de placement (> 1 % de {n_total}) "
                f"— suspension trop dense ou bug, on n'exporte pas")
        if n_failed > 0:
            print(f"P4 : {n_failed} échecs de placement (<= 1 %, tolérés) "
                  f"— n_actual={n_actual}")
        # apesanteur de bout en bout : placement, pré-relaxation, presse.
        # Restaurer la gravité ici = remettre la pluie balistique du run
        # seed 1 (vmax 29 mm/ms, débris à 347 mm).
        p.setGravity(0, 0, 0.0, physicsClientId=client)
        print(f"pré-relaxation (sans gravité) : {PRELAX_STEPS} pas")
        t0 = time.perf_counter()
        for _ in range(PRELAX_STEPS):
            p.stepSimulation(physicsClientId=client)
        print(f"  pré-relaxation : {PRELAX_STEPS} pas "
              f"({time.perf_counter()-t0:.0f} s)")
        zero_g = True
        seed_out = args.seed
        out = args.out or f"pressed_random_{n_total}_{args.seed}.npz"
    else:
        # ---------------- route A : tas chargé vivant ----------------
        print(f"route A : chargement {args.input}")
        fibers, meta = load_npz_live(client, args.input, args.frict0,
                                     ROLL_START, ROLL_START)
        n_total = n_actual = len(fibers)
        seed_out = int(meta["seed"]) if "seed" in meta else None
        source = args.input
        print(f"pré-relaxation (sous gravité) : {PRELAX_STEPS} pas")
        t0 = time.perf_counter()
        for _ in range(PRELAX_STEPS):
            p.stepSimulation(physicsClientId=client)
        print(f"  pré-relaxation : {PRELAX_STEPS} pas "
              f"({time.perf_counter()-t0:.0f} s)")
        out = args.out or f"pressed_{args.input.rsplit('.', 1)[0]}.npz"

    # ---- cage À LA TAILLE DU TAS (v2.3), APRÈS chargement : rien ne
    # dépasse de l'étendue mesurée, aucune fibre ne peut être entaillée
    walls_built = not args.no_walls
    if walls_built:
        ext_x, ext_y = pile_extent(fibers)
        build_wells(client, ext_x, ext_y)

    # ---- relevé post-pré-relaxation : une éjection (cage, placement)
    # se voit ICI, avant le premier incrément — pas à l'incrément 20
    vmax0, z_top0 = _survey(client, fibers)
    print(f"état après pré-relaxation : sommet {z_top0:.3f} mm, "
          f"vmax {vmax0:.3f} mm/ms")

    # ----------------------------- presse -----------------------------
    run_press(client, fibers, args.target, args.frict0, args.frict1,
              zero_g=zero_g, damping=LIQUID_DAMPING if zero_g else DAMPING)

    # ----------------------------- export -----------------------------
    _survey(client, fibers)          # dernier état : pos/dir à jour
    half_m = MEASURE_SIDE / 2.0
    centers, directions, lengths, radii, core = [], [], [], [], []
    escaped = 0
    for f in fibers:
        pos, d = f["pos"], f["dir"]
        centers.append(pos); directions.append(d)
        lengths.append(f["length"]); radii.append(f["radius"])
        ck = (abs(pos[0]) <= half_m) and (abs(pos[1]) <= half_m)
        core.append(ck)
        if (abs(pos[0]) > half + BUFFER_MARGIN
                or abs(pos[1]) > half + BUFFER_MARGIN):
            escaped += 1
    save_pile(out, centers, directions, lengths, radii, core,
              route="press", code_version=CODE_VERSION,
              seed=seed_out, n_total=n_total, n_actual=n_actual,
              target=args.target, frict0=args.frict0, frict1=args.frict1,
              walls=walls_built, escaped=escaped, source=source)

    z_top_core = max((top_exact(c, d, L, r)
                      for c, d, L, r, ck in zip(centers, directions,
                                                lengths, radii, core)
                      if ck), default=0.0)
    print(f"=== {len(fibers)} fibres pressées, {escaped} hors domaine, "
          f"hauteur finale {z_top_core:.3f} mm -> {out} ===")
    if z_top_core > args.target * 1.05:
        print(f"AVERTISSEMENT : sommet {z_top_core:.3f} mm > cible "
              f"{args.target} mm — tas NON consolidé (presse arrêtée avant "
              f"consolidation, ou matériau trop raide pour la rampe)")
    p.disconnect(physicsClientId=client)


if __name__ == "__main__":
    main()
