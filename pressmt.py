"""
pressMT.py — presse MULTI-PROCESSUS par tuiles (v1.5, route B).

v1.5 : ISOLEMENT DES TUILLES + COUCHAGE (press v2.7). Incident du
3/11 : campagne 16x4750 box 6 tol 0 — à l'incr 200 (plaque ~1,02 mm,
phi ~5,8 %) les 16 tuiles dérivent SIMULTANÉMENT : le réseau de
fibres percole (vmax 0,4-4,2, zéro fibre au-dessus de la plaque).
La garde fait son travail, mais UNE tuile restée à 0,205 mm/ms
après 10 000 pas plaque-fixe lève RuntimeError -> pool.imap_unordered
-> le with-block TUE les 15 tuiles saines encore en pause : ~20 h de
presse, ZÉRO export. Dorénavant le worker capture l'échec et rend
un compte-rendu : les survivantes s'exportent, l'exit code porte
l'échec (tuiles échouées à relancer une à une après diagnostic —
elles sont déterministes). Plomberie --pause-blocks / --max-pauses
(press v2.7) jusqu'au worker ; couchage (gravité dès la première
dérive) et couchage_z (hauteur de percolation) stockés dans le npz
de tuile.

v1.4 : plomberie --overlap-tol (press v2.6) — tolérance de
chevauchement initial au placement, bornée (cap 6 µm = moitié du
plus petit rayon), résolue par la pré-relaxation. Contexte : le
sweep 16x4750 a montré le mur du RSA sans chevauchement (3,4 %
d'échecs à box 4, ~0,4 % projeté à box 6) ; avec 3 µm de
tolérance, box 2-3 redevient viable et la plaque balaie 2-3 fois
moins loin. tol=0 par défaut = régression exacte de la v1.3.
overlap_tol stocké dans le npz de tuile.

v1.3 : plomberie --nips — l'argument traverse jusqu'à
press.run_press (v2.5). nips=1 par défaut = descente continue
(régression exacte de v1.2).

v1.2 : placement par RE-TIRAGE COMPLET (cellule + xy + azimut + z à
chaque essai). v1.1 fixait la cellule de grille d'une fibre et ne
retentait que z : au-delà de ~15 fibres/mm³ la cellule est condamnée
à tous les z et la fibre est perdue (10000 fibres en 4 tuiles =
67 fibres/mm³ -> 358 échecs sur 2500, garde déclenchée). Le
re-tirage consomme une cellule seulement si un essai réussit : les
cellules bouchées sont contournées et la limite devient le vrai
jamming RSA de la suspension. Ajouts : --place-only (sonde de
faisabilité en secondes, sans presse ni export), densité de
brouillard affichée, avertissement au-delà du validé (6,7-12,5
fibres/mm³).

v1.1 : quiet_consolidate passe à 2000 pas de fenêtre de silence
(une chute libre partant du repos franchit vquiet en ~510 pas —
une fenêtre de 500 pouvait se faire escroquer par l'accélération).

Architecture : un worker par tuile, un monde pybullet par worker.
Bullet 2 (impulsions séquentielles) ne parallélise pas un seul monde :
la seule vraie parallélisation est des mondes INDÉPENDANTS = des
processus séparés (multiprocessing, contexte « spawn » pour ne pas
hériter de l'état C++ de pybullet). Chaque worker refait tout seul :
monde, murs de tuile, suspension, pré-relaxation, presse, silence,
export.

Vues d'ensemble :
  - Le domaine n'est PAS réduit : les tuiles restent de vrais volumes
    de filtration, au plancher RVE de 5×L de côté (L moyen 1 mm ->
    5×5 mm minimum, jamais 0,5 mm : plus petit qu'une fibre, les murs
    mangent tout).
  - Chaque tuile exporte ses fibres en COORDONNÉES LOCALES ; le flag
    core marque la ZONE UTILE = tuile moins la bande de rebut
    BURN_IN par côté (recommandé >= L_max/2 + marge = 1,25 mm) :
    la bande proche des murs est artificielle (écoulement latéral
    bloqué, densité biaisée par le clip anti-mur au placement).
  - HONNÊTETÉ : l'ensemble des cœurs n'est PAS un mat continu, c'est
    un échantillonnage de RVE. Porosité, profils, élévation, voxels :
    statistiques légitimes. Percolation, chemins d'écoulement,
    intégrité du mat assemblé : INTERDITS — le montage colle des
    pièces sans fibres traversantes aux coutures.
  - Anti-narcolepsie intégré : sleepThreshold=0 sur chaque fibre
    (un corps endormi ignore la gravité et escroque les critères
    d'immobilité — dossier des fibres fantômes). Si ce pybullet
    rejette le paramètre, le worker prévient ; vérifier d'abord le
    micro-test du sommeil avant de faire foi.
  - press.run_press est réutilisé TEL QUEL : la plaque qu'il
    construit est dimensionnée sur le domaine 10×10 (12×12 mm) et
    déborde sur les murs de la tuile — sans conséquence (statique vs
    cinématique : aucun contact résolu). Les gardes, la rampe P1 et
    le test de consolidation restent ceux du press v2.7 (v1.3 : le
    --nips est passé tel quel au worker ; v1.5 : --pause-blocks,
    --max-pauses et le dict couchage aussi).

Coûts attendus (à 10 000 fibres) : la phase chère est la compaction
dense, dont le coût par pas suit ~le nombre de contacts -> une tuile
de 2 500 coûte ~4× moins par pas ; 4 tuiles en parallèle -> ~3-4× en
temps mur. Avec --box-z 1.5 (le press balaie 7,7 mm de brouillard
pour rien en boîte de 8) : ~2-3× de plus. Total ~6-10× : le run de
dix heures devient 1-2 h.

Usage :
    python pressMT.py --random 10000 --seed 2 --target 0.3
    python pressMT.py --random 10000 --seed 2 --nx 2 --ny 2 \
        --box-z 1.5 --burn-in 1.25 --workers 4
    # pilote phase dense (v1.5 + press v2.7) : la tuile (0,0) partage
    # le flux RNG (seed, 0, 0) avec la tuile (0,0) de la campagne de
    # même seed — placement identique, presse comparable à l'incr 200 :
    python pressMT.py --random 4750 --seed 2 --nx 1 --ny 1 --box-z 6 \
        --pause-blocks 60 --max-pauses 50 --workers 1 --outdir pilote_c569

Validation croisées (AVANT de croire l'ensemble) :
    1. une tuile seule vs le monolith : profil de porosité du cœur
       de tuile vs profil de la zone centrale du run monolith de
       même seed — à quelques % près, l'ensemble est certifié ;
    2. le diagnostic flottants sur chaque sortie de tuile (des
       dorants suspendus passeraient tous les gates actuels) ;
    3. validate : la fraction core attendue d'une tuile n'est PAS
       49 % mais ((tile_side - 2*burn_in)/tile_side)^2 = 25 % à
       5 mm / 1,25 mm — la porte uniformité doit connaître la
       géométrie de la tuile (métadonnées tile_* dans le npz).

Système d'unités : (mm, mg, ms). Nécessite : pybullet, numpy, scipy.
REQUIERT le press.py déployé (v2.7) importable sans effet de bord
(garde __main__ ; vérifier : grep -n "__main__" press.py ; la v1.5
exige run_press(pause_blocks=..., max_pauses=..., info=...), qui
n'existe qu'à partir de la v2.7).
"""

import argparse
import math
import multiprocessing
import os
import time

import numpy as np
import pybullet as p

import press                    # monde, plaque, échantillonnage, gardes (v2.4)
from structure import CELLULOSE_DENSITY, save_pile


# ------------------------------------------------------------ monde tuile

def build_tile_walls(client, half):
    """Cage aux dimensions de la tuile — mêmes conventions que
    press.build_wells (épaisseur 0,05, faces internes à +/- half,
    hauteur 0..16 mm, friction fixe)."""
    bodies = []
    specs = [([half, 0.05, press.WALL_H], [0.0, half + 0.05, press.WALL_H]),
             ([half, 0.05, press.WALL_H], [0.0, -(half + 0.05), press.WALL_H]),
             ([0.05, half, press.WALL_H], [half + 0.05, 0.0, press.WALL_H]),
             ([0.05, half, press.WALL_H], [-(half + 0.05), 0.0, press.WALL_H])]
    for half_ext, pos in specs:
        col = p.createCollisionShape(p.GEOM_BOX, halfExtents=half_ext,
                                     physicsClientId=client)
        wall = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                                 basePosition=pos, physicsClientId=client)
        p.changeDynamics(wall, -1, lateralFriction=press.WALL_FRICTION,
                         rollingFriction=0.0, spinningFriction=0.0,
                         restitution=0.0, physicsClientId=client)
        bodies.append(wall)
    return bodies


def _no_sleep(client, body):
    """Anti-narcolepsie : seuil de sommeil Bullet = 0,4 en unités
    brutes (soit 0,4 mm/ms ici, au-dessus des vitesses réelles du
    dépôt) et temps de désactivation 2,0 s brutes (= 2 ms ici) : tout
    corps laissé tranquille 1000 pas s'endort et IGNORE la gravité.
    sleepThreshold=0 ferme la porte."""
    try:
        p.changeDynamics(body, -1, sleepThreshold=0.0, physicsClientId=client)
    except Exception:
        print("ATTENTION : sleepThreshold non supporté par ce pybullet "
              "— risque de fibres endormies en vol, vérifier avant usage")


def tile_suspension(client, n, rng, tile_side, box_z, lat, roll, spin, damping,
                    tol=0.0):
    """Suspension homogène dans une tuile : grille 2D brassée + jitter,
    azimut uniforme, z anti-chevauchement. Spécificité tuile : clip
    anti-mur du centre sur la projection de la demi-longueur (biais
    de densité près des murs, mangé par la bande de rebut).

    v1.2 — re-tirage complet : à chaque essai, la fibre tire une
    cellule LIBRE au hasard dans la piscine + jitter xy + azimut + z.
    v1.1 fixait la cellule et ne retentait que z : au-delà de ~15
    fibres/mm³, la cellule est condamnée à tous les z et la fibre
    échoue (358/2500 observés à 67 fibres/mm³). Ici la cellule n'est
    consommée que si un essai réussit ; la vraie limite devient le
    jamming RSA de la suspension, à mesurer avec --place-only.
    v1.4 : tol borne le chevauchement initial autorisé (0 = régression
    v1.3 ; cap press.OVERLAP_TOL_MAX = 6 µm) ; résolu par la
    pré-relaxation — usage DEM standard."""
    half = tile_side / 2.0
    diameters = press.sample_diameters(rng, n)
    lengths = press.sample_lengths(rng, n)
    g = int(math.ceil(math.sqrt(n)))
    cellsize = tile_side / g
    order = np.arange(g * g)
    rng.shuffle(order)
    grid = [(k // g, k % g) for k in order]
    pool = list(range(len(grid)))          # cellules libres

    placed, hashcells, n_failed = [], {}, 0
    for k in range(n):
        r = float(diameters[k]) / 2.0
        L = float(lengths[k])
        center, direction = None, None
        for _ in range(press.MAX_PLACE_TRIES):
            if not pool:
                break
            j = int(rng.integers(len(pool)))
            gi, gj = grid[pool[j]]
            cx = -half + (gi + rng.uniform(0.15, 0.85)) * cellsize
            cy = -half + (gj + rng.uniform(0.15, 0.85)) * cellsize
            phi = rng.uniform(0.0, 2.0 * np.pi)
            direction = np.array([np.cos(phi), np.sin(phi), 0.0])
            ex = 0.5 * L * abs(direction[0])
            ey = 0.5 * L * abs(direction[1])
            limx = half - press.CLEAR - ex
            limy = half - press.CLEAR - ey
            cx = min(max(cx, -limx), limx)
            cy = min(max(cy, -limy), limy)
            z = rng.uniform(r, box_z - r)
            cand = np.array([cx, cy, z])
            if not press._touches_any(placed, hashcells, cand,
                                      direction, L, r, tol=tol):
                center = cand
                pool.pop(j)                # cellule consommée
                break
        if center is None:
            n_failed += 1
            continue
        mass = CELLULOSE_DENSITY * math.pi * r ** 2 * L
        body = press.make_fiber_body(client, center, direction, L, r,
                                    mass, lat, roll, spin, damping)
        _no_sleep(client, body)
        placed.append(dict(body=body, center=center, direction=direction,
                           length=L, radius=r))
        press._hash_add(hashcells, len(placed) - 1, center, L, r)
    return placed, n_failed


def quiet_consolidate(client, fibers, min_steps=1000, max_steps=10000,
                      vquiet=0.01, quiet_steps=2000):
    """Consolidation sous gravité jusqu'au silence (discipline v2.4).
    Rétablit la gravité (inoffensif si run_press l'a déjà restaurée)
    puis boucle jusqu'à `quiet_steps` pas consécutifs sous VQUIET.

    quiet_steps = 2000 (v1.1 : était 500) : une fibre en chute libre
    PARTANT DU REPOS atteint vquiet = 0,01 mm/ms après
    0,01/9,81e-3 ≈ 1,02 ms soit ~510 pas — une fenêtre de 500 pas
    consécutifs peut être traversée ENTIÈREMENT pendant
    l'accélération, et le silence est déclaré alors qu'une fibre
    tombe encore (signature observée : tuiles sorties à 1001 pas,
    flottants à ~44 µm d'altitude). 2000 pas ne peut pas être
    traversé par un seul épisode de chute ; chaque rebond remet le
    compteur à zéro."""
    p.setGravity(0, 0, press.GRAVITY, physicsClientId=client)
    quiet = 0
    for step in range(max_steps):
        p.stepSimulation(physicsClientId=client)
        if step % 50 == 0:
            vmax = 0.0
            for f in fibers:
                lin, _ = p.getBaseVelocity(f["body"], physicsClientId=client)
                vmax = max(vmax, float(np.linalg.norm(lin)))
            quiet = quiet + 50 if vmax < vquiet else 0
            if (step + 1) >= min_steps and quiet >= quiet_steps:
                print(f"    repos atteint à {step + 1} pas "
                      f"(vmax={vmax:.4f} mm/ms)")
                return
    print("    AVERTISSEMENT : silence non atteint dans le budget de pas")


def export_tile(client, fibers, job, info=None):
    """Export au contrat unifié : route=press + géométrie de la tuile.
    Coordonnées LOCALES à la tuile ; core = zone utile après rebut.
    v1.5 : couchage et couchage_z (hauteur de percolation — grandeur
    physique du mat) passés au npz via le dict `info` de run_press."""
    press._survey(client, fibers)          # dernier état : pos/dir à jour
    half = job["tile_side"] / 2.0
    half_core = half - job["burn_in"]
    centers, directions, lengths, radii, core = [], [], [], [], []
    for f in fibers:
        pos, d = f["pos"], f["dir"]
        centers.append(pos)
        directions.append(d)
        lengths.append(f["length"])
        radii.append(f["radius"])
        core.append(abs(pos[0]) <= half_core and abs(pos[1]) <= half_core)
    save_pile(job["out"], centers, directions, lengths, radii, core,
              route="press", code_version=press.CODE_VERSION,
              seed=job["seed"], n_total=job["n_tile"], n_actual=len(fibers),
              target=job["target"], frict0=job["frict0"], frict1=job["frict1"],
              overlap_tol=job.get("overlap_tol", 0.0),
              pause_blocks=job.get("pause_blocks", press.PAUSE_MAX_BLOCKS),
              max_pauses=job.get("max_pauses", press.MAX_PAUSES),
              couchage=bool(info.get("couchage", False)) if info else False,
              couchage_z=info.get("couchage_z") if info else None,
              walls=True, escaped=0, pressmt=1,
              tile_ix=job["ix"], tile_iy=job["iy"],
              tile_side=job["tile_side"], burn_in=job["burn_in"])
    n_core = sum(core)
    if info and info.get("couchage"):
        print(f"    tuile ({job['ix']},{job['iy']}) : couchage à plaque "
              f"{info['couchage_z']:.3f} mm, "
              f"{info.get('total_pauses', 0)} pauses au total")
    print(f"=== tuile ({job['ix']},{job['iy']}) : {len(fibers)} pressées, "
          f"{n_core} en zone utile ({2 * half_core:.2f} mm de côté) "
          f"-> {job['out']} ===")


# ------------------------------------------------------------ worker

def worker_tile(job):
    """Un monde pybullet complet pour une tuile. job = dict simple
    (picklable). Semence par tuile : (seed, ix, iy) -> SeedSequence.
    v1.5 : ISOLEMENT — une tuile qui lève (garde de placement,
    budget de pause épuisé, éclat de fin de presse) rend un
    compte-rendu d'échec au lieu de tuer le pool : les autres
    tuiles continuent et s'exportent (l'incident du 3/11 a perdu
    ~20 h de presse sur 15 tuiles saines à cause d'UNE tuile)."""
    ix, iy = job["ix"], job["iy"]
    rng = np.random.default_rng((job["seed"], ix, iy))
    client = p.connect(p.DIRECT)
    try:
        return _worker_run(job, rng, client, ix, iy)
    except Exception as e:               # v1.5 : isolation de tuile
        return dict(ix=ix, iy=iy, status="failed",
                    error=f"{type(e).__name__}: {e}", out=job["out"])
    finally:
        p.disconnect(physicsClientId=client)


def _worker_run(job, rng, client, ix, iy):
    """Corps du worker (v1.5) — séparé de worker_tile pour que
    l'isolation couvre TOUT : placement, pré-relaxation, presse,
    consolidation, export."""
    p.setGravity(0, 0, 0.0, physicsClientId=client)     # apesanteur route B
    p.setPhysicsEngineParameter(fixedTimeStep=press.DT,
                                numSubSteps=press.NUM_SUBSTEPS,
                                numSolverIterations=50,
                                physicsClientId=client)

    half = job["tile_side"] / 2.0
    col = p.createCollisionShape(p.GEOM_BOX,
                                 halfExtents=[2 * half, 2 * half, 1e-2],
                                 physicsClientId=client)
    floor = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                              basePosition=[0, 0, -1e-2],
                              physicsClientId=client)
    p.changeDynamics(floor, -1, lateralFriction=0.1, rollingFriction=0.0,
                     spinningFriction=0.0, restitution=0.0,
                     physicsClientId=client)
    build_tile_walls(client, half)

    fibers, n_failed = tile_suspension(client, job["n_tile"], rng,
                                       job["tile_side"], job["box_z"],
                                       job["frict0"], press.ROLL_START,
                                       press.ROLL_START, press.LIQUID_DAMPING,
                                       tol=job.get("overlap_tol", 0.0))
    if job.get("place_only"):
        dens = job["n_tile"] / (job["tile_side"] ** 2 * job["box_z"])
        print(f"tuile ({ix},{iy}) : PLACEMENT SEUL — {len(fibers)} placées, "
              f"{n_failed} échecs "
              f"({100.0 * n_failed / max(1, job['n_tile']):.1f} %), "
              f"brouillard {dens:.1f} fibres/mm³")
        return dict(ix=ix, iy=iy, status="ok", n=len(fibers),
                    failed=n_failed, place_only=True)
    if n_failed > max(1, int(0.01 * job["n_tile"])):
        raise RuntimeError(f"tuile ({ix},{iy}) : {n_failed} échecs de "
                           f"placement (> 1 %) — brouillard trop dense : "
                           f"sonder --place-only ou augmenter --box-z")

    t0 = time.perf_counter()
    for _ in range(press.PRELAX_STEPS):
        p.stepSimulation(physicsClientId=client)
    print(f"tuile ({ix},{iy}) : {len(fibers)} fibres, pré-relaxation "
          f"{press.PRELAX_STEPS} pas ({time.perf_counter() - t0:.0f} s)")

    info = {}                 # v1.5 : couchage + hauteur de percolation
    press.run_press(client, fibers, job["target"], job["frict0"],
                    job["frict1"], verbose=True, zero_g=True,
                    damping=press.LIQUID_DAMPING, nips=job.get("nips", 1),
                    pause_blocks=job.get("pause_blocks",
                                        press.PAUSE_MAX_BLOCKS),
                    max_pauses=job.get("max_pauses", press.MAX_PAUSES),
                    info=info)
    quiet_consolidate(client, fibers)
    export_tile(client, fibers, job, info)
    return dict(ix=ix, iy=iy, status="ok", n=len(fibers),
                failed=n_failed, out=job["out"])


# ------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Presse multi-processus par tuiles (route B).")
    ap.add_argument("--random", type=int, required=True, metavar="N",
                    help="fibres TOTALES sur le domaine (réparties par tuile)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--target", type=float, default=0.3)
    ap.add_argument("--nips", type=int, default=1,
                    help="presse en N nips (v2.5 : setpoints géométriques "
                         "+ consolidation inter-nips ; 1 = descente "
                         "continue). Si la tuile pilote éclate en fin de "
                         "course, réessayer avec --nips 2 ou 3.")
    ap.add_argument("--overlap-tol", type=float, default=0.0,
                    help="v1.4 : interpénétration initiale autorisée au "
                         "placement, en mm (0 = régression v1.3 ; "
                         "recommandé 0.003 = ~0,3 rayon ; cap "
                         f"{press.OVERLAP_TOL_MAX}). Résolue par la "
                         "pré-relaxation ; calibrer au --place-only.")
    ap.add_argument("--pause-blocks", type=int,
                    default=press.PAUSE_MAX_BLOCKS,
                    help=f"v2.7 : budget de pause plaque-fixe en blocs "
                         f"de {press.PAUSE_STEPS} pas (défaut "
                         f"{press.PAUSE_MAX_BLOCKS} ; recommandé 60 en "
                         "phase dense)")
    ap.add_argument("--max-pauses", type=int, default=press.MAX_PAUSES,
                    help=f"v2.7 : épisodes de dérive consécutifs avant "
                         f"abort (défaut {press.MAX_PAUSES}) ; après "
                         "couchage une avalanche = une pause, c'est "
                         "l'opération normale : 50 en campagne")
    ap.add_argument("--frict0", type=float, default=press.F0_DEFAULT,
                    help="début de rampe (défauts du press, NE PAS écraser)")
    ap.add_argument("--frict1", type=float, default=press.F1_DEFAULT,
                    help="fin de rampe (= attachement du settle)")
    ap.add_argument("--nx", type=int, default=2)
    ap.add_argument("--ny", type=int, default=2)
    ap.add_argument("--tile-side", type=float, default=5.0,
                    help="plancher RVE : >= 5x la longueur de fibre")
    ap.add_argument("--box-z", type=float, default=1.5,
                    help="hauteur de suspension par tuile (cheat : 1,5 mm "
                         "et non 8 — la colonne haute ne sert à rien)")
    ap.add_argument("--burn-in", type=float, default=1.25,
                    help="bande de rebut par côté (>= L_max/2 + marge)")
    ap.add_argument("--workers", type=int, default=0,
                    help="0 = min(tuiles, cpu_count)")
    ap.add_argument("--place-only", action="store_true",
                    help="placer la suspension et sortir (sonde de "
                         "faisabilité : pas de pré-relaxation, presse, "
                         "ni export — quelques dizaines de secondes)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    if not (0.0 <= args.overlap_tol <= press.OVERLAP_TOL_MAX):
        raise SystemExit(f"--overlap-tol dans [0, {press.OVERLAP_TOL_MAX}] "
                         f"mm (recommandé : 0.003)")
    if args.pause_blocks < 1:
        raise SystemExit("--pause-blocks >= 1")
    if args.max_pauses < 1:
        raise SystemExit("--max-pauses >= 1")

    n_tiles = args.nx * args.ny
    n_tile = args.random // n_tiles
    rem = args.random - n_tile * n_tiles
    jobs, k = [], 0
    for ix in range(args.nx):
        for iy in range(args.ny):
            out = os.path.join(
                args.outdir,
                f"pressedMT_{args.random}_{args.seed}_t{ix}{iy}.npz")
            jobs.append(dict(ix=ix, iy=iy, n_tile=n_tile + (1 if k < rem else 0),
                             seed=args.seed, tile_side=args.tile_side,
                             box_z=args.box_z, burn_in=args.burn_in,
                             target=args.target, frict0=args.frict0,
                             frict1=args.frict1, nips=args.nips,
                             overlap_tol=args.overlap_tol,
                             pause_blocks=args.pause_blocks,
                             max_pauses=args.max_pauses, out=out,
                             place_only=args.place_only))
            k += 1

    w = args.workers or min(n_tiles, os.cpu_count() or 1)
    fog = (args.random / n_tiles) / (args.tile_side ** 2 * args.box_z)
    print(f"pressMT : {args.random} fibres en {n_tiles} tuiles de "
          f"{args.tile_side} mm ({args.nx}x{args.ny}), suspension "
          f"{args.tile_side}x{args.tile_side}x{args.box_z} mm "
          f"({fog:.1f} fibres/mm³), tolérance placement "
          f"{args.overlap_tol*1e3:.1f} µm, {w} workers, "
          f"rampe {args.frict0}->{args.frict1}, pauses "
          f"{args.pause_blocks}x{press.PAUSE_STEPS} pas "
          f"(max {args.max_pauses} épisodes)")
    if fog > 15.0:
        print("ATTENTION : brouillard au-delà des densités validées "
              "(6,7-12,5 fibres/mm³) — jamming RSA possible. Sonder "
              "avec --place-only ; sinon augmenter --box-z.")
    t0 = time.perf_counter()
    ok, failed = [], []                        # v1.5 : tuiles isolées
    ctx = multiprocessing.get_context("spawn")   # jamais fork avec pybullet
    with ctx.Pool(processes=w) as pool:
        for res in pool.imap_unordered(worker_tile, jobs):
            t = time.perf_counter() - t0
            if res.get("place_only"):
                print(f"[{t:7.0f} s] sonde : {res}")
                continue
            if res.get("status") == "failed":
                failed.append(res)
                print(f"[{t:7.0f} s] ÉCHEC tuile ({res['ix']},{res['iy']}) : "
                      f"{res['error']} — les autres tuiles continuent")
            else:
                ok.append(res)
                print(f"[{t:7.0f} s] tuile finie : {res}")
    print(f"=== pressMT terminé en {time.perf_counter() - t0:.0f} s : "
          f"{len(ok)} tuile(s) exportée(s), {len(failed)} échec(s) ===")
    if failed:
        raise SystemExit(
            f"pressMT : {len(failed)} tuile(s) en échec — les survivantes "
            "sont exportées ; une tuile échouée est DÉTERMINISTE (même "
            "commande = même échec) : diagnostiquer son log, ajuster "
            "--pause-blocks/--max-pauses ou la physique, puis la relancer "
            "seule")


if __name__ == "__main__":
    main()
