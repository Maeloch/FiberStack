"""
pressMT.py — presse MULTI-PROCESSUS par tuiles (v1.1, route B).

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
    cinématique : aucun contact résolu). Les gardes P2, la rampe P1
    et le test de consolidation restent ceux du press v2.4.

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
REQUIERT le press.py déployé (v2.4) importable sans effet de bord
(garde __main__ ; vérifier : grep -n "__main__" press.py).
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


def tile_suspension(client, n, rng, tile_side, box_z, lat, roll, spin, damping):
    """Suspension homogène dans une tuile : grille 2D brassée + jitter
    (mêmes conventions que press.random_suspension), azimut uniforme,
    z anti-chevauchement. Spécificité tuile : un candidat dont la
    capsule déborderait d'une paroi est RAPPROCHÉ vers l'intérieur
    (clip sur la projection de la demi-longueur) — léger biais de
    densité près des murs, mangé par la bande de rebut."""
    half = tile_side / 2.0
    diameters = press.sample_diameters(rng, n)
    lengths = press.sample_lengths(rng, n)
    g = int(math.ceil(math.sqrt(n)))
    cellsize = tile_side / g
    order = np.arange(g * g)
    rng.shuffle(order)
    grid = [(k // g, k % g) for k in order]

    placed, hashcells, n_failed = [], {}, 0
    for k in range(n):
        r = float(diameters[k]) / 2.0
        L = float(lengths[k])
        gi, gj = grid[k % len(grid)]
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

        center = None
        for _ in range(press.MAX_PLACE_TRIES):
            z = rng.uniform(r, box_z - r)
            cand = np.array([cx, cy, z])
            if not press._touches_any(placed, hashcells, cand,
                                      direction, L, r):
                center = cand
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


def export_tile(client, fibers, job):
    """Export au contrat unifié : route=press + géométrie de la tuile.
    Coordonnées LOCALES à la tuile ; core = zone utile après rebut."""
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
              walls=True, escaped=0, pressmt=1,
              tile_ix=job["ix"], tile_iy=job["iy"],
              tile_side=job["tile_side"], burn_in=job["burn_in"])
    n_core = sum(core)
    print(f"=== tuile ({job['ix']},{job['iy']}) : {len(fibers)} pressées, "
          f"{n_core} en zone utile ({2 * half_core:.2f} mm de côté) "
          f"-> {job['out']} ===")


# ------------------------------------------------------------ worker

def worker_tile(job):
    """Un monde pybullet complet pour une tuile. job = dict simple
    (picklable). Semence par tuile : (seed, ix, iy) -> SeedSequence."""
    ix, iy = job["ix"], job["iy"]
    rng = np.random.default_rng((job["seed"], ix, iy))
    client = p.connect(p.DIRECT)
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
                                       press.ROLL_START, press.LIQUID_DAMPING)
    if n_failed > max(1, int(0.01 * job["n_tile"])):
        raise RuntimeError(f"tuile ({ix},{iy}) : {n_failed} échecs de "
                           f"placement (> 1 %) — suspension trop dense ou bug")

    t0 = time.perf_counter()
    for _ in range(press.PRELAX_STEPS):
        p.stepSimulation(physicsClientId=client)
    print(f"tuile ({ix},{iy}) : {len(fibers)} fibres, pré-relaxation "
          f"{press.PRELAX_STEPS} pas ({time.perf_counter() - t0:.0f} s)")

    press.run_press(client, fibers, job["target"], job["frict0"],
                    job["frict1"], verbose=True, zero_g=True,
                    damping=press.LIQUID_DAMPING)
    quiet_consolidate(client, fibers)
    export_tile(client, fibers, job)
    p.disconnect(physicsClientId=client)
    return dict(ix=ix, iy=iy, n=len(fibers), failed=n_failed, out=job["out"])


# ------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Presse multi-processus par tuiles (route B).")
    ap.add_argument("--random", type=int, required=True, metavar="N",
                    help="fibres TOTALES sur le domaine (réparties par tuile)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--target", type=float, default=0.3)
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
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

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
                             frict1=args.frict1, out=out))
            k += 1

    w = args.workers or min(n_tiles, os.cpu_count() or 1)
    print(f"pressMT : {args.random} fibres en {n_tiles} tuiles de "
          f"{args.tile_side} mm ({args.nx}x{args.ny}), suspension "
          f"{args.tile_side}x{args.tile_side}x{args.box_z} mm, {w} workers, "
          f"rampe {args.frict0}->{args.frict1}")
    t0 = time.perf_counter()
    ctx = multiprocessing.get_context("spawn")   # jamais fork avec pybullet
    with ctx.Pool(processes=w) as pool:
        for res in pool.imap_unordered(worker_tile, jobs):
            print(f"[{time.perf_counter() - t0:7.0f} s] tuile finie : {res}")
    print(f"=== pressMT terminé en {time.perf_counter() - t0:.0f} s ===")


if __name__ == "__main__":
    main()