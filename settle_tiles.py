"""
settle_tiles.py — campagne par tuiles avec recouvrement, pour dépasser
le coût O(N^2) du broadphase (profilé : ~3,7-11 µs par corps et par pas
dans stepSimulation).

Domaine 10x10 mm découpé en TILE_GRID x TILE_GRID tuules (2x2 -> 5 mm
de côté). Deux vagues en damier : vague 1 sans contexte, vague 2 avec
les fibres de la vague 1 situées dans la marge reconstruites en corps
gelés. Chaque tuile est un monde PyBullet indépendant (processus
séparé) : le coût par pas ne voit que les corps de la tuule + marge.

APPROXIMATION DE COUTURE (documentée, à valider) : près des frontières
internes, les fibres de la vague 1 reposent un peu plus bas (voisin
vide au moment de leur dépôt), celles de la vague 2 un peu plus haut.
Protocole de validation : python settle_tiles.py 1000 0, comparer le
tas fusionné à pile_1000_0.npz / pile_1000_1.npz (élévation, porosité
par tranches, occupation, audit de support).

Usage :
    python settle_tiles.py 30000 0     # campagne
    python settle_tiles.py 1000 0      # validation vs séquentiel
"""

import sys
import time
import numpy as np
import pybullet as p
from multiprocessing import Pool

from settle import (DOMAIN_SIDE, BUFFER_MARGIN, MEASURE_SIDE,
                    FRICTION, ROLLING_FRICTION, SPINNING_FRICTION,
                    DAMPING, DT, NUM_SUBSTEPS, FALL_MAX_STEPS,
                    STICK_MAX_STEPS,
                    sample_diameters, sample_lengths,
                    _quat_align_z, _wait_until_settled,
                    place_fiber, drop_fiber)

TILE_CODE_VERSION = 6     # génération "tuilée" : npz non échangeables
                           # avec la ligne séquentielle (CODE_VERSION=5)
TILE_GRID = 2             # 2x2 tuules de 5 mm (fibres de 1 mm : la
                          # couture doit rester << la tuile)
WAVES = [[(0, 0), (1, 1)], [(0, 1), (1, 0)]]   # damier, 2 vagues


def tile_region(tx, ty):
    """Limites (x0, x1, y0, y1) du cœur de la tuile (pas de la marge)."""
    s = DOMAIN_SIDE / TILE_GRID
    x0 = -DOMAIN_SIDE / 2 + tx * s
    y0 = -DOMAIN_SIDE / 2 + ty * s
    return x0, x0 + s, y0, y0 + s


def load_margin_context(client, npz_paths, region, margin=BUFFER_MARGIN):
    """Reconstruit en corps gelés (masse 0) les fibres des exports voisins
    dont l'enveloppe intersecte [région ± margin]. Ne renvoie rien : les
    corps ne servent que de support figé, ils ne sont pas re-déposés."""
    x0, x1, y0, y1 = region
    n_ctx = 0
    for path in npz_paths:
        data = np.load(path)
        for c, u, L, r in zip(data["centers"], data["directions"],
                              data["lengths"], data["radii"]):
            reach = float(L) / 2 + float(r)
            if (c[0] + reach < x0 - margin or c[0] - reach > x1 + margin or
                    c[1] + reach < y0 - margin or c[1] - reach > y1 + margin):
                continue                      # hors de portée d'interaction
            quat = _quat_align_z(u)
            col = p.createCollisionShape(p.GEOM_CAPSULE, radius=float(r),
                                         height=max(float(L) - 2 * float(r),
                                                    float(r) * 0.1),
                                         physicsClientId=client)
            body = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                                     basePosition=[float(c[0]), float(c[1]),
                                                   float(c[2])],
                                     baseOrientation=quat, physicsClientId=client)
            p.changeDynamics(body, -1, lateralFriction=FRICTION,
                             rollingFriction=ROLLING_FRICTION,
                             spinningFriction=SPINNING_FRICTION,
                             restitution=0.0, physicsClientId=client)
            n_ctx += 1
    return n_ctx


def run_tile(args):
    """Un processus = une tuule = un monde PyBullet complet.
    Dépose n_tile fibres dont le xy est dans le cœur de la tuule,
    avec la même logique de relance que build_pile."""
    tx, ty, n_tile, seed, context_paths = args
    t_start = time.perf_counter()

    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, -9.81e-3, physicsClientId=client)
    p.setPhysicsEngineParameter(fixedTimeStep=DT, numSubSteps=NUM_SUBSTEPS,
                                numSolverIterations=50, physicsClientId=client)
    half = DOMAIN_SIDE / 2.0
    floor_col = p.createCollisionShape(p.GEOM_BOX,
                                       halfExtents=[2 * half, 2 * half, 1e-2],
                                       physicsClientId=client)
    floor = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=floor_col,
                              basePosition=[0, 0, -1e-2], physicsClientId=client)
    p.changeDynamics(floor, -1, lateralFriction=FRICTION,
                     rollingFriction=ROLLING_FRICTION,
                     spinningFriction=SPINNING_FRICTION, restitution=0.0,
                     physicsClientId=client)

    x0, x1, y0, y1 = tile_region(tx, ty)
    n_ctx = 0
    if context_paths:
        n_ctx = load_margin_context(client, context_paths,
                                    (x0, x1, y0, y1))

    rng = np.random.default_rng(seed)
    diameters = sample_diameters(rng, n_tile)
    lengths = sample_lengths(rng, n_tile)

    fibers = []
    max_h = 0.0
    n_ejected = n_retried = n_abandoned = 0
    for i in range(n_tile):
        radius = diameters[i] / 2.0
        length = lengths[i]
        xy = rng.uniform(x0, x1, size=2)
        phi = rng.uniform(0, 2 * np.pi)
        direction = np.array([np.cos(phi), np.sin(phi), 0.0])

        attempt = 0
        while True:
            z_start = place_fiber(client, xy, direction, radius, length, max_h)
            result = drop_fiber(client, xy, z_start, direction,
                                radius=radius, length=length,
                                fall_steps=FALL_MAX_STEPS)
            attempt += 1
            pathological = (result is None or result[0] is None)
            if not pathological:
                body, pos, dirf, found, used = result
                if not found or used >= STICK_MAX_STEPS:
                    pathological = True
            if pathological:
                if result is not None and result[0] is not None:
                    p.removeBody(result[0], physicsClientId=client)
                if attempt == 1:
                    n_retried += 1
                    xy = rng.uniform(x0, x1, size=2)
                    continue
                if result is None or result[0] is None:
                    n_ejected += 1
                else:
                    n_abandoned += 1
                break
            top = pos[2] + length / 2 + radius
            max_h = max(max_h, top)
            fibers.append(dict(center=pos, direction=dirf, length=length,
                               radius=radius))
            break

    out = f"tile_{tx}{ty}.npz"
    np.savez(out,
             centers=np.array([f["center"] for f in fibers]),
             directions=np.array([f["direction"] for f in fibers]),
             lengths=np.array([f["length"] for f in fibers]),
             radii=np.array([f["radius"] for f in fibers]),
             tile=(tx, ty), seed=seed, n_ctx=n_ctx)
    print(f"[{tx},{ty}] {len(fibers)} fibres ({n_retried} relancées, "
          f"{n_ejected} perdues, {n_abandoned} abandonnées), "
          f"ctx={n_ctx}, max_h={max_h:.4f} mm, "
          f"{time.perf_counter()-t_start:.1f} s", flush=True)
    p.disconnect(physicsClientId=client)
    return out


def merge_tiles(tile_paths, n_total, seed, path_out="pile_tiles.npz"):
    """Fusionne les exports : chaque fibre a été déposée une seule fois
    (par la tuule de son xy). Le flag core est recalculé globalement
    (±MEASURE_SIDE/2), pas hérité des tuiles."""
    centers, directions, lengths, radii = [], [], [], [], 
    for tp in tile_paths:
        data = np.load(tp)
        centers.append(data["centers"])
        directions.append(data["directions"])
        lengths.append(data["lengths"])
        radii.append(data["radii"])
    C = np.concatenate(centers)
    U = np.concatenate(directions)
    L = np.concatenate(lengths)
    R = np.concatenate(radii)
    hm = MEASURE_SIDE / 2
    core = (np.abs(C[:, 0]) <= hm) & (np.abs(C[:, 1]) <= hm)
    np.savez(path_out, centers=C, directions=U, lengths=L, radii=R,
             core=core, code_version=TILE_CODE_VERSION, seed=seed,
             n_total=n_total)
    print(f"fusion : {len(C)} fibres, {core.sum()} en zone de mesure "
          f"({MEASURE_SIDE:.2f} mm de côté), hauteur "
          f"{(C[:,2]+L/2+R).max():.4f} mm -> {path_out}")
    return path_out


def main():
    n_total = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    n_per_tile = n_total // (TILE_GRID * TILE_GRID)
    t_start = time.perf_counter()

    tile_paths = []
    for wave, tiles in enumerate(WAVES):
        jobs = []
        for tx, ty in tiles:
            # contexte = exports des vagues précédentes (toutes les tuules
            # déjà faites ; load_margin_context filtre par la marge)
            jobs.append((tx, ty, n_per_tile, seed + 100 * tx + 10 * ty + wave,
                         list(tile_paths)))
        print(f"--- vague {wave} : {len(jobs)} tuules en parallèle ---",
              flush=True)
        with Pool(len(jobs)) as pool:
            tile_paths += pool.map(run_tile, jobs)

    merge_tiles(tile_paths, n_total, seed,
                path_out=f"pile_tiles_{n_total}_{seed}.npz")
    print(f"total : {time.perf_counter()-t_start:.1f} s")


if __name__ == "__main__":
    main()
