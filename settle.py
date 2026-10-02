"""
settle_si4.py — dépôt séquentiel de fibres de cellulose, PyBullet.
Version "campagne" : la chute libre est remplacée par un raycast de
l'ombre de la fibre (place_fiber) + chute courte + micro-relaxation.
La physique de calage est intégralement conservée.

Système d'unités : (mm, mg, ms). Conversions :
  1 mm = 1e-3 m ; 1 mg = 1e-9 kg ; 1 mm/ms = 1 m/s ; g = 9.81e-3 mm/ms².
Cellulose : rho = 1.6 mg/mm³ (= 1600 kg/m³).

Historique des correctifs :
  - damping converti en unités locales
  - axe des fibres HORIZONTAL à la génération (azimut uniforme)
  - plancher 2x le domaine + garde anti-chute sous le plancher
  - garde anti-éjection + détection d'explosion (vel_limit=0.5,
    soit 2.5x V_START — un seuil de 5 mm/ms était incohérent avec
    les unités : une éjection à 9 mm ne demande que 0.42 mm/ms)
  - repos mesuré par déplacement (insensible au micro-jitter du solveur)
  - frictionAnchor désactivé (impulsions parasites sur contacts empilés)
  - reprise de tas depuis un export .npz, contrôle code_version
  - chronométrage par fibre
  - place_fiber (raycast) : chute courte depuis l'ombre de la fibre,
    O(N) pas de chute économisés — coût par fibre quasi constant
  - N_RAYCAST_PROBES=51 : espacement ~20 µm, détecte tout support
    perpendiculaire (9 sondes laissaient 4 chances sur 5 de rater un
    support de 20 µm de large -> fibre créée en chevauchement ->
    éjection par le solveur)

Nécessite : pip install pybullet scipy matplotlib numpy
"""

import math
import time
import numpy as np
from scipy.stats import truncnorm
import pybullet as p

# ------------------------------------------------------------- paramètres

CODE_VERSION = 6          # entier obligatoire : int(code_version) comparé
                          # au chargement ; une décimale casserait le contrôle

N_FIBERS = 1000
SEED = 1

FILTER_THICKNESS = 0.3        # mm (300 µm, cible finale)

# Diamètre fibre : lognormale tronquée 15-25 µm
FIBER_DIAMETER_LOW = 15e-3    # mm
FIBER_DIAMETER_HIGH = 25e-3   # mm
_d_mu = 0.5 * (math.log(FIBER_DIAMETER_LOW) + math.log(FIBER_DIAMETER_HIGH))
_d_sigma = (math.log(FIBER_DIAMETER_HIGH) - math.log(FIBER_DIAMETER_LOW)) / 4.0

# Longueur : normale tronquée à 2σ autour de 1 mm (σ = 15 %)
FIBER_LENGTH_MEAN = 1.0       # mm
FIBER_LENGTH_STD = 0.15       # mm

CELLULOSE_DENSITY = 1.6       # mg/mm³

# Domaine ~10x la longueur de fibre ; marge tampon = approximation des
# conditions périodiques (fibres de la marge simulées, exclues des stats)
DOMAIN_SIDE = 10.0            # mm
BUFFER_MARGIN = 1.5 * FIBER_LENGTH_MEAN
MEASURE_SIDE = DOMAIN_SIDE - 2 * BUFFER_MARGIN

# Dépôt
V_START = 0.2                # mm/ms (= 0.2 m/s, vitesse terminale)
CONTACT_CLEARANCE = 0.05      # mm : fibre posée à cette distance de son appui
N_RAYCAST_PROBES = 51        # espacement ~20 µm le long de l'axe

# Temps
DT = 2e-3                     # ms (= 2e-6 s)
NUM_SUBSTEPS = 4
FALL_MAX_STEPS = 800          # chute courte : clearance à V_START = 250
                              # pas, marge x3 (raycast approximatif)
STICK_MAX_STEPS = 30000

# Frottements "à fond" : à cette échelle un contact est un attachement.
FRICTION = 1.0
ROLLING_FRICTION = 0.3
SPINNING_FRICTION = 0.3
FRICTION_ANCHOR = False

NORMALIZE_MASSES = False      # True = masse 1 partout (test artefact solveur)

# damping : 90 % de dissipation par seconde réelle, converti en unités
# locales — Bullet traite le damping par "seconde" simulée.
DAMPING = 1.0 - 0.1 ** (1.0 / 1000.0)   # ≈ 0.00230


# ------------------------------------------------------------ échantillonnage

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


# ------------------------------------------------------------ quaternions

def _quat_align_z(direction):
    """Quaternion (x,y,z,w) alignant l'axe local z de la capsule sur
    `direction` (l'axe local d'une GEOM_CAPSULE est z)."""
    z_axis = np.array([0.0, 0.0, 1.0])
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    axis = np.cross(z_axis, d)
    axis_norm = np.linalg.norm(axis)
    angle = np.arccos(np.clip(z_axis @ d, -1.0, 1.0))
    if axis_norm < 1e-12:
        return [0.0, 0.0, 0.0, 1.0]
    axis = axis / axis_norm
    return list(np.sin(angle / 2.0) * axis) + [float(np.cos(angle / 2.0))]


def _quat_to_direction(quat):
    x, y, z, w = quat
    qv = np.array([x, y, z])
    v = np.array([0.0, 0.0, 1.0])
    t = 2.0 * np.cross(qv, v)
    v_rot = v + w * t + np.cross(qv, t)
    return v_rot / np.linalg.norm(v_rot)


# ------------------------------------------------------------ stabilisation

def _wait_until_settled(body, client, max_steps, move_tol=5e-4,
                        check_every=100, required_consecutive=5,
                        vel_limit=0.5):
    """Repos = immobilité entre contrôles, 5 fois d'affilée. Retourne le
    pas de repos, max_steps si plafond, ou -1 si explosion détectée
    (vitesse > vel_limit = 2.5x V_START après contact : artefact du
    solveur, pas de la physique)."""
    prev_pos = None
    consecutive = 0
    for step in range(max_steps):
        p.stepSimulation(physicsClientId=client)
        if step % check_every == 0:
            lin, ang = p.getBaseVelocity(body, physicsClientId=client)
            if np.linalg.norm(lin) > vel_limit:
                return -1    # explosion : l'appelant retire le corps
            pos, quat = p.getBasePositionAndOrientation(body, physicsClientId=client)
            if prev_pos is not None:
                moved = max(abs(pos[k] - prev_pos[k]) for k in range(3))
                if moved < move_tol:
                    consecutive += 1
                    if consecutive >= required_consecutive:
                        return step
                else:
                    consecutive = 0
            prev_pos = pos
    return max_steps


# ------------------------------------------------------------ placement raycast

def place_fiber(client, xy, direction, radius, length, max_h):
    """Sonde l'ombre verticale de la fibre (raycasts le long de l'axe,
    N_RAYCAST_PROBES points, du sommet du tas vers le sol) et renvoie le
    z de départ tel que la fibre arrive à CONTACT_CLEARANCE de son appui
    le plus haut. Si l'ombre ne touche rien, l'appui est le sol (z=0)."""
    u = np.asarray(direction, dtype=float)
    z_top = max_h + length          # les rayons partent au-dessus de tout
    ts = np.linspace(-0.5, 0.5, N_RAYCAST_PROBES)
    starts = []
    for t in ts:
        px = xy[0] + t * length * u[0]
        py = xy[1] + t * length * u[1]
        starts.append([px, py, z_top])
    ends = [[s[0], s[1], -0.02] for s in starts]
    hits = p.rayTestBatch(starts, ends, physicsClientId=client)
    # hit = (bodyId, linkIndex, hitFraction, position, normal)
    z_support = 0.0
    for h in hits:
        if h[0] != -1:
            z_impact = z_top + h[2] * (-0.02 - z_top)
            z_support = max(z_support, z_impact)
    # le point bas de la fibre (axe horizontal : centre - rayon)
    return z_support + radius + CONTACT_CLEARANCE


# ------------------------------------------------------------ dépôt

def drop_fiber(client, xy, z_start, direction, radius, length,
               fall_steps=None):
    """Crée la fibre (axe aligné sur direction), part à -V_START, premier
    contact, stabilisation, gel (masse -> 0). Renvoie
    (body, center, direction_finale, contact_trouvé, pas_utilisés),
    ou (None, ...) si fibre perdue (explosion/éjection/chute)."""
    quat = _quat_align_z(direction)
    col = p.createCollisionShape(p.GEOM_CAPSULE, radius=radius,
                                 height=max(length - 2 * radius, radius * 0.1),
                                 physicsClientId=client)
    mass = (CELLULOSE_DENSITY * math.pi * radius ** 2 * length
            if not NORMALIZE_MASSES else 1.0)
    body = p.createMultiBody(baseMass=mass, baseCollisionShapeIndex=col,
                             basePosition=[xy[0], xy[1], z_start],
                             baseOrientation=quat, physicsClientId=client)
    p.changeDynamics(body, -1, lateralFriction=FRICTION,
                     rollingFriction=ROLLING_FRICTION,
                     spinningFriction=SPINNING_FRICTION,
                     frictionAnchor=FRICTION_ANCHOR,
                     restitution=0.0,
                     linearDamping=DAMPING, angularDamping=DAMPING,
                     physicsClientId=client)
    p.resetBaseVelocity(body, [0.0, 0.0, -V_START], [0.0, 0.0, 0.0],
                        physicsClientId=client)
    if fall_steps is None:
        fall_steps = FALL_MAX_STEPS

    found_contact = False
    for _ in range(fall_steps):
        p.stepSimulation(physicsClientId=client)
        if p.getContactPoints(bodyA=body, physicsClientId=client):
            found_contact = True
            break

    used = 0
    if found_contact:
        used = _wait_until_settled(body, client, STICK_MAX_STEPS)
        if used == -1:
            print("  <-- EXPLOSION DU SOLVEUR : fibre ignorée")
            p.removeBody(body, physicsClientId=client)
            return None, None, None, False, 0

    pos, quat_final = p.getBasePositionAndOrientation(body, physicsClientId=client)
    # garde-fou anti-éjection : jamais plus haut que son départ
    if found_contact and pos[2] > z_start + 0.5 * length:
        print(f"  <-- ÉJECTION SUSPECTE (z={pos[2]:.3f} > départ) : fibre ignorée")
        p.removeBody(body, physicsClientId=client)
        return None, None, None, False, used
    # garde anti-chute hors domaine : sous le plancher = perdue
    if found_contact and pos[2] < -radius:
        print(f"  <-- PERDUE SOUS LE PLANCHER (z={pos[2]:.3f}) : fibre ignorée")
        p.removeBody(body, physicsClientId=client)
        return None, None, None, False, used
    p.changeDynamics(body, -1, mass=0, physicsClientId=client)  # figé
    return body, np.array(pos), _quat_to_direction(quat_final), found_contact, used


# ------------------------------------------------------------ reprise .npz

def resume_pile(client, npz_path):
    """Reconstruit un tas depuis un export .npz : chaque fibre devient une
    capsule gelée (masse 0) à sa position/orientation. À appeler après
    connect/gravity/plancher, avant la boucle de dépôt.
    NB : pas de relaxation du tas reconstruit. Le seed du tirage n'est
    PAS dans le npz : utilise un seed différent par lot."""
    data = np.load(npz_path)
    if int(data["code_version"]) != CODE_VERSION:
        raise RuntimeError("npz d'une autre génération de code — ne pas reprendre")

    centers, directions = data["centers"], data["directions"]
    lengths, radii, core = data["lengths"], data["radii"], data["core"]

    fibers = []
    max_h = 0.0
    for c, u, L, r, ck in zip(centers, directions, lengths, radii, core):
        quat = _quat_align_z(u)
        col = p.createCollisionShape(p.GEOM_CAPSULE, radius=float(r),
                                     height=max(float(L) - 2 * float(r),
                                                float(r) * 0.1),
                                     physicsClientId=client)
        body = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col,
                                 basePosition=[float(c[0]), float(c[1]), float(c[2])],
                                 baseOrientation=quat, physicsClientId=client)
        p.changeDynamics(body, -1, lateralFriction=FRICTION,
                         rollingFriction=ROLLING_FRICTION,
                         spinningFriction=SPINNING_FRICTION, restitution=0.0,
                         physicsClientId=client)
        fibers.append(dict(center=np.array(c), direction=np.array(u),
                           length=float(L), radius=float(r),
                           core=bool(ck)))
        max_h = max(max_h, float(c[2]) + float(L) / 2 + float(r))
    return fibers, max_h


# ------------------------------------------------------------ boucle principale

def build_pile(n_fibers, seed=0, verbose=True, resume_npz=None):
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, -9.81e-3, physicsClientId=client)
    p.setPhysicsEngineParameter(fixedTimeStep=DT, numSubSteps=NUM_SUBSTEPS,
                                numSolverIterations=50, physicsClientId=client)

    half = DOMAIN_SIDE / 2.0
    half_floor = 2.0 * half
    floor_col = p.createCollisionShape(p.GEOM_BOX,
                                       halfExtents=[half_floor, half_floor, 1e-2],
                                       physicsClientId=client)
    floor = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=floor_col,
                              basePosition=[0, 0, -1e-2], physicsClientId=client)
    p.changeDynamics(floor, -1, lateralFriction=FRICTION,
                     rollingFriction=ROLLING_FRICTION,
                     spinningFriction=SPINNING_FRICTION, restitution=0.0,
                     physicsClientId=client)

    rng = np.random.default_rng(seed)
    diameters = sample_diameters(rng, n_fibers)
    lengths = sample_lengths(rng, n_fibers)

    if resume_npz is not None:
        fibers, max_h = resume_pile(client, resume_npz)
        n_core = sum(f["core"] for f in fibers)
        print(f"reprise : {len(fibers)} fibres reconstruites, "
              f"max_h={max_h:.4f} mm")
    else:
        fibers, max_h, n_core = [], 0.0, 0

    n_ejected = 0      # perdues après 2 tentatives (éjection/explosion/plancher)
    n_retried = 0      # tentatives relancées (1re tentative ratée)
    n_abandoned = 0    # fibres abandonnées après 2 tentatives
    t_start = time.perf_counter()
    t_mark = t_start
    for i in range(n_fibers):
        radius = diameters[i] / 2.0
        length = lengths[i]
        xy = rng.uniform(-half, half, size=2)
        phi = rng.uniform(0, 2 * np.pi)
        direction = np.array([np.cos(phi), np.sin(phi), 0.0])

        attempt = 0
        while True:
            z_start = place_fiber(client, xy, direction, radius, length, max_h)
            result = drop_fiber(client, xy, z_start, direction,
                                radius=radius, length=length,
                                fall_steps=FALL_MAX_STEPS)
            attempt += 1

            # --- cas pathologique : perdue OU aucun contact OU pas au repos
            pathological = (result is None or result[0] is None)
            if not pathological:
                body, pos, dirf, found, used = result
                if not found or used >= STICK_MAX_STEPS:
                    pathological = True

            if pathological:
                # retirer le corps : drop_fiber l'a déjà fait pour les
                # pertes ; sinon (gelé en vol), on le fait ici
                if result is not None and result[0] is not None:
                    p.removeBody(result[0], physicsClientId=client)
                if attempt == 1:
                    n_retried += 1
                    if verbose and not found:
                        print(f"fibre {i:4d}  <-- AUCUN CONTACT : redéposée")
                    elif verbose and result is not None and result[0] is not None:
                        print(f"fibre {i:4d}  <-- PAS AU REPOS : redéposée")
                    xy = rng.uniform(-half, half, size=2)
                    continue
                # 2e échec : abandonnée, PAS ajoutée au tas
                if result is not None and result[0] is None:
                    n_ejected += 1      # éjection/explosion/plancher
                else:
                    n_abandoned += 1    # non-repos ou sans contact
                if verbose:
                    print(f"fibre {i:4d}  d={diameters[i]*1e3:5.2f}um  "
                          f"L={length:5.3f}mm  <-- ABANDONNÉE (2 échecs)")
                break

            # --- dépôt réussi
            top = pos[2] + length / 2 + radius
            max_h = max(max_h, top)
            core = (abs(pos[0]) <= MEASURE_SIDE / 2) and (abs(pos[1]) <= MEASURE_SIDE / 2)
            n_core += core
            fibers.append(dict(center=pos, direction=dirf, length=length,
                               radius=radius, core=core))
            if verbose and (i % 20 == 0):
                print(f"fibre {i:4d}  d={diameters[i]*1e3:5.2f}um  "
                      f"L={length:5.3f}mm  pas_stab={used:5d}  "
                      f"pos_z={pos[2]:7.4f}mm  max_h={max_h:7.4f}mm  core={core}")
            break

        if verbose and (i + 1) % 100 == 0:
            t_now = time.perf_counter()
            dt_chunk = t_now - t_mark
            t_mark = t_now
            print(f"    [chrono] fibres {i-98}-{i+1} : "
                  f"{dt_chunk/100*1000:.1f} ms/fibre "
                  f"(total {t_now-t_start:.1f} s)")

    if verbose:
        print(f"=== {len(fibers)} fibres déposées "
              f"({n_retried} relancées, {n_ejected} perdues, "
              f"{n_abandoned} abandonnées), "
              f"{n_core} en zone de mesure ({MEASURE_SIDE:.2f} mm de côté), "
              f"hauteur {max_h:.4f} mm, "
              f"{time.perf_counter()-t_start:.1f} s ===")
    p.disconnect(physicsClientId=client)
    return fibers

# ------------------------------------------------------------ export & contrôle

def export_npz(fibers, path="pile.npz"):
    np.savez(path,
             centers=np.array([f["center"] for f in fibers]),
             directions=np.array([f["direction"] for f in fibers]),
             lengths=np.array([f["length"] for f in fibers]),
             radii=np.array([f["radius"] for f in fibers]),
             core=np.array([f["core"] for f in fibers]),
             code_version=CODE_VERSION, seed=SEED, n_total=N_FIBERS)


def check_penetration(fibers, cell=0.5):
    """Pénétration maximale entre fibres voisines (distance axe-axe -
    somme des rayons). Doit rester << rayon (~µm). NB : une pénétration
    résiduelle de ~15-25 µm est SYSTÉMIQUE (repos du solveur Bullet sur
    masses 1e-7 mg) : constante sur 3 runs indépendants, à documenter
    comme biais connu plutôt qu'à corriger."""
    cells = {}
    for idx, f in enumerate(fibers):
        c = f["center"]
        reach = f["length"] / 2 + f["radius"]
        ix0 = int((c[0] - reach) // cell); ix1 = int((c[0] + reach) // cell)
        iy0 = int((c[1] - reach) // cell); iy1 = int((c[1] + reach) // cell)
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                cells.setdefault((ix, iy), []).append(idx)

    def seg_dist(f1, f2):
        # points les plus proches entre deux segments (Ericson RTCD §5.1.9)
        p1 = f1["center"] - 0.5 * f1["length"] * f1["direction"]
        q1 = f1["center"] + 0.5 * f1["length"] * f1["direction"]
        p2 = f2["center"] - 0.5 * f2["length"] * f2["direction"]
        q2 = f2["center"] + 0.5 * f2["length"] * f2["direction"]
        d1, d2, r = q1 - p1, q2 - p2, p1 - p2
        a, e, f_ = float(d1 @ d1), float(d2 @ d2), float(d2 @ r)
        EPS = 1e-12
        if a <= EPS and e <= EPS:
            return float(np.linalg.norm(p1 - p2))
        if a <= EPS:
            s, t = 0.0, np.clip(f_ / e, 0.0, 1.0)
        else:
            c = float(d1 @ r)
            if e <= EPS:
                t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
            else:
                b = float(d1 @ d2)
                denom = a * e - b * b
                s = np.clip((b * f_ - c * e) / denom, 0.0, 1.0) if denom > EPS else 0.0
                t = (b * s + f_) / e
                if t < 0.0:
                    t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
                elif t > 1.0:
                    t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)
        c1 = p1 + s * d1
        c2 = p2 + t * d2
        return float(np.linalg.norm(c1 - c2))

    max_pen, n_pairs = 0.0, 0
    seen = set()
    for members in cells.values():
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                i, j = members[a], members[b]
                if i >= j:
                    continue
                key = (i, j)
                if key in seen:
                    continue
                seen.add(key)
                f1, f2 = fibers[i], fibers[j]
                if abs(f1["center"][2] - f2["center"][2]) > (0.5 * (f1["length"] + f2["length"])
                                                             + f1["radius"] + f2["radius"]):
                    continue
                d = seg_dist(f1, f2)
                pen = (f1["radius"] + f2["radius"]) - d
                n_pairs += 1
                max_pen = max(max_pen, pen)
    print(f"pénétration max : {max_pen*1e3:.4f} µm sur {n_pairs} paires "
          f"(attendu : << rayon ≈ {min(f['radius'] for f in fibers)*1e3:.1f} µm)")
    return max_pen


# ------------------------------------------------------------ visualisation

def plot_fibers(fibers, save_path=None, n_sides=8, core_only=False):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if core_only:
        fibers = [f for f in fibers if f["core"]]

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("viridis")

    heights = np.array([f["center"][2] for f in fibers])
    h_min, h_max = heights.min(), heights.max()
    span_h = max(h_max - h_min, 1e-12)

    all_pts = []
    for f, h in zip(fibers, heights):
        color = cmap((h - h_min) / span_h)
        alpha = 0.95 if f["core"] else 0.25
        u = f["direction"] / np.linalg.norm(f["direction"])
        arbitrary = np.array([0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        v1 = np.cross(u, arbitrary); v1 /= np.linalg.norm(v1)
        v2 = np.cross(u, v1)
        theta = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
        circle = np.outer(np.cos(theta), v1) + np.outer(np.sin(theta), v2)
        e1 = f["center"] - 0.5 * f["length"] * u
        e2 = f["center"] + 0.5 * f["length"] * u
        ring1 = e1 + f["radius"] * circle
        ring2 = e2 + f["radius"] * circle
        faces = [[ring1[i], ring1[(i + 1) % n_sides],
                  ring2[(i + 1) % n_sides], ring2[i]] for i in range(n_sides)]
        faces += [list(ring1), list(ring2)]
        ax.add_collection3d(Poly3DCollection(faces, facecolor=color,
                                            edgecolor="none", alpha=alpha))
        all_pts.extend([e1, e2])

    m = MEASURE_SIDE / 2.0
    ax.plot([-m, m, m, -m, -m], [-m, -m, m, m, -m], [0, 0, 0, 0, 0],
            color="crimson", linewidth=1.2, linestyle="--")

    all_pts = np.array(all_pts)
    lo, hi = all_pts.min(axis=0), all_pts.max(axis=0)
    lo[2] = 0.0
    pad = 0.05 * max(hi - lo)
    lo, hi = lo - pad, hi + pad

    xx, yy = np.meshgrid([lo[0], hi[0]], [lo[1], hi[1]])
    ax.plot_surface(xx, yy, np.zeros_like(xx), color="0.85", alpha=0.3)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(0, hi[2])
    try:
        ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], (hi[2]) * 5))
    except AttributeError:
        pass
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("z (mm)")
    ax.set_title(f"{len(fibers)} fibres ({sum(f['core'] for f in fibers)} en zone de mesure)")
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig, ax


# ------------------------------------------------------------ main

if __name__ == "__main__":
    # 1) Validation : 1000 fibres, comparer à pile_0.npz / pile.npz.
    # 2) Campagne par lots : seed différent à chaque lot, reprise npz.
    fibers = build_pile(n_fibers=N_FIBERS, seed=SEED)
    export_npz(fibers, f"pile_{N_FIBERS}_{SEED}.npz")
    check_penetration(fibers, cell=0.5)
    plot_fibers(fibers, save_path=f"pile_{N_FIBERS}_{SEED}.png")
    print("=== terminé ===")
