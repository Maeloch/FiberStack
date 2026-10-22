"""
settle.py — dépôt séquentiel de fibres de cellulose (PyBullet), v6.1.

Le dépôt reste celui validé par les runs 300/1000 (v5) : raycast de
l'ombre de la fibre (place_fiber) + chute courte + stabilisation,
relance des cas pathologiques (retrait + redépôt, abandon au 2e échec).

v6 — refactor de la revue de code :
  - géométrie, constantes et contrat npz importés de structure.py
  - garde de VITESSE ANGULAIRE pendant la stabilisation (v5 ne testait
    que la vitesse linéaire : une fibre en toupie à > 1 rad/ms est un
    artefact du solveur, pas de la physique)
  - hauteur convention EXACTE (top_exact) : le tas ne gonfle plus de
    L/2 par fibre verticale. Les max_h v6 ne sont PAS comparables aux
    max_h v5 (compte ~0,2-0,6 mm plus bas) — les AUTRES colonnes du
    log, elles, doivent être identiques
  - « aucun contact » traité comme une PERTE : v5 gelait le corps en
    plein air et le comptait comme déposé (l'origine des fibres
    « accrochées en l'air »)
  - contrat npz unifié : route="settle", code_version=6, reprise
    limitée aux routes settle/tiles ; les npz v5 sont refusés
    (attendu : on ne reprend pas un tas d'une autre génération)
  - N_RAYCAST_PROBES maintenu à 51 pour préserver la régression
    bit-exacte v5 -> v6 ; le passage à 64 (détecter un appui de 15 µm
    partout) se fera dans un commit séparé avec son propre banc

v6.1 — anti-narcolepsie (le correctif qui manquait à la route A) :
  - sleepThreshold=0 sur chaque fibre DÈS la création (drop_fiber).
    Le seuil de sommeil Bullet vaut 0,4 en unités brutes (0,4 mm/ms
    ici : au-dessus de V_START) et le corps s'endort après ~1000 pas
    tranquilles : une fibre en stabilisation lente s'endort, IGNORE
    la gravité et escroque le test d'immobilité de _wait_until_settled
    — dossier des fibres fantômes (241/999 flottants). Même correctif
    que press.py v2.5 et pressMT.py ; paramètre sondé une fois (si ce
    pybullet le rejette : avertissement unique puis silence).
  - CODE_VERSION 6 -> 7 : la physique change, donc la génération
    change. Conséquence du contrôle strict à la reprise : tous les
    npz v6 (routes settle/tiles) sont refusés. C'est VOULU — les
    baselines route A (pile_300_0, pile_1000_1) contiennent des
    dormeurs potentiels : à régénérer en v6.1, pas à reprendre.
  - resume_pile n'a pas besoin du correctif : les corps reconstruits
    sont gelés (masse 0), un corps statique ne s'endort pas.
  - DIVERGENCE ATTENDUE vs v6 : un run v6.1 (même seed, même n) peut
    différer d'un run v6 exactement là où v6 laissait une fibre
    s'endormir en cours de stabilisation (pas_stab plus longs, fibres
    qui finissent plus bas, moins de flottants). Toute AUTRE
    divergence est un bug.

NEUTRALITÉ COMPORTEMENTALE : un run v6 (même seed, même n) doit
reproduire le run v5 colonne par colonne (d, L, pas_stab, pos_z, core) ;
seule max_h change d'échelle. Toute divergence démarre soit sur une
ligne EXPLOSION (la garde angulaire attrape un cas que v5 laissait
passer — le relancer, ne pas l'ignorer), soit sur une fibre « sans
contact » que v5 gelait en l'air.

Système d'unités : (mm, mg, ms). Conversions :
  1 mm = 1e-3 m ; 1 mg = 1e-9 kg ; 1 mm/ms = 1 m/s ; g = 9.81e-3 mm/ms².
Cellulose : rho = 1.6 mg/mm³ (= 1600 kg/m³).

Nécessite : pip install pybullet scipy matplotlib numpy
"""

import math
import sys
import time

import numpy as np
from scipy.stats import truncnorm
import pybullet as p

from structure import (DOMAIN_SIDE, MEASURE_SIDE, CELLULOSE_DENSITY,
                       FIBER_DIAMETER_LOW, FIBER_DIAMETER_HIGH,
                       FIBER_LENGTH_MEAN, FIBER_LENGTH_STD,
                       top_exact, quat_align_z, quat_to_direction,
                       seg_dist, save_pile, load_pile)

# ------------------------------------------------------------- paramètres

CODE_VERSION = 7          # contrôle strict à la reprise (int obligatoire)
                          # v6.1 : physique nouvelle (plus de sommeil) ->
                          # génération nouvelle ; npz v6 refusés,
                          # baselines route A à régénérer

N_FIBERS = 2000
SEED = 0

# Dépôt
V_START = 0.2                 # mm/ms (= 0.2 m/s, vitesse terminale)
CONTACT_CLEARANCE = 0.05      # mm : posée à cette distance de son appui
N_RAYCAST_PROBES = 51         # MAINTENU à 51 (régression bit-exacte) ;
                              # passage à 64 = commit séparé + banc dédié

# Temps
DT = 2e-3                     # ms (= 2e-6 s)
NUM_SUBSTEPS = 4
FALL_MAX_STEPS = 800          # clearance à V_START = 125 pas, marge x6
STICK_MAX_STEPS = 30000

# Frottements « à fond » : à cette échelle un contact est un attachement.
FRICTION = 1.0
ROLLING_FRICTION = 0.3
SPINNING_FRICTION = 0.3
FRICTION_ANCHOR = False

NORMALIZE_MASSES = False      # True = masse 1 partout (test artefact solveur)

# 90 % de dissipation par seconde réelle, en unités locales.
DAMPING = 1.0 - 0.1 ** (1.0 / 1000.0)   # ≈ 0.00230

# Garde angulaire (v6) : > 1 rad/ms = 1000 rad/s — artefact du solveur.
ANG_VEL_LIMIT = 1.0           # rad/ms

_d_mu = 0.5 * (math.log(FIBER_DIAMETER_LOW) + math.log(FIBER_DIAMETER_HIGH))
_d_sigma = (math.log(FIBER_DIAMETER_HIGH) - math.log(FIBER_DIAMETER_LOW)) / 4.0


# ------------------------------------------------------------ échantillonnage

def sample_diameters(rng, n):
    """Lognormale troncée 15-25 µm. NE PAS MODIFIER le schéma de tirage :
    la régression bit-exacte v5 -> v6 dépend de la consommation exacte
    du RNG (ordre et tailles des tirages)."""
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
    """Normale tronquée à 2σ autour de 1 mm (σ = 15 %)."""
    return truncnorm.rvs(-2.0, 2.0, loc=FIBER_LENGTH_MEAN,
                         scale=FIBER_LENGTH_STD, size=n, random_state=rng)


# ------------------------------------------------------------ stabilisation

def _wait_until_settled(body, client, max_steps, move_tol=5e-4,
                        check_every=100, required_consecutive=5,
                        vel_limit=0.5, ang_vel_limit=ANG_VEL_LIMIT):
    """Repos = immobilité entre contrôles, 5 fois d'affilée. Retourne le
    pas de repos, max_steps si plafond, ou -1 si explosion détectée
    (vitesse linéaire > 2.5x V_START, ou vitesse angulaire > 1 rad/ms
    depuis la v6 : artefact du solveur, pas de la physique)."""
    prev_pos = None
    consecutive = 0
    for step in range(max_steps):
        p.stepSimulation(physicsClientId=client)
        if step % check_every == 0:
            lin, ang = p.getBaseVelocity(body, physicsClientId=client)
            if np.linalg.norm(lin) > vel_limit:
                return -1    # explosion linéaire : l'appelant retire le corps
            if np.linalg.norm(ang) > ang_vel_limit:
                return -1    # explosion angulaire (v6)
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
    le plus haut. Si l'ombre ne touche rien, l'appui est le sol (z=0).
    NB : le point de départ ne dépend PAS de max_h (le hitFraction est
    relatif au segment) — changer de convention de hauteur ne change pas
    la physique du dépôt."""
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

_SLEEP_SUPPORTED = None    # sonde v6.1 : None = pas encore testé


def _no_sleep(client, body):
    """Anti-narcolepsie (v6.1). Même correctif que press v2.5 :
    sleepThreshold=0 dès la création, paramètre sondé une fois.
    Un corps endormi pendant la stabilisation escroque
    _wait_until_settled (l'immobilité est mesurée sur un corps que
    la gravité n'atteint plus) — dossier 241/999 flottants."""
    global _SLEEP_SUPPORTED
    if _SLEEP_SUPPORTED is False:
        return
    try:
        p.changeDynamics(body, -1, sleepThreshold=0.0,
                         physicsClientId=client)
        _SLEEP_SUPPORTED = True
    except Exception:
        _SLEEP_SUPPORTED = False
        print("ATTENTION : sleepThreshold non supporté par ce pybullet "
              "— risque de fibres endormies en stabilisation "
              "(dossier narcolepsie)")


def drop_fiber(client, xy, z_start, direction, radius, length,
               fall_steps=None):
    """Crée la fibre (axe aligné sur direction), part à -V_START, premier
    contact, stabilisation, gel (masse -> 0). Renvoie
    (body, center, direction_finale, contact_trouvé, pas_utilisés),
    ou (None, ...) si fibre perdue (explosion/éjection/chute)."""
    quat = quat_align_z(direction)
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
    _no_sleep(client, body)     # v6.1 : jamais de mise en sommeil
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
    # PAS de gel d'un corps sans contact : v5 le figeait en plein air
    # (fibres « accrochées en l'air »). Le corps est laissé au monde,
    # c'est l'appelant qui décide (retrait + relance en v6).
    if not found_contact:
        return body, np.array(pos), _quat_to_direction(quat_final), False, 0
    p.changeDynamics(body, -1, mass=0, physicsClientId=client)  # figé
    return body, np.array(pos), _quat_to_direction(quat_final), found_contact, used


def _quat_to_direction(quat):
    """Alias local : garder le nom historique tout en utilisant la copie
    unique de structure.py."""
    return quat_to_direction(quat)


# ------------------------------------------------------------ reprise .npz

def resume_pile(client, npz_path):
    """Reconstruit un tas depuis un export .npz : chaque fibre devient une
    capsule gelée (masse 0) à sa position/orientation. À appeler après
    connect/gravity/plancher, avant la boucle de dépôt.
    v6 : réservé aux routes settle/tiles, code_version strict — les npz
    v5 sont refusés (regénérer, ou reprendre un npz v6).
    v6.1 : les npz v6 sont refusés aussi (génération nouvelle :
    anti-narcolepsie) — régénérer les baselines route A. Les corps
    reconstruits sont gelés (masse 0) : un corps statique ne dort
    pas, aucun correctif nécessaire ici.
    NB : pas de relaxation du tas reconstruit."""
    data, meta = load_pile(npz_path, routes=("settle", "tiles"))
    if "code_version" not in meta or int(meta["code_version"]) != CODE_VERSION:
        raise RuntimeError(
            f"npz code_version={meta.get('code_version')} != {CODE_VERSION} "
            f"— tas d'une autre génération de code, ne pas reprendre")

    centers, directions = data["centers"], data["directions"]
    lengths, radii, core = data["lengths"], data["radii"], data["core"]

    fibers = []
    max_h = 0.0
    for c, u, L, r, ck in zip(centers, directions, lengths, radii, core):
        quat = quat_align_z(u)
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
        # v6 : hauteur EXACTE (v5 : pire-cas vertical)
        max_h = max(max_h, top_exact(c, u, float(L), float(r)))
    return fibers, max_h


# ------------------------------------------------------------ boucle principale

def build_pile(n_fibers, seed=0, verbose=True, resume_npz=None):
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, -9.81e-3, physicsClientId=client)  # mm/ms²
    p.setPhysicsEngineParameter(fixedTimeStep=DT, numSubSteps=NUM_SUBSTEPS,
                                numSolverIterations=50, physicsClientId=client)

    half = DOMAIN_SIDE / 2.0
    half_floor = 2.0 * half            # plancher 2x le domaine
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

    # compteurs de la relance : initialisés HORS de la boucle fibre
    # (le bug v5 de les réinitialiser à chaque fibre masquait les totaux)
    n_retried = 0
    n_lost = 0
    n_abandoned = 0
    n_new = 0
    t_start = time.perf_counter()
    t_mark = t_start
    for i in range(n_fibers):
        radius = diameters[i] / 2.0
        length = lengths[i]
        xy = rng.uniform(-half, half, size=2)
        # axe horizontal, azimut uniforme : une fibre qui tombe dans une
        # suspension descend à plat, pas bout en premier.
        phi = rng.uniform(0, 2 * np.pi)
        direction = np.array([np.cos(phi), np.sin(phi), 0.0])

        # ---------- dépôt, avec relance des cas pathologiques ----------
        attempt = 0
        while True:
            attempt += 1
            z_start = place_fiber(client, xy, direction, radius, length, max_h)
            result = drop_fiber(client, xy, z_start, direction,
                                radius=radius, length=length,
                                fall_steps=FALL_MAX_STEPS)

            # cas 1 : perte (explosion / éjection / plancher / AUCUN
            # CONTACT — v6 : un corps sans contact est retiré et relancé,
            # v5 le gelait en l'air et le comptait déposé)
            lost = (result is None or result[0] is None or not result[3])
            if lost:
                if result is not None and result[0] is not None:
                    p.removeBody(result[0], physicsClientId=client)
                if attempt == 1:
                    n_retried += 1
                    xy = rng.uniform(-half, half, size=2)   # nouveau point
                    continue
                n_lost += 1
                if verbose:
                    print(f"fibre {i:4d}  d={diameters[i]*1e3:5.2f}um  "
                          f"L={length:5.3f}mm  <-- ABANDONNÉE (2 pertes)")
                break

            body, pos, dirf, found, used = result

            # cas 2 : pas au repos -> le corps est gelé en pleine
            # bascule : on le retire du monde AVANT de redéposer
            if used >= STICK_MAX_STEPS:
                p.removeBody(body, physicsClientId=client)
                if attempt == 1:
                    n_retried += 1
                    if verbose:
                        print(f"fibre {i:4d}  <-- PAS AU REPOS : redéposée")
                    xy = rng.uniform(-half, half, size=2)
                    continue
                n_abandoned += 1
                if verbose:
                    print(f"fibre {i:4d}  <-- ABANDONNÉE (2x pas au repos)")
                break

            # cas 3 : dépôt réussi
            top = top_exact(pos, dirf, length, radius)   # v6 : exacte
            max_h = max(max_h, top)
            core = (abs(pos[0]) <= MEASURE_SIDE / 2) and (abs(pos[1]) <= MEASURE_SIDE / 2)
            n_core += core
            n_new += 1
            fibers.append(dict(center=pos, direction=dirf, length=length,
                               radius=radius, core=core))
            if verbose and (i % 20 == 0):
                print(f"fibre {i:4d}  d={diameters[i]*1e3:5.2f}um  "
                      f"L={length:5.3f}mm  pas_stab={used:5d}  "
                      f"pos_z={pos[2]:7.4f}mm  max_h={max_h:7.4f}mm  core={core}")
            break

        # chrono : le coût par fibre doit rester ~constant
        if verbose and (i + 1) % 100 == 0:
            t_now = time.perf_counter()
            dt_chunk = t_now - t_mark
            t_mark = t_now
            print(f"    [chrono] fibres {i-98}-{i+1} : "
                  f"{dt_chunk/100*1000:.1f} ms/fibre "
                  f"(total {t_now-t_start:.1f} s)")

    if verbose:
        print(f"=== {len(fibers)} fibres déposées "
              f"({n_retried} relancées, {n_lost} perdues, "
              f"{n_abandoned} abandonnées), "
              f"{n_core} en zone de mesure ({MEASURE_SIDE:.2f} mm de côté), "
              f"hauteur {max_h:.4f} mm, "
              f"{time.perf_counter()-t_start:.1f} s ===")
        # invariant du banc : ce qui sort + ce qui échoue = ce qui entre
        if n_new + n_lost + n_abandoned == n_fibers:
            print(f"invariant : {n_new} déposées + {n_lost} perdues + "
                  f"{n_abandoned} abandonnées = {n_fibers} OK")
        else:
            print(f"*** INVARIANT VIOLÉ : {n_new} + {n_lost} + "
                  f"{n_abandoned} != {n_fibers} ***")
    p.disconnect(physicsClientId=client)
    return fibers


# ------------------------------------------------------------ export & contrôle

def export_npz(fibers, path, seed, n_total):
    """Délègue au contrat unifié de structure.py : route obligatoire,
    valeurs None omises, code_version vérifié à la reprise."""
    save_pile(path,
              centers=[f["center"] for f in fibers],
              directions=[f["direction"] for f in fibers],
              lengths=[f["length"] for f in fibers],
              radii=[f["radius"] for f in fibers],
              core=[f["core"] for f in fibers],
              route="settle", code_version=CODE_VERSION,
              seed=seed, n_total=n_total)


def check_penetration(fibers, cell=0.5):
    """Pénétration maximale entre fibres voisines (distance axe-axe -
    somme des rayons). Doit rester << rayon (~µm). NB : une pénétration
    résiduelle de ~15-25 µm est SYSTÉMIQUE (repos du solveur Bullet sur
    masses 1e-7 mg) : constante sur 3 runs indépendants, à documenter
    comme biais connu plutôt qu'à corriger.
    v6 : seg_dist vient de structure.py (copie unique)."""
    cells = {}
    for idx, f in enumerate(fibers):
        c = f["center"]
        reach = f["length"] / 2 + f["radius"]
        ix0 = int((c[0] - reach) // cell); ix1 = int((c[0] + reach) // cell)
        iy0 = int((c[1] - reach) // cell); iy1 = int((c[1] + reach) // cell)
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                cells.setdefault((ix, iy), []).append(idx)

    def ends(f):
        return (f["center"] - 0.5 * f["length"] * f["direction"],
                f["center"] + 0.5 * f["length"] * f["direction"])

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
                d = seg_dist(*ends(f1), *ends(f2))
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
    # Banc de validation : 300 fibres seed 0 (régression colonne par
    # colonne vs log v5), puis 1000 (comparaison pile_1000_0/1).
    # Arguments optionnels : python settle.py [n] [seed]
    n = int(sys.argv[1]) if len(sys.argv) > 1 else N_FIBERS
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else SEED

    fibers = build_pile(n_fibers=n, seed=seed)
    export_npz(fibers, f"pile_{n}_{seed}.npz", seed=seed, n_total=n)
    check_penetration(fibers, cell=0.5)
    plot_fibers(fibers, save_path=f"pile_{n}_{seed}.png")
    print("=== terminé ===")
