"""
structure.py — module partagé FiberStack : constantes du domaine,
géométrie des capsules, entrée/sortie des .npz.

Source unique pour ce qui était dupliqué ou incohérent entre
settle / analyse / press / crossply (revue de code v6) :
  - DOMAIN_SIDE / MEASURE_SIDE (fini le littéral 7.0 dans analyse.py)
  - hauteur EXACTE d'une capsule : top_exact / bottom_exact
    (l'ancienne convention « pire-cas vertical » c_z + L/2 gonflait
    de ~L/2 toute fibre couchée — les hauteurs des logs v5 et les
    hauteurs v6 ne sont PAS comparables directement)
  - quaternion d'alignement (correctif antiparallèle : direction ~ -z
    ne renvoie plus l'identité)
  - seg_dist : distance segment-segment (Ericson RTCD 5.1.9),
    copie unique (fini les copies divergentes dans settle/press)
  - save_pile / load_pile : contrat npz unifié, `route` obligatoire
    à l'écriture, métadonnées affichées au chargement, refus des
    routes non autorisées quand une liste est passée

Unités : (mm, mg, ms).
"""

import numpy as np

# ------------------------------------------------------------- domaine

DOMAIN_SIDE = 10.0                                  # mm
BUFFER_MARGIN = 1.5                                  # mm (approx. périodique)
MEASURE_SIDE = DOMAIN_SIDE - 2 * BUFFER_MARGIN       # 7.0 mm

CELLULOSE_DENSITY = 1.6                             # mg/mm³ (1600 kg/m³)
FIBER_DIAMETER_LOW = 15e-3                          # mm
FIBER_DIAMETER_HIGH = 25e-3                         # mm
FIBER_LENGTH_MEAN = 1.0                            # mm
FIBER_LENGTH_STD = 0.15                            # mm
FILTER_THICKNESS = 0.3                              # mm (cible C569)

# ------------------------------------------------------------ géométrie


def top_exact(center, direction, length, radius):
    """Sommet RÉEL d'une capsule : extension verticale exacte.
    Convention unique du projet (S1 de la revue)."""
    return float(np.asarray(center)[2]
                 + 0.5 * float(length) * abs(np.asarray(direction)[2])
                 + float(radius))


def bottom_exact(center, direction, length, radius):
    """Base réelle d'une capsule (extension verticale exacte)."""
    return float(np.asarray(center)[2]
                 - 0.5 * float(length) * abs(np.asarray(direction)[2])
                 - float(radius))


def quat_align_z(direction):
    """Quaternion (x,y,z,w) alignant l'axe local z d'une GEOM_CAPSULE
    sur `direction`. Correctif antiparallèle : direction ~ -z renvoie
    la rotation de 180° autour de x (v5 renvoyait l'identité, ce qui
    alignait sur +z — latent, déclenché par tout tilt >= 90°)."""
    z_axis = np.array([0.0, 0.0, 1.0])
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    cos_a = float(np.clip(z_axis @ d, -1.0, 1.0))
    if cos_a < -1.0 + 1e-9:
        return [1.0, 0.0, 0.0, 0.0]         # rotation de pi autour de x
    axis = np.cross(z_axis, d)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-12:
        return [0.0, 0.0, 0.0, 1.0]         # déjà alignée sur +z
    axis = axis / axis_norm
    angle = np.arccos(cos_a)
    return list(np.sin(angle / 2.0) * axis) + [float(np.cos(angle / 2.0))]


def quat_to_direction(quat):
    """Direction de l'axe local z d'un corps, depuis son quaternion."""
    x, y, z, w = quat
    qv = np.array([x, y, z])
    v = np.array([0.0, 0.0, 1.0])
    t = 2.0 * np.cross(qv, v)
    v_rot = v + w * t + np.cross(qv, t)
    return v_rot / np.linalg.norm(v_rot)


def seg_dist(p1, q1, p2, q2):
    """Distance minimale entre les segments [p1,q1] et [p2,q2]
    (Ericson, Real-Time Collision Detection 5.1.9). Copie unique :
    toute correction se fait ICI et nulle part ailleurs."""
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
            s = (np.clip((b * f_ - c * e) / denom, 0.0, 1.0)
                 if denom > EPS else 0.0)
            t = (b * s + f_) / e
            if t < 0.0:
                t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)
    return float(np.linalg.norm(p1 + s * d1 - (p2 + t * d2)))

# ------------------------------------------------------------ npz


def save_pile(path, centers, directions, lengths, radii, core,
              route, code_version, **params):
    """Export unifié. `route` OBLIGATOIRE : settle | tiles | press |
    crossply | routed. `params` : tout ce qui décrit la génération (seed,
    n_total, target, spacing...) — les valeurs None sont omises
    (une métadonnée absente vaut mieux qu'une métadonnée fausse)."""
    scalars = {k: v for k, v in params.items() if v is not None}
    np.savez(path,
             centers=np.asarray(centers, dtype=float),
             directions=np.asarray(directions, dtype=float),
             lengths=np.asarray(lengths, dtype=float),
             radii=np.asarray(radii, dtype=float),
             core=np.asarray(core, dtype=bool),
             route=route, code_version=code_version, **scalars)


def load_pile(path, routes=None):
    """Charge un npz + affiche ses métadonnées. Si `routes` est donné,
    refuse les fichiers dont la route n'est pas dans la liste (évite
    de reprendre un tas pressé comme un tas de dépôt). Les fichiers
    sans route (ancien format, ex. settle_tiles v1) passent avec un
    avertissement."""
    data = np.load(path, allow_pickle=False)
    meta = {k: data[k].item() for k in data.files if data[k].ndim == 0}
    route = meta.get("route", "inconnue")
    print(f"[npz] {path} : route={route}"
          + (", " + ", ".join(f"{k}={v}" for k, v in meta.items()
                              if k != "route") if meta else ""))
    if route == "inconnue":
        print("[npz] AVERTISSEMENT : fichier sans route (ancien format)")
    elif routes is not None and route not in routes:
        raise RuntimeError(f"route « {route} » non autorisée ici "
                           f"(attendu : {routes}) — {path}")
    return data, meta
