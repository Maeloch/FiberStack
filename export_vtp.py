"""
export_vtp.py — conversion d'une pile FiberStack (.npz) en VTK PolyData (.vtp)
pour exploration interactive 3D dans ParaView.

Principe : chaque fibre devient une ligne de 2 points (extrémités de l'axe),
avec le rayon porté en donnée par point. C'est le filtre "Tube" de ParaView qui
gonfle chaque ligne en cylindre, en lisant le rayon dans le tableau "radius".
On n'écrit donc PAS la géométrie des cylindres : fichier minuscule (quelques
centaines de ko pour 4000 fibres), tubes calculés à la volée par ParaView.

Usage :
    python export_vtp.py pile_300_0.npz
    python export_vtp.py pile_1000_0.npz -o tas1000.vtp --scale 3
    python export_vtp.py crossply_62um_12c.npz --core
    python export_vtp.py pressed_random_1000_1.npz --scale 5

Options :
    -o / --output   nom du .vtp (défaut : <nom_du_npz>.vtp)
    --scale S       multiplicateur du rayon (défaut 1.0 = rayons vrais).
                    ATTENTION : à l'échelle vraie, un rayon de 10 µm dans un
                    domaine de 10 mm = 0,1 % de la boîte — les tubes sont
                    fins mais visibles ; S=2..5 aide à l'inspection. Le
                    facteur ne doit JAMAIS servir à une métrique, uniquement
                    au confort visuel (les volumes affichés seraient faux).
    --core          n'exporter que les fibres de la zone de mesure (7×7 mm).

Dans ParaView (>= 5.6) :
    1. Open → sélectionner le .vtp → Apply.
    2. Sélectionner l'objet puis Filters → Alphabetical → Tube → Apply,
       avec dans le panneau Properties :
         - "Vary Radius" = "By Absolute Scalar"
         - "Radius Array" = "radius"
         - "Number of sides" = 12 (6 par défaut suffit en vue lointaine)
         - "Capping" = coché (extrémités plates fermées)
       (le rayon du filtre lui-même est ignoré en mode By Absolute Scalar)
    3. Coloration : dans la toolbar, le sélecteur de couleur → "core"
       (1 = zone de mesure, 0 = marges) ; ou une coordonnée (PointsZarray)
       pour voir la stratification verticale.
    4. Exploration : Slice (coupe à z fixé,Filter → Common), Clip (moitié
       du tas), ou Spin View pour tourner. Rendu : Surface avec éclairage
       par défaut suffit.
    5. Capture : File → Save Screenshot (résolution libre).

Notes :
    - unités : tout est en mm dans les npz ; ParaView affichera "mm" partout.
    - les fibres de settle descendent parfois à z ≈ -0,02 mm (biais de
      pénétration Bullet documenté) — c'est physique de solveur, pas un bug
      d'export.
    - les directions des npz sont déjà unitaires ; on n'y touche pas.
    - ce script ne dépend que de numpy : pas d'import de structure.py, il
      accepte donc aussi les vieux fichiers sans route (il n'utilise que
      centers/directions/lengths/radii/core).
"""

import argparse
import os
import sys

import numpy as np


def load_fibers(path):
    """Lit un npz FiberStack ; tolère l'absence de 'core' (vieux formats)."""
    data = np.load(path, allow_pickle=False)

    for key in data.files:
        # métadonnées scalaires : on les affiche, exactement comme load_pile
        val = data[key]
        if val.shape == ():
            print(f"[npz] {os.path.basename(path)} : {key} = {val.item()}")

    needed = ("centers", "directions", "lengths", "radii")
    missing = [k for k in needed if k not in data.files]
    if missing:
        sys.exit(f"erreur : {path} ne contient pas {missing} — npz FiberStack attendu")

    centers = np.asarray(data["centers"], dtype=np.float64)
    directions = np.asarray(data["directions"], dtype=np.float64)
    lengths = np.asarray(data["lengths"], dtype=np.float64).reshape(-1)
    radii = np.asarray(data["radii"], dtype=np.float64).reshape(-1)
    core = (np.asarray(data["core"]).astype(bool).reshape(-1)
            if "core" in data.files else np.ones(len(lengths), dtype=bool))

    n = len(lengths)
    if not (len(centers) == len(directions) == len(radii) == len(core) == n):
        sys.exit("erreur : tableaux de longueurs incohérentes dans le npz")
    return centers, directions, lengths, radii, core


def write_vtp(path, centers, directions, lengths, radii, core, scale):
    """Écrit un PolyData XML ascii : 1 ligne (2 points) par fibre.

    Données par point (dupliquées aux 2 extrémités, mêmes valeurs) :
      - radius : rayon déjà multiplié par --scale (lu par Tube/By Absolute
        Scalar)
      - core   : 0/1, pour coloration
    """
    n = len(lengths)
    half = (lengths * 0.5)[:, None]
    ends_a = centers - directions * half
    ends_b = centers + directions * half

    points = np.empty((2 * n, 3), dtype=np.float64)
    points[0::2] = ends_a
    points[1::2] = ends_b

    r_pts = np.repeat(radii * scale, 2)
    core_pts = np.repeat(core.astype(np.uint8), 2)

    offsets = 2 * np.arange(1, n + 1)
    connectivity = np.arange(2 * n)

    def fmt_block(arr, fmt):
        return " ".join(format(v, fmt) for v in arr.reshape(-1))

    with open(path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="PolyData" version="0.1" '
                'byte_order="LittleEndian">\n')
        f.write(f'<PolyData>\n')
        f.write(f'  <Piece NumberOfPoints="{2 * n}" NumberOfLines="{n}" '
                'NumberOfVerts="0" NumberOfStrips="0" NumberOfPolys="0">\n')
        f.write('    <Points>\n')
        f.write('      <DataArray type="Float32" NumberOfComponents="3" '
                f'format="ascii">\n        {fmt_block(points, ".6g")}\n'
                '      </DataArray>\n')
        f.write('    </Points>\n')
        f.write('    <Lines>\n')
        f.write('      <DataArray type="Int64" Name="connectivity" '
                f'format="ascii">\n        {fmt_block(connectivity, "d")}\n'
                '      </DataArray>\n')
        f.write('      <DataArray type="Int64" Name="offsets" '
                f'format="ascii">\n        {fmt_block(offsets, "d")}\n'
                '      </DataArray>\n')
        f.write('    </Lines>\n')
        f.write('    <PointData Scalars="radius">\n')
        f.write('      <DataArray type="Float32" Name="radius" '
                f'format="ascii">\n        {fmt_block(r_pts, ".6g")}\n'
                '      </DataArray>\n')
        f.write('      <DataArray type="UInt8" Name="core" '
                f'format="ascii">\n        {fmt_block(core_pts, "d")}\n'
                '      </DataArray>\n')
        f.write('    </PointData>\n')
        f.write('  </Piece>\n')
        f.write('</PolyData>\n')
        f.write('</VTKFile>\n')


def main():
    parser = argparse.ArgumentParser(
        description="Exporte une pile npz FiberStack en .vtp pour ParaView.")
    parser.add_argument("npz", help="fichier pile .npz")
    parser.add_argument("-o", "--output", default=None,
                        help="fichier .vtp de sortie (défaut : <npz>.vtp)")
    parser.add_argument("--scale", type=float, default=1.0,
                        help="multiplicateur de rayon, confort visuel uniquement "
                             "(défaut 1.0 = rayons vrais)")
    parser.add_argument("--core", action="store_true",
                        help="n'exporter que les fibres de la zone de mesure")
    args = parser.parse_args()

    centers, directions, lengths, radii, core = load_fibers(args.npz)
    n_all = len(lengths)

    if args.core:
        mask = core
        centers = centers[mask]
        directions = directions[mask]
        lengths = lengths[mask]
        radii = radii[mask]
        core = core[mask]

    n = len(lengths)
    out = args.output
    if out is None:
        out = os.path.splitext(args.npz)[0] + (".core.vtp" if args.core else ".vtp")

    write_vtp(out, centers, directions, lengths, radii, core, args.scale)

    r_um = radii.min() * 1e3 * args.scale
    r_um_max = radii.max() * 1e3 * args.scale
    print(f"=== {n} fibres exportées ({n_all} dans le npz"
          f"{' , --core' if args.core else ''}) -> {out} ===")
    print(f"    rayons tubes : {r_um:.1f}..{r_um_max:.1f} µm "
          f"(échelle ×{args.scale:g}) — Tube/By Absolute Scalar/array 'radius'")


if __name__ == "__main__":
    main()
