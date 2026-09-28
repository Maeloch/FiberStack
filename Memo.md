# Filtres C569 — mémo d'état du projet

## Objectif

Générer une structure 3D réaliste d'un lit de fibres de cellulose (filtre type C569) par dépôt séquentiel, pour en extraire porosité locale et, à terme, simuler l'écoulement. Référence : Whatman C569, \~0,3 g/cm³ (fraction solide \~19 %).

## Physique et conventions

- Unités (mm, mg, ms) : g = 9,81e-3 mm/ms², forces en µN. Zone de confort float32 de Bullet.
- Fibres : capsules, diamètre lognormal 15-25 µm, longueur 1 mm ± 15 %, densité 1,6 mg/mm³, masse \~1e-7 mg.
- Dépôt : axe horizontal à la génération (azimut uniforme), vitesse d'arrivée 0,2 mm/ms (vitesse terminale), frottements élevés (contact = collage, approximation des forces van der Waals).
- Domaine 10×10 mm, zone de mesure centrale 7×7 mm, marge tampon 1,5 mm (approximation conditions périodiques).

## Méthode actuelle : settle\_si4 (dépôt par raycast)

1. `place_fiber` : 51 raycasts verticaux sondent l'ombre de la fibre depuis le sommet du tas → position de départ à 50 µm de l'appui le plus haut. Élimine la chute libre.
2. Chute courte (≤ 800 pas) puis micro-relaxation dynamique (repos = immobilité &lt; 0,5 µm sur 5 contrôles de 100 pas).
3. Gel (masse 0), puis fibre suivante.
4. Gardes d'intégrité : éjection (z &gt; départ), chute sous plancher, explosion (v &gt; 0,5 mm/ms après contact).
5. Export .npz avec `code_version=4` ; reprise par `resume_pile` (contrôle de version au chargement).

## Diagnostics établis (à ne pas re-dérouler)

- **Pénétration résiduelle \~15-25 µm systémique** : profondeur de repos du solveur Bullet sur masses minuscules. Constante sur 3 runs indépendants. Biais connu, pas un bug. Éliminable en post-traitement si nécessaire.
- **frictionAnchor=True = éjections en cascade** : impulsions parasites sur contacts empilés. Toujours False.
- **Biais vertical initial (axe tombant bout en premier) = forêt de fibres verticales** : artefact corrigé par l'axe horizontal initial.
- **Coût de la chute libre en O(N²)** : si3 terminait à 40 s/fibre (projeté 200 jours pour 30 000). Le raycast (si4) donne \~160 ms/fibre constant.
- **Sondage de l'ombre insuffisant = éjections** : 9 sondes laissaient \~80 % de chances de rater un support perpendiculaire de 20 µm → fibre créée en chevauchement → catapultée. 51 sondes = espacement 20 µm.
- **Un tas ne se reprend qu'avec du code de même génération** (`code_version` dans le npz).

## Résultats de référence

- si3 (chute complète, 1000 fibres, seed 0) : 0 perdue, 0 non-repos, 494 core, hauteur 1,95 mm, 4h15. Porosité globale \~99,85 % (attendu : clairsemé à 1000 fibres).
- Élévation médiane \~0,3°, queue 90 % \~6° — possiblement trop plat vs C569 réel (5-20° typiques). À trancher avant campagne.
- si4 (raycast, 100 fibres) : 1 éjection (bug sondage, corrigé), 162 ms/fibre plat.

## Prochaines étapes

1. **Validation si4** : 1000 fibres seed 0, comparer à pile\_0.npz (élévation, profil de porosité, occupation surfacique). Le chrono doit rester plat.
2. **Trancher la raideur d'orientation** (dispersion initiale `direction[2] = rng.uniform(-0.1, 0.1)`) selon comparaison micrographie C569.
3. **Campagne \~30 000 fibres** par lots de 10 000 avec reprise, seeds distincts (fraction solide C569 \~19 % visée).
4. **Analyse voxel** (grille 3D solide/fluide, rayon de pore via distance transform) : porosité hydrodynamique comparable au C569.
5. Comparaison image à image avec micrographie si disponible (élévation, occupation surfacique).

## Dossier Git

À faire : `git init` dans le dossier de travail — deux pertes de fichiers déjà occasionnées par des changements de machine sans synchronisation.
