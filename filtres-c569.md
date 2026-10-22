# Filtres C569 — mémo d'état du projet

## Objectif

Générer une structure 3D réaliste d'un lit de fibres de cellulose (filtre type C569) par dépôt séquentiel, pour en extraire porosité locale et, à terme, simuler l'écoulement. Référence : Whatman C569, ~0,3 g/cm³ (fraction solide ~19 %).

## Physique et conventions

- Unités (mm, mg, ms) : g = 9,81e-3 mm/ms², forces en µN. Zone de confort float32 de Bullet.
- Fibres : capsules, diamètre lognormal 15-25 µm, longueur 1 mm ± 15 %, densité 1,6 mg/mm³, masse ~1e-7 mg.
- Dépôt : axe horizontal à la génération (azimut uniforme), vitesse d'arrivée 0,2 mm/ms (vitesse terminale), frottements élevés (contact = collage, approximation des forces van der Waals).
- Domaine 10×10 mm, zone de mesure centrale 7×7 mm, marge tampon 1,5 mm (approximation conditions périodiques).

## Historique : route A, settle_si4 (dépôt séquentiel par raycast)

1. `place_fiber` : 51 raycasts verticaux sondent l'ombre de la fibre depuis le sommet du tas → position de départ à 50 µm de l'appui le plus haut. Élimine la chute libre.
2. Chute courte (≤ 800 pas) puis micro-relaxation dynamique (repos = immobilité < 0,5 µm sur 5 contrôles de 100 pas).
3. Gel (masse 0), puis fibre suivante.
4. Gardes d'intégrité : éjection (z > départ), chute sous plancher, explosion (v > 0,5 mm/ms après contact).
5. Export .npz avec `code_version=4` ; reprise par `resume_pile` (contrôle de version au chargement).

## Diagnostics établis (à ne pas re-dérouler)

- **Pénétration résiduelle ~15-25 µm systémique** : profondeur de repos du solveur Bullet sur masses minuscules. Constante sur 3 runs indépendants. Biais connu, pas un bug. Éliminable en post-traitement si nécessaire.
- **frictionAnchor=True = éjections en cascade** : impulsions parasites sur contacts empilés. Toujours False.
- **Biais vertical initial (axe tombant bout en premier) = forêt de fibres verticales** : artefact corrigé par l'axe horizontal initial.
- **Coût de la chute libre en O(N²)** : si3 terminait à 40 s/fibre (projeté 200 jours pour 30 000). Le raycast (si4) donne ~160 ms/fibre constant.
- **Sondage de l'ombre insuffisant = éjections** : 9 sondes laissaient ~80 % de chances de rater un support perpendiculaire de 20 µm → fibre créée en chevauchement → catapultée. 51 sondes = espacement 20 µm.
- **Un tas ne se reprend qu'avec du code de même génération** (`code_version` dans le npz).

## Résultats de référence

- si3 (chute complète, 1000 fibres, seed 0) : 0 perdue, 0 non-repos, 494 core, hauteur 1,95 mm, 4h15. Porosité globale ~99,85 % (attendu : clairsemé à 1000 fibres).
- Élévation médiane ~0,3°, queue 90 % ~6° — possiblement trop plat vs C569 réel (5-20° typiques). À trancher avant campagne.
- si4 (raycast, 100 fibres) : 1 éjection (bug sondage, corrigé), 162 ms/fibre plat.

## Méthode actuelle : presse route B (press.py v2.4, pressMT v1.1)

- press.py v2.4 : suspension homogène (grille brassée + jitter), pré-relaxation 1500 pas en apesanteur, plaque par incréments de 25 µm × 200 pas, rampe de frottement 0,15 → 1,0 (lubrifié → attachement, la physique de presse humide), consolidation sous gravité, export à l'état reposé. PLATE_FRICTION = 0,05 constante (plaque glissante, v3).
- pressMT v1.1 : même physique par tuiles indépendantes en multi-processus (spawn), tuile 5×5×1,5 mm, bande de rebut burn-in 1,25 mm par côté, core = zone utile 2,5×2,5 (25 % des fibres attendu), fenêtre de silence 2000 pas (une chute libre partant du repos franchit 0,01 mm/ms en ~510 pas : 500 était escroquable).
- Frottements : Bullet multiplie les coefficients par paire — les défauts 0,15/1,0 sont la bonne physique ; ne jamais écraser par CLI (0,02/0,05 = glace, cause de l'exode du run 10000_2 : 915 core au lieu de ~4900).
- Narcolepsie Bullet : seuil de sommeil 0,4 en unités brutes (0,4 mm/ms) — un corps tranquille 1000 pas s'endort et IGNORE la gravité → 24-26 % de flottants, cause du tas hanté. sleepThreshold=0 par fibre → 0-3/250 en production. Casse réglée ; reste à l'appliquer à settle.py et regénérer les baselines route A.
- L'ensemble des cœurs de tuiles est un échantillonnage de RVE : porosité, profils, élévation, voxels = légitimes ; percolation et chemins d'écoulement continus = interdits (coutures sans fibres traversantes).

## Résultats de référence — presse (22/10/2026, Koios, repo 6a48e47)

- Monolithe 1000 seed 0 (box 8,0, 309 incréments) : 85 min, VERT, 477 core (47,7 %), hauteur 0,221 mm, φ ≈ 1,1-1,4 %, élévation médiane 1,45° / q90 4,04°, pénétration max 11,4 µm. Un mat à ~1 % de solide tient ~0,22 mm : la cible C569 (0,3 mm à 19 %) est une question de chargement (4750 fibres par tuile 5×5, ou 19000 en monolithe 10×10), pas de presse.
- pressMT 1000 seed 2 (4×250, box 1,5) : 343 s, 4/4 VERT, core 62-64 (25,1 %), z_max cœur 0,130-0,161 mm, φ_tuile ≈ 2,1 %. Sorties de silence à 1951-2451 pas = fenêtre 2000 vérifiée.
- Écart ouvert cœur-tuile vs monolithe : à chargement surfacique égal (10 fibres/mm²), la tuile finit ×1,5-2 plus dense. Candidats : densité volumique du brouillard (6,7 fibres/mm³ à box 1,5 vs 1,25 à box 8), nombre d'incréments (49 vs 309), sanglage par murs (±2,5 vs ±5,63). Discrimination : monolithe seed 2 box 1,5 (même brouillard, 49 incréments) puis compare_tile.py ; si l'écart persiste, balayage tile_side 5 → 7,5 pour mesurer la longueur de sanglage.
- Hyperion : anomalie L ≈ 0,5 mm vue sur les tuiles analysées là-bas. Le code Koios est propre partout (sondage sample_lengths 0,70-1,30 ; monolithe et tuiles VERT). Clone Hyperion à resynchroniser (git pull), re-sonder, re-valider ; ses npz en quarantaine jusqu'ici. Collision de noms : pressed_random_1000_0.npz existe sur les deux machines avec des physiques différentes (glace 0,02/0,05 côté Hyperion, h 0,047 = témoin négatif ; propre 0,15/1,0 côté Koios, h 0,221) — les npz restent intradatables tant que le nom ne code pas la physique.

## Incidents du 23/10/2026 — le mur de la phase dense

- Jamming de placement : pressMT 10000 (4×2500, box 1,5) = 67 fibres/mm³ contre 6,7-12,5 validées → 358 échecs (14,3 %) sur une tuile, garde déclenchée. Cause : v1.1 fixait la cellule de grille et ne retentait que z — au-delà de ~15 fibres/mm³ la cellule est condamnée à tous les z. pressMT v1.2 : re-tirage complet (cellule + xy + azimut + z), --place-only pour sonder la faisabilité en secondes, densité de brouillard affichée.
- Explosion de presse : monolithe 10000 seed 2 box 8,0 — 45,6 h (2,66 s/pas ; scaling ~n^1,6 : ×40 temps pour ×10 fibres), puis éclat en fin de presse : vmax 0,785 à incr 300/309, vmax 16,07 mm/ms à incr 309, sommet 34 µm AU-DESSUS de la plaque (= un pas à 16 mm/ms : fibre catapultée à travers la plaque). φ n'était que ~4-10 % à l'éclat : ce n'est pas la densité globale, c'est un bourrage local (fibre coincée + friction ~1,0 atteinte par la rampe + plaque qui pousse). À 1000 fibres, vmax est resté à 0,07 de bout en bout.
- Lectures : (1) le bloqueur de la campagne C569 (190 fibres/mm², φ 19 % à 0,3 mm) n'est ni la parallélisation ni le placement — c'est la stabilité de la presse en phase dense ; (2) les deux runs du jour étaient la MÊME densité aréique (100/mm², φ 10 %) — le monolithe montre ce que les tuiles auraient trouvé en fin de presse ; (3) npz du monolithe explosé : quarantaine, jamais une baseline (validate ne voit pas l'histoire — un mat ressayé après un éclat à 16 mm/ms peut être VERT et rester compromis).
- Fix press.py v2.5 à livrer : presse multi-nips (presser par étapes avec consolidation intermédiaire — une vraie calandre presse en nips multiples), incréments fins (5 µm) sur le dernier demi-mm, garde DÉRIVE durcie (pause jusqu'au vrai silence, pas 400 pas fixes ; refus d'export si éclat non résorbé).

## Prochaines étapes

1. Resync Hyperion : git pull dans ~/Git/fiberstack, sondage sample_lengths, validate de ses npz — quarantaine maintenue jusqu'à re-validation.
2. Certification tuile/monolithe : press.py --random 1000 --seed 2 --target 0.3 --box-z 1.5 (~30 min), puis compare_tile.py (cœur de tuile vs zone centrale 2,5×2,5 du monolithe). Rapport φ ≈ 1 = certifié.
3. press.py v2.5 — stabilité phase dense (multi-nips, incréments fins finaux, garde DÉRIVE durcie) : condition d'entrée de toute campagne ; le fichier v2.4 déployé (commit e5b8fee) doit servir de base.
4. Sonde de placement : pressmt.py --random 76000 --nx 4 --ny 4 --box-z 2/3/4 --place-only (4750/tuile, retenir le plus petit box-z sous 1 % d'échecs), puis une tuile 4750 complète pilote — exigence : zéro DÉRIVE.
5. Campagne : 16 tuiles × 4750 fibres (une tuile = un C569 complet), workers 16, une à deux nuits de Koios ; statistiques d'ensemble (3 seeds × 4 tuiles).
6. Vérifier la fenêtre de silence interne de press.py (« repos atteint à 6000 pas » = plafond du budget ? grep -n "quiet\\|QUIET" press.py) ; appliquer le fix sleep à settle.py et regénérer les baselines route A (pile_300_0, pile_1000_1).
7. Analyse voxel (grille 3D solide/fluide, rayon de pore par distance transform) : porosité hydrodynamique comparable au C569 ; comparaison micrographie si disponible.

## Dossier Git

- Repo propre et poussé (origin/main = HEAD 6a48e47) : 1d895c5 gitignore pycache, e5b8fee press v2.4 + pressMT, 6a48e47 pressMT. Les scripts patchés (validate bornes L ±2σ, analyse zone utile tuile, pressmt v1.1) sont déployés ET committés — la règle commit-avant-run est enfin tenue.
- Reste : committer les mouvements de données (deletions + dossiers non suivis) ; renommer « 1000 0/ », « 300 0/ », « 1000 1 1/ » (les espaces cassent les scripts un jour ou l'autre) vers old/ ; gitignore les sorties (*.npz, *.vtp) — le code se versionne, les données se comptent ; pile_1000_1.npz modifié : committer délibérément si c'est la régénération settle v6, sinon old/.