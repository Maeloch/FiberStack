Filtres C569 — mémo d'état du projet

Objectif

Générer une structure 3D réaliste d'un lit de fibres de cellulose (filtre type C569) par dépôt séquentiel, pour en extraire porosité locale et, à terme, simuler l'écoulement. Référence : Whatman C569, ~0,3 g/cm³ (fraction solide ~19 %).

Physique et conventions

Unités (mm, mg, ms) : g = 9,81e-3 mm/ms², forces en µN. Zone de confort float32 de Bullet.

Fibres : capsules, diamètre lognormal 15-25 µm, longueur 1 mm ± 15 %, densité 1,6 mg/mm³, masse ~1e-7 mg.

Dépôt : axe horizontal à la génération (azimut uniforme), vitesse d'arrivée 0,2 mm/ms (vitesse terminale), frottements élevés (contact = collage, approximation des forces van der Waals).

Domaine 10×10 mm, zone de mesure centrale 7×7 mm, marge tampon 1,5 mm (approximation conditions périodiques).

Historique : route A, settle_si4 (dépôt séquentiel par raycast)

place_fiber : 51 raycasts verticaux sondent l'ombre de la fibre depuis le sommet du tas → position de départ à 50 µm de l'appui le plus haut. Élimine la chute libre.

Chute courte (≤ 800 pas) puis micro-relaxation dynamique (repos = immobilité < 0,5 µm sur 5 contrôles de 100 pas).

Gel (masse 0), puis fibre suivante.

Gardes d'intégrité : éjection (z > départ), chute sous plancher, explosion (v > 0,5 mm/ms après contact).

Export .npz avec code_version=4 ; reprise par resume_pile (contrôle de version au chargement).

Diagnostics établis (à ne pas re-dérouler)

Pénétration résiduelle ~15-25 µm systémique : profondeur de repos du solveur Bullet sur masses minuscules. Constante sur 3 runs indépendants. Biais connu, pas un bug. Éliminable en post-traitement si nécessaire.

frictionAnchor=True = éjections en cascade : impulsions parasites sur contacts empilés. Toujours False.

Biais vertical initial (axe tombant bout en premier) = forêt de fibres verticales : artefact corrigé par l'axe horizontal initial.

Coût de la chute libre en O(N²) : si3 terminait à 40 s/fibre (projeté 200 jours pour 30 000). Le raycast (si4) donne ~160 ms/fibre constant.

Sondage de l'ombre insuffisant = éjections : 9 sondes laissaient ~80 % de chances de rater un support perpendiculaire de 20 µm → fibre créée en chevauchement → catapultée. 51 sondes = espacement 20 µm.

Un tas ne se reprend qu'avec du code de même génération (code_version dans le npz).

Résultats de référence

si3 (chute complète, 1000 fibres, seed 0) : 0 perdue, 0 non-repos, 494 core, hauteur 1,95 mm, 4h15. Porosité globale ~99,85 % (attendu : clairsemé à 1000 fibres).

Élévation médiane ~0,3°, queue 90 % ~6° — possiblement trop plat vs C569 réel (5-20° typiques). À trancher avant campagne.

si4 (raycast, 100 fibres) : 1 éjection (bug sondage, corrigé), 162 ms/fibre plat.

Méthode actuelle : presse route B (press.py v2.4, pressMT v1.1)

press.py v2.4 : suspension homogène (grille brassée + jitter), pré-relaxation 1500 pas en apesanteur, plaque par incréments de 25 µm × 200 pas, rampe de frottement 0,15 → 1,0 (lubrifié → attachement, la physique de presse humide), consolidation sous gravité, export à l'état reposé. PLATE_FRICTION = 0,05 constante (plaque glissante, v3).

pressMT v1.1 : même physique par tuiles indépendantes en multi-processus (spawn), tuile 5×5×1,5 mm, bande de rebut burn-in 1,25 mm par côté, core = zone utile 2,5×2,5 (25 % des fibres attendu), fenêtre de silence 2000 pas (une chute libre partant du repos franchit 0,01 mm/ms en ~510 pas : 500 était escroquable).

Frottements : Bullet multiplie les coefficients par paire — les défauts 0,15/1,0 sont la bonne physique ; ne jamais écraser par CLI (0,02/0,05 = glace, cause de l'exode du run 10000_2 : 915 core au lieu de ~4900).

Narcolepsie Bullet : seuil de sommeil 0,4 en unités brutes (0,4 mm/ms) — un corps tranquille 1000 pas s'endort et IGNORE la gravité → 24-26 % de flottants, cause du tas hanté. sleepThreshold=0 par fibre → 0-3/250 en production. Casse réglée ; reste à l'appliquer à settle.py et regénérer les baselines route A.

L'ensemble des cœurs de tuiles est un échantillonnage de RVE : porosité, profils, élévation, voxels = légitimes ; percolation et chemins d'écoulement continus = interdits (coutures sans fibres traversantes).

Résultats de référence — presse (22/10/2026, Koios, repo 6a48e47)

Monolithe 1000 seed 0 (box 8,0, 309 incréments) : 85 min, VERT, 477 core (47,7 %), hauteur 0,221 mm, φ ≈ 1,1-1,4 %, élévation médiane 1,45° / q90 4,04°, pénétration max 11,4 µm. Un mat à ~1 % de solide tient ~0,22 mm : la cible C569 (0,3 mm à 19 %) est une question de chargement (4750 fibres par tuile 5×5, ou 19000 en monolithe 10×10), pas de presse.

pressMT 1000 seed 2 (4×250, box 1,5) : 343 s, 4/4 VERT, core 62-64 (25,1 %), z_max cœur 0,130-0,161 mm, φ_tuile ≈ 2,1 %. Sorties de silence à 1951-2451 pas = fenêtre 2000 vérifiée.

Certification tuile/monolithe RÉUSSIE (23/10 soir, Hyperion) : monolithe 1000 seed 2 box 1,5 (25 min, VERT, 482 core, h 0,266, φ zone 1,78 %) vs cœurs de tuiles même seed : φ 2,03 ± 0,30 % → rapport 1,14, monolithe à ~1,7σ de la moyenne tuiles (non significatif à N=4) ; densité aréique en zone 2,5×2,5 : 63 fibres (mono) vs 62-63 (tuiles) — identique au fibre près ; élévations (1,0-1,8° médiane) et profils en profondeur relative superposés ; occupation surfacique 42-45 % partout. L'écart ×1,5-2 craint le 22/10 était un artefact de comparaisons confondues (box 8 vs 1,5, zone 7×7 vs cœur, aire de cage étalée). Route tuiles certifiée à l'échelle 1000 ; le vrai test du sanglage par murs sera la tuile pilote 4750 (dispersion inter-tuiles ~1/√19 plus resserrée).

RECERTIFICATION v2.5 (24/10, Hyperion, tout en 7cdf4e6) : monolithe 1000 seed 2 box 1,5 relancé en v2.5 (30m44s, VERT, 482 core, z_max core 0,293 mm, repos à 7000 pas / 4 contrôles, vmax 0,0012) + 4 tuiles v1.3 (445 s, 4/4 VERT, 0 échec de placement à 6,7 fibres/mm³, cores 62-63, pénétrations 7,5-9,6 µm, élévations médianes 1,2-2,1°). compare_tile : φ mono 1,79 vs tuiles 1,57±0,31 → rapport 0,88 (le v2.4 donnait 1,14 : le signe de l'écart a changé, les deux à <1σ à N=4 — non significatif, et φ reste une statistique dominée par z_max, un max sur une seule fibre : t01 z_max 0,148 → φ 1,97). L'INVARIANT FORT est ailleurs et identique : densité aréique 63 fibres (mono, zone 2,5) vs 62-63 (cœurs), profils en profondeur relative superposés, occupations surfaciques 42-45 %. Certification tenue, cette fois avec les DEUX côtés en v2.5. Coût mesuré : monolithe 1000 = 30 min ; extrapolation tuile 4750 (scaling n^1,6) ≈ 14 h/tuile → 16 tuiles en parallèle sur Koios ≈ une nuit, conforme au plan. NB : une tuile a fini sommet 0,256 sous plaque 0,300 (mat décroché de la plaque en fin de course — bénin à ~1,5 % de solide, à surveiller au pilote 4750) ; en-têtes de colonnes du profil compare_tile collés (cosmétique, sans urgence).

Physique des profils (commune aux deux routes) : base dense (φ local 3-7 % dans le premier dixième) puis moutonnement vers 100 % en surface — la « hauteur » est portée par les fibres dressées (z_max = max, pas médiane). Profil en beignet : le monolithe a un z_max central 0,177 pour 0,266 en zone ±3,5 — la périphérie est plus épaisse que le centre (extrusion latérale douce vers la cage, même signature que l'exode glace, en mode bénin) ; les tuiles montrent la même asymétrie cœur/tuile. Pour toute comparaison µCT : échantillonner la zone centrale, jamais le bord.

Déterminisme pressmt vérifié (bonus du 23/10) : deux runs même seed/même v1.2 sur Hyperion → fichiers identiques au dernier chiffre (φ, h, core). Reproductibilité bit-à-bit. Les tuiles Koios du 22/10 différaient car générées en v1.1 (placement à cellule fixe).

Hyperion resynchronisé et propre (auth réparée, pull, v1.2 + validate/analyse patchés déployés) : run 1000 seed 2 → 4/4 VERT, L 0,705-1,299, 0 échec de placement à 6,7 fibres/mm³. Anomalie L≈0,5 : clone fossile, close. Quarantaine levée pour les sorties post-resync ; anciens npz (monolithe glace, tuiles L≈0,5) dans old/ comme pièces à conviction. NB : le repo GitHub est PUBLIC (lecture anonyme = pourquoi le pull marchait) — visibilité à confirmer comme voulue.

Incidents du 23/10/2026 — le mur de la phase dense

Jamming de placement : pressMT 10000 (4×2500, box 1,5) = 67 fibres/mm³ contre 6,7-12,5 validées → 358 échecs (14,3 %) sur une tuile, garde déclenchée. Cause : v1.1 fixait la cellule de grille et ne retentait que z — au-delà de ~15 fibres/mm³ la cellule est condamnée à tous les z. pressMT v1.2 : re-tirage complet (cellule + xy + azimut + z), --place-only pour sonder la faisabilité en secondes, densité de brouillard affichée.

Explosion de presse : monolithe 10000 seed 2 box 8,0 — 45,6 h (2,66 s/pas ; scaling ~n^1,6 : ×40 temps pour ×10 fibres), puis éclat en fin de presse : vmax 0,785 à incr 300/309, vmax 16,07 mm/ms à incr 309, sommet 34 µm AU-DESSUS de la plaque (= un pas à 16 mm/ms : fibre catapultée à travers la plaque). φ n'était que ~4-10 % à l'éclat : ce n'est pas la densité globale, c'est un bourrage local (fibre coincée + friction ~1,0 atteinte par la rampe + plaque qui pousse). À 1000 fibres, vmax est resté à 0,07 de bout en bout.

Lectures : (1) le bloqueur de la campagne C569 (190 fibres/mm², φ 19 % à 0,3 mm) n'est ni la parallélisation ni le placement — c'est la stabilité de la presse en phase dense ; (2) les deux runs du jour étaient la MÊME densité aréique (100/mm², φ 10 %) — le monolithe montre ce que les tuiles auraient trouvé en fin de presse ; (3) npz du monolithe explosé : quarantaine, jamais une baseline (validate ne voit pas l'histoire — un mat ressayé après un éclat à 16 mm/ms peut être VERT et rester compromis).

Fix press.py v2.5 LIVRÉ (canvas, 23/10) : presse multi-nips (presser par étapes avec consolidation intermédiaire — une vraie calandre presse en nips multiples), incréments fins (5 µm) sur le dernier 0,1 mm, garde DÉRIVE durcie (pauses jusqu'au vrai silence, comptées par blocs, plafond 25 ; refus d'export si éclat non résorbé), anti-narcolepsie et silence 4-contrôles en prime. À déployer et committer avant tout run de campagne.

Prochaines étapes

FAIT (23/10) : resync Hyperion + certification tuile/monolithe (rapport 1,14) + route D (routed.py, contrôle négatif booléen).

press.py v2.5 LIVRÉ, COMMITTÉ et DÉPLOYÉ (7cdf4e6, 23/10 soir, vérifié patch par patch sur GitHub ; coquille dans le message : « press 2.4 » — le contenu est v2.5) — anti-narcolepsie (sleepThreshold=0 par fibre), silence = 4 contrôles consécutifs sous 0,01 mm/ms espacés de 500 pas, garde DÉRIVE par blocs re-contrôlés + porte de fin de presse (refus d'export d'un éclat), --nips multi-nips (setpoints géométriques, consolidation inter-nips), incréments fins 5 µm sous 0,1 mm de la cible, CODE_VERSION=2.5, nips dans les métadonnées. pressMT v1.3 avec (plomberie --nips ; le re-tirage v1.2 rattrapé au passage). Hyperion : à jour, monolithe + 4 tuiles relancés en v2.5 le 24/10 — RECERTIFICATION RÉUSSIE (voir Résultats de référence). Koios : pull, puis grep -n "CODE_VERSION" press.py (attendu 2.5) avant toute campagne. RESTE À PASSER : test-2 (flottants) sur pressed_random_1000_2.npz v2.5 (attendu ~0 — c'est LE test de l'anti-narcolepsie sur la sortie) ; les anciens npz v2.4 n'ont plus de rôle de baseline. Vérifier aussi settle.py : v6.1 COMMITTÉ (7cdf4e6) — sleepThreshold=0 à la création de chaque fibre + CODE_VERSION 6→7 (la reprise refuse donc tout npz v6 : les baselines route A pile_300_0 et pile_1000_1 sont à régénérer, elles contiennent des dormeurs potentiels) ; à puller sur Koios demain. Divergence attendue vs v6 : uniquement là où v6 laissait une fibre s'endormir (pas_stab plus longs, moins de flottants) ; toute autre divergence est un bug.

Sonde de placement FAITE (24/10, Hyperion) : échecs 24,9 % (box 2, 95 f/mm³) / 10,4 % (box 3, 63,3) / 3,4 % (box 4, 47,5) — décroissance ~÷2,7 par mm de box-z : extrapolation box 5 ≈ 1,1 %, box 6 ≈ 0,4 %. La porte <1 % exigerait donc box-z 6, soit un brouillard à φ ≈ 1 % — c'est le mur du placement séquentiel SANS chevauchement (RSA) pour des fibres d'élancement 50, cohérent avec la littérature empilements aléatoires de bâtonnets (fraction décroissante avec l'élancement). PAS un bug. CONTREPARTIE : la plaque balaie box-z − 0,3 mm → box 6 = 244 incréments (×3,8 vs box 1,5) ; à 4750 fibres et coût du pas ~n^1,6, l'estimation « une nuit de Koios » (calée sur box 1,5) est morte : box 4 ≈ ×2,6 → ~25 h/tuile, box 6 ≈ ×3,8 → ~40 h/tuile. CONTESTATION DE LA PORTE : le critère <1 % date de la pathologie v1.1 (cellules condamnées) ; avec le re-tirage v1.2 un échec est du jamming RSA pur — fibres manquantes aléatoires, brouillard uniforme, zéro chevauchement garanti ; et à box 4 le brouillard (φ 1,5 %) est exactement la densité d'une vraie suspension de pâte humide (0,5-2 %). FOURCHE À TRANCHER : (a) box 6, zéro échec, zéro code, ~2 jours de Koios ; (b) box 4 + compensation --random ~4920/tuile pour placer 4750 (garde à élargir de 1 % à ~5 %, dispersion de chargement ±0,4 % par tuile, pilote à recalibrer la densité du cœur) ; (c) RECOMMANDÉ : tolérance de chevauchement bornée au placement (≤ ~0,3 rayon au lieu de CLEAR=+10 µm), résolue par la pré-relaxation — box 2-3 suffit, presse ×1,3-1,9, la nuit de Koios revient, et la configuration relaxée est plus proche d'une suspension brownienne que le RSA ; MAIS press/pressmt v1.4 + commit + pilote de validation AVANT campagne (risque : budget de pré-relaxation 1500 pas à réviser, pénétrations initiales à surveiller via validate). DÉCISION (24/10 soir) : (a) RETENU pour la campagne — pilote box 6 puis 16 tuiles box 6 sur Koios, zéro code, ~12-16 h/tuile (244 incréments à ~1 s/pas à 4750) ; (b) ÉCARTÉ (aucun avantage) ; (c) v1.4 LIVRÉ (canvas press v2.6 + pressmt v1.4) : --overlap-tol borné par OVERLAP_TOL_MAX = 6 µm (moitié du plus petit rayon), recommandé 0,003 mm (~0,3 rayon), défaut 0 = régression exacte v2.5/v1.3, garde vmax > 0,01 mm/ms après pré-relaxation si tol > 0, overlap_tol stocké dans le npz. HONNÊTETÉ sur l'argument « fibres aplatties aux croisements » (retenu par l'utilisateur pour l'interprétation) : le mécanisme du tol est la résolution par la pré-relaxation, pas une déformation simulée ; la fiction d'aplatissement vit dans la pénétration résiduelle 15-25 µm du solveur, documentée depuis le début. Et une configuration relaxée est plus proche d'une suspension brownienne que le RSA (corrélatiions parasites). CALIBRAGE v1.4 FAIT (3/11, Koios, --place-only 16×4750, tol 3 µm) : box 2 : 20,3 % d'échecs en moyenne (24,9 % à tol 0), box 3 : 7,0 % (10,4 %). VERDICT : le tol fait exactement son travail mais reste MARGINAL — 3 µm = −10 % de volume exclu (distance d'exclusion ~30 → 27 µm) = l'équivalent de ~0,2-0,3 mm de box-z. Le jamming des bâtonnets est gouverné par le volume exclu ~L²d : dominé par la LONGUEUR (1 mm), pas par le rayon (10 µm) — une concession de 3 µm ne peut pas défaire un problème d'élancement 50. Rendre box 2-3 viables exigerait 10-15 µm de tol = interpénétrations profondes que la pré-relaxation ne résout pas proprement. Le pronostic initial « box 2-3 redevient viable avec le tol » est RETIRÉ (erreur du 24/10). Projections à tol 3 µm : box 4 ≈ 2 % (porte échouée), box 5 ≈ 0,6 % (porte passée, 40 incréments économisés vs box 6), box 6 ≈ 0,2 %. DÉCISION MAINTENUE : campagne box 6 tol 0 — provenance la plus propre pour la métrologie (zéro chevauchement initial, brouillard dilué réaliste) ; le tol reste un outil marginal pour un box-z contrait. L'A/B tol-vs-RSA n'a plus d'enjeu de campagne ; à ne rouvrir que si la question « le RSA biaise-t-il les statistiques du mat final ? » devient pertinente (pilote box 5 + tol vs box 6 tol 0, même seed, à nombre placé égal). NB : MAX_PLACE_TRIES = 20 dans press.py (partagé routes A/B) — le monter à 50-60 réduit la marge d'échecs à box constant.

Tuile pilote 4750 complète, exigence zéro DÉRIVE ; y vérifier φ cœur ≈ 19 % à 0,3 mm et le sanglage par murs (la certification 1000 ne préjuge pas de la densité C569).

Campagne : 16 tuiles × 4750 fibres, workers 16, une à deux nuits de Koios ; statistiques d'ensemble (3 seeds × 4 tuiles).

Fenêtre de silence press.py : VÉRIFIÉE en v2.5 (repos monolithe à 7000 pas, tuiles 4000-4500, puis quiet pressmt 1951). Reste : régénérer les baselines route A en settle v6.1 (pile_300_0, pile_1000_1 — leurs npz v6 sont refusés par le contrôle strict).

Analyse voxel (grille 3D solide/fluide, rayon de pore par distance transform) puis comparaison routed vs press à φ égal — que vaut la physique ? ; comparaison micrographie si disponible.

Dossier Git

Git au 3/11 (vérifié marqueur par marqueur) : origin/main = f78d3b6 « maj » (2/11, 23:01) = déploiement press v2.6 + pressMT v1.4 (--overlap-tol, port fidèle des canvases : OVERLAP_TOL_MAX, tol dans _touches_any, garde pré-relax, plomberie pressmt, tol défaut 0 = régression exacte) ; 2aa480e (24/10) = Memo.md (l'utilisateur tient désormais son mémo dans le repo) + validate/analyse/structure ; 7cdf4e6 (23/10) = settle v6.1 + press v2.5 + pressMT v1.3. Koios resynchronisé le 3/11 par git reset --hard origin/main après conflit de pull (il était resté à 6a48e47, trois commits de retard) : les modifications locales écrasées (pressmt.py patch main des 22/10 + pile_1000_1.npz/png) sont sans objet — le v1.4 du repo contient strictement plus que le patch local, et la baseline v6 est refusée par le settle v6.1 de toute façon. Enseignement : conflit de pull = git stash d'abord, jamais reset --hard en premier réflexe. Piège historique à ne pas refer : c1b9635 « press v2.5, pressMT v1.3 » était VIDE (fichiers édités non stagés) — rituel obligatoire : git show --stat HEAD AVANT tout push. ÉCART À CORRIGER dans press.py : la route B n'écrit pas overlap_tol dans le npz (pressMT, lui, l'écrit) — un monolithe généré avec tol > 0 ne serait pas traçable dans son fichier ; ajouter overlap_tol=(args.overlap_tol if args.random is not None else None) au save_pile et committer.

Reste : committer les mouvements de données (deletions + dossiers non suivis) ; renommer « 1000 0/ », « 300 0/ », « 1000 1 1/ » (les espaces cassent les scripts un jour ou l'autre) vers old/ ; gitignore les sorties (*.npz, *.vtp) — le code se versionne, les données se comptent ; pile_1000_1.npz modifié : committer délibérément si c'est la régénération settle v6, sinon old/.
