# FiberStack

Génération de la microstructure 3D d'un lit de fibres de cellulose par dépôt séquentiel simulé (PyBullet), pour reproduire la structure d'un filtre de type Whatman C569 (\~19 % de fraction solide). Objectif final : porosité locale, distribution des orientations, et à terme simulation d'écoulement à travers le milieu.

## Modèle physique

- Unités : (mm, mg, ms) — g = 9,81e-3 mm/ms², zone de confort float32 de Bullet.
- Fibres : capsules rigides, diamètre lognormale tronquée 15-25 µm, longueur 1 mm ± 15 %, densité 1,6 g/cm³ (masse \~1e-7 mg).
- Dépôt : vitesse terminale 0,2 m/s, axe horizontal initial (azimut uniforme) — une fibre tombe à plat, pas bout en premier.
- Adhésion : frottements latéral/roulement/rotation élevés — à 10 µm, van der Waals ≫ poids, un contact est un attachement. `frictionAnchor=False` obligatoire (impulsions parasites sinon).
- Fibre posée = gelée (masse 0) : dépôt séquentiel, le tas n'évolue jamais rétroactivement.
- Domaine 10×10 mm, zone de mesure centrale 7×7 mm, marge tampon 1,5 mm (approximation des conditions périodiques).

## Méthode de dépôt (settle.py)

1. `place_fiber` : 51 raycasts verticaux sondent l'ombre de la fibre sur le tas → position de départ à 50 µm de l'appui le plus haut (élimine la chute libre, coût O(N²) → quasi constant par fibre).
2. Chute courte (≤ 800 pas) puis micro-relaxation : repos = immobilité &lt; 0,5 µm sur 5 contrôles de 100 pas ; garde d'explosion (v &gt; 2,5× la vitesse terminale après contact).
3. Cas pathologique (éjection, explosion, chute sous plancher, non-repos, aucun contact) : retrait du corps, redépôt une fois, abandon au deuxième échec — jamais de fibre gelée en cours de bascule.
4. Export `.npz` avec `code_version`, `seed`, `n_total` ; reprise par `resume_pile` avec contrôle de génération.

## Campagne (settle\_tiles.py)

Découpage 2×2 en tuiles de 5 mm, deux vagues en damier (`multiprocessing`) : chaque tuile est un monde PyBullet indépendant, la marge des tuiles voisines est reconstruite en corps gelés. Approximation de couture aux frontières internes, validée par comparaison aux runs séquentiels. Usage :

```bash
python settle_tiles.py 1000 0    # validation vs runs séquentiels
python settle_tiles.py 30000 0   # campagne C569
```

## Installation et usage

```bash
pip install pybullet scipy matplotlib numpy
python settle.py        # constantes N_FIBERS / SEED en tête de fichier
python analyse.py pile_1000_0.npz
```

`analyse.py` sort : porosité globale et par tranches de 1/20e de hauteur (zone de mesure), élévation des axes par rapport à l'horizontale (médiane, queue 90 %), occupation surfacique du plan (complément du « voile » vu de dessus).

## Protocole de validation (obligatoire avant tout commit)

1. Banc 300 fibres, seed 0 — attendu : hauteur \~0,70 mm, 0 abandonnée, invariant `déposées + perdues + abandonnées = N`. Référence : `pile_300_0.npz`.
2. Run 1000 fibres, comparé aux références `pile_1000_0.npz` / `pile_1000_1.npz` : élévation médiane \~0,31°, queue 3-6°, occupation \~42-43 %, porosité \~99,8 %.
3. Toute modification de constante (temps, frottements, critère de repos) passe par le banc avant d'entrer en production. Deux runs historiques (check\_every=25 ; 2000 fibres non validées) ont coûté 27 h de calcul en ignorant ce protocole.

## Biais et limites connus (documentés, pas des bugs)

- Pénétration résiduelle \~20-23 µm entre fibres en contact : profondeur de repos du solveur sur masses minuscules. Constante sur 6 runs indépendants, \~1 rayon de fibre, &lt; 0,1 % du volume.
- Gel sans réarrangement : une fibre perchée sur un seul appui ne se réencastrera jamais (dans la réalité, les suivantes la pousseraient). Audit de support dans `analyse.py` pour quantifier ; l'arche entre deux appuis est par ailleurs une structure réelle de filtre.
- Reproductibilité : même tas exige même seed ET même n\_fibers ET même version de code (le flux RNG diverge sinon) — d'où `seed` et `n_total` dans les exports.

## Fichiers


| Fichier           | Rôle                                                    |
| ----------------- | ------------------------------------------------------- |
| `settle.py`       | Dépôt séquentiel, monde unique — référence validée      |
| `settle_tiles.py` | Campagne par tuules, multiprocessing                    |
| `analyse.py`      | Porosité, orientations, occupation surfacique           |
| `Memo.md`         | État du projet, diagnostics établis, prochaines étapes  |
| `pile_*.npz`      | Tas de référence (calibration des tests non-régression) |


## Références

- Filtres de cellulose : Whatman C569 (grade 40 homologue), \~0,3 g/cm³.
- Historique des essais aérosols (essais A-D, SPM) : G. Hoarau, thèse 2020, tableaux 5 et 25-28.
