# Corrélation entre la composition du GC de l'ARN ribosomique et la température de croissance dans Archaea.

## 1. Introduction

### 1.1 Contexte

Les archées sont des micro-organismes qui vivent à des températures très variées, des milieux tempérés jusqu'à des sources chaudes (de mésophiles ~20°C à hyperthermophiles ~100°C).

Deux signaux moléculaires sont connus pour être corrélés à la température de croissance optimale (OGT) : le %GC de l'ARN ribosomique et la composition en acides aminés des protéines, ce sont des **thermomètres moléculaires**.
On peut donc, en principe, utiliser ces signaux pour reconstruire la température des ancêtres (des organismes qu'on ne peut évidemment plus observer directement).

### 1.2 L'étude de référence : Groussin et Gouy (2011)

- Leur méthode : deux étapes séparées

Les auteurs ont étudié 35 espèces d'archées dont le génome est complètement séquencé, à partir de deux sources d'information : l'ARNr et 72 familles de protéines.

1. **Reconstruire la composition en GC des ancêtres:**
Construisent une phylogénie des Archées (35 génomes), utilisent des modèles non-homogènes d'évolution moléculaire pour reconstruire les séquences ancestrales (ARNr + protéines) à chaque nœud de l'arbre.
2. **Traduire ces compositions en températures:** À partir de ces séquences ancestrales reconstruites, ils calculent le %GC (ou la composition en acides aminés) à chaque nœud.
Chez les espèces actuelles, ils mesurent la relation entre composition et OGT (r = 0,95 pour le GC de l'ARNr, r = 0,84 pour un indice de composition des protéines). Cette  droite de régression GC% ↔ OGT sert de thermomètre : ils y placent les compositions ancestrales pour en déduire la température des ancêtres.

Pour vérifier que ces corrélations ne viennent pas simplement de la parenté entre espèces, ils les recalculent avec des contrastes phylogénétiquement indépendants (PIC), qui confirment le lien.

***Resultat***: 
- Ils trouvent que l'ancêtre commun des archées était **hyperthermophile** : environ 90 °C d'après l'ARNr et 82 °C d'après les protéines avec un intervalle de confiance de 74-89 °C
- Plusieurs lignées, notamment chez les Euryarchées, se sont ensuite adaptées progressivement à des milieux plus froids.
- Les lignées des milieux tempérés évoluent plus vite (branches plus longues) : la température apparaît comme un déterminant majeur de la vitesse d'évolution moléculaire chez les archées.



## 2 Problematique

- Leur approche se fait en deux étapes séparées : (1) reconstruction du GC% ancestral par un modèle moléculaire, puis (2) régression post-hoc GC%→OGT.
- Conséquence, l'incertitude de l'étape 1 et celle de l'étape 2 ne sont jamais combinées correctement. On sous-estime probablement l'incertitude réelle sur les températures ancestrales,on obtient donc probablement des intervalles de confiance trop optimistes sur les températures ancestrales.
- De plus, leur régression GC↔OGT est calibrée sur les espèces actuelles uniquement (corrélation "statique"), alors que ce qu'on veut vraiment tester, c'est si GC et température évoluent ensemble le long des branches de l'arbre (corrélation "dynamique"/évolutive). Ce n'est pas rigoureusement la même chose (c'est d'ailleurs pour cette raison qu'ils utilisent en renfort les contrastes indépendants de Felsenstein, PIC, mais toujours dans une étape séparée).



## 3 Objectifs

Ce projet a deux objectifs, qui s'enchaînent :

1. **Modéliser l'évolution corrélée de la teneur en GC de l'ARNr et de la température de croissance à travers Archaea, et d'utiliser ce modèle pour estimer la corrélation entre le taux de GC et la température le long de l'arbre des archées.**
   Quand la température optimale de croissance (OGT) d'une lignée augmente au cours de l'évolution, son taux de GC augmente-t-il aussi ? Cette corrélation est mesurée par un paramètre du modèle, `beta`.

2. **Utiliser cette corrélation pour inférer les températures ancestrales aux nœuds internes de l'arbre, notamment à l'ancêtre commun des Archées.**
   La racine de l'arbre représente l'ancêtre commun des espèces étudiées : vivait-il dans un milieu chaud ?

> **Question de recherche :** peut-on estimer, dans un seul modèle phylogénétique bayésien, à la fois la corrélation évolutive entre GC et température et la température des ancêtres ?

Cette question a déjà été abordée par Groussin et Gouy (2011) avec une autre méthode. Nous présentons d'abord leur travail et ses limites, puis la façon dont notre modèle y répond.



## 4. **Materiel et méthodes**

***Données***

L'analyse porte sur 33 espèces d'archées. Le dossier `data/` contient :

| Fichier | Contenu |
| --- | --- |
| `data/archaea.nex` | Alignement de l'ARNr : 1 801 positions des régions en double brin, dont la composition est la plus liée à la température |
| `data/archaea.tree` | Arbre phylogénétique des 33 espèces (topologie et longueurs de branches, considérées comme fixes) |
| `data/archaea_temp_gc.nex` | Température optimale de croissance (OGT) et taux de GC de l'ARNr de chaque espèce actuelle |
| `data/rna.itgc` | Fichier auxiliaire (85 taxons), non utilisé dans l'analyse actuelle |


### 3.1 Ce que nous changeons

Nous reprenons la même question avec l'ARNr, mais en remplaçant les deux étapes par un seul modèle bayésien, implémenté dans RevBayes.

| Limite de l'approche en deux étapes | Réponse de notre modèle |
| --- | --- |
| L'incertitude du thermomètre n'est pas propagée | La corrélation (`beta`), les GC ancestraux et les températures ancestrales sont estimés ensemble : les intervalles de crédibilité intègrent toutes ces incertitudes |
| Le thermomètre traite les espèces comme indépendantes | La corrélation porte sur les changements le long des branches : la parenté entre espèces fait partie du modèle |
| La température n'évolue pas le long de l'arbre | La température évolue elle-même le long des branches : les températures actuelles renseignent directement celles des ancêtres |

Notre modèle ne règle pas en revanche la dernière limite : il n'utilise que l'ARNr (voir la section 8).

### 3.2 Comment fonctionne le modèle

Le modèle décrit ce qui se passe le long de chaque branche de l'arbre :

- **La température varie un peu, au hasard.** Plus la branche est longue, plus elle peut changer. L'ampleur de ces variations est réglée par `sigma_T`.
- **Le GC varie aussi, en partie à cause de la température :**

  ```
  changement de GC ≈ beta × changement de température + variation propre au GC
  ```

  La part que la température n'explique pas est réglée par `sigma_GC`.
- **Le GC de chaque branche oriente l'évolution des séquences.** Sur une branche riche en GC, les substitutions vers G et C sont plus fréquentes. L'alignement d'ARNr apporte donc lui aussi de l'information sur le GC des ancêtres.

Aux extrémités de l'arbre, la température et le GC des espèces actuelles sont connus. À partir de ces observations et des séquences, le modèle estime `beta` et les valeurs ancestrales, chacune avec son incertitude (distribution a posteriori).

<details>
<summary><b>Détails techniques du modèle</b></summary>

Pour la branche qui va du nœud parent `pa(i)` au nœud `i`, de longueur `bl(i)`, avec Normale(moyenne, écart-type) :

```
log T(i)    ~ Normale( log T(pa(i)),  sigma_T × √bl(i) )
logit GC(i) ~ Normale( logit GC(pa(i)) + beta × [log T(i) - log T(pa(i))],  sigma_GC × √bl(i) )
```

- La température est modélisée en logarithme (elle reste positive) et le GC en logit (il reste entre 0 et 1) : ce sont deux mouvements browniens le long de l'arbre, couplés par `beta`.
- Le GC d'une branche est la moyenne du GC à ses deux extrémités. Il définit une matrice de substitution T92 propre à la branche, avec un rapport transitions/transversions `kappa` commun à tout l'arbre. Les fréquences des bases à la racine sont déduites du GC de la racine.
- La vraisemblance de l'alignement est calculée avec `dnPhyloCTMC` sur l'arbre fixé.
- Lois a priori : `sigma_T` et `sigma_GC` suivent une exponentielle de taux 1 ; `kappa` une exponentielle de taux 0,1 ; `beta` ~ Normale(0, 10) ; log T à la racine ~ Normale(ln 75, 10) ; logit GC à la racine ~ Normale(0, 10).
- Inférence par MCMC : 100 000 générations, un échantillon toutes les 10 générations.

</details>

## 5. Analyses réalisées

Les analyses ont été menées dans cet ordre :

1. **Vérification sans arbre.** Corrélation de Pearson et régression linéaire entre GC et OGT chez les 33 espèces actuelles, pour vérifier que nos données contiennent le signal attendu.
2. **Modèle complet** (`archaea_simple.Rev`). Tout est estimé : `beta` (objectif 1) et la température de la racine (objectif 2).
3. **Deux analyses de contrôle**, pour comprendre comment l'information sur le GC atteint la température :
   - GC de la racine fixé à environ 86 % (`archaea_simpleRoot.Rev`) ;
   - GC de la racine fixé et `beta = 0`, c'est-à-dire sans lien entre GC et température (`archaea_GCRoot_β0.Rev`).

Chaque analyse MCMC produit 10 001 échantillons (fichiers `analyses/*.log`). Les 10 % premiers sont écartés (burn-in). Les résultats sont résumés par la moyenne a posteriori et l'intervalle de crédibilité à 95 % (quantiles 2,5 % et 97,5 %).

## 6. Résultats

### 6.1 Objectif 1 : la corrélation entre GC et température

**Sans tenir compte de l'arbre**, la relation est très forte :

![Relation entre le taux de GC et la température optimale de croissance](figures/gc_vs_temperature.png)

- Corrélation de Pearson : **r = 0,949** (p < 2,2 × 10⁻¹⁶)
- R² de la régression : **0,90**
- Nombre d'espèces : **33**

Cette valeur est très proche de celle de Groussin et Gouy (r = 0,95) : nos données reproduisent bien le signal attendu. Mais cette analyse traite les espèces comme indépendantes, alors que des espèces proches se ressemblent en partie par héritage.

**Le long de l'arbre**, avec le modèle complet :

- `beta` : moyenne **0,96**, intervalle de crédibilité à 95 % **[0,66 ; 1,41]** ;
- probabilité a posteriori que `beta` soit positif : **> 0,999** (tous les échantillons retenus sont positifs).

Quand la température d'une lignée augmente au cours de l'évolution, son taux de GC tend donc lui aussi à augmenter. Ce résultat est stable : dans le contrôle 1 (GC de la racine fixé), `beta` vaut 0,98 [0,71 ; 1,42].

Les deux résultats ne disent pas la même chose : le r de Pearson décrit une ressemblance entre espèces actuelles, alors que `beta` décrit une association entre les changements des deux caractères au cours de l'évolution, en tenant compte de la parenté.

### 6.2 Objectif 2 : la température de la racine

| Analyse | GC de la racine | `beta` | Température de la racine (IC 95 %) |
| --- | --- | --- | ---: |
| Modèle complet | estimé : 82,4 % [80,5 ; 84,2] | estimé | **80,9 °C** [58,8 ; 109,4] |
| Contrôle 1 | fixé à 85,9 % | estimé | **105,6 °C** [75,3 ; 141,0] |
| Contrôle 2 | fixé à 85,9 % | fixé à 0 | **72,2 °C** [47,9 ; 103,7] |

**Résultat principal.** Avec le modèle complet, la température de l'ancêtre commun est estimée à environ 81 °C, avec une incertitude large : entre 59 et 109 °C (intervalle à 95 %). Cet intervalle exclut un ancêtre mésophile et contient les estimations de Groussin et Gouy (90 °C avec l'ARNr, 82 °C avec les protéines) : les deux approches concluent à un ancêtre chaud. Il ne permet pas en revanche de trancher entre thermophile et hyperthermophile (seuil de 80 °C dans l'article).

**Contrôles.** Fixer le GC de la racine à 85,9 %, une valeur plus élevée que celle estimée par le modèle (82,4 %), fait monter la température de la racine à environ 106 °C : l'information sur le GC passe donc vers la température. Si l'on supprime en plus le lien (`beta = 0`), le GC ne renseigne plus sur la température, qui n'est plus estimée qu'à partir des températures actuelles : environ 72 °C. `beta` est donc bien le canal par lequel le GC informe la température.

Ces écarts sont à lire avec prudence : les intervalles des trois analyses se chevauchent largement. Ils montrent surtout que la reconstruction est sensible à ce que l'on impose à la racine.

## 7. Conclusion

Oui, il est possible d'estimer dans un seul modèle à la fois la corrélation évolutive entre GC et température et la température des ancêtres :

1. **Le GC de l'ARNr et la température ont évolué ensemble chez les archées** : `beta` est positif (0,96, intervalle de crédibilité à 95 % [0,66 ; 1,41]).
2. **Cette corrélation permet d'estimer la température de l'ancêtre commun** : environ 81 °C (entre 59 et 109 °C), ce qui rejoint la conclusion de Groussin et Gouy d'une origine chaude des archées, obtenue ici par une méthode différente.
3. **`beta` est le mécanisme qui relie les deux caractères** : sans lui, le GC n'apporte plus d'information sur la température.

## 8. Limites et perspectives

### Limites

- **Données restreintes.** 33 espèces et un seul signal moléculaire, l'ARNr, alors que Groussin et Gouy utilisaient aussi les protéines.
- **Un seul lien pour tout l'arbre.** Le modèle suppose la même relation entre GC et température (un seul `beta`) dans toutes les lignées.
- **Arbre fixé.** L'incertitude sur la topologie et les longueurs de branches n'est pas prise en compte.
- **Modèle de substitution simplifié.** Le modèle T92 est plus simple que le modèle HKY85 retenu par l'article pour l'ARNr.
- **Sensibilité à la racine.** L'estimation de 105,6 °C obtenue en fixant le GC de la racine dépasse nettement les valeurs de l'article.
- **Incertitude large sur la racine.** Les intervalles de crédibilité de la température de la racine couvrent environ 50 °C : les valeurs centrales doivent toujours être lues avec leur intervalle.
- **Convergence.** Dans les modèles où `beta` est estimé, `beta`, `sigma_T` et `sigma_GC` ont des tailles d'échantillon effectives (ESS) d'environ 100 à 160, en dessous du seuil de 200 habituellement recommandé. Les chaînes sont stables (les deux moitiés de chaque chaîne donnent des moyennes proches), mais des analyses plus longues rendraient ces estimations plus précises.

### Perspectives

- Ajouter la composition des protéines comme second thermomètre moléculaire.
- Permettre à la relation entre GC et température de varier entre les grands groupes d'archées.
- Tester d'autres topologies pour mesurer la sensibilité des reconstructions à l'arbre.
- Comparer formellement les modèles avec et sans lien (contrôles 1 et 2) par un facteur de Bayes.
- Allonger les analyses MCMC (ou ajuster les propositions) pour atteindre un ESS d'au moins 200 pour tous les paramètres.
- Comparer notre reconstruction à une reconstruction en deux étapes, comme celle de l'article, sur les mêmes données.

## Références

Groussin M., Gouy M. (2011). Adaptation to Environmental Temperature Is a Major Determinant of Molecular Evolutionary Rates in Archaea. *Molecular Biology and Evolution*, 28(9), 2661-2674. https://doi.org/10.1093/molbev/msr098

Boussau B., Blanquart S., Necsulea A., Lartillot N., Gouy M. (2008). Parallel adaptations to high temperatures in the Archaean eon. *Nature*, 456, 942-945.
