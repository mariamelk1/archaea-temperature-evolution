# Température et évolution moléculaire chez les archées

## 1. Introduction

### 1.1 Contexte

Les archées sont des micro-organismes qui vivent à des températures très variées, des milieux tempérés jusqu'à des sources chaudes proches de 100 °C.

La température laisse une trace dans leurs molécules. L'ARN ribosomique (ARNr) contient de nombreuses régions en double brin, stabilisées par l'appariement des bases. Les paires G-C, liées par trois liaisons hydrogène, sont plus stables que les paires A-U, qui n'en ont que deux. Les espèces des milieux chauds ont donc un ARNr plus riche en GC.

Le taux de GC de l'ARNr peut ainsi servir de **thermomètre moléculaire**. C'est particulièrement utile pour les ancêtres : on ne peut pas mesurer leur température, mais on peut reconstruire leurs séquences.

### 1.2 Objectifs

Ce projet a deux objectifs, qui s'enchaînent :

1. **Estimer la corrélation entre le taux de GC et la température le long de l'arbre des archées.**
   Quand la température optimale de croissance (OGT) d'une lignée augmente au cours de l'évolution, son taux de GC augmente-t-il aussi ? Cette corrélation est mesurée par un paramètre du modèle, `beta`.

2. **Utiliser cette corrélation pour estimer la température des ancêtres, en particulier celle de la racine.**
   La racine de l'arbre représente l'ancêtre commun des espèces étudiées : vivait-il dans un milieu chaud ?

> **Question de recherche :** peut-on estimer, dans un seul modèle phylogénétique bayésien, à la fois la corrélation évolutive entre GC et température et la température des ancêtres ?

Cette question a déjà été abordée par Groussin et Gouy (2011) avec une autre méthode. Nous présentons d'abord leur travail et ses limites, puis la façon dont notre modèle y répond.

## 2. L'étude de référence : Groussin et Gouy (2011)

### 2.1 Leur méthode : deux étapes séparées

Les auteurs ont étudié 35 espèces d'archées dont le génome est complètement séquencé, à partir de deux sources d'information : l'ARNr et 72 familles de protéines.

1. **Reconstruire la composition des ancêtres.** Sur un arbre phylogénétique fixé, ils estiment à chaque nœud ancestral le taux de GC de l'ARNr et la composition en acides aminés des protéines. Ils utilisent des modèles d'évolution dits « non homogènes », où la composition peut changer d'une branche à l'autre.
2. **Traduire ces compositions en températures.** Chez les espèces actuelles, ils mesurent la relation entre composition et OGT (r = 0,95 pour le GC de l'ARNr, r = 0,84 pour un indice de composition des protéines). Cette droite sert de thermomètre : ils y placent les compositions ancestrales pour en déduire la température des ancêtres.

Pour vérifier que ces corrélations ne viennent pas simplement de la parenté entre espèces, ils les recalculent avec des contrastes phylogénétiquement indépendants (PIC), qui confirment le lien.

### 2.2 Leurs résultats

- L'ancêtre commun des archées était **hyperthermophile** : environ 90 °C d'après l'ARNr et 82 °C d'après les protéines.
- Plusieurs lignées, notamment chez les Euryarchées, se sont ensuite adaptées progressivement à des milieux plus froids.
- Les lignées des milieux tempérés évoluent plus vite (branches plus longues) : la température apparaît comme un déterminant majeur de la vitesse d'évolution moléculaire chez les archées.

### 2.3 Les limites de leur approche

- **L'incertitude se perd entre les deux étapes.** La droite GC-température est utilisée comme si elle était exacte. Les auteurs le reconnaissent : cette incertitude n'est pas prise en compte, et les intervalles réels seraient plus larges.
- **Le thermomètre traite les espèces comme indépendantes.** Les PIC vérifient la corrélation, mais la droite qui convertit le GC en température est ajustée sur les valeurs brutes des espèces.
- **La température n'évolue pas dans leur modèle.** Les températures des espèces actuelles servent seulement à calibrer le thermomètre ; elles n'informent pas directement les ancêtres à travers l'arbre.
- **Les résultats dépendent des données choisies.** Leur estimation pour l'ancêtre commun (de 74 à 89 °C) est nettement plus élevée que celle d'une étude antérieure (Boussau et al., 2008 : de 59 à 73 °C), écart qu'ils attribuent au choix des gènes et à l'incertitude des thermomètres.

## 3. Notre approche : un seul modèle

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

## 4. Données

L'analyse porte sur 33 espèces d'archées. Le dossier `data/` contient :

| Fichier | Contenu |
| --- | --- |
| `data/archaea.nex` | Alignement de l'ARNr : 1 801 positions des régions en double brin, dont la composition est la plus liée à la température |
| `data/archaea.tree` | Arbre phylogénétique des 33 espèces (topologie et longueurs de branches, considérées comme fixes) |
| `data/archaea_temp_gc.nex` | Température optimale de croissance (OGT) et taux de GC de l'ARNr de chaque espèce actuelle |
| `data/rna.itgc` | Fichier auxiliaire (85 taxons), non utilisé dans l'analyse actuelle |

## 5. Analyses réalisées

Les analyses ont été menées dans cet ordre :

1. **Vérification sans arbre.** Corrélation de Pearson et régression linéaire entre GC et OGT chez les 33 espèces actuelles, pour vérifier que nos données contiennent le signal attendu.
2. **Modèle complet** (`archaea_simple.Rev`). Tout est estimé : `beta` (objectif 1) et la température de la racine (objectif 2).
3. **Deux analyses de contrôle**, pour comprendre comment l'information sur le GC atteint la température :
   - GC de la racine fixé à environ 86 % (`archaea_simpleRoot.Rev`) ;
   - GC de la racine fixé et `beta = 0`, c'est-à-dire sans lien entre GC et température (`archaea_GCRoot_β0.Rev`).

## 6. Résultats

### 6.1 Objectif 1 : la corrélation entre GC et température

**Sans tenir compte de l'arbre**, la relation est très forte :

![Relation entre le taux de GC et la température optimale de croissance](figures/gc_vs_temperature.png)

- Corrélation de Pearson : **r = 0,949** (p < 2,2 × 10⁻¹⁶)
- R² de la régression : **0,90**
- Nombre d'espèces : **33**

Cette valeur est très proche de celle de Groussin et Gouy (r = 0,95) : nos données reproduisent bien le signal attendu. Mais cette analyse traite les espèces comme indépendantes, alors que des espèces proches se ressemblent en partie par héritage.

**Le long de l'arbre**, avec le modèle complet, l'intervalle de crédibilité à 95 % de `beta` est **[0,5 ; 2,8]**. Il est entièrement positif : quand la température d'une lignée augmente au cours de l'évolution, son taux de GC tend lui aussi à augmenter.

Les deux résultats ne disent pas la même chose : le r de Pearson décrit une ressemblance entre espèces actuelles, alors que `beta` décrit une association entre les changements des deux caractères au cours de l'évolution, en tenant compte de la parenté.

### 6.2 Objectif 2 : la température de la racine

| Analyse | GC de la racine | `beta` | Température de la racine |
| --- | --- | --- | ---: |
| Modèle complet | estimé | estimé | **80,91 °C** |
| Contrôle 1 | fixé (≈ 86 %) | estimé | **105,61 °C** |
| Contrôle 2 | fixé (≈ 86 %) | fixé à 0 | **72,15 °C** |

**Résultat principal.** Avec le modèle complet, l'ancêtre commun est estimé à environ 81 °C, juste au-dessus du seuil de 80 °C qui sépare thermophiles et hyperthermophiles dans l'article. C'est un peu moins que l'estimation de Groussin et Gouy avec l'ARNr (90 °C) et proche de celle obtenue avec les protéines (82 °C) : les deux approches concluent à un ancêtre chaud.

**Contrôles.** Fixer un GC élevé à la racine fait monter sa température à près de 106 °C : l'information sur le GC passe donc bien vers la température. Si l'on supprime en plus le lien (`beta = 0`), le GC ne renseigne plus sur la température, qui n'est plus estimée qu'à partir des températures actuelles : environ 72 °C. `beta` est donc bien le canal par lequel le GC informe la température. Ces contrôles montrent aussi que la reconstruction est sensible à ce que l'on impose à la racine.

## 7. Conclusion

Oui, il est possible d'estimer dans un seul modèle à la fois la corrélation évolutive entre GC et température et la température des ancêtres :

1. **Le GC de l'ARNr et la température ont évolué ensemble chez les archées** : `beta` est positif, avec un intervalle de crédibilité à 95 % de [0,5 ; 2,8].
2. **Cette corrélation permet d'estimer la température de l'ancêtre commun** : environ 81 °C, ce qui rejoint la conclusion de Groussin et Gouy d'une origine chaude des archées, obtenue ici par une méthode différente.
3. **`beta` est le mécanisme qui relie les deux caractères** : sans lui, le GC n'apporte plus d'information sur la température.

## 8. Limites et perspectives

### Limites

- **Données restreintes.** 33 espèces et un seul signal moléculaire, l'ARNr, alors que Groussin et Gouy utilisaient aussi les protéines.
- **Un seul lien pour tout l'arbre.** Le modèle suppose la même relation entre GC et température (un seul `beta`) dans toutes les lignées.
- **Arbre fixé.** L'incertitude sur la topologie et les longueurs de branches n'est pas prise en compte.
- **Modèle de substitution simplifié.** Le modèle T92 est plus simple que le modèle HKY85 retenu par l'article pour l'ARNr.
- **Sensibilité à la racine.** L'estimation de 105,61 °C obtenue en fixant le GC de la racine dépasse nettement les valeurs de l'article.
- **Convergence.** Les paramètres des nœuds profonds peuvent être moins bien échantillonnés par le MCMC ; la convergence et les tailles d'échantillon effectives (ESS) doivent être vérifiées dans Tracer.

### Perspectives

- Ajouter la composition des protéines comme second thermomètre moléculaire.
- Permettre à la relation entre GC et température de varier entre les grands groupes d'archées.
- Tester d'autres topologies pour mesurer la sensibilité des reconstructions à l'arbre.
- Comparer formellement les modèles avec et sans lien (contrôles 1 et 2) par un facteur de Bayes.
- Comparer notre reconstruction à une reconstruction en deux étapes, comme celle de l'article, sur les mêmes données.

## Références

Groussin M., Gouy M. (2011). Adaptation to Environmental Temperature Is a Major Determinant of Molecular Evolutionary Rates in Archaea. *Molecular Biology and Evolution*, 28(9), 2661-2674. https://doi.org/10.1093/molbev/msr098

Boussau B., Blanquart S., Necsulea A., Lartillot N., Gouy M. (2008). Parallel adaptations to high temperatures in the Archaean eon. *Nature*, 456, 942-945.