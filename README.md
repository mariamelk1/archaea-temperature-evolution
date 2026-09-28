# Temperature and Molecular Evolution in Archaea

## Biological background

Archaea are microorganisms capable of living across a very wide range of environmental temperatures, including extremely hot environments. Environmental temperature can influence molecular evolution and leave a detectable signal in the molecular composition of organisms.

In particular, the nucleotide composition of ribosomal RNA (rRNA) in Archaea is strongly associated with temperature. Archaeal species adapted to higher temperatures tend to have a higher GC content in their rRNA than species adapted to lower temperatures.

This project is based on the study by Groussin and Gouy (2011), *Adaptation to Environmental Temperature Is a Major Determinant of Molecular Evolutionary Rates in Archaea*. The authors investigated the evolution of optimal growth temperature (OGT) across the archaeal domain and its relationship with molecular evolution. Their results support an important role of environmental temperature in archaeal molecular evolution and suggest that the common ancestor of the studied Archaea was adapted to high temperatures.

## Project objectives

The general objective of this project is to investigate the evolutionary relationship between optimal growth temperature (OGT) and rRNA GC content in Archaea while taking their phylogenetic history into account.

The project has two main objectives.

1. **Test whether changes in optimal growth temperature and rRNA GC content are associated during archaeal evolution.**

   We use a phylogenetic model in RevBayes in which temperature and GC content evolve along the branches of the archaeal phylogeny. The parameter `beta` describes the association between evolutionary changes in temperature and changes in GC content.

   A positive value of `beta` indicates that increases in temperature tend to be associated with increases in GC content, whereas a value close to zero indicates little evidence for correlated evolutionary change.

2. **Reconstruct ancestral optimal growth temperatures, particularly the temperature at the root of the archaeal phylogeny.**

   Using the observed temperatures of extant species and the phylogenetic tree, the model estimates ancestral temperatures at internal nodes. Of particular interest is the posterior estimate of the root temperature, which allows us to investigate whether the common ancestor represented by our phylogeny was adapted to high temperatures.

Thus, the project moves from an observed correlation among present-day species to a phylogenetic model that investigates how temperature and GC content may have evolved together and reconstructs ancestral environmental adaptations.

## Data

The dataset contains 33 archaeal species.

- `data/archaea.nex`: rRNA sequence alignment containing 1801 nucleotide sites.
- `data/archaea.temp`: optimal growth temperature (OGT) for each species.
- `data/archaea.tree`: phylogenetic tree describing the evolutionary relationships among the archaeal species.
- `data/archaea_traits.nex`: continuous dataset containing GC content and OGT for each species.
- `data/archaea_temp_gc.nex`: continuous data used in the RevBayes phylogenetic model, containing temperature and GC content for the terminal species.

The 1801 nucleotide positions in the rRNA alignment correspond to double-stranded regions of the rRNA, which are particularly informative for studying the relationship between nucleotide composition and environmental temperature.

## Exploratory analysis: relationship between GC content and temperature

As an initial exploratory analysis, the relationship between GC content and optimal growth temperature was examined without accounting for phylogenetic relationships among species.

Both variables are continuous, so a Pearson correlation was used to evaluate the strength and direction of their linear association. A simple linear regression was also fitted to quantify the relationship between OGT and GC content.

![Relationship between GC content and optimal growth temperature](figures/gc_vs_temperature.png)

A strong positive association was observed between optimal growth temperature and GC content across the 33 archaeal species (Pearson's r = 0.949, p < 2.2 × 10^-16).

The linear regression explained approximately 90% of the observed variation in GC content (R² = 0.900). Therefore, species adapted to higher optimal growth temperatures tend to have higher GC content in their rRNA sequences.

However, this exploratory analysis treats the 33 species as statistically independent observations. This assumption is problematic because species share evolutionary history. Closely related species may have similar temperatures and GC contents partly because they inherited characteristics from common ancestors.

## Phylogenetic analysis

To account for shared evolutionary history, a phylogenetic model is implemented in RevBayes.

The model describes the evolution of optimal growth temperature along the phylogenetic tree using a Brownian-motion process on log-transformed temperature.

GC content is also allowed to evolve along the tree. Its evolutionary change is related to the change in temperature through the parameter `beta`.

The relationship can be summarized conceptually as:

change in GC content ~ beta × change in temperature + evolutionary variation

More precisely, the model operates on log-transformed temperature and logit-transformed GC content.

The main parameter of interest is therefore `beta`.

- `beta > 0`: increases in temperature tend to be associated with increases in GC content.
- `beta ≈ 0`: little evidence for an evolutionary association.
- `beta < 0`: increases in temperature tend to be associated with decreases in GC content.

The model also reconstructs ancestral values of temperature and GC content at the internal nodes of the phylogeny.

One of the main quantities of interest is the reconstructed temperature at the root of the tree. This estimate can be compared with the hypothesis proposed by Groussin and Gouy (2011) that the ancestral archaeal lineage was adapted to high temperatures.

Posterior distributions obtained from the MCMC analysis will be used to estimate `beta`, ancestral temperatures, and their associated uncertainty.

## Reference

Groussin, M. & Gouy, M. (2011). *Adaptation to Environmental Temperature Is a Major Determinant of Molecular Evolutionary Rates in Archaea*. Molecular Biology and Evolution, 28(9), 2661–2674. https://doi.org/10.1093/molbev/msr098