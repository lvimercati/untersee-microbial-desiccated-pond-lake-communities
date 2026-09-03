# Lake Untersee Cross-Habitat Microbial Communities — Analysis Code

Code accompanying: *Environmental filtering shapes divergent bacterial and
eukaryotic community structure and connectivity across desiccated, pond, and
lake mat habitats in the Untersee Oasis, Antarctica* (Vimercati et al.,
submitted to *Applied and Environmental Microbiology*).

All analyses exclude the lake filament tuft samples (see manuscript Methods).

## Pipeline order

Scripts are numbered in the order they are typically run. Each is
self-contained given its stated inputs, but later scripts generally depend
on the rarefied ASV tables and metadata objects produced early in
`01_16S_processing_and_analysis.R` / `02_18S_processing_and_analysis.R`.

| # | Script | Description |
|---|--------|-------------|
| 01 | `01_16S_processing_and_analysis.R` | DADA2 processing, ASV inference, taxonomy assignment, and core composition/diversity analyses for 16S rRNA gene data |
| 02 | `02_18S_processing_and_analysis.R` | Same pipeline for 18S rRNA gene data |
| 03 | `03_core_microbiome.R` | Core microbiome analysis (prevalence across detection thresholds) |
| 04 | `04_indicator_species.R` | Indicator species analysis (bacteria, cyanobacteria, eukarya) |
| 05 | `05_betaNTI_assembly.R` | β-nearest taxon index (βNTI) community assembly analysis |
| 06 | `06_feast_sourcetracking.R` | FEAST source-tracking analysis, all six sink/domain combinations |
| 07 | `07_feast_figure.R` | Builds the FEAST source-contribution pie chart figure from FEAST output |
| 08 | `08_cooccurrence_networks.R` | Cross-domain (16S–18S) co-occurrence network construction and figures |
| 09 | `09_faprotax_functional_analysis.R` | FAPROTAX functional group prediction and all associated figures/statistics |
| 10 | `10_16S_taxa_stats_allLevels.R` | Full taxa-level summary statistics (Kruskal-Wallis + Dunn) across all taxonomic levels |
| 11 | `11_composition_figures_allLevels.R` | Stacked composition bar charts at Class/Order/Family/Genus levels, both domains |
| 12 | `12_mantel_16S.R` | Mantel test / distance-decay analysis, 16S rRNA gene data |
| 13 | `13_mantel_18S.R` | Mantel test / distance-decay analysis, 18S rRNA gene data |
| 14 | `14_mantel_master_combined.R` | Sources scripts 12 and 13 in isolated environments and builds the combined Figure 7 (16S/18S side by side, shared legend) |
| 15 | `15_combined_diversity_and_composition_figures.R` | Combined 16S+18S figures: the alpha/beta diversity 4-panel figure (Figure 5) and the phylum-level composition figure (Figure 2). Requires objects from scripts 01 and 02. |

## Requirements

R (≥ 4.0) with the following packages: dada2, ShortRead, dplyr, tidyr,
Hmisc, ggplot2, vegan, phyloseq, mctoolsr, microbiome, iCAMP, picante,
ape, igraph, compositions, FEAST, patchwork, cowplot, RColorBrewer,
FSA, multcompView, scales.

## Data availability

Raw sequence data are deposited under BioProject PRJNA1366640 (this
study) and PRJNA638357 (previously published benthic lake mat data,
Greco et al. 2020).
