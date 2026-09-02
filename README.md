# Evolutionarily constrained positioning of low-adaptation synonymous codon segments relative to transmembrane-domain boundaries in Enterobacterales

This repository contains the analysis code, frozen parameters, processed data, statistical outputs, and figure-generation scripts associated with the study:

**Evolutionarily constrained positioning of low-adaptation synonymous codon segments relative to transmembrane-domain boundaries in Enterobacterales**

## Overview

This study tests whether locally low-adaptation synonymous codon segments show evolutionarily conserved positioning relative to homologous transmembrane-domain (TMD) boundaries across Enterobacterales membrane-protein families.

The analysis uses:

- strict one-to-one orthologous groups across 10 representative Enterobacterales genomes;
- reviewed membrane-protein topology information supplemented by DeepTMHMM;
- species-specific synonymous codon adaptation metrics;
- homologous TMD starts and ends as primary positional anchors;
- within-protein synonymous-codon permutation null models;
- sensitivity and fixed-presence positional analyses;
- homologous non-TMD pseudo-anchor controls;
- one-to-one matched soluble-protein family controls.

The analysis plan and primary parameter choices were frozen before outcome-based interpretation.

## Repository structure

```text
config/
    Frozen analysis parameters and protocol settings.

docs/
    Analysis specifications and manuscript-related documentation.

scripts/
    Python and shell scripts used for data processing, statistical analysis,
    null-model generation, matched-control construction, and figure generation.

data/processed/
    Processed datasets required for reproducing the principal analyses.

results/
    Final statistical outputs and analysis summaries.

figures/
    Main and supplementary figures.

supplementary/
    Supplementary tables and supplementary figure legends.

## License

Analysis code is released under the MIT License. See `LICENSE` for details.

## Reproducibility

`core_scripts/` contains portable scripts that reproduce the released supplementary figures from the processed datasets and final result tables included in this repository.

`workflow_scripts/` contains the original analysis scripts used during the study. These scripts preserve the project-relative paths and intermediate-file structure of the original working directory and are provided for workflow transparency rather than as a fully self-contained raw-data pipeline.

The processed data and result tables included here are sufficient to inspect the principal analyses and reproduce the released summary figures without the original large intermediate working files.

## Repository and archive

GitHub repository:
https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning

Archived release:
https://doi.org/10.5281/zenodo.22254908
