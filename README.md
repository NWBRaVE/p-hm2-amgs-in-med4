# P-HM2 AMGs in MED4
Flux analysis of cyanophage P-HM2 auxilliary metabolic gene (AMG) impacts in *Prochlorococcus marinus* MED4.

This repository contains various scripts and figures used to produce the following manuscript:


> Rozum, Sineath, Kim, Bohutsky, Johnson, Evans, Pollock, Qian, Cheung, Wu, Feng. (2024).<br>
*Synergy and antagonism in a genome-scale model of metabolic hijacking by bacteriophage.*

This README will be updated with publication information when it becomes available. A preprint is available on [bioarxiv](https://doi.org/10.1101/2024.12.11.628001).

# Requirements and Installation

The scripts and notebooks were run using Python version 3.12.8. See `requirements.txt` for Python packages installed during execution.

The baseline model, [iSO595](https://github.com/segrelab/Prochlorococcus_Model), is imported and modified (e.g., to incorporate phage biomass) in `generate_biomass_functions.ipynb`.

The remaining notebooks produce tables and figures presented in the manuscript. The Python files contain various utility functions for the analysis.

# Acknowledgments

This work is supported by the NW-BRaVE for Biopreparedness project funded by the U. S. Department of Energy (DOE), Office of Science, Office of Biological and Environmental Research, under FWP 81832. A portion of this research was performed on a project award (https://doi.org/10.46936/staf.proj.2023.61054/60012367) from the Environmental Molecular Sciences Laboratory, a DOE Office of Science User Facility sponsored by the Biological and Environmental Research program under Contract No. DE-AC05-76RL01830. Pacific Northwest National Laboratory is a multi-program national laboratory operated by Battelle for the DOE under Contract DE-AC05-76RLO 1830.

The following people contributed to this work:

Jordan C. Rozum, William Sineath, Doo Nam Kim, Pavlo Bohutskyi, Connah Johnson, James Evans, David Pollock, Wei-Jun Qian, Margaret S. Cheung, Ruonan Wu, Song Feng.
