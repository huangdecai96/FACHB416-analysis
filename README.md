# FACHB-416 analysis code

Code supporting *Morphological plasticity as a potential buffer for gene-family
shifts in a green alga*. Prepared from the final Science Advances manuscript
inputs on 12 November 2026. Draft repository version: 1.0.0.

## Get the matching data

Extract `FACHB416-source-data-v1.0.0-draft04.zip` beside this repository:

```
working-folder/
  FACHB416-analysis/
  FACHB416-source-data/
```

Repository URLs and DOIs will be added after deposits are created. The data
archive supplies a README, field inventory, figure mapping and checksum index.
Raw sequence deposits remain at NGDC: GWHJXYB00000000.1 and CRA053458,
BioProject PRJCA074263. No reviewer-access credentials are included here.

## Reproduce the final statistics

Python 3.12 was used for verification. Install dependencies in your own environment:

```
python -m pip install -r requirements.txt
python reproduce_statistics.py --data ../FACHB416-source-data
Rscript reproduce_growth.R ../FACHB416-source-data
```

Python rebuilds set scores from the full VST matrix and membership, Welch/BH24,
interaction F-tests/BH4, exact correlations/BH28, phenotype ANOVA, CPP, six DEG
summaries and family-level CAFE consensus counts. R repeats all 18 final 17-day
logistic fits with 125 starting points, then ANOVA and Tukey comparisons.

For the raw-count differential-expression analysis, install DESeq2 through
Bioconductor separately, then run:

```
Rscript reproduce_deseq2.R ../FACHB416-source-data
```

Optional output argument: `python reproduce_statistics.py --data DATA --out NEW_DIR`,
or `Rscript reproduce_growth.R DATA NEW_DIR` / `Rscript reproduce_deseq2.R DATA NEW_DIR`.
Outputs go into separate results subdirectories by default. Nonempty output folders
are refused, so an old run is not overwritten. Input data are never updated.

## Versions and scope

See environment_verified.txt. The original DESeq2 analysis used R 4.5.2 / DESeq2
1.49.9; replay with R 4.5.1 / DESeq2 1.50.2 retained all six significant-gene sets
and directions. The accepted original tables remain the submitted values.
The `figures` folder contains canonical Figure1.R–Figure5.R and final-numbered supplementary scripts. Shared functions actually used by those scripts are in figures/helpers. Upstream scientific methods are in methods/. Old plotting wrappers, duplicate style files, copied test outputs and server packaging utilities are excluded.

The main pairwise family is **24 tests**, not six separate tests per set.
Interaction and correlation families remain BH4 and BH28. Historical 15-day fits
are not an alternative to the final 17-day analysis. Final figure numbers are
defined in the data archive (the corrected CPP time series is Fig. S6).

Individual-cell dimension records (1,274 observations in 18 cultures) and both
organelle-map annotation sources are supplied. Dimension sampling and colony
distribution counting were independent. The dimension table contains sample/cell
identifiers, lengths, widths and volumes, without colony-category metadata.

The archived plotted volume means are not identical to direct arithmetic means
of the measured cells. The audit script reports this distinction without changing
the manuscript inputs. This draft does not assert numerical agreement between
the two. No manuscript figures or results were changed in this preparation.

Pixel-identical reproduction of every final panel is not claimed. Figure 6 is the final PowerPoint illustration; no R drawing script is included. See Organelles/README.md for the recovered web-map inputs.

## Citation and license

Cite the fixed Figshare code version and associated data DOI after they are
assigned; a mutable GitHub branch is not the archival citation. The manuscript
has no publication DOI assigned in this preparation.

Code author and copyright holder: Decai Huang. Copyright (c) 2026 Decai Huang.
Original project code is licensed under the MIT License; see LICENSE and
AUTHORS.md. Third-party software and separately identified third-party material
retain their own licenses. The separate source-data archive uses CC BY 4.0.


## Selected final drawing replay

The final Fig. 4 and BH24 Fig. 5 drawing definitions can also be replayed to PDF
and PNG, without changing the supplied tables or manuscript figures:

```
python reproduce_figures.py --data ../FACHB416-source-data --out ../drawing-check
```

Individual figures can be run directly:

```
Rscript figures/Figure4.R ../FACHB416-source-data ../Figure4-output
Rscript figures/Figure5.R ../FACHB416-source-data ../Figure5-output
```

Figure1.R generates the scientific panels before manual PowerPoint edits. Figure2.R exports the retained approved composition rather than substituting an unused redraw; it writes PNG only.

The optional runner accepts --figures followed by canonical figure names. Each output directory must be new and separate from source data. The default replay selects Figure4 and Figure5. Shared R packages are ggplot2, patchwork, ragg and base grid; Figure1 also uses circlize; Figure2 uses png. Fonts use Arial. No installed environment is changed by these scripts.

See figures/Non_R_figures.md for the PowerPoint, Hi-C image and web-generated organelle figures. Repackaging keeps the final numerical inputs and selected drawing definitions; final manual layout edits are not claimed to be regenerated automatically.

## Audit the recovered cell dimensions

```
python audit_cell_dimensions.py --data ../FACHB416-source-data
```

This verifies dimensions and actual counts, computes direct cell means, and
compares them with archived plotting values without changing either input.
