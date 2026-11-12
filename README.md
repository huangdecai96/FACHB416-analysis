# FACHB-416 analysis code

Analysis and figure-generation scripts supporting *Morphological plasticity as a potential buffer for gene-family shifts in a green alga*.

## Data and setup

Source dataset: *Source data for morphological plasticity and gene-family variation in Tetradesmus obliquus FACHB-416*.

Reserved Figshare DOI: [10.6084/m9.figshare.34995331](https://doi.org/10.6084/m9.figshare.34995331).
The dataset is currently unpublished; the DOI will become active when the dataset is published. Private access for editors and reviewers is supplied separately with the submission.

Extract `FACHB416_Figure_Data.zip` beside the code folder. If downloading a GitHub source archive, rename the extracted code folder to `FACHB416-analysis` to use the layout below:

```text
working-folder/
  FACHB416-analysis/
  FACHB416-source-data/
```

The data archive contains a short `README.md` with units and data-specific notes, and `figure_data_map.tsv` linking inputs to manuscript figures. Genome assembly and RNA-seq deposits remain at NGDC under GWHJXYB00000000.1 and CRA053458, respectively (BioProject PRJCA074263).

Run the commands below from the `FACHB416-analysis` directory. Python package requirements are in [requirements.txt](requirements.txt); tested versions are listed in [environment_verified.txt](environment_verified.txt). R and the packages used by each R script must be installed separately.

## Reproduce statistics

```sh
python -m pip install -r requirements.txt
python reproduce_statistics.py --data ../FACHB416-source-data
Rscript reproduce_growth.R ../FACHB416-source-data
```

The Python script rebuilds gene-set scores from the VST matrix and membership table, Welch tests with BH24 adjustment, interaction tests with BH4 adjustment, exact Spearman tests with BH28 adjustment, phenotype ANOVA, CPP, DEG counts and CAFE consensus counts. The R script repeats all 18 final 17-day logistic fits, followed by ANOVA and Tukey comparisons.

For differential-expression analysis from raw gene counts, install DESeq2 through Bioconductor and run:

```sh
Rscript reproduce_deseq2.R ../FACHB416-source-data
```

The original DESeq2 analysis used R 4.5.2 / DESeq2 1.49.9; the verified replay used R 4.5.1 / DESeq2 1.50.2 and retained all six significant-gene sets and directions. The archived original tables remain the manuscript values.

Optional output arguments are `--out NEW_DIR` for the Python statistics script, or a second positional argument `NEW_DIR` for the growth and DESeq2 scripts. Use new output directories; nonempty output directories are refused. Input data are not modified.

## Reproduce figures

The default runner renders Figures 4 and 5:

```sh
python reproduce_figures.py --data ../FACHB416-source-data --out ../drawing-check
```

Individual figures can be run directly:

```sh
Rscript figures/Figure4.R ../FACHB416-source-data ../Figure4-output
Rscript figures/Figure5.R ../FACHB416-source-data ../Figure5-output
```

The `figures/` directory contains `Figure1.R` through `Figure5.R` and scripts named by supplementary figure number. Shared functions are in `figures/helpers/`. The optional runner accepts `--figures` followed by these figure names. Output directories must be new and separate from source data.

Figure 1 scripts generate the scientific panels before manual PowerPoint edits. Figure 2 exports and crops the retained approved composition and writes PNG only. Figure 6 was prepared in PowerPoint and has no R drawing script. See [figures/Non_R_figures.md](figures/Non_R_figures.md) for manual composition, Hi-C and web-generated organelle figures. Organelle inputs are under `Organelles/` in the data archive; its root `README.md` and `figure_data_map.tsv` describe their use.

R plotting dependencies include ggplot2, patchwork, ragg and base grid; Figure 1 also uses circlize and Figure 2 uses png. Figures use Arial. Final manual layout edits are not regenerated automatically.

## Scope and data notes

- Pairwise comparisons use one BH family of 24 tests. Interaction and correlation analyses use separate families of 4 and 28 tests. DESeq2 adjustment is within each contrast.
- Upstream methods are retained in `methods/`; see [methods/README.md](methods/README.md) for their external toolchain and input requirements. The Figshare figure-data archive is not a complete copy of every upstream intermediate.
- Cell dimensions and colony distributions were measured independently. The dimension table contains 1,274 individual-cell records across 18 cultures, without colony-category metadata. Archived plotted culture-volume means differ from direct arithmetic means of those measurements; this discrepancy remains unresolved. Both tables are retained.

The following command checks measurement counts and dimensions and compares direct cell means with archived plotting values without changing the inputs:

```sh
python audit_cell_dimensions.py --data ../FACHB416-source-data
```

## Citation and license

Code repository: [huangdecai96/FACHB416-analysis](https://github.com/huangdecai96/FACHB416-analysis). Cite Decai Huang, the release tag and the commit SHA used. The GitHub release identifies a code version; the Figshare DOI above identifies the separate dataset.

Code author and copyright holder: **Decai Huang**. Copyright (c) 2026 Decai Huang. Original project code is licensed under the **MIT License**; see [LICENSE](LICENSE) and [AUTHORS.md](AUTHORS.md). The source-data archive uses **CC BY 4.0**. Third-party software and separately identified third-party material retain their own licenses.
