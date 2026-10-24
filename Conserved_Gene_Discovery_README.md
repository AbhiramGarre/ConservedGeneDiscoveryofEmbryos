# Conserved Gene Discovery — Preliminary Analysis Package

## Research question

Which orthologous gene-expression programs are shared across early zebrafish embryonic spatial patterning and mouse gastrulation/early mesoderm differentiation?

## What was actually analyzed

Two refine.bio ZIP archives were inspected and the full expression matrices were loaded. The results here are descriptive quality checks and within-species gene variability rankings. **No mouse–zebrafish ortholog mapping has been completed, and no conserved genes have been identified yet.**

| Data | Zebrafish | Mouse |
| --- | ---: | ---: |
| Refine.bio dataset | `f2f28d77-0c5c-40aa-aa62-509e40294760` | `5ffeadd4-da71-4a27-a895-13caea41400a` |
| Experiment | SRP055996 | ERP013381 |
| Genes/features | 31,882 | 41,249 |
| Samples | 1,138 | 1,006 |
| Missing or nonfinite expression values | 0 | 0 |
| Exact zero values | 0 | 0 |
| Duplicate gene IDs | 0 | 0 |
| Minimum expression value | 0.008134 | 0.357073 |
| Maximum expression value | 4.976743 | 182.460678 |
| Quantile normalized | yes | yes |
| Scale by | NONE | NONE |
| Refine.bio KS warning | yes | yes |

### Warnings and interpretation

- Samples are columns, genes/features are rows. The archive values are **processed quantile-normalized expression values**, not raw UMI count matrices. Their absence of zeroes cannot be used to estimate dropout/detection or RNA library size. Direct count-based models (e.g., standard count assumptions in differential expression) are inappropriate without original raw counts.
- The original publication reports 851 *filtered* zebrafish cells; these archives contain 1,138 zebrafish samples. The original mouse study reports 1,205 cells; these archives contain 1,006 mouse samples. Reconcile sample sets with original publications and author metadata before biological interpretation.
- The zebrafish sample metadata is mostly unlabeled beyond the uniform 50%-epiboly stage; mouse samples lack cell-state/stage labels in the bundled TSV. Author-supplied cell annotations and stage labels must be obtained externally if we want state-specific interpretation.
- Cell types and developmental stages differ between the organisms. Direct expression correlations between zebrafish and mouse cells are biologically unjustified unless comparable cell states are defined first.
- Ortholog mapping is mandatory. Zebrafish gene duplication frequently produces one-to-many homology relationships. Start with high-confidence Ensembl Compara **one-to-one** orthologs, and analyze duplicated gene families separately.
- Even when a gene has high variability in both species, this is only an exploratory signal, not proof that its function or regulatory control is conserved.

## Files in the package

- `embryo_dataset_qc.json`: summary of exact numbers computed from both matrices.
- `mouse_gene_variability.csv`, `zebrafish_gene_variability.csv`: processed gene means and log1p variance across cells (not raw-count HVGs).
- `mouse_sample_summary.csv`, `zebrafish_sample_summary.csv`: diagnostic means and sums of processed values; **these are not raw library sizes**.
- `rank_conserved_candidates.py`: optional next-stage script to rank jointly variable, verified one-to-one orthologs; requires a mapping file.
- `inspect_embryos.py`: the code used for the full-matrix QC calculations (requires ZIP archives in the same directory).

## How to create a verified ortholog mapping

Retrieve a mouse-to-zebrafish Ensembl BioMart/Compara ortholog table from <https://www.ensembl.org/biomart/martview>. Choose the *Mouse genes* dataset. Select mouse Ensembl Gene ID and zebrafish homolog Ensembl Gene ID, along with the zebrafish homolog orthology type. Rename the exported columns to:

```
mouse_gene_id,zebrafish_gene_id,orthology_type
```

Retain relationships explicitly marked `ortholog_one2one`. Include gene names/symbols if available. Check for ID version mismatches and report the proportion of genes matching each input expression matrix. Keep track of Ensembl release and archive/build to ensure reproducibility.

## Run a first exploratory ortholog screen

```
python rank_conserved_candidates.py \
  --mouse mouse_gene_variability.csv \
  --zebrafish zebrafish_gene_variability.csv \
  --orthologs verified_mouse_zebrafish_orthologs.csv \
  --output conserved_gene_results
```

The script computes species-specific ranks of *log1p-expression variability* over the observable one-to-one ortholog universe. It reports intersection of top 10% variable orthologs and a hypergeometric overlap p-value. **It is a screening method, not a discovery validation.** It does not compare absolute expression levels across species.

## Recommended full scientific project workflow

1. Obtain the paper's original QC criteria, cell-state/stage metadata, and, if feasible, non-quantile-normalized expression/counts for both studies.
2. Perform cell QC and within-species clustering using an approach appropriate for the underlying assay values. Annotate zebrafish margin/mesendoderm and mouse epiblast/nascent mesoderm populations with author labels and known marker sets.
3. Get versioned, one-to-one Ensembl orthologs and harmonize identifiers; quantify mapping coverage.
4. Within each species, identify cell-state-associated gene signatures in biologically analogous subpopulations using differential expression or rank-based approaches with suitable controls.
5. Compare mapped ortholog gene signatures with overlap/enrichment tests (e.g., Fisher/hypergeometric with a defensible tested-gene universe); test directionality consistency only for biologically comparable contrasts.
6. Identify shared programs using GO/Reactome enrichment with ortholog-aware background and, when feasible, gene-module analysis; evaluate stability under downsampling.
7. Report which candidates reproduce across alternative filtering thresholds. Explicitly distinguish sequence orthology, shared expression signatures, and functional conservation.

## Original studies and sources

- Satija et al., *Spatial reconstruction of single-cell gene expression*, Nature Biotechnology (2015), DOI 10.1038/nbt.3192. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4430369/>
- Scialdone et al., *Resolving Early Mesoderm Diversification through Single Cell Expression Profiling*, Nature (2016), DOI 10.1038/nature18633. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4947525/>
- Ensembl BioMart documentation: <https://www.ensembl.org/info/data/biomart/how_to_use_biomart.html>

The project is an analysis of public data and cannot, by itself, prove causal or functional gene conservation.
