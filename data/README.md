# data/

This folder is **not versioned**. The count matrices, metadata and fusion calls come from patient-derived and cell line models, so only this file is tracked. It documents what the scripts expect to find here.

To run the pipeline on your own experiment, place the files below in this folder, following the naming and column conventions.

## Layout

```
data/
├── counts_<DATASET>.csv        # one count matrix per dataset
├── metadata_<DATASET>.csv      # one sample sheet per dataset
└── fusions/                    # optional, only for the fusion scripts
    ├── <Sample>_star-fusion.fusion_predictions.abridged.coding_effect.tsv
    └── <Sample>.xlsx
```

`<DATASET>` is the value set in the `DATASET` variable of the CONFIG block at the top of each script.

## Count matrices: `counts_<DATASET>.csv`

Raw, **unnormalized** gene-level counts. One row per gene, one column per sample.

| Column | Content |
|--------|---------|
| `id_gene` | Ensembl gene ID with version, e.g. `ENSG00000000003.10` |
| `gene_name` | Gene symbol |
| `gene_type` | Gene biotype |
| one column per sample | Integer counts |

Sample column names must match the `Sample` column of the metadata file. The separator (`;` or `,`) is detected automatically.

## Sample sheets: `metadata_<DATASET>.csv`

Semicolon-separated, one row per sample.

| Column | Content |
|--------|---------|
| `Sample` | Sample name, identical to the column name in the count matrix |
| `Model` | Experimental model |
| `Cell_Line` | Cell line or patient-derived culture |
| `Condition` | `Parental` or `Resistance` |
| `Treatment` | Treatment applied to the sample |
| `Drug` | Optional. Drug name, needed when a dataset has more than one drug |

## Fusion calls: `fusions/`

Used only by the fusion scripts (`rna_08.1_fusions_ret.R` and the `ALK_08.1_*` scripts).

| File | Format | Read by |
|------|--------|---------|
| `<Sample>_star-fusion.fusion_predictions.abridged.coding_effect.tsv` | Tab-delimited STAR-Fusion output (`#FusionName`, `JunctionReadCount`, `SpanningFragCount`, ..., `FFPM`) | `scripts/specific/ALK/` |
| `<Sample>.xlsx` | STAR-Fusion plus FusionInspector output, one workbook per sample | `scripts/specific/RET/` |

## Data availability

The sequencing data are not included in this repository. Accession details will be added after publication.
