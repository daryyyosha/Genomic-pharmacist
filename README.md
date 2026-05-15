# Genomic Pharmacist

`Genomic Pharmacist` is a custom pharmacogenomic interpretation pipeline for VCF-based drug response analysis.

The project implements a lightweight Python workflow that analyzes pharmacogenetic markers from a VCF file and generates drug-specific recommendations without using PharmCAT. The pipeline extracts target variants, classifies SNP genotypes, calls star-allele diplotypes, applies manually curated pharmacogenomic rules, and generates both TSV output tables and a human-readable HTML report.

---

## Project task

The goal of this project was to analyze genetic markers that influence individual response to selected drugs.

These markers may include:

- single nucleotide polymorphisms (SNPs);
- star alleles;
- diplotypes;
- HLA-associated markers;
- combined pharmacogenetic rules.

The task required implementing a custom interpretation pipeline instead of using PharmCAT, a widely used pharmacogenomic clinical annotation tool [1].

The pipeline was developed to process pharmacogenomic markers from a reduced 1000 Genomes VCF dataset and produce interpretable drug recommendations based on manually curated pharmacogenomic rules.

---

## Input data

The project was designed for 1000 Genomes samples aligned to the **GRCh38 / hg38 human genome build**.

Two datasets were provided for the task:

| Dataset | Description |
|---|---|
| `Pharma_subset.vcf.gz` | Reduced VCF with 10 samples and approximately 38,000 variants relevant for pharmacological profiling |
| `1kG_Full_dataset.vcf.gz` | Full WGS VCF with 110 samples and approximately 70.9 million variants, including the same 10 pharmacogenomic samples and 100 additional European samples |

The current implementation focuses on the reduced pharmacogenomic subset VCF.

In the final local run, the input file was used as:

```text
input/Pharma_subset.vcf
```

The input VCF is not included in this repository because genomic input files can be large and may contain individual-level genetic data.

The target coordinates in `pgx_variants.tsv` are intended for the **GRCh38 / hg38** human genome build. The GRCh38 reference FASTA file is not included in the repository because it is a large public reference file and is not required for running the current pipeline.

---

## Repository contents

The repository contains the pipeline script, manually curated rule tables, and output files generated from the final run.

```text
Genomic-pharmacist/
│
├── README.md
├── pgx_pipeline.py
│
├── pgx_variants.tsv
├── star_alleles.tsv
├── rulessnp_rules.tsv
├── rulesstar_alleles.tsv
├── rules_comb_star_alleles_snp.tsv
├── rules_comb_star_alleles.tsv
│
├── extracted_variants.tsv
├── coverage_report.tsv
├── gene_calls.tsv
├── drug_report.tsv
├── summary_by_sample.tsv
└── sample_report.html
```

The repository includes the following files:

| File | Description |
|---|---|
| `pgx_pipeline.py` | Main Python pipeline script |
| `pgx_variants.tsv` | Coordinates of target pharmacogenetic markers |
| `star_alleles.tsv` | Star-allele definitions |
| `rulessnp_rules.tsv` | Direct SNP-based drug recommendation rules |
| `rulesstar_alleles.tsv` | Diplotype-to-phenotype and drug recommendation rules |
| `extracted_variants.tsv` | Variants extracted from the input VCF |
| `coverage_report.tsv` | Quality-control report for marker calling |
| `gene_calls.tsv` | SNP genotype calls and star-allele diplotype calls |
| `drug_report.tsv` | Final pharmacogenomic recommendation table |
| `summary_by_sample.tsv` | Compact summary of recommendations per sample |
| `sample_report.html` | Human-readable pharmacogenetic report |

---

## Pipeline overview

The implemented workflow converts a reduced VCF file into a structured pharmacogenetic report.

```text
Pharma_subset.vcf
        |
        |  VCFExtractor
        |  target coordinate extraction using pgx_variants.tsv
        v
extracted_variants.tsv
        |
        |  PGxCaller
        |  SNP genotype classification and star-allele diplotype calling
        v
gene_calls.tsv + coverage_report.tsv
        |
        |  RulesEngine
        |  SNP rules + star rules + combined rules
        v
drug_report.tsv + summary_by_sample.tsv
        |
        |  Reporter
        |  HTML rendering
        v
sample_report.html
```

The main pipeline script is:

```text
pgx_pipeline.py
```

It consists of four logical modules:

| Module | Purpose |
|---|---|
| `VCFExtractor` | Reads the input VCF and extracts only target pharmacogenetic variants |
| `PGxCaller` | Converts extracted VCF genotypes into SNP calls and star-allele diplotypes |
| `RulesEngine` | Applies manually curated pharmacogenomic interpretation rules |
| `Reporter` | Generates the final HTML report and summary tables |

---

## Rule tables

The interpretation is based on manually curated TSV rule tables.

| Rule table | Purpose |
|---|---|
| `pgx_variants.tsv` | Coordinates of target pharmacogenetic markers for GRCh38 / hg38 |
| `star_alleles.tsv` | Definitions of star alleles through marker IDs |
| `rulessnp_rules.tsv` | Direct SNP-based drug recommendation rules |
| `rulesstar_alleles.tsv` | Diplotype-to-phenotype and drug recommendation rules |

---

## Star allele curation

During the project, the `star_alleles.tsv` table was additionally curated.

A key issue was that some star-allele definitions cannot be represented only by a marker ID. For several markers, the same rsID can correspond to more than one possible alternative allele. Therefore, the table was extended with the `required_allele` column.

Example:

```text
CYP2D6    *50    rs267608302    G
```

This means that the `G` allele at `rs267608302` is required for this star-allele definition.

This step was necessary because some pharmacogenetic markers are multi-allelic. In such cases, the marker ID alone is not sufficient for reliable interpretation.

Ambiguous cases were manually checked using pharmacogenomic sources. Rows without reliable allele assignment were excluded from the final working table to avoid unsupported star-allele calls.

---

## HLA-B*58:01 and allopurinol

Allopurinol was considered through the clinically important `HLA-B*58:01` association, which is related to the risk of severe cutaneous adverse reactions [3].

However, the reduced VCF used in this project did not contain:

- direct HLA typing data;
- the proxy marker `rs9263726`;
- another validated HLA-B*58:01 marker available for this workflow.

Therefore, the allopurinol / `HLA-B*58:01` rule was excluded from the active interpretation workflow.

This decision was made to avoid false or unsupported clinical interpretation when the required HLA information is absent. Future support for this drug would require direct HLA typing from FASTQ/BAM data, a full VCF with validated HLA proxy markers, or integration of a dedicated HLA typing tool.

---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/daryyosha/Genomic-pharmacist.git
cd Genomic-pharmacist
```

### 2. Install dependencies

The pipeline requires Python 3 and `pandas`.

```bash
pip install pandas
```

### 3. Prepare the input VCF

Place the input VCF locally, for example:

```text
input/Pharma_subset.vcf
```

The input VCF is not included in the repository.

### 4. Run the pipeline

Since the rule tables are stored in the repository root, use the current directory as `--rules_dir`:

```bash
python pgx_pipeline.py \
  --vcf input/Pharma_subset.vcf \
  --rules_dir . \
  --out_dir output/
```

### 5. Check generated files

The pipeline produces intermediate and final output files:

```text
output/extracted_variants.tsv
output/coverage_report.tsv
output/gene_calls.tsv
output/reports/drug_report.tsv
output/reports/summary_by_sample.tsv
output/reports/sample_report.html
```

In this repository, the generated output files from the final run were uploaded to the repository root for demonstration.

---

## Output files

### `extracted_variants.tsv`

This file contains target variants extracted from the input VCF.

Columns:

```text
sample
gene
marker_id
marker_type
chrom
pos
ref
alt
raw_gt
genotype
```

This file is the direct result of target marker extraction.

---

### `coverage_report.tsv`

This is a quality-control file showing how many target markers were successfully called for each sample.

Columns:

```text
sample
total_targeted_markers
successfully_called
call_rate
```

In the final run, all 10 samples had:

```text
195 / 195 successfully called markers
call_rate = 1.0
```

---

### `gene_calls.tsv`

This is the main intermediate file with SNP and star-allele calls.

Columns:

```text
sample
gene
marker_id
call_type
result
confidence
warnings
```

It includes:

- `SNP` calls for individual markers;
- `STAR` calls for gene-level diplotypes.

Example call types and results:

```text
SNP     HOM_REF
SNP     HET
SNP     HOM_ALT
STAR    *1/*1
STAR    *1/*3
STAR    *3/*3
```

---

### `drug_report.tsv`

This is the main interpretation output.

Columns:

```text
sample
drug
rule_type
target
patient_result
phenotype
effect
recommendation
```

It contains the final matched pharmacogenomic rules and drug-specific recommendations.

---

### `summary_by_sample.tsv`

This file summarizes how many recommendation rules were matched for each sample and drug.

Columns:

```text
sample
drug
rule_matches
```

---

### `sample_report.html`

This is a human-readable pharmacogenetic report generated from `drug_report.tsv`.

The report groups results by patient and displays:

- drug;
- target gene or marker;
- patient genotype or diplotype;
- phenotype;
- expected effect;
- clinical recommendation.

---

## Final run summary

The final run was performed on 10 samples:

```text
HG00367
HG01061
HG01892
HG03069
HG04227
NA18939
NA19007
NA20299
NA21091
NA21120
```

Main output statistics:

| Output | Count |
|---|---:|
| Samples analyzed | 10 |
| Target markers per sample | 195 |
| Extracted variant records | 1950 |
| Successfully called markers per sample | 195 |
| Call rate | 1.0 |
| Gene/SNP call records | 420 |
| Final drug recommendation records | 54 |
| Drugs with generated recommendations | 6 |

The final report includes recommendations for the following drugs:

| Drug | Main pharmacogenetic target |
|---|---|
| Sertraline | `CYP2C19`, `CYP2B6` |
| Clopidogrel | `CYP2C19` |
| Efavirenz | `CYP2B6` |
| Voriconazole | `CYP2C19` |
| Tacrolimus | `CYP3A5` |
| Rosuvastatin | `ABCG2` / SNP-based rule |

---

## Method notes

### VCF parsing

The pipeline reads the VCF line by line and extracts only variants whose coordinates are present in `pgx_variants.tsv`.

Both chromosome formats are supported:

```text
chr6
6
```

The parser supports both plain VCF and gzip-compressed VCF.

---

### Genotype classification

Patient genotypes are converted into simplified classes:

| VCF genotype | Pipeline class |
|---|---|
| `0/0`, `0|0` | `HOM_REF` |
| `0/1`, `1/0`, `0|1`, `1|0` | `HET` |
| `1/1`, `1|1` | `HOM_ALT` |
| `./.`, `.`, missing genotype | `NO_CALL` |

---

### Star-allele calling

Star alleles are inferred from the presence of defining markers in the patient genotype.

If no defining variants are detected, the gene is interpreted as:

```text
*1/*1
```

If one star allele is detected, the result is reported as:

```text
*1/*X
```

If the defining markers are homozygous, the result can be reported as:

```text
*X/*X
```

If multiple possible star alleles are detected without phasing, the result is reported with lower confidence.

---

### Rule priority

The rule engine applies interpretation rules in the following order:

1. combined star + star rules;
2. combined star + SNP rules;
3. single star-allele rules;
4. single SNP rules.

This prevents less specific single-gene rules from overriding more specific combined rules.

---

## Limitations

This project is a research prototype and has several limitations.

1. HLA alleles are not directly inferred.

   `HLA-B*58:01` was not interpreted because the available reduced VCF did not contain HLA typing data or the proxy marker `rs9263726`.

2. CYP2D6 copy number variation is not fully resolved.

   The current pipeline does not perform CNV calling or hybrid allele detection.

3. Star-allele phasing is not performed.

   If multiple star alleles are compatible with the observed markers, the result is reported with lower confidence.

4. Interpretation depends on the reduced VCF content.

   If a marker is absent from the subset VCF, it cannot be interpreted unless the full VCF or sequencing reads are available.

5. The pipeline is not a clinical decision support system.

   The output should be interpreted as a research prototype and not as a standalone clinical recommendation.

---

## Requirements

```text
Python 3.10+
pandas
```

Install dependencies:

```bash
pip install pandas
```

---

## Files not included

The following files are intentionally not included in the public repository:

```text
*.vcf
*.vcf.gz
*.tbi
*.fa
*.fna
*.fasta
*.fai
*.fastq
*.bam
*.cram
```

These files are large and may contain genomic or reference data. The GRCh38 reference FASTA is not needed in the repository because the rule table already stores target coordinates for GRCh38 / hg38.

---

## References

[1] Klein, T. E., & Ritchie, M. D. (2018). PharmCAT: A Pharmacogenomics Clinical Annotation Tool. *Clinical Pharmacology & Therapeutics*, 104(1), 19–22. https://doi.org/10.1002/cpt.928

[2] Clinical Pharmacogenetics Implementation Consortium (CPIC). https://www.clinpgx.org/cpic

[3] Dean L, Kane M. Allopurinol Therapy and HLA-B*58:01 Genotype. 2013 Mar 26 [Updated 2020 Dec 9]. In: Pratt VM, Scott SA, Pirmohamed M, et al., editors. *Medical Genetics Summaries* [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2012-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK127547/

---

## Summary

`Genomic Pharmacist` implements a custom pharmacogenomic interpretation workflow from VCF to drug recommendations.

The project demonstrates how manually curated pharmacogenetic knowledge can be integrated with targeted variant extraction, SNP genotype classification, star-allele interpretation, combined rule matching, and automated report generation without relying on PharmCAT.

Project developed by Daria Makarova and Anastasia Solovyeva (@nika-abl) 
