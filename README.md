# Genomic Pharmacist

Custom pharmacogenomic interpretation pipeline for VCF-based analysis of drug response markers.

## Project task

The goal of this project was to analyze genetic markers that influence individual response to selected drugs. These markers may include single nucleotide polymorphisms (SNPs), star alleles, HLA haplotypes, and combined pharmacogenetic rules.

The task required implementing a custom interpretation pipeline without using PharmCAT, a commonly used pharmacogenomic interpretation tool. The pipeline takes a VCF file as input, extracts pharmacogenetically relevant variants, determines SNP genotypes and star-allele diplotypes, maps them to metabolizer phenotypes, and generates drug-specific recommendations.

Pharmacogenomic rules were manually curated using publicly available resources such as ClinPGx, CPIC/PharmGKB-related annotations, and drug-specific pharmacogenetic guidelines.

## Input data

The project was designed for 1000 Genomes samples aligned to the hg38 genome build.

Two datasets were provided:

1. `Pharma_subset.vcf.gz`  
   A reduced VCF containing 10 samples and approximately 38,000 variants relevant for pharmacological profiling.


The current implementation focuses on the reduced pharmacogenomic subset VCF.

## Supported drug examples

The pipeline was developed for pharmacogenomic interpretation of drugs such as:

- Sertraline
- Rosuvastatin
- Efavirenz
- Clopidogrel
- Voriconazole
- Simvastatin
- Tacrolimus
- Warfarin
- Isoniazid
- Ivacaftor

Allopurinol was considered through the `HLA-B*58:01` association. However, since the reduced VCF did not contain direct HLA typing data or the proxy marker `rs9263726`, this rule was excluded from the active interpretation workflow.

## Pipeline architecture

The pipeline is divided into five logical modules:

```text
VCFExtractor
    ↓
GenotypeCaller
    ↓
DiplotypeCaller
    ↓
RulesEngine
    ↓
Reporter
