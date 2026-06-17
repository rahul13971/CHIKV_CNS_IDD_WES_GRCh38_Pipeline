# Whole-Exome Variant Calling and Annotation Pipeline for a Neuroimmunology WES Sample

## Project overview

I built a complete whole-exome sequencing (WES) variant calling workflow using public sequencing data related to neuroimmunology. The sample used in this version of the project is SRR13106578, which comes from a study focused on central nervous system inflammatory demyelinating disease following Chikungunya virus infection.

The main goal was not to make clinical claims, but to build a clear and reproducible pipeline starting from raw FASTQ files and ending with annotated and prioritized variants. I used this project to understand how raw sequencing reads are processed, aligned to a human reference genome, converted into variant calls, filtered, and then annotated for biological interpretation.

## Dataset

- Sample analyzed: SRR13106578
- Data type: Whole-exome sequencing
- Reference genome: GRCh38 / Homo_sapiens_assembly38
- Disease context: CNS inflammatory demyelinating disease associated with Chikungunya virus infection

Only one sample was processed in this version of the project. Because of that, the results are interpreted as variant annotation and prioritization, not as disease association or biomarker discovery.

## Tools used

The main tools used in this workflow were:

- FastQC for raw read quality assessment
- Trimmomatic for adapter and quality trimming
- BWA-MEM for read alignment
- Samtools for BAM sorting, indexing, and alignment statistics
- GATK MarkDuplicates for duplicate marking
- GATK HaplotypeCaller for variant calling
- GATK VariantFiltration for hard filtering
- Ensembl VEP for variant annotation
- Snakemake for workflow organization
- Conda and WSL Ubuntu for the computing environment

## Workflow summary

The pipeline follows these major steps:

1. Raw FASTQ quality check
2. Read trimming
3. Alignment to the GRCh38 human reference genome
4. BAM sorting and indexing
5. Duplicate marking
6. Variant calling
7. Variant filtering
8. PASS-only VCF generation
9. VEP annotation
10. Functional variant prioritization

A Snakemake workflow was also created to organize the core steps of the pipeline and make the analysis easier to reproduce.

## Main results

The workflow successfully processed the WES sample from raw reads to annotated variants.

Key results from the pipeline:

- Trimmed paired reads: 68,413,384 read pairs
- Mapping rate: 99.20%
- Properly paired reads: 96.83%
- Duplicates marked: 17,000,416
- Raw variants: 174,933
- PASS variants: 104,972
- VEP picked annotation records: 98,546

## Variant annotation

Variants that passed filtering were annotated using Ensembl VEP with the GRCh38 cache. I first generated a full VEP annotation file, then used the VEP `--pick` option to keep one representative annotation per variant. This made the output easier to summarize and interpret.

The VEP `--pick` annotation produced the following impact summary:

- HIGH impact: 334
- MODERATE impact: 9,715
- LOW impact: 15,551
- MODIFIER impact: 72,946

For downstream interpretation, I focused mainly on HIGH and MODERATE impact variants, protein-coding variants, and variant consequences that are more likely to affect gene or protein function.

## Priority variant filtering

I created a priority variant table by focusing on protein-coding variants with potentially important functional effects. This included:

- HIGH impact variants
- frameshift variants
- stop-gained variants
- splice donor variants
- splice acceptor variants
- missense variants predicted to be damaging by SIFT or PolyPhen

The final priority variant table contained 2,649 variants.

Major consequence categories in the priority table included:

- missense_variant: 2,314
- frameshift_variant: 132
- stop_gained: 55
- splice_donor_variant: 21
- splice_acceptor_variant: 20
- start_lost: 11
- stop_lost: 7

## Biological interpretation

Most of the prioritized variants were missense variants, which means they may change amino acids in proteins. A smaller number of variants had stronger predicted effects, such as frameshift, stop-gained, and splice-site consequences. These classes are important because they can disrupt protein structure, shorten proteins, or affect RNA splicing.

Some genes appeared multiple times in the priority table, including PABPC3, FCGBP, AHNAK2, RBMX, PABPC1, IGFN1, FLG, SLC35G4, MUC16, ENTPD2, and TYRO3. These genes were treated as candidates for preliminary review only.

Since this version of the project uses one WES sample, I did not treat these findings as disease-causing or disease-associated variants. A larger case-control analysis would be needed to make stronger biological or clinical conclusions.

## Snakemake workflow

A Snakemake workflow was created for the core variant calling steps. This helped organize the analysis and made the pipeline easier to rerun or extend later.

The workflow currently tracks:

- FastQC output
- trimmed FASTQ files
- aligned and sorted BAM files
- duplicate-marked BAM files
- raw VCF files
- filtered VCF files
- PASS-only VCF files

## Limitations

This project has a few important limitations:

- Only one WES sample was processed.
- No case-control comparison was performed.
- The interpretation is based on predicted variant effects, not experimental validation.
- The results should not be used for clinical decision-making.
- Additional samples would be needed for disease association, variant burden testing, or pathway-level analysis.

## Future improvements

Possible next steps include:

- Processing additional case and control samples
- Adding VEP annotation directly into Snakemake
- Filtering against immune and neuroinflammatory gene lists
- Performing pathway enrichment analysis
- Creating summary plots for variant consequences and impact categories
- Comparing variant patterns across multiple samples

## Project status

Version 1 of this project is complete. It includes a working WES variant calling pipeline, GRCh38-based alignment, GATK variant calling, VEP annotation, variant prioritization, and a Snakemake workflow for reproducibility.

This project demonstrates my ability to build and organize a bioinformatics pipeline using commonly used tools in genomics and variant analysis.
