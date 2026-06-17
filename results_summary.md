# Results Summary

## Sample Information

* Sample analyzed: SRR13106578
* Data type: Whole-exome sequencing
* Reference genome: GRCh38 / Homo_sapiens_assembly38
* Project focus: CNS inflammatory demyelinating disease related to Chikungunya virus infection

## Read Processing

Raw paired-end reads were processed using Trimmomatic to remove adapters and low-quality bases.

* Trimmed paired reads: 68,413,384 read pairs
* Forward-only surviving reads: 89,050
* Reverse-only surviving reads: 2,113
* Dropped reads: 11,588

## Alignment Summary

The trimmed reads were aligned to the GRCh38 human reference genome using BWA-MEM.

* Total reads: 136,876,156
* Mapped reads: 135,774,740
* Mapping rate: 99.20%
* Properly paired reads: 96.83%

## Duplicate Marking

Duplicate reads were marked using GATK MarkDuplicates.

* Duplicates marked: 17,000,416

## Variant Calling and Filtering

Variants were called using GATK HaplotypeCaller and filtered using GATK VariantFiltration.

* Raw variants: 174,933
* PASS variants: 104,972

GATK VariantFiltration labels variants based on filter criteria. A separate PASS-only VCF file was created to keep variants that passed the selected filters.

## VEP Annotation

PASS-filtered variants were annotated using Ensembl VEP with the GRCh38 offline cache. The VEP `--pick` option was used to keep one representative annotation per variant.

* VEP picked annotation records: 98,546

## VEP Impact Summary

* HIGH impact: 334
* MODERATE impact: 9,715
* LOW impact: 15,551
* MODIFIER impact: 72,946

## Priority Variant Summary

Priority variants were selected from protein-coding variants with predicted functional relevance.

The priority variant table contained 2,649 variants.

Main consequence categories:

* missense_variant: 2,314
* frameshift_variant: 132
* stop_gained: 55
* splice_donor_variant: 21
* splice_acceptor_variant: 20
* start_lost: 11
* stop_lost: 7

## Interpretation Note

This project was completed using one WES sample. Therefore, the results should be viewed as variant annotation and prioritization, not as disease association or clinical interpretation.

A larger case-control analysis would be needed to make stronger biological conclusions.
