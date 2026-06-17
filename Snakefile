import csv

configfile: "config/config.yaml"

with open(config["sample_table"]) as f:
    reader = csv.DictReader(f, delimiter="\t")
    samples = {row["sample"]: row for row in reader}

SAMPLES = list(samples.keys())
THREADS = config["threads"]
REF = config["reference"]
VEP_CACHE = config.get("vep_cache", "vep_cache")

rule all:
    input:
        expand("qc_results/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("qc_results/{sample}_2_fastqc.html", sample=SAMPLES),
        expand("variants/{sample}_pass_only.vcf", sample=SAMPLES),
        expand("annotation/{sample}_PASS_vep_pick.tsv", sample=SAMPLES),
        expand("interpretation/{sample}_priority_variants.tsv", sample=SAMPLES)

rule fastqc:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"]
    output:
        html1="qc_results/{sample}_1_fastqc.html",
        html2="qc_results/{sample}_2_fastqc.html"
    shell:
        """
        mkdir -p qc_results
        fastqc {input.r1} {input.r2} -o qc_results
        """

rule trim_reads:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"]
    output:
        r1_paired="trimmed_reads/{sample}_R1_paired.fastq.gz",
        r1_unpaired="trimmed_reads/{sample}_R1_unpaired.fastq.gz",
        r2_paired="trimmed_reads/{sample}_R2_paired.fastq.gz",
        r2_unpaired="trimmed_reads/{sample}_R2_unpaired.fastq.gz"
    threads: THREADS
    shell:
        """
        mkdir -p trimmed_reads logs

        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.r1_paired} {output.r1_unpaired} \
          {output.r2_paired} {output.r2_unpaired} \
          ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
          LEADING:3 TRAILING:3 MINLEN:36 \
          2>&1 | tee logs/trimmomatic_{wildcards.sample}.log
        """

rule align_reads:
    input:
        r1="trimmed_reads/{sample}_R1_paired.fastq.gz",
        r2="trimmed_reads/{sample}_R2_paired.fastq.gz",
        ref=REF
    output:
        bam="alignment/{sample}_aligned_sorted.bam",
        bai="alignment/{sample}_aligned_sorted.bam.bai"
    threads: THREADS
    shell:
        """
        mkdir -p alignment logs tmp

        bwa mem -t {threads} \
          -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1\\tPU:{wildcards.sample}" \
          {input.ref} {input.r1} {input.r2} \
          2> logs/bwa_mem_{wildcards.sample}.log \
          | samtools sort -@ {threads} -m 1G \
              -T tmp/{wildcards.sample}_sorttmp \
              -o {output.bam}

        samtools index {output.bam}
        """

rule mark_duplicates:
    input:
        bam="alignment/{sample}_aligned_sorted.bam"
    output:
        bam="alignment/{sample}_deduped.bam",
        bai="alignment/{sample}_deduped.bam.bai",
        metrics="alignment/{sample}_duplicate_metrics.txt"
    shell:
        """
        mkdir -p alignment logs

        gatk MarkDuplicates \
          -I {input.bam} \
          -O {output.bam} \
          -M {output.metrics} \
          2>&1 | tee logs/markduplicates_{wildcards.sample}.log

        samtools index {output.bam}
        """

rule haplotypecaller:
    input:
        bam="alignment/{sample}_deduped.bam",
        bai="alignment/{sample}_deduped.bam.bai",
        ref=REF
    output:
        vcf="variants/{sample}_raw_variants.vcf"
    threads: THREADS
    shell:
        """
        mkdir -p variants logs

        gatk HaplotypeCaller \
          -R {input.ref} \
          -I {input.bam} \
          -O {output.vcf} \
          --native-pair-hmm-threads {threads} \
          2>&1 | tee logs/haplotypecaller_{wildcards.sample}.log
        """

rule filter_variants:
    input:
        vcf="variants/{sample}_raw_variants.vcf",
        ref=REF
    output:
        vcf="variants/{sample}_filtered_variants.vcf"
    shell:
        """
        mkdir -p variants logs

        gatk VariantFiltration \
          -R {input.ref} \
          -V {input.vcf} \
          -O {output.vcf} \
          --filter-expression "QD < 2.0" --filter-name "QD2" \
          --filter-expression "DP < 10" --filter-name "LowDepth" \
          --filter-expression "MQ < 40.0" --filter-name "MQ40" \
          2>&1 | tee logs/variantfiltration_{wildcards.sample}.log
        """

rule pass_only:
    input:
        vcf="variants/{sample}_filtered_variants.vcf"
    output:
        vcf="variants/{sample}_pass_only.vcf"
    shell:
        """
        grep "^#" {input.vcf} > {output.vcf}
        grep -v "^#" {input.vcf} | awk '$7=="PASS"' >> {output.vcf}
        """

rule vep_annotation_pick:
    input:
        vcf="variants/{sample}_pass_only.vcf",
        ref=REF
    output:
        tsv="annotation/{sample}_PASS_vep_pick.tsv",
        html="annotation/{sample}_PASS_vep_pick_summary.html"
    threads: THREADS
    shell:
        """
        mkdir -p annotation logs

        vep \
          --input_file {input.vcf} \
          --output_file {output.tsv} \
          --tab \
          --cache \
          --offline \
          --dir_cache {VEP_CACHE} \
          --assembly GRCh38 \
          --species homo_sapiens \
          --fasta {input.ref} \
          --symbol \
          --canonical \
          --biotype \
          --hgvs \
          --sift b \
          --polyphen b \
          --check_existing \
          --pick \
          --pick_order canonical,appris,tsl,biotype,rank,ccds,length \
          --force_overwrite \
          --fork {threads} \
          --stats_file {output.html} \
          2>&1 | tee logs/vep_PASS_pick_{wildcards.sample}.log
        """

rule priority_variants:
    input:
        tsv="annotation/{sample}_PASS_vep_pick.tsv"
    output:
        tsv="interpretation/{sample}_priority_variants.tsv"
    shell:
        """
        mkdir -p interpretation

        {{
        echo -e "Variant\tLocation\tGene\tEnsembl_Gene\tConsequence\tImpact\tBiotype\tCanonical\tSIFT\tPolyPhen\tHGVSc\tHGVSp\tKnown_ID\tClinical_Significance"

        awk -F'\t' 'BEGIN{{OFS="\t"}}
        /^#/ {{next}}
        $21=="protein_coding" && ($14=="HIGH" || $7 ~ /frameshift_variant/ || $7 ~ /stop_gained/ || $7 ~ /splice_acceptor_variant/ || $7 ~ /splice_donor_variant/ || ($7 ~ /missense_variant/ && ($23 ~ /deleterious/ || $24 ~ /damaging/))) {{
          print $1,$2,$18,$4,$7,$14,$21,$22,$23,$24,$25,$26,$13,$28
        }}' {input.tsv}
        }} > {output.tsv}
        """

