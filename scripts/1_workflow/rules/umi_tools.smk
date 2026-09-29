
# ================
# sort + run bedtools
# ================

rule per_cell_sort:
    input:
        sam=ALIGN_DIR + "/{sample}/aligned_original.sam"
    output:
        bam=UMITOOLS_DIR + "/{sample}/per_cell.sorted.bam"
    threads: config["threads"]
    resources:
        mem_mb=32000,
        runtime=60
    benchmark:
        BENCH_DIR + "/{sample}/sort_per_cell.benchmark.txt"
    shell:
        """
        module load samtools
        samtools sort {input.sam} -o {output.bam} -@ {threads}
        samtools index {output.bam}
        """

rule bedtools_intersect:
    input:
        sam=UMITOOLS_DIR + "/{sample}/per_cell.sorted.bam"
    output:
        bed=UMITOOLS_DIR + "/{sample}/overlaps.bed",
        bam=UMITOOLS_DIR + "/{sample}/per_gene.bedtools.sorted.bam"
    resources:
        mem_mb=64000,
        runtime=240
    benchmark:
        BENCH_DIR + "/{sample}/bedtools_intersect.benchmark.txt"
    params:
        gtf=lambda wildcards: SAMPLES[wildcards.sample]["gene_gtf"]
    shell:
        """
        module load bedtools
        module load samtools

        bedtools intersect -a {input.sam} -b {params.gtf} -wo -sorted -bed \
            > {output.bed}
            
        source /vast/projects/davidson_longread/cheng.o/proj/25.06.17-np-benchmark/snakemake/venv/bin/activate
        python ./scripts/add_gene_after_bamtools.py {output.bed} {input.sam} {output.bam}

        samtools index {output.bam}
        """

# ================
# run UMI-tools
# Note: we run both dedup and group
# Performance benchmarking is based on the `dedup` performance
# Both commands are run with the same input parameters
#
# Reasoning:
# umi_tools dedup will produce a representative read per group but
# has no options to 1) show/assign the discarded reads and 2) tag
# whether a representative is from a singleton or duplicate group.
# this limits comparisons for benchmarking. To fix this, we run
# umi-tools group as well, which performs no deduplication, to see
# each read's assigned group.
# ================

rule run_umi_tools_bedtools_dedup:
    input:
        mapped=UMITOOLS_DIR + "/{sample}/per_gene.bedtools.sorted.bam"
    output:
        bam=UMITOOLS_DIR + "/{sample}/umitools_bedtools.bam"
    threads: 1
    resources:
        mem_mb=64000,
        runtime=60*24
    log:
        LOG_DIR + "/umi_tools/{sample}_bedtools_dedup.log"
    benchmark:
        BENCH_DIR + "/{sample}/umi_tools_bedtools_dedup.benchmark.txt"
    shell:
        """
        source /vast/projects/davidson_longread/cheng.o/proj/25.06.17-np-benchmark/snakemake/venv/bin/activate

        umi_tools dedup \
            --stdin={input.mapped} \
            --stdout={output.bam} \
            --temp-dir="/vast/scratch/users/cheng.o/tmp" \
            --extract-umi-method=tag \
            --umi-tag=UB \
            --cell-tag=CB \
            --per-cell \
            --per-gene \
            --gene-tag=XT \
            --assigned-status-tag=XS \
            --no-sort-output 2>&1 | tee {log}
        """

rule run_umi_tools_bedtools_group:
    input:
        mapped=UMITOOLS_DIR + "/{sample}/per_gene.bedtools.sorted.bam"
    output:
        tsv=UMITOOLS_DIR + "/{sample}/umitools_bedtools.tsv"
    threads: 1
    resources:
        mem_mb=64000,
        runtime=60*24
    log:
        LOG_DIR + "/umi_tools/{sample}_bedtools_group.log"
    benchmark:
        BENCH_DIR + "/{sample}/umi_tools_bedtools_group.benchmark.txt"
    shell:
        """
        source /vast/projects/davidson_longread/cheng.o/proj/25.06.17-np-benchmark/snakemake/venv/bin/activate

        umi_tools group \
            -I {input.mapped} \
            --group-out={output.tsv} \
            --temp-dir="/vast/scratch/users/cheng.o/tmp" \
            --extract-umi-method=tag \
            --umi-tag=UB \
            --cell-tag=CB \
            --per-cell \
            --per-gene \
            --gene-tag=XT \
            --assigned-status-tag=XS \
            --no-sort-output 2>&1 | tee {log}
        """

# Assign the duplicate group status to GG by using the output of
# umi_tools dedup and umi_tools group
rule prepare_umitools_result:
    input:
        bam=UMITOOLS_DIR + "/{sample}/umitools_{t}.bam",
        tsv=UMITOOLS_DIR + "/{sample}/umitools_{t}.tsv"
    output: RESULTS_DIR + "/{sample}/umitools_{t}.sam"
    threads: 1
    resources:
        mem_mb=32000,
        runtime=90
    shell:
        """
        source ./venv/bin/activate

        python ./scripts/prepare_umitools_result.py {input.bam} {input.tsv} -o {output}
        """