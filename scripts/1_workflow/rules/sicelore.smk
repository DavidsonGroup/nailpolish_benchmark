# prepare a mapped .sam file for SiCeLoRe by adding
# tags (like CS) which would normally be inserted
# upstream by preprocessing steps
rule sicelore_prepare:
    input:
        sam=ALIGN_DIR + "/{sample}/aligned_original.sam"
    output:
        bam=SICELORE_DIR + "/{sample}/with_cs_tag.bam"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120
    benchmark:
        BENCH_DIR + "/{sample}/sicelore_prepare.benchmark.txt"
    shell:
        """
        ./scripts/prepare_for_sicelore.sh {input.sam} {output.bam}
        """

# run SiCeLoRe
rule sicelore_compute_consensus:
    input:
        bam=SICELORE_DIR + "/{sample}/with_cs_tag.bam"
    output:
        fastq=SICELORE_DIR + "/{sample}/sicelore.fastq"
    threads: config["threads"]
    resources:
        mem_mb=460000,
        runtime=1440*7,
        # temporary to run some datasets
        slurm_partition="long"
    log: LOG_DIR + "/{sample}/sicelore_compute_consensus.log"
    benchmark:
        BENCH_DIR + "/{sample}/sicelore_compute_consensus.benchmark.txt"
    shell:
        """
        # add spoa to PATH
        export PATH=./binary/spoa/build/bin:$PATH

        java -jar -Xmx400g binary/sicelore-2.1/Jar/Sicelore-2.1.jar ComputeConsensus \
            I={input.bam} \
            O={output.fastq} \
            T={threads} \
            CELLTAG=CB \
            UMITAG=UB \
            CDNATAG=CS \
            TMPDIR=/vast/scratch/users/cheng.o/tmp/   2>&1 | tee {log}
        """

# add the duplicate group status into the format expected by the
# benchmarking scripts (e.g. in the GG tag)
rule sicelore_assign_gg:
    input: SICELORE_DIR + "/{sample}/sicelore.fastq"
    output: SICELORE_DIR + "/{sample}/sicelore_with_gg_tag.fasta"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    benchmark:
        BENCH_DIR + "/{sample}/sicelore_assign_gg.benchmark.txt"
    shell:
        """
        ./scripts/sicelore_add_gg_tag.sh {input} {output}
        """

# map the consensus called reads for benchmarking
rule minimap2_align_sicelore:
    input:
        consensus=SICELORE_DIR + "/{sample}/sicelore_with_gg_tag.fasta",
        ref=lambda wildcards: ACTIVE_SAMPLES[wildcards.sample]["genome"]
    output:
        sam=RESULTS_DIR + "/{sample}/sicelore.sam"
    threads: config["minimap_threads"]
    resources:
        mem_mb=64000,
        runtime=480
    benchmark:
        BENCH_DIR + "/{sample}/minimap2_align_sicelore.benchmark.txt"
    shell:
        """
        module load minimap2
        minimap2 -ax splice \
            -y \
            --secondary=no \
            -c \
            --MD \
            --eqx \
            -t {threads} {input.ref} {input.consensus} > {output.sam}
        """