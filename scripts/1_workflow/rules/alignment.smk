rule minimap2_align:
    input:
        consensus=NP_DIR + "/{sample}/consensus.fastq",
        ref=lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        sam=ALIGN_DIR + "/{sample}/called.sam"
    threads: config["minimap_threads"]
    resources:
        mem_mb=64000,
        runtime=1440
    benchmark:
        BENCH_DIR + "/{sample}/minimap2_align.benchmark.txt"
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


rule minimap2_align_slim:
    input:
        consensus_slim=NP_DIR + "/{sample}/consensus_slim.fastq",
        ref=lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        sam=ALIGN_DIR + "/{sample}/aligned_slim.sam"
    threads: config["minimap_threads"]
    resources:
        mem_mb=64000,
        runtime=1440
    benchmark:
        BENCH_DIR + "/{sample}/minimap2_align_slim.benchmark.txt"
    shell:
        """
        module load minimap2
        minimap2 -ax splice \
            -y \
            --secondary=no \
            -c \
            --MD \
            --eqx \
            -t {threads} {input.ref} {input.consensus_slim} > {output.sam}
        """

rule minimap2_align_original:
    input:
        fastq=lambda wildcards: SAMPLES[wildcards.sample]["path"],
        ref=lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        sam=ALIGN_DIR + "/{sample}/aligned_original.sam"
    threads: config["minimap_threads"]
    resources:
        mem_mb=64000,
        runtime=1440
    benchmark:
        BENCH_DIR + "/{sample}/minimap2_align_original.benchmark.txt"
    shell:
        """
        module load minimap2
        minimap2 -ax splice \
            -y \
            --secondary=no \
            -c \
            --MD \
            --eqx \
            -t {threads} {input.ref} {input.fastq} > {output.sam}
        """