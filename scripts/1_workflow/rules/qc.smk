rule fastqc:
    input:
        fastq=lambda wildcards: ACTIVE_SAMPLES[wildcards.sample]["path"]
    output:
        directory(QC_DIR + "/{sample}")
    resources:
        mem_mb=32000,
        runtime=120
    threads: 1
    shell:
        """
        module load fastqc

        mkdir -p {output}
        fastqc --memory 9999 {input.fastq} -o {output}
        """
