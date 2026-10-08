rule np_index:
    input:
        fastq=lambda wildcards: SAMPLES[wildcards.sample]["path"]
    output:
        fastq=NP_DIR + "/{sample}/file.fastq",
        index=NP_DIR + "/{sample}/file.fastq.nailpolish.idx",
        log=NP_DIR + "/{sample}/np_index.log"
    params:
        nailpolish=config["nailpolish"],
        barcode_regex=lambda wildcards: SAMPLES[wildcards.sample]["barcode_regex"]
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60
    benchmark:
        BENCH_DIR + "/{sample}/np_index.benchmark.txt"
    shell:
        """
        ln -s {input.fastq} {output.fastq} 

        {params.nailpolish} index {output.fastq} --barcode-regex '{params.barcode_regex}' 2>&1 | tee {output.log}

         echo $SLURM_JOB_ID >> {output.log}
        """


rule np_summary:
    input:
        fastq=NP_DIR + "/{sample}/file.fastq",
        index=NP_DIR + "/{sample}/file.fastq.nailpolish.idx"
    output:
        summary=NP_DIR + "/{sample}/summary.html",
        stdout=NP_DIR + "/{sample}/summary_stdout.txt"
    params:
        nailpolish=config["nailpolish"]
    threads: 1
    resources:
        mem_mb=32000,
        runtime=15
    benchmark:
        BENCH_DIR + "/{sample}/np_summary.benchmark.txt"
    shell:
        """
        {params.nailpolish} summary {input.fastq} -o {output.summary} 2>&1 | tee {output.stdout}
        """


rule np_consensus:
    input:
        fastq=NP_DIR + "/{sample}/file.fastq",
        index=NP_DIR + "/{sample}/file.fastq.nailpolish.idx"
    output:
        consensus=NP_DIR + "/{sample}/consensus.fastq",
        log=NP_DIR + "/{sample}/np.log"
    params:
        nailpolish=config["nailpolish"],
        custom_params=lambda wildcards: SAMPLES[wildcards.sample].get("nailpolish_consensus_params", "")
    threads: config["threads"]
    resources:
        mem_mb=64000,
        runtime=2160,
        slurm_extra="'--nodelist=il-n[01-20]'"
    benchmark:
        BENCH_DIR + "/{sample}/np_consensus.benchmark.txt"
    shell:
        """
        {params.nailpolish} consensus {input.fastq} -t {threads} {params.custom_params} --report-original-reads --extra-stats -o {output.consensus} 2> >(tee {output.log} >&2) 
        echo $SLURM_JOB_ID >> {output.log}
        """


# Run Nailpolish without outputting the original reads
# Used for performance benchmarking
rule np_consensus_slim:
    input:
        fastq=NP_DIR + "/{sample}/file.fastq",
        index=NP_DIR + "/{sample}/file.fastq.nailpolish.idx"
    output:
        consensus=NP_DIR + "/{sample}/consensus_slim.fastq",
        log=NP_DIR + "/{sample}/np2.log"
    params:
        nailpolish=config["nailpolish"],
        custom_params=lambda wildcards: SAMPLES[wildcards.sample].get("nailpolish_consensus_params", "")
    threads: config["threads"]
    resources:
        mem_mb=32000,
        runtime=180,
        # run on il-n[01-20] nodes for benchmarking
        slurm_extra="'--nodelist=il-n[01-20]'"
    benchmark:
        BENCH_DIR + "/{sample}/np_consensus_slim.benchmark.txt"
    shell:
        """
        {params.nailpolish} consensus {input.fastq} -t {threads} {params.custom_params} -o {output.consensus} 2> >(tee {output.log} >&2) 
        echo $SLURM_JOB_ID >> {output.log}
        """

rule prepare_nailpolish_result:
    input: ALIGN_DIR + "/{sample}/called.sam"
    output: RESULTS_DIR + "/{sample}/nailpolish.sam"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=90
    shell:
        """
        source ./venv/bin/activate  

        python ./scripts/prepare_nailpolish_result.py {input} --output {output}
        """


# no clustering
rule np_consensus_no_cluster:
    input:
        fastq=NP_DIR + "/{sample}/file.fastq",
        index=NP_DIR + "/{sample}/file.fastq.nailpolish.idx"
    output:
        consensus=NP_DIR + "/{sample}/consensus_no_cluster.fastq",
        log=NP_DIR + "/{sample}/np_no_cluster.log"
    params:
        nailpolish=config["nailpolish"],
        custom_params=lambda wildcards: SAMPLES[wildcards.sample].get("nailpolish_consensus_params", "")
    threads: config["threads"]
    resources:
        mem_mb=64000,
        runtime=2160,
        slurm_extra="'--nodelist=il-n[01-20]'"
    benchmark:
        BENCH_DIR + "/{sample}/np_consensus_no_cluster.benchmark.txt"
    shell:
        """
        {params.nailpolish} consensus {input.fastq} -t {threads} {params.custom_params} --report-original-reads --extra-stats --no-clustering -o {output.consensus} 2> >(tee {output.log} >&2) 
        echo $SLURM_JOB_ID >> {output.log}
        """

rule minimap2_align_no_cluster:
    input:
        consensus=NP_DIR + "/{sample}/consensus_no_cluster.fastq",
        ref=lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        sam=ALIGN_DIR + "/{sample}/called_no_cluster.sam"
    threads: config["minimap_threads"]
    resources:
        mem_mb=64000,
        runtime=1440
    benchmark:
        BENCH_DIR + "/{sample}/minimap2_align_no_cluster.benchmark.txt"
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

rule prepare_nailpolish_result_no_cluster:
    input: ALIGN_DIR + "/{sample}/called_no_cluster.sam"
    output: RESULTS_DIR + "/{sample}/nailpolish_no_cluster.sam"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=90
    shell:
        """
        source ./venv/bin/activate  

        python ./scripts/prepare_nailpolish_result.py {input} --output {output}
        """