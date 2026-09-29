rule summarise_sam:
    input: RESULTS_DIR + "/{sample}/{process}.sam"
    output:
        csv=RESULTS_DIR + "/{sample}/{process}.csv",
        parquet=RESULTS_DIR + "/{sample}/{process}.parquet"
    threads: 1
    benchmark:
        BENCH_DIR + "/{sample}/summarise_{process}.benchmark.txt"
    resources:
        mem_mb=96000,
        runtime=120
    shell:
        """
        source ./venv/bin/activate 

        python ./scripts/summarise_sam.py {input} -o {output.csv} --parquet {output.parquet}
        """