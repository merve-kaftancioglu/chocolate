ALL.extend([expand('{directory}{sample}{extension}',
                  directory = TRIM_DIR,
                  sample = SAMPLE_NAMES,
                  extension = '.trimmed.fastq'),])
                    
rule trim:
    input:
        FASTQ_DIR + '{sample}.fastq',
    output:
        TRIM_DIR + '{sample}.trimmed.fastq',
    params:
        trimq = config['trim']['trimq'],
        minlen = config['trim']['minlen'],
        adaptor = config['trim']['adaptor']
    conda:
        'envs/bbmap_env.yml'
    threads: config['thread_info']['trim']
    log:
        LOG_DIR + 'trim/{sample}.trim.log'
    benchmark:
        BENCH_DIR + 'trim/{sample}.trimBenchmark.txt'
    shell:
        """
        bbduk.sh \
        in={input} \
        out={output} \
        ref={params.adaptor} \
        threads={threads} \
        k=19 \
        mink=5 \
        hdist=1 \
        hdist2=0 \
        ktrim=r \
        qtrim=r \
        minlength={params.minlen} \
        trimq={params.trimq} \
        2> {log}
        """
