ALL.extend([expand('{directory}{sample}{extension}',
                  directory = ALN_DIR,
                  sample = SAMPLE_NAMES,
                  extension = '.sam')])
                    
rule bwa_align:
    input:
        index = REF_DIR + ORG + '.bwt',
        R1 = TRIM_DIR + '{sample}_1.fastq',
        R2 = TRIM_DIR + '{sample}_2.fastq',
    output:
        sam = ALN_DIR + '{sample}.sam',
    params:
        org = REF_DIR + ORG
    conda:
        'envs/bwa_env.yml'
    threads: config['thread_info']['alignment']
    log:
        bwa = LOG_DIR + 'bwaAlign/{sample}.bwaAlign.log',
    benchmark:
        BENCH_DIR + 'bwaAlign/{sample}.bwaAlignBenchmark.txt'
    shell:
        """
        bwa mem \
        -M \
        -t {threads} \
        {params.org} \
        {input.R1} \
        {input.R2} \
        > {output.sam} \
        2> {log.bwa}
        """
