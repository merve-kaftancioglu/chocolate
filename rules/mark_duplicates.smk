ALL.extend([expand('{directory}{sample}{extension}',
                  directory = MD_DIR,
                  sample = SAMPLE_NAMES,
                  extension = ['.mkdp.bam','.mkdp.bai', '_dup_metrics.txt']),])
                    
rule mark_duplicates:
    input:
        ALN_DIR + '{sample}.sorted.bam',
    output:
        bam = MD_DIR + '{sample}.mkdp.bam',
        bai = MD_DIR + '{sample}.mkdp.bai',
        metric = MD_DIR + '{sample}_dup_metrics.txt'
    params:
        assumesorted = config['mark_duplicates']['assumesorted'],
        validationstringency = config['mark_duplicates']['validationstringency'],
        createindex = config['mark_duplicates']['createindex'],
        tmp = MD_DIR + '{sample}_mkdp'
    conda:
        'envs/picard_env.yml'
    benchmark:
        BENCH_DIR + 'mark_duplicates/{sample}.mark_duplicatesBenchmark.txt',
    shell:
        """
        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        ASSUME_SORTED={params.assumesorted} \
        METRICS_FILE={output.metric} \
        CREATE_INDEX={params.createindex} \
        VALIDATION_STRINGENCY={params.validationstringency} \
        MAX_FILE_HANDLES=1000 \
        TMP_DIR={params.tmp}
        rm -rf {params.tmp}
        """
