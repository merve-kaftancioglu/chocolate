ALL.extend([expand('{directory}{sample}{extension}',
                  directory = MD_DIR,
                  sample = SAMPLE_NAMES,
                  extension = ['.filtered.bam','.filtered.bam.bai']),])
                    
rule filter_bam:
    input:
        bam = MD_DIR + '{sample}.mkdp.bam',
        bai = MD_DIR + '{sample}.mkdp.bai',
    output:
        bam = MD_DIR + '{sample}.filtered.bam',
        bai = MD_DIR + '{sample}.filtered.bam.bai',
    params:
        mapquality = config['filter_bam']['mapquality'],
        tmp = MD_DIR + '{sample}_filter'
    conda:
        'envs/samtools_env.yml'
    threads: config['thread_info']['index']
    benchmark:
        BENCH_DIR + 'filter_bam/{sample}.filter_bamBenchmark.txt',
    shell:
        """
        samtools view \
        -@ {threads} \
        -b \
        -h \
        -f 3 \
        -F 4 \
        -F 8 \
        -F 256 \
        -F 1024 \
        -F 2048 \
        -q {params.mapquality} \
        {input.bam} \
        > {output.bam}
        rm -rf {params.tmp}
        samtools index -@ {threads} -b {output.bam}
        """
