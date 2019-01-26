ALL.extend([expand('{directory}{sample}{extension}',
                  directory = MD_DIR,
                  sample = SAMPLE_NAMES,
                  extension = ['.filtered.sorted.bam','.filtered.sorted.bam.bai']),])
                    
rule filter_bam:
    input:
        bam = MD_DIR + '{sample}.mkdp.bam',
        bai = MD_DIR + '{sample}.mkdp.bai',
    output:
        bam = MD_DIR + '{sample}.filtered.bam',
        sorted = ALN_DIR + '{sample}.filtered.sorted.bam',
        bai = MD_DIR + '{sample}.filtered.sorted.bam.bai',
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
        -q 1 \
        -F 4 \
        -F 256 \
        -F 1024 \
        -q {params.mapquality} \
        {input.bam} \
        > {output.bam}
        rm -rf {params.tmp}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        samtools index -@ {threads} -b {output.sorted}
        """
