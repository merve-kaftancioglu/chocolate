ALL.extend([expand('{directory}{sample}{extension}',
                  directory = ALN_DIR,
                  sample = SAMPLE_NAMES,
                  extension =['.bam',
                              '.sorted.bam',
                              '.sorted.bam.bai']),])

rule sort_index_bam:
    input:
        ALN_DIR + '{sample}.sam'
    output:
        bam = ALN_DIR + '{sample}.bam',
        sorted = ALN_DIR + '{sample}.sorted.bam',
        index = ALN_DIR + '{sample}.sorted.bam.bai',
    conda:
        'envs/samtools_env.yml'
    threads: config['thread_info']['index']
    benchmark:
        BENCH_DIR + 'sortIndexBam/{sample}.bamSortIndexBenchmark.txt'
    shell:
        """
        samtools view -@ {threads} -bS {input} > {output.bam}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        samtools index -@ {threads} -b {output.sorted}
        """
