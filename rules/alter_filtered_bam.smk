ALL.extend([expand('{directory}{sample}{extension}',
                  directory = ALTERED_DIR,
                  sample = SAMPLE_NAMES,
                  extension = ['.filtered.bam','.filtered.sorted.bam','.filtered.sorted.bam.bai']),])
                    
rule alter_filtered_bam:
    input:
        MD_DIR + '{sample}.filtered.sorted.bam',
    output:
        bam = ALTERED_DIR + '{sample}.filtered.bam',
        sorted = ALTERED_DIR + '{sample}.filtered.sorted.bam',
        bai = ALTERED_DIR + '{sample}.filtered.sorted.bam.bai',
    params:
        tmp = ALTERED_DIR + '{sample}.tmp.bam'
    conda:
        'envs/samtools_env.yml'
    threads: config['thread_info']['index']
    shell:
        """
        samtools view \
        -@ {threads} \
        -h {input} | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > {params.tmp}
        mv {params.tmp} {output.bam}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        samtools index -@ {threads} -b {output.sorted}
        """
