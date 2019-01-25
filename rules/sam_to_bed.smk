ALL.extend([expand('{directory}{sample}{extension}',
                  directory = MD_DIR,
                  sample = SAMPLE_NAMES,
                  extension ='.filtered.bed'),])

rule sam_to_bed:
    input:
        MD_DIR + '{sample}.filtered.bam',
    output:
        MD_DIR + '{sample}.filtered.bed',
    params:
        script = 'scripts/SAMtoBED2.py'
    conda:
        'envs/samtools_env.yml'
    threads: config['thread_info']['index']
    shell:
        """
        samtools view -@ {threads} -h {input} | python {params.script} -i - -o {output} -x -v -e
        """
