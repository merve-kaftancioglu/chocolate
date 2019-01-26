ALL.extend([expand('{directory}{sample}{extension}',
                  directory = FASTQ_DIR,
                  sample = SAMPLE_NAMES,
                  extension = '.fastq')])       

rule get_fastq: 
    output:
        FASTQ_DIR + 'Sample_{sample}.fastq',
    params:
        fastq_dir= config['input_dir']['fastq']
    shell:
        """
        zcat {params.fastq_dir}Sample_{wildcards.sample}/{wildcards.sample}_*_001.fastq.gz > {output}
        """
