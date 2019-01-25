ALL.extend([expand('{directory}{sample}_{read}{extension}',
                  directory = FASTQ_DIR,
                  sample = SAMPLE_NAMES,
                  read = READS,
                  extension = '.fastq')])       

rule get_fastq: 
    output:
        FASTQ_DIR + 'Sample_{sample}_{read}.fastq',
    params:
        fastq_dir= config['input_dir']['fastq']
    shell:
        """
        zcat {params.fastq_dir}Sample_{wildcards.sample}/{wildcards.sample}_*_{wildcards.read}_001.fastq.gz > {output}
        """
