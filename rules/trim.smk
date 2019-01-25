ALL.extend([expand('{directory}{sample}_{read}{extension}',
                  directory = TRIM_DIR,
                  sample = SAMPLE_NAMES,
                  read = ['1','2'],
                  extension = '.fastq'),])
                    
rule trim:
    input:
        R1 = FASTQ_DIR + '{sample}_R1.fastq',
        R2 = FASTQ_DIR + '{sample}_R2.fastq'
    output:
        R1 = TRIM_DIR + '{sample}_1.fastq',
        R2 = TRIM_DIR + '{sample}_2.fastq'
    params:
        min_overlap = config['trim_info']['min_overlap'],
        sample = TRIM_DIR + '{sample}',
        max_quality = config['trim_info']['max_quality']
    conda:
        'envs/ngmerge_env.yml'
    threads: config['thread_info']['trim']
    shell:
        """
        NGmerge \
        -1 {input.R1} \
        -2 {input.R2} \
        -n {threads} \
        -a \
        -e {params.min_overlap} \
        -o {params.sample} \
        -u {params.max_quality} \
        -v
        """
