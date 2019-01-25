ALL.extend([expand('{directory}{sample}_{read}{extension}',
                  directory = MULTI_DIR,
                  sample = SAMPLE_NAMES,
                  read = READS,
                  extension = ['_screen.html','_screen.txt'])])

rule fastq_screen_multi:
    input:
        fastq = FASTQ_DIR + '{sample}_{read}.fastq',
    output:
        html = MULTI_DIR + '{sample}_{read}_screen.html',
        txt = MULTI_DIR + '{sample}_{read}_screen.txt',
    params:
        output_dir = MULTI_DIR,
        subset = config['fastq_screen']['subset'],
        aligner = config['fastq_screen']['aligner'],
        config_file = config['fastq_screen']['multi_species']
    conda:
        'envs/fastqscreen_env.yml'
    threads: config['thread_info']['fastqc']
    shell:
        """
        fastq_screen \
        --aligner {params.aligner} \
        --subset {params.subset} \
        --conf {params.config_file} \
        --outdir {params.output_dir} \
        --threads {threads} \
        {input}
        """
