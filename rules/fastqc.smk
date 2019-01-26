ALL.extend([expand('{directory}{sample}{extension}',
                  directory = QC_DIR,
                  sample = SAMPLE_NAMES,
                  extension = '_fastqc.zip')])

rule fastqc:
    input:
        fastq = FASTQ_DIR + '{sample}.fastq',
    output:
        fastq = QC_DIR + '{sample}_fastqc.zip',
    params:
        qc_dir = QC_DIR,
        format = 'fastq',
    conda:
        'envs/fastqc_env.yml'
    threads: config['thread_info']['fastqc']
    shell:
        """
        mkdir -p {params.qc_dir}
        fastqc \
        --noextract \
        -f {params.format} \
        -o {params.qc_dir} \
        -t {threads} \
        {input}
        """
