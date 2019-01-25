ALL.extend([QC_DIR + 'multiqc_report.html',])

rule multiqc:
    input:
        expand('{directory}{sample}_{read}{extension}', #fastqc raw
              directory = QC_DIR,
              sample = SAMPLE_NAMES,
              read = READS,
              extension = '_fastqc.zip'),
        expand('{directory}{sample}_{read}{extension}', #fastq screen biotype
              directory = BIOTYPE_DIR,
              sample = SAMPLE_NAMES,
              read = READS,
              extension = ['_screen.html','_screen.txt']),
        expand('{directory}{sample}_{read}{extension}', #fastqc screen multispecies
              directory = MULTI_DIR,
              sample = SAMPLE_NAMES,
              read = READS,
              extension = ['_screen.html','_screen.txt']),
        expand('{directory}{sample}{extension}', #macs2
              directory = PEAKS_DIR,
              sample = PE_BED_NAMES,
              extension = '_peaks.broadPeak'),
        expand('{directory}{sample}{extension}', #bam 
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.bam')
    output:
        QC_DIR + 'multiqc_report.html',
    params:
        QC_DIR,
        BIOTYPE_DIR,
        MULTI_DIR,
        PEAKS_DIR,
        MD_DIR
    conda:
        'envs/multiqc_env.yml'
    shell:
        """
        multiqc \
        --force \
        -s \
        --interactive \
        -n {output} \
        {params}
        """
