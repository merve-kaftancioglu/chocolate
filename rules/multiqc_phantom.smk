ALL.extend([QC_DIR + 'multiqc_report.html',])

rule multiqc_phantom:
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
              extension = '.filtered.bam'),
        expand('{directory}{sample}{extension}', #phantompeakqual
              directory = PHANTOM_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.phantom.txt')
    output:
        QC_DIR + 'multiqc_report.html',
    params:
        QC_DIR,
        BIOTYPE_DIR,
        MULTI_DIR,
        PEAKS_DIR,
        MD_DIR,
        PHANTOM_DIR
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
