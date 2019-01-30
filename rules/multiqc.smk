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
        expand('{directory}{sample}{extension}', #bam 
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.bam'),
        expand('{directory}{file}', #deeptools
              directory = DEEP_DIR,
              file = ['pcaplot_multiBamSummary.png', 
              'heatmap_SpearmanCorr_multiBamSummary.png', 
              'scatterplot_PearsonCorr_multiBamSummary.png']), 
        expand('{directory}{name}{extension}',
                          directory = DEEP_DIR,
                          name = 'fingerprints',
                          extension = ['.pdf','.txt']),
        expand('{directory}{sample}{extension}', #phantompeakqual
              directory = PHANTOM_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.phantom.txt') \
              if RUN_PHANTOM.upper() == 'Y' \
              else \
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
        expand('{directory}{sample}{extension}', #bam 
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.bam'),
        expand('{directory}{file}', #deeptools
              directory = DEEP_DIR,
              file = ['pcaplot_multiBamSummary.png', 
              'heatmap_SpearmanCorr_multiBamSummary.png', 
              'scatterplot_PearsonCorr_multiBamSummary.png']), 
        expand('{directory}{name}{extension}',
                          directory = DEEP_DIR,
                          name = 'fingerprints',
                          extension = ['.pdf','.txt']),
         
    output:
        QC_DIR + 'multiqc_report.html',
    params:
        QC_DIR,
        BIOTYPE_DIR,
        MULTI_DIR,
        MD_DIR,
        DEEP_DIR,
        PHANTOM_DIR \
        if RUN_PHANTOM.upper() == 'Y' \
        else \
        QC_DIR,
        BIOTYPE_DIR,
        MULTI_DIR,
        MD_DIR,
        DEEP_DIR,
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
