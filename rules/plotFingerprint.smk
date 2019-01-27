ALL.extend([expand('{directory}{name}{extension}',
                  directory = DEEP_DIR,
                  name = 'fingerprints',
                  extension = ['.pdf','.txt']),])
                    
rule plotFingerprint:
    input:
        expand('{directory}{sample}{extension}',
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam')
    output:
        pdf = DEEP_DIR + 'fingerprints.pdf',
        counts = DEEP_DIR + 'fingerprints.txt'
    params:
        extendReads = config['plotFingerprint']['extendReads'],
        numberOfSamples = config['plotFingerprint']['numberOfSamples'],
        binSize = config['plotFingerprint']['binSize'],
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        plotFingerprint \
        --bamfiles {input} \
        --plotFile {output.pdf} \
        --outRawCounts {output.counts} \
        --extendReads {params.extendReads} \
        --skipZeros \
        --ignoreDuplicates \
        --numberOfSamples {params.numberOfSamples} \
        --binSize {params.binSize} \
        --plotFileFormat pdf \
        --numberOfProcessors {threads} \
        --smartLabels \
        --plotTitle "Fingerprints"
        """
