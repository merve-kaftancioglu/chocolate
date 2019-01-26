ALL.extend([DEEP_DIR + 'multiBamSummary.npz',])
                    
rule multiBamSummary:
    input:
        expand('{directory}{sample}{extension}',
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam')
    output:
        DEEP_DIR + 'multiBamSummary.npz',
    params:
        extendReads = config['multiBamSummary']['extendReads'],
        binSize = config['multiBamSummary']['binSize']
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        multiBamSummary \
            bins \
            --bamfiles {input} \
            --extendReads {params.extendReads} \
            --binSize {params.binSize} \
            --numberOfProcessors {threads} \
            --ignoreDuplicates \
            --centerReads \
            -o {output}
        """
