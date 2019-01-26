ALL.extend([expand('{directory}{name}{extension}',
                  directory = DEEP_DIR,
                  name = SAMPLE_NAMES,
                  extension = '.bw'),])
                    
rule bamCoverage:
    input:
        MD_DIR + '{sample}.filtered.sorted.bam'
    output:
        DEEP_DIR + '{sample}.bw',
    params:
        extendReads = config['bamCoverage']['extendReads'],
        normalizeUsing = config['bamCoverage']['normalizeUsing'],
        binSize = config['bamCoverage']['binSize']
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        bamCoverage \
            --bam {input} \
            --normalizeUsing {params.normalizeUsing} \
            --extendReads {params.extendReads} \
            --binSize {params.binSize} \
            --numberOfProcessors {threads} \
            --outFileName {output}
        """
