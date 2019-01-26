ALL.extend([PHANTOM_DIR + 'cross_correlation.txt',
            PHANTOM_DIR + 'cross_correlation_processed.txt'])
                    
rule calculateNSCRSC:
    input:
        expand('{directory}{sample}{extension}',
              directory = PHANTOM_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.phantom.txt')
    output:
        cc = PHANTOM_DIR + 'cross_correlation.txt',
        ccp = PHANTOM_DIR + 'cross_correlation_processed.txt',
    params:
        script = 'scripts/calculateNSCRSC.r',
        outdir = PHANTOM_DIR
    threads: config['thread_info']['phantompeakqual']
    shell:
        """
        cat {input} > {output.cc}
        {params.script} {output.cc} {params.outdir}
        """
