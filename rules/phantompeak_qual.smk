ALL.extend([expand('{directory}{sample}{extension}',
                  directory = PHANTOM_DIR,
                  sample = SAMPLE_NAMES,
                  extension = '.filtered.sorted.phantom.txt'),])
                    
rule phantompeak_qual:
    input:
        bam = MD_DIR + '{sample}.filtered.sorted.bam',
        bai = MD_DIR + '{sample}.filtered.sorted.bam.bai',
    output:
        PHANTOM_DIR + '{sample}.filtered.sorted.phantom.txt',
    params:
        outdir = PHANTOM_DIR,
        tmp = PHANTOM_DIR + '{sample}.tmp'
    conda:
        'envs/phantompeakqualtools_env.yml'
    threads: config['thread_info']['phantompeakqual']
    shell:
        """
        run_spp.R \
        -c={input.bam} \
        -savp \
        -rf \
        -p={threads} \
        -odir={params.outdir} \
        -out={output}
        """
