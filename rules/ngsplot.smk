ALL.extend([NGS_DIR + 'Genebody.pdf',
            NGS_DIR + 'TSS.pdf',])
                    
rule ngsplot:
    input:
        expand('{directory}{sample}{extension}',
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam')
    output:
        ngs_conf = NGS_DIR + 'ngsplot_config',
        gb = NGS_DIR + 'Genebody.pdf',
        tss = NGS_DIR + 'TSS.pdf',
    params:
        generate = 'scripts/ngsplot_config.r',
        outdir = NGS_DIR,
        gb = NGS_DIR + 'Genebody',
        tss = NGS_DIR + 'TSS',
        genome = config['ngsplot']['genome']
    conda:
        'envs/ngsplot_env.yml'
    shell:
        """
        ngs_config_generate.r {params.outdir} {input}
        ngs.plot.r \
        -G {params.genome} \
        -R genebody \
        -C {output.ngs_conf} \
        -O {params.gb} \
        -D ensembl \
        -FL 300
        ngs.plot.r \
        -G {params.genome} \
        -R tss \
        -C {output.ngs_conf} \
        -O {params.tss} \
        -FL 300
        """
