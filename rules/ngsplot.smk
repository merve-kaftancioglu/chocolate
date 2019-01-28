ALL.extend([expand('{directory}{region}{type}{extension}',
                  directory = NGS_DIR,
                  region = ['Genebody','TSS'],
                  type = ['.avgprof','.heatmap'],
                  extension = '.pdf'),])
                    
rule ngsplot:
    input:
        expand('{directory}{sample}{extension}',
              directory = ALTERED_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam') \
        if ENSEMBL.upper() == 'Y' \
        else \
        expand('{directory}{sample}{extension}',
              directory = MD_DIR,
              sample = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam')
    output:
        ngs_conf = NGS_DIR + 'ngsplot_config',
        file = expand('{directory}{region}{type}{extension}',
                     directory = NGS_DIR,
                     region = ['Genebody','TSS'],
                     type = ['.avgprof','.heatmap'],
                     extension = '.pdf')
    params:
        generate = 'scripts/ngs_config_generate.r',
        outdir = NGS_DIR,
        gb = NGS_DIR + 'Genebody',
        tss = NGS_DIR + 'TSS',
        genome = config['ngsplot']['genome']
    shell:
        """
        source activate ngsplot_env
        {params.generate} {params.outdir} {input}
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
