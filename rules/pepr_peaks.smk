ALL.extend(expand('{directory}{filename}{extension}',
                   directory=PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/',
                   filename=PEPR_GROUP_NAMES,
                   extension=['__PePr_peaks.bed',
                              '__PePr_parameters.txt']))

rule pepr_peaks:
    input:
        expand('{directory}{filename}{extension}',
              directory = MD_DIR,
              filename = SAMPLE_NAMES,
              extension = '.filtered.sorted.bam'),
    output:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks.bed',
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_parameters.txt',
    params:
        format = config['pepr_info']['format'],
        peak_type = config['pepr_info']['peak_type'],
        name = lambda wildcards: '{chip}'.format(**wildcards),
        chip = lambda wildcards: ','.join(expand('{directory}{sample}{extension}', 
                                          directory = MD_DIR,
                                          sample = PEPR_GROUPS_DICT[wildcards.chip]['chip'],
                                          extension = '.filtered.sorted.bam')),
        input = lambda wildcards:','.join(expand('{directory}{sample}{extension}', 
                                          directory = MD_DIR,
                                          sample = PEPR_GROUPS_DICT[wildcards.chip]['input'],
                                          extension = '.filtered.sorted.bam')),
        outdir = PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/',
        threshold = config['pepr_info']['threshold'],
        window = config['pepr_info']['window'],
        shift = config['pepr_info']['shiftsize']
    conda:
        'envs/pepr_env.yml'
    threads: config['thread_info']['pepr']
    shell:
         """
         PePr -c {params.chip} \
              -i {params.input} \
              -f {params.format} \
              --peaktype {params.peak_type} \
              --output-directory {params.outdir} \
              --num-processors {threads} \
              --normalization intra-group \
              --name {params.name} \
              --windowsize {params.window} \
              --shiftsize {params.shift}
         """
