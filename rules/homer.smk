ALL.extend([expand('{directory}{filename}{extension}',
                    directory = HOMER_DIR + config['pepr_info']['peak_type'] + '/',
                    filename = PEPR_GROUP_NAMES,
                    extension = '_homer.done')])

rule homer:
    input:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_fixed_peaks.bed',
    output:
        touch(HOMER_DIR + config['pepr_info']['peak_type'] + '/{chip}_homer.done')
    params:
         ref = REF_DIR + ORG + '.fa',
         base = '{chip}',
         out_base = HOMER_DIR + config['pepr_info']['peak_type'] + '/{chip}'
    conda:
        'envs/homer_env.yml'
    threads:
        config['thread_info']['homer']
    shell:
         """
         findMotifsGenome.pl \
         {input} \
         {params.ref}\
         {params.out_base} \
         -p {threads}
         """
