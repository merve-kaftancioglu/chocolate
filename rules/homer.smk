ALL.extend(expand([HOMER_DIR + 'homer.done']))

rule homer:
    input:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_fixed_peaks.bed',
    output:
        touch(HOMER_DIR + 'homer.done')
    params:
         ref = REF_DIR + ORG + '.fa',
         base = '{chip}',
         out_base = PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}'
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
