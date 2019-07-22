ALL.extend(expand([HOMER_DIR + 'homer.done']))

rule homer:
    input:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_fixed_peaks.bed',
    output:
        touch(HOMER_DIR + 'homer.done')
    params:
         ref = config['']
    conda: 'envs/homer_env.yml'
    threads: config['thread_info']['homer']
    shell:
         """
         findMotifsGenome.pl \
         {input} \
         ${REF}/Xentr41.fa \
         ${HOMER_OUTPUT}/${BASE} \
         -p {threads}
         """
