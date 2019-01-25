ALL.extend([HTML_DIR + 'ataqv_html.done',])
                    
rule ataqv_html:
    input:
        expand('{directory}{sample}{extension}',
                directory = ATAQV_DIR,
                sample = PE_BED_NAMES,
                extension = '.ataqv.json')
    output:
        touch(HTML_DIR + 'ataqv_html.done')
    params:
        output_dir = HTML_DIR
    shell:
        """
        source activate ataqv_env
        mkarv \
        -f \
        {params.output_dir} \
        {input}
        """
