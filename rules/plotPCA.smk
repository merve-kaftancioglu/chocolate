ALL.extend([DEEP_DIR + 'pcaplot_multiBamSummary.png',
            DEEP_DIR + 'pcaplot_multiBamSummary.txt'])
                    
rule plotPCA:
    input:
        DEEP_DIR + 'multiBamSummary.npz',
    output:
        png = DEEP_DIR + 'pcaplot_multiBamSummary.png',
        data = DEEP_DIR + 'pcaplot_multiBamSummary.txt'
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        plotPCA \
        -in {input} \
        -o {output.png} \
        --plotTitle "Principal Component Analysis Plot" \
        --outFileNameData {output.data}
        """
