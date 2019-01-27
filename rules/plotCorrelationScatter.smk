ALL.extend([DEEP_DIR + 'scatterplot_PearsonCorr_multiBamSummary.png',
            DEEP_DIR + 'scatterplot_PearsonCorr_multiBamSummary_matrix.txt'])
                    
rule plotCorrelationScatter:
    input:
        DEEP_DIR + 'multiBamSummary.npz',
    output:
        png = DEEP_DIR + 'scatterplot_PearsonCorr_multiBamSummary.png',
        matrix = DEEP_DIR + 'scatterplot_PearsonCorr_multiBamSummary_matrix.txt'
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        plotCorrelation \
        -in {input} \
        -o {output.png} \
        --outFileCorMatrix {output.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson Correlation of Read Counts" \
        --whatToPlot scatterplot
        """
