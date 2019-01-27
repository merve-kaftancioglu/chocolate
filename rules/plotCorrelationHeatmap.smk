ALL.extend([DEEP_DIR + 'heatmap_SpearmanCorr_multiBamSummary.png',
            DEEP_DIR + 'heatmap_SpearmanCorr_multiBamSummary_matrix.txt'])
                    
rule plotCorrelationHeatmap:
    input:
        DEEP_DIR + 'multiBamSummary.npz',
    output:
        png = DEEP_DIR + 'heatmap_SpearmanCorr_multiBamSummary.png',
        matrix = DEEP_DIR + 'heatmap_SpearmanCorr_multiBamSummary_matrix.txt'
    conda:
        'envs/deeptools_env.yml'
    threads: config['thread_info']['deeptools']
    shell:
        """
        plotCorrelation \
        -in {input} \
        -o {output.png} \
        --outFileCorMatrix {output.matrix} \
        --corMethod spearman \
        --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap \
        --colorMap RdYlBu \
        --plotNumbers
        """
