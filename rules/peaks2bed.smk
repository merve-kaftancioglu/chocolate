ALL.extend([expand('{directory}{sample}{extension}',
                  directory = NUC_DIR,
                  sample = PE_BED_NAMES,
                  extension = '_peaks.broadPeak.bed'),])
                    
rule peaks2bed:
    input:
        PEAKS_DIR + '{sample}_peaks.broadPeak',
    output:
        NUC_DIR + '{sample}_peaks.broadPeak.bed',
    params:
        genome_size = config['peaks2bed']['sizes'],
    conda:
        'envs/bedtools_env.yml'
    shell:
        """
        bedtools slop \
        -b 200 \
        -i {input} \
        -g {params.genome_size} \
        | sort -k1,1 -k2,2n | bedtools merge > {output}
        """
