ALL.extend(expand('{directory}{filename}{extension}',
                   directory=PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/',
                   filename=PEPR_GROUP_NAMES,
                   extension= '__PePr_peaks_fixed_chr_annot.bed'))

rule pepr_chipseeker:
    input:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks_fixed_chr.bed',
    output:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks_fixed_chr_annot.bed',
    params:
        txdb = config['pepr_chipseeker']['txdb'],
        annodb = config['pepr_chipseeker']['annodb'],
        outDir = PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/',
        script = 'scripts/pepr_chipseeker.R'
    shell:
         """
         Rscript {params.script} \
         --peaks {input} \
         --txdb {params.txdb} \
         --annodb {params.annodb} \
         --outDir {params.outDir}
         """
