ALL.extend(expand('{directory}{filename}{extension}',
                   directory=PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/',
                   filename=PEPR_GROUP_NAMES,
                   extension=['__PePr_peaks_fixed.bed',
                              '__PePr_peaks_fixed_chr.bed']))

rule pepr_peaks_parse:
    input:
        PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks.bed',
    output:
        fix = PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks_fixed.bed',
        fix_chr = PEPR_DIR + config['pepr_info']['peak_type'] + '/peaks/{chip}__PePr_peaks_fixed_chr.bed',
    shell:
         """
         cat {input} \
         | awk 'BEGIN{{FS=OFS="\t"}} {{gsub(/.0$/,"",$2)}};1 {{gsub(/.0$/,"",$3)}};1' \
         > {output.fix}
         cat {output.fix} \
         | sed 's/^/chr/g' > {output.fix_chr}
         """
