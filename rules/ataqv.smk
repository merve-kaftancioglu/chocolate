ALL.extend([expand('{directory}{sample}{extension}',
                  directory = ATAQV_DIR,
                  sample = PE_BED_NAMES,
                  extension = ['.ataqv.json','.ataqv.out']),])
                    
rule ataqv:
    input:
        list = MD_DIR + 'Collect_beds.done',
        bed = PEAKS_DIR + '{sample}_peaks.broadPeak',
        bam = MD_DIR + '{sample}.filtered.bam'
    output:
        json = ATAQV_DIR + '{sample}.ataqv.json',
        txt = ATAQV_DIR + '{sample}.ataqv.out'
    params:
        tss = config['ataqv']['tss'],
        autosomalreference = config['ataqv']['autosomalreference'],
        extension = config['ataqv']['extension'],
        organism = config['ataqv']['organism'],
        mitochondrial = config['ataqv']['mitochondrial'],
    threads: config['thread_info']['ataqv']

    shell:
        """
        source activate ataqv_env
        ataqv \
        --peak-file {input.bed} \
        --tss-file {params.tss} \
        --tss-extension {params.extension} \
        --threads {threads} \
        --autosomal-reference-file {params.autosomalreference} \
        --mitochondrial-reference-name {params.mitochondrial} \
        --metrics-file {output.json} \
        {params.organism} \
        {input.bam} > {output.txt}
        """
