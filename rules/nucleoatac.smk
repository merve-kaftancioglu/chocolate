ALL.extend([expand('{directory}{sample}{extension}',
                  directory = NUC_DIR,
                  sample = PE_BED_NAMES,
                  extension = '.nucleoatac.done'),])
                    
rule nucleoatac:
    input:
        bed = NUC_DIR + '{sample}_peaks.broadPeak.bed',
        bam = MD_DIR + '{sample}.filtered.bam'
    output:
        touch(NUC_DIR + '{sample}.nucleoatac.done'),
    params:
        genome_size = config['peaks2bed']['sizes'],
        name = NUC_DIR + '{sample}',
        fa = REF_DIR + ORG + '.fa'
    threads: config['thread_info']['nucleoatac'],
    conda:
        'envs/nucleoatac_env.yml'
    shell:
        """
        nucleoatac run \
        --bed {input.bed} \
        --bam {input.bam} \
        --fasta {params.fa} \
        --cores {threads} \
        --out {params.name}
        """
