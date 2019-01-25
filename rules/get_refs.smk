ALL.extend([expand('{directory}{sample}{extension}',
                  directory = REF_DIR,
                  sample = ORG,
                  extension = ['.fa','.gtf'])])       

rule get_refs: 
    input:
        gtf = config['input_dir']['genome_gtf'],
        fa = config['input_dir']['genome_fasta']
    output:
        gtf = REF_DIR + ORG + '.gtf',
        fa = REF_DIR + ORG + '.fa'
    shell:
        """
        cp {input.gtf} {output.gtf}
        cp {input.fa} {output.fa}
        """
