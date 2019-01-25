ALL.extend([expand('{directory}{org}{extension}',
                  directory = REF_DIR,
                  org = ORG,
                  extension = ['.amb','.ann','.bwt','.pac','.sa'])])

rule bwa_index:
    input:
        REF_DIR + ORG + '.fa',
    output: 
        expand('{directory}{org}{extension}',
              directory = REF_DIR,
              org = ORG,
              extension = ['.amb','.ann','.bwt','.pac','.sa'])
    params:
        prefix = ORG,
        tmp = expand('{org}{extension}',
                    org = ORG,
                    extension = ['.amb','.ann','.bwt','.pac','.sa'])
    conda:
        'envs/bwa_env.yml'
    shell:
        """
        bwa index \
        -p {params.prefix} \
        {input}
        mv {params.tmp} {output}
        """
