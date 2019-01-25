ALL.extend([MD_DIR + 'Collect_beds_list_populated.done'])

rule populate_collected_beds:
    input:
        MD_DIR + 'Collect_beds.done'
    output:
        touch(MD_DIR + 'Collect_beds_list_populated.done')
    shell:
        """
        cat {input}
        """
