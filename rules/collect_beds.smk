ALL.extend([MD_DIR + 'Collect_beds.done'])

rule collect_beds:
    input:
        expand('{directory}{sample}{extension}',
                directory = MD_DIR,
                sample = SAMPLE_NAMES,
                extension ='.filtered.bed')
    output:
        beds = MD_DIR + 'Collect_beds.done'
    params:
        fileList = glob.glob(MD_DIR + '*.filtered.bed'),
    run:
       non_zero = [os.stat(x).st_size > 0 for x in params.fileList]
       bed_list = [i for (i, v) in zip(params.fileList, non_zero) if v]
       with open(output.beds, 'w', newline='') as csvfile:
                lwriter = csv.writer(csvfile, delimiter=',', 
                                   quoting=csv.QUOTE_NONE)
                lwriter.writerow(bed_list)
      
