ALL.extend([COUNTS_DIR + 'merged_counts.txt'])       

rule merge_counts:
    input:
        expand('{directory}{sample}{extension}',
                directory = COUNTS_DIR,
                sample = SAMPLE_NAMES,
                extension = '_starReadsPerGene.out.tab'),
    output:
        counts = COUNTS_DIR + 'merged_counts.txt'
    params:
        fileList = glob.glob(COUNTS_DIR + '*.tab'),
        strandedness = config['library_info']['strandedness']
    run:
        dfList = []
        names = [(os.path.basename(file).split('_star')[0]) for file in params.fileList]
        names.insert(0, 'GeneId')
        for file in params.fileList:
            dfList.append(pd.read_csv(file, 
                                     header=3, 
                                     delim_whitespace=True,
                                     usecols=[0, params.strandedness]))
        df = reduce(lambda x, y: pd.merge(x, y, on='N_ambiguous'), dfList)
        df.columns = names
        df.sort_index(axis=1, inplace=True)
        df.to_csv(output.counts,
                  sep='\t',
                  encoding='utf-8',
                  header=True,
                  index=False)
