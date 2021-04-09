#!/usr/bin/env python3
"""
July 16th
Hillary Elrick

Using the annotated peaks tsv file from homer, for every gene in the study, check if its
ensembl id is in the peaks. Then, if its enrichment was annotated as a TSS-promoter region,
record the peak value
"""

import argparse
import pandas as pd

def addAnnotation(clean_data, annotationFile, newColumnName):
    """
    Given the data and annotated peaks, match them up
    """
    avgPeakScore = []
    for index, row in clean_data.iterrows():
        # get all the annotation data associated with this gene's ensid
        gene_ensid = str(row['ENSID'])
        peaks = annotationFile.loc[(annotationFile['Nearest Ensembl'] == gene_ensid) & (annotationFile['Annotation'].str.contains('promoter-TSS'))]

        if peaks.shape[0] == 0:
            avgPeakScore.append('') # no peaks associated with this gene
        elif peaks.shape[0] == 1:
            avgPeakScore.append(str(float(peaks.iloc[0]['Peak Score'])))
        elif peaks.shape[0] > 1:
            # take the avg
            peakscores = []
            for subindex, subrow in peaks.iterrows():
                peakscores.append(float(subrow['Peak Score']))
            avg = sum(peakscores) / float(len(peakscores))
            #print("multiple, had to average: " + str(avg))
            avgPeakScore.append(avg)
    
    clean_data[newColumnName] = avgPeakScore
    
    return

def addIntronAnnotation(clean_data, annotationFile, newColumnName):
    """
    Given the data and annotated peaks, match them up
    """
    avgPeakScore = []
    for index, row in clean_data.iterrows():
        # get all the annotation data associated with this gene's ensid
        gene_ensid = str(row['ENSID'])
        peaks = annotationFile.loc[(annotationFile['Nearest Ensembl'] == gene_ensid) & (annotationFile['Annotation'].str.contains('intron'))]

        if peaks.shape[0] == 0:
            avgPeakScore.append('') # no peaks associated with this gene
        elif peaks.shape[0] == 1:
            avgPeakScore.append(str(float(peaks.iloc[0]['Peak Score'])))
        elif peaks.shape[0] > 1:
            # take the avg
            peakscores = []
            for subindex, subrow in peaks.iterrows():
                peakscores.append(float(subrow['Peak Score']))
            avg = sum(peakscores) / float(len(peakscores))
            #print("multiple, had to average: " + str(avg))
            avgPeakScore.append(avg)
    
    clean_data[newColumnName] = avgPeakScore
    
    return

PERCENTILE_RANK = False
def map_probes(probeset, entrez_ids):
  """ parse out the probe id to entrezID mapping """  
  entrez_idx = None
  mapping = {}
  with open(probeset) as probes:
    for line in probes:
      if line.startswith('ID'):
        entrez_idx = line.split('\t').index('ENTREZ_GENE_ID')
      elif entrez_idx:
        # if the index has been defined then we're past the header
        row = [x.strip() for x in line.split('\t')]
        # if we're doing percentile rank, we need all the mappings, otherwise can just track the mappings of interest
        if PERCENTILE_RANK:
          if '///' in row[entrez_idx]:
            # multile genes add an entry for every gene overlapped by the probe
            # TODO: FIX; THIS IS A MANY TO MANY MAPPING ISSUE 
            # since this only happens once in this dataset, I'm just using the first one but can also use last (or develop a solution that works for all cases...)
            mapping[row[0]] = row[entrez_idx].split(' /// ')[0]
            """ # option to use the last one 
            for entrez_id in [x for x in row[entrez_idx].split(' /// ')]:
              print('Entrez ID:'+str(entrez_id)+' in probe that maps to multiple genes')
               mapping[row[0]] = entrez_id[0]           
            """
            print('MANY TO MANY: '+str(row[0])+"->"+str(row[entrez_idx]))
          else:
              mapping[row[0]] = row[entrez_idx]
        elif row[entrez_idx] in entrez_ids:
          mapping[row[0]] = row[entrez_idx]

  return mapping

def get_expression(data_series, probes_to_genes):
  """ using the geo series, annotate a gene_list excel """
  with open(data_series, 'r') as mtx:
    stage_columns = {'all_stages': {'sample_ids': []}} # will always need an average, other stages are determined by the file
    sample_ids = None
    for line in mtx:
      if line.startswith('!Sample_title'):
        sample_stages = [x.strip().replace('"','').split(",")[0] for x in line.split("\t")[1:]] # this line is likely dataset specific.
      elif line.startswith('"ID_REF"'): # this comes after the sample titles
        sample_ids = [x.strip().replace('"','') for x in line.split("\t")[1:]]
        # now have the ids and their stages, convert to dict
        """
        if named differently, may need to modify this.
        ultimately, stage_columns should be a dictionary with the following properties:
        - the keys are the stage names. 
        - each 'stage' dict should have a key 'sample_ids' that has a list the sample_ids belonging to that stage.
        {
          'stage1': {
            'sample_ids': ['sample_id1','sample_id2', ..., 'sample_idn']
          },
          'stage2': {
            'sample_ids': ['sample_idn+1', ...]
          },
          ...
        }
        """
        for i in range(0, len(sample_stages)):
          if sample_stages[i] not in stage_columns:
            stage_columns[sample_stages[i]] = {'sample_ids': []}
          stage_columns[sample_stages[i]]['sample_ids'].append(sample_ids[i])
          stage_columns['all_stages']['sample_ids'].append(sample_ids[i]) # add every sample to this
      elif sample_ids is not None:
        row = [x.strip().replace('"','') for x in line.split('\t')]
        """
        here, the stage_columns dictionary is being updated with the expression data for each gene.
        {
          'stage1': {
            'sample_ids': ['sample_id1','sample_id2', ..., 'sample_idn'],
            'genes': { <- **NEW KEY**
              'entrezID-1': ['sample_id1ExpLevel', 'sample_id2ExpLevel', ..., 'sample_idnExpLevel'],
              'entrezID-2': ['sample_id1ExpLevel', 'sample_id2ExpLevel', ..., 'sample_idnExpLevel'],
              ... (if PERCENTILE_RANK is True, all in dataset are recorded otherwise, just the genes of interest )
            }
          },
          ...
        }
        """
        if row[0] in probes_to_genes:
          # get gene from probe
          entrez_id = probes_to_genes[row[0]]
          # add the average expression for all the samples in a stage for the gene
          for stage, stage_data in stage_columns.items():
            stage_data['genes'] = {} if 'genes' not in stage_data else stage_data['genes'] # initialize
            for sample_id in stage_data['sample_ids']:
              # get the index of the sample_id in the row
              sample_idx = sample_ids.index(sample_id) + 1
              if entrez_id not in stage_data['genes']:
                stage_data['genes'][entrez_id] = [float(row[sample_idx])]
              else:
                stage_data['genes'][entrez_id].append(float(row[sample_idx]))

  return stage_columns

def get_percentile_rank(stage_columns):
  """ rank of a gene's expression at a given stage for the entire dataset """
  percentile_ranks = {}
  for stage, stage_data in stage_columns.items():
    percentile_ranks[stage] = {}
    avg_exps = [(entrez_id, sum(exp)/len(exp)) for entrez_id, exp in stage_data['genes'].items()]
    # sort by avg expression
    avg_exps.sort(key = lambda x: x[1])
    # calculate percentile ranks
    for rank in range(len(avg_exps)):
      percentile_ranks[stage][avg_exps[rank][0]] = rank/len(avg_exps)

  return percentile_ranks

def addExpression(genelist, stage_expression, percentile_rank=None):
  """ average the expression at every stage for each available Entrez ID"""
  # add each stage to the excel
  excel_entrez_ids = genelist['Entrez_ID']
  for stage, exp_data in stage_expression.items():
    stage_col = []
    percentile_col = []
    for entrez_id in genelist['Entrez_ID']:
      if not pd.isnull(entrez_id) and str(int(entrez_id)) in exp_data['genes']:
        avg = round(sum(exp_data['genes'][str(int(entrez_id))]) / len(exp_data['genes'][str(int(entrez_id))]),3)
        stage_col.append(avg)
        if percentile_rank:
          percentile_col.append(percentile_rank[stage][str(int(entrez_id))])
      else:
        stage_col.append('')
        if percentile_rank:
          percentile_col.append('')
    
    genelist[stage] = stage_col
    if percentile_rank:
        genelist[stage+"_PercentileRank"] = percentile_col 
     
  return

def get_staining_overlap(chrom, guide_cutsites, cytogenic_bands):
  """ for a given guide, return a comma separated list of the staining regions it overlaps """

  guide_stains_overlap = []
  if str(chrom) in cytogenic_bands:
    for location in cytogenic_bands[chrom]:
      if guide_cutsites[0] >= location[0] and guide_cutsites[0] < location[1]:
        guide_stains_overlap.append(location[2])
      elif guide_cutsites[-1] >= location[0] and guide_cutsites[-1] < location[1]:
        guide_stains_overlap.append(location[2])
  else:
    print(chrom + " not found in cytogenic bands file")
    print(cytogenic_bands.keys())

  if guide_stains_overlap:
    return ",".join(guide_stains_overlap)
  else:
    return ''

def addCytogenicBands(merged_data, cytogenicBandingFile):
  """ intersect the guide coordinates with the cytogenic bands """

  # make dictionary of the cytogenic bands
  cytogenic_bands = {}
  for line in open(cytogenicBandingFile):
    chrom, start, end, _, stain = [x.strip() for x in line.split("\t")]
    chrom = chrom.replace('chr','').lower()
    if chrom not in cytogenic_bands:
      cytogenic_bands[chrom] = [[int(start),int(end),stain]]
    else:
      cytogenic_bands[chrom].append([int(start),int(end),stain])

  import ast
  staining_column = [] 
  for index, row in merged_data.iterrows():
    guide_cutsites = ast.literal_eval(row['Sorted Cut-Sites'])
    #guide_cutsites = [x.strip() for x in ast.literal_eval(row['Sorted Cut-Sites'])]
    chrom = row['Chromosome (+ strand)']
    staining_column.append(get_staining_overlap(chrom, guide_cutsites, cytogenic_bands))
  
  merged_data['Staining Overlap'] = staining_column

def addCutPositionInformation(merged_data, chromSizes):
  """ for every experiment's guide, take the average of the cutsites and divide by the chromosome length """
  # make dictionary of chrom sizes
  sizesDict = {}
  for line in open(chromSizes):
    chrom, size = [x.strip() for x in line.split("\t")]
    chrom = chrom.replace('chr','').lower()
    if chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 'x', 'y']:
      # don't bother using the non-standard chromosomes since there aren't any guides there
      sizesDict[chrom] = int(size)

  import ast
  position_column = [] 
  for index, row in merged_data.iterrows():
    guide_cutsites = ast.literal_eval(row['Sorted Cut-Sites'])
    #guide_cutsites = [x.strip() for x in ast.literal_eval(row['Sorted Cut-Sites'])]
    chrom = row['Chromosome (+ strand)']
    avg_cutsite = sum(guide_cutsites)/len(guide_cutsites)
    relativePosition =  avg_cutsite/sizesDict[chrom]
    position_column.append(relativePosition) 

  merged_data['Relative Chromosomal Position'] = position_column
  return

def addOmimAnnotation(merged_data, OmimAnnotationFile):
  """ for every experiment, if the human ortholog has an omim annotation, record it """
  omim_genes = dict.fromkeys(list(OmimAnnotationFile['ENSID']))
  has_omim = []
  for index, row in merged_data.iterrows():
    human_ensid = str(row['Human ENSID'])
    if human_ensid in omim_genes:
      has_omim.append('t')
    else:
      has_omim.append('f')

  merged_data['Has Omim Annotation'] = has_omim
  return

def main():

    # parse command line args
    parser = argparse.ArgumentParser(description="annotate data prior to secondary analysis")
    parser.add_argument('input', type=argparse.FileType('r'))
    parser.add_argument('--verbose', '-v', action='store_true', help='print verbose cleaning steps')
    parser.add_argument('output', type=argparse.FileType('w'))
    args = parser.parse_args()

    # import in the exported dataset from the clean_data.py script
    clean_data = pd.read_excel(args.input.name, engine='openpyxl')
  
    # import the mgi to ensid mapping and entrez mapping
    mapping_data = pd.read_excel(pd.ExcelFile('annotation/genelists/MGI_ENSID_Mapping.xlsx', engine='openpyxl'))
    entrez_mapping = pd.read_excel(pd.ExcelFile('annotation/genelists/MGI_Entrez_Mapping.xlsx', engine='openpyxl'))

    # do a left join on the mgi id and entrez id
    merged_data = pd.merge(clean_data, mapping_data, how='left', left_on='Gene MGI Accession ID', right_on='MGI_ID')
    merged_data.to_excel(args.output.name, index=False)
    merged_data = pd.merge(merged_data, entrez_mapping, how='left', left_on='Gene MGI Accession ID', right_on='MGI_ID')

    # Add Human Ortholog info (OMIM)
    OmimAnnotationFile = pd.read_csv('annotation/human_orthologs/omim_to_ensembl.tsv', sep='\t')
    addOmimAnnotation(merged_data, OmimAnnotationFile)

    # Add positional information about cutsites (relative to chromosome)
    chromSizes = 'annotation/mm10.chrom.sizes'
    addCutPositionInformation(merged_data, chromSizes)

    # Add Cytogenic Banding information
    cytogenicBandingFile = 'annotation/cytoBand.txt'
    addCytogenicBands(merged_data, cytogenicBandingFile)
  
    # Add ChIP Seq Methylation data
    MethylationAnnotationFile = pd.read_csv('annotation/histone_modification/mm10_H3K27me3_annotatedpeaks.tsv', sep='\t')
    addIntronAnnotation(merged_data, MethylationAnnotationFile, 'AverageH3K27me3_Intron_PeakScore')
    # Add ChIP Seq Acetylation data
    AcetylationAnnotationFile = pd.read_csv('annotation/histone_modification/mm10_H3K27ac_annotatedpeaks.tsv', sep='\t')
    addAnnotation(merged_data, AcetylationAnnotationFile, 'AverageH3K27ac_PeakScore')
    
    # Add Embryonic Expression data
    # first, map the entrez ids to individual probes
    entrez_ids = [str(int(x)) for x in merged_data['Entrez_ID'].tolist() if not pd.isnull(x)]
    probeset = 'annotation/GEO/GPL1261-56135.txt'
    probes_to_genes = map_probes(probeset, entrez_ids)
    data_series = 'annotation/GEO/GSE11224_series_matrix.txt'
    stage_expression = get_expression(data_series, probes_to_genes)
    percentile_rank = get_percentile_rank(stage_expression)
    addExpression(merged_data, stage_expression, percentile_rank)

    # Export
    merged_data.to_excel(args.output.name, index=False)

if __name__ == "__main__":
    main()


