#!/usr/bin/env python
"""
PYTHON 3.8.5

January 16th, 2019
Hillary Elrick, adapted from Greg Clark's code

First, perform steps to clean IMITS data including:
- Joining the Model Production & Colony Information
- Filtering out Controls, Not Reporting to Public,
Experiments from Centres with < 50 completed experiments,
D10A, non-standard PAMs

Then, calculate the cut sizes based on the guides

Finally, add the viability information

"""
import os
import pandas as pd
pd.options.mode.chained_assignment = None
import random
import numpy as np
from datetime import datetime

# indicates whether functions should report information about their cleaning steps (set in main function)
VERBOSE = True

def getSuccesses(successList):
    """
    Accepts a list of all GLT experiments from iMits. Records the Production Centre, MGI Accession ID,
    Gene Marker Symbol, and Method. This is used later in two places. First, to mark experiments where a successful
    attempt exists. It's also used when filtering the repeated dataset to remove repeat subsets for which a successful
    attempt exists at a different centre.
    """
    successGenes = []
    successAttempts = []
    for index, row in successList.iterrows():
        successGenes.append(str(row['MGI Accession ID']))
        attemptRow = [str(row['MGI Accession ID'])]
        attemptRow.append(str(row['Marker Symbol']))
        attemptRow.append(str(row['Production Centre']))
        attemptRow.append(str(row['Method']))
        successAttempts.append(attemptRow)

    successGenes = list(set(successGenes))
    return successGenes, successAttempts


def filterData(joinedData, viabilityReport):
    """
    Accepts the dataframe that has the merged Mouse Production and F1 Colony information.
    Removes experiments that do not fit the criteria or are missing essential information.
    Also removes centres below the cutoff from the viability report dataframe.
    """
    # make a deep copy of the original so can compare later
    fd = joinedData.copy(deep=True)

    removedCount = 0
    for index, row in fd.iterrows():
        if str(row['Report Micro Injection Progress To Public']) == 'f':
            fd.drop(index, inplace=True)
            removedCount += 1
        elif str(row['Experimental (exclude from production reports)']) == 't':
            fd.drop(index, inplace=True)
            removedCount += 1
        elif str(row['mRNA Nuclease']) == 'D10A' or str(row['Protein Nuclease']) == 'D10A':
            fd.drop(index, inplace=True)
            removedCount += 1
        elif str(row['Report F1 Colony To Public']) == 'f':
            fd.drop(index, inplace=True)
            removedCount += 1
        elif str(row['Status Name']) not in ["Genotype confirmed", "Micro-injection aborted", "Founder obtained"]:
            # this will remove unset or Micro-Injection in progress experiments
            fd.drop(index, inplace=True)
            removedCount += 1
        elif not pd.isnull(row['Embryo Transfer Day']) and str(row['Embryo Transfer Day']) == 'Next Day':
            # there are not very many of these so removing to keep the dataset consistent
            fd.drop(index, inplace=True)
            removedCount += 1
        elif str(row['Genotype Confirmed']) == 't' and pd.isnull(row['Mutant Fasta Sequence']):
            # if genotype is confirmed, allele QC should be present
            fd.drop(index, inplace=True)
            removedCount += 1
        elif not pd.isnull(row['Allele Type']):
            '''
            N.B. this must be the last criteria checked, otherwise elif statements after it
            won't be applied if they also have their Allele Type set
            '''
            # if the allele type is not null we want to ensure it's a deletion
            if str(row['Allele Type']) != 'Deletion':
                fd.drop(index, inplace=True)
                removedCount += 1
            elif not pd.isnull(row['Allele Subtype']):
                # if it is a deletion but its subtype is Indel, exclude
                if str(row['Allele Subtype']) == 'Indel':
                    fd.drop(index, inplace=True)
                    removedCount += 1

    validRemoved = removedCount
    excludingString = f'Removed {validRemoved} experiments that don\'t match criteria for analysis:'
    excludingString += "\nReport Micro Injection Progress = 'f'"
    excludingString += "\nExperimental = 't'"
    excludingString += "\nmRNA or Protein Nucelase = 'D10'"
    excludingString += "\nReport F1 Colony to Public = 'f'"
    excludingString += "\nStatus Name = 'Micro-injection in progress'"
    excludingString += "\nEmbryo Transfer Day is 'Next Day'"
    excludingString += "\nGenotype is confirmed but Allele QC data is missing"
    excludingString += "\nAllele Type is set and is not Deletion"
    excludingString += "\nAllele Subtype is Indel"
    excludingString += "\n"
    print(excludingString)

    for index, row in fd.iterrows():
        if str(row['Delivery Method']) not in ["Cytoplasmic Injection", "Pronuclear Injection", "Electroporation"]:
            fd.drop(index, inplace=True)
            removedCount += 1
        elif not pd.isnull(row['#Embryos Injected']) and int(row['#Embryos Injected']) > 1000:
            fd.drop(index, inplace=True)
            removedCount += 1
        elif pd.isnull(row['#Founder Pups Born']):
            if str(row['Status Name']) != "Micro-injection aborted":
                # number of pups should be set, remove record
                fd.drop(index, inplace=True)
                removedCount += 1

    invalidRemoved = removedCount - validRemoved
    invalidString = f'Removed {invalidRemoved} experiments with data issues:'
    invalidString += "\nDelivery Method"
    invalidString += "\n#Founder Pups Born is not set (for non MI aborted experiments)"
    invalidString += "\n>1000 Embryos Injected"
    invalidString += "\n"
    print(invalidString)

    # remove experiments with HR, HDR, donor insertions, or subset of donor insertions detected
    # additionally, if there are experiments with one of HR/HDR/donor explicitly set to 0, > 0
    # deletion events detected, and no founders selected for breeding, we can surmise that this
    # experiment was never intended to be a deletion and should be removed
    wrongEvent = ['#G0 HR event detected',
                  '#G0 HDR event detected',
                  '#G0 all donor insertions detected',
                  '#G0 subset of donors inserted detected']
    g0RemovedCount = 0
    for index, row in fd.iterrows():
        removed = False
        for event in wrongEvent:
            if not pd.isnull(row[event]) and int(row[event]) > 0:
                fd.drop(index, inplace=True)
                removed = True
                break  # don't want to remove twice
        if not removed and columnsSet(wrongEvent, row):
            # if it wasn't previously removed and one of the wrong columns is set,
            # check if there were deletion events detected but no founders selected for breeding (either unset or 0)
            if not pd.isnull(row['#G0 deletion event detected']) and (int(row['#G0 deletion event detected']) > 0):
                if pd.isnull(row['#Founders Selected For Breeding']) or (int(row['#Founders Selected For Breeding']) == 0):
                    fd.drop(index, inplace=True)
                    removed = True
        if removed:
            g0RemovedCount += 1

    print(f'Removed {g0RemovedCount} experiments with HR, HDR, or donor insertions detected or 0 HR, HDR')
    print(f'Or donor insertions detected, > 0 deletion events detected, but no founders selected for breeding')

    # determine centres above the cutoff
    completedCounts = {}
    for index, row in fd.iterrows():
        # counting the number of successful deletion attempts per centre
        if row['Status Name'] == 'Genotype confirmed' and row['Allele Type'] == 'Deletion':
            centre = str(row['Production Centre'])
            if centre in completedCounts:
                completedCounts[centre] += 1
            else:
                completedCounts[centre] = 1

    if VERBOSE:
        print(completedCounts)

    # remove centres below the threshold from the dict
    completedCounts = {centre: count for centre, count in completedCounts.items() if count > 50}

    if VERBOSE:
        print(completedCounts)

    # remove rows from centres that don't have > 50
    centreRemovedCount = 0
    for index, row in fd.iterrows():
        if str(row['Production Centre']) not in completedCounts:
            fd.drop(index, inplace=True)
            centreRemovedCount += 1

    print(f'Removed {centreRemovedCount} experiments from centres with < 50 Genotype Confirmed Deletion experiments')

    '''
    We also want to remove these rows from the viability data. However, the centre
    names don't exactly match up so the below_threshold list was manually generated 
    after looking at the counts
    '''
    below_threshold = ['HMGU', 'RBRC', 'MARC']
    viabRemovedCount = 0
    for index, row in viabilityReport.iterrows():
        if str(row['phenotyping_center']) in below_threshold:
            viabilityReport.drop(index, inplace=True)
            viabRemovedCount += 1

    print(f'Removed {viabRemovedCount} viability calls from centres with < 50 Genotype Confirmed Deletion experiments')

    return fd


def removeNonProteinCoding(filteredData, proteinCodingGenes):
    """
    Only want to keep the experiments on genes that are protein-coding
    """
    # make a deep copy of the original so can compare later
    fd = filteredData.copy(deep=True)

    # create list of all protein coding genes' mgi ids
    proteinCoding_MGI_IDs = proteinCodingGenes['MGI_ID'].tolist()

    # see which ones aren't protein coding
    remove = fd[~fd['Gene MGI Accession ID'].isin(proteinCoding_MGI_IDs)]['Gene Marker Symbol'].tolist()
    if VERBOSE:
        print("Removed " + str(len(remove)) + " non protein-coding genes listed below:")
        # print out those removed for being non-coding
        for gene in remove:
            print(gene)
    else:
        print("Removed " + str(len(remove)) + " non protein-coding genes")

    # remove those not in the protein coding list
    fd = fd[fd['Gene MGI Accession ID'].isin(proteinCoding_MGI_IDs)]

    return fd


def columnsSet(subset, row):
    """
    Given a row, checks if even one of the columns in subset are not null.
    Used primarily when seeing if any HR/HDR/Donor Insertion values are set as
    this gets ugly when done on one line
    """

    for column in subset:
        if not pd.isnull(row[column]):
            return True

    # none of columns in the subset were set
    return False


def addCas9Concentration(fd):
    """
    Removes experiments that have issues with their mRNA or Protein Concentration in
    terms of missing data, conflicting data, or outliers. Then adds a 'Concentration
    Category' field to the remaining records. Currently:
    mRNA < 50 : Low
    mRNA >= 50 and < 100 : Medium
    mRNA >= 100 : High
    Protein Concentration < 500 : Low
    Protein Concentration >= 500 and < 1000 : Medium
    Protein Concentration >= 1000 : High
    """
    conc_mRNA = []
    conc_protein = []

    removedCount = 0
    for index, row in fd.iterrows():
        if not pd.isnull(row['mRNA Concentration']) and not pd.isnull(row['Protein Concentration']):
            removedCount += 1
            if VERBOSE: print("Both mRNA and Protein set: " + str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
        elif pd.isnull(row['mRNA Concentration']) and pd.isnull(row['Protein Concentration']):
            removedCount += 1
            if VERBOSE: print("Neither mRNA or Protein concentration set: " + str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
        elif not pd.isnull(row['mRNA Concentration']):
            if int(row['mRNA Concentration']) <= 5:
                fd.drop(index, inplace=True)
                removedCount += 1
                if VERBOSE: print("Outlier: " + str(row['Mi Attempt URL']) + " [mRNA] = " + str(row['mRNA Concentration']))
            elif int(row['mRNA Concentration']) >= 1344:
                fd.drop(index, inplace=True)
                removedCount += 1
                if VERBOSE: print("Outlier: " + str(row['Mi Attempt URL']) + " [mRNA] = " + str(row['mRNA Concentration']))
            else:
                conc_mRNA.append(int(row['mRNA Concentration']))
        elif not pd.isnull(row['Protein Concentration']):
            if int(row['Protein Concentration']) < 10:
                removedCount += 1
                fd.drop(index, inplace=True)
                if VERBOSE: print("Outlier: " + str(row['Mi Attempt URL']) + " [Protein] " + str(row['Protein Concentration']))
            else:
                conc_protein.append(int(row['Protein Concentration']))

    print(f'Removed {removedCount} experiments for concentration outliers or conflicting values (both mRNA and Protein set)')

    # set variable indicating whether mRNA or protein was used
    conditions = [(fd['mRNA Concentration'].notnull()) & (fd['Protein Concentration'].isnull()),
                  (fd['mRNA Concentration'].isnull()) & (fd['Protein Concentration'].notnull())]
    cas9Type = ["mRNA", "Protein"]
    fd['Cas9 Type'] = np.select(conditions, cas9Type)

    # set the levels
    thresholds = [(fd['mRNA Concentration'] < 50) & (fd['Protein Concentration'].isnull()),
             (fd['mRNA Concentration'] >= 50) & (fd['mRNA Concentration'] < 100) & (fd['Protein Concentration'].isnull()),
             (fd['mRNA Concentration'] >= 100) & (fd['Protein Concentration'].isnull()),
             (fd['Protein Concentration'] < 500) & (fd['mRNA Concentration'].isnull()),
             (fd['Protein Concentration'] >= 500) & (fd['Protein Concentration'] < 1000) & (fd['mRNA Concentration'].isnull()),
             (fd['Protein Concentration'] >= 1000) & (fd['mRNA Concentration'].isnull())
             ]
    categories = ["Low", 'Medium', 'High', 'Low', 'Medium', 'High']

    # update 03/27/2019: Not reporting the concentration level categories
    #fd['Cas9 Concentration Level'] = np.select(thresholds, categories)

    if VERBOSE:
        # sake of curiosity
        """
        print("# mRNA: " + str(len(conc_mRNA)))
        print("Max mRNA: " + str(max(conc_mRNA)))
        print("Min mRNA: " + str(min(conc_mRNA)))
        print("Mean mRNA: " + str(np.mean(conc_mRNA)))
        print("mRNA Q1: " + str(np.percentile(conc_mRNA, 25)))
        print("mRNA Q2 (Median): " + str(np.percentile(conc_mRNA, 50)))
        print("mRNA Q3: " + str(np.percentile(conc_mRNA, 75)))
        print("# protein: " + str(len(conc_protein)))
        print("Max protein: " + str(max(conc_protein)))
        print("Min protein: " + str(min(conc_protein)))
        print("Mean protein: " + str(np.mean(conc_protein)))
        print("Protein Q1: " + str(np.percentile(conc_protein, 25)))
        print("Protein Q2 (Median): " + str(np.percentile(conc_protein, 50)))
        print("Protein Q3: " + str(np.percentile(conc_protein, 75)))
        """


def addCutSizes(fd):
    """
    Casts the chr fields to lowercase (so x is consistent), add cuts for the guides,
    calculate cutsizes, and remove data that has unsupported PAMs
    """

    chrColumns = ['Chromosome (+ strand)',
                  'Chromosome (+ strand).1',
                  'Chromosome (+ strand).2',
                  'Chromosome (+ strand).3']
    for col in chrColumns:
        fd[col] = fd[col].str.lower()

    guideColumns = ['gRNA Sequence (+ strand)',
                    'gRNA Sequence (+ strand).1',
                    'gRNA Sequence (+ strand).2',
                    'gRNA Sequence (+ strand).3']
    allCuts = []
    allGuideNums = []
    allCutsizes = []
    removedCount = 0
    for index, row in fd.iterrows():
        guideCuts = calculateCutSize(row, guideColumns)
        if guideCuts:  # we don't add if the PAM isn't supported
            maxCut = abs(guideCuts[-1]-guideCuts[0])
            # anything < 33bp or >10kbp is likely error, experimental, or non-deletion; remove it
            if maxCut > 10000 or maxCut < 33:
                fd.drop(index, inplace=True)
                removedCount += 1
                if VERBOSE:
                    print("Removing " + row['Mi Attempt URL'] + ", Calculated Cut: " + str(maxCut))
            else:
                allCuts.append(guideCuts)
                allGuideNums.append(len(guideCuts))
                allCutsizes.append(maxCut)  # subtract first from last
        else:  # unsupported PAM, remove record
            fd.drop(index, inplace=True)
            removedCount += 1
            if VERBOSE:
                print("Removing " + row['Mi Attempt URL'] + ", Unsupported PAM")

    fd['Sorted Cut-Sites'] = allCuts
    fd['Num Guides'] = allGuideNums
    fd['Max Cut-size'] = allCutsizes

    print(f'Removed {removedCount} experiments with an unsupported PAM, or Max Cut Size < 33bp or > 10kbp')


def markRepeated(fd):
    geneSymbols = fd["Gene Marker Symbol"]
    # find all the experiments using a repeated gene
    experimentsRepeatedGenes = fd[geneSymbols.isin(geneSymbols[geneSymbols.duplicated()])]
    # get just the genes
    repeatedGenes = experimentsRepeatedGenes["Gene Marker Symbol"].tolist()

    flagColumn = []
    for index, row in fd.iterrows():
        if not pd.isnull(row['Gene Marker Symbol']) and str(row['Gene Marker Symbol']) in repeatedGenes:
            flagColumn.append('t')
        else:
            flagColumn.append('f')

    fd['RepeatedGene'] = flagColumn

    return


def addEssentiality(essentialityReport, fd):
    '''
    Add the cell essentiality from Human genes mapped to mouse orthologs. Match by MGI ID of the Mouse Gene
    '''

    essential_column = []
    for index, row in fd.iterrows():
        mgi_id = str(row['Gene MGI Accession ID'])
        # select the row(s) (there should only be one) in the essentialityReport that match this attempt's MGI ID
        subset = essentialityReport.loc[essentialityReport['MGI_ID'] == mgi_id]
        if subset.shape[0] == 1:
            essential_column.append(str(subset.iloc[0]['Cellular_Essential']))
        elif subset.shape[0] > 1:
            print("ERROR: multiple essentiality results for " + str(row['Gene Marker Symbol']))
            essential_column.append('')
        else:
            essential_column.append('')

    fd['Cellular Essential'] = essential_column


def addCalculatedColumns(filteredData, viabilityReport, essentialityReport):
    """
    Accepts the data that's been filtered add adds columns that are derived from the existing data.
    Returns the modified dataframe
    """
    # make a deep copy of original to compare later
    fd = filteredData.copy(deep=True)

    '''ADDING COLUMNS WITH MORE EXTENSIVE CALCULATIONS/LOGIC FIRST'''
    # (records also filtered during G0, Cas9 concentration, and CutSize calculations)
    addCas9Concentration(fd)
    addCutSizes(fd)
    addViability(viabilityReport, fd)
    addEssentiality(essentialityReport, fd)

    # set nan (blank) to 0 for columns that are required for calculations
    required = ['#Embryos Injected',
                '#Embryos Transfered',
                '#Founder Pups Born',
                '#G0 deletion event detected',
                '#Founders Selected For Breeding']
    fd[required] = fd[required].fillna(value=0)


    '''DOING LOGICAL CHECKS ON DATA'''
    qc_issues = {}
    deletionSubtypes = ['Exon Deletion', 'Inter-exdel deletion', 'Intra-exdel deletion', 'Whole-gene deletion']
    logicalIssues = 0
    for index, row in fd.iterrows():
        if row['#Embryos Transfered'] > row['#Embryos Injected']:
            if VERBOSE:
                print("Logical Issue, more embryos transferred than injected")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#Founder Pups Born'] > row['#Embryos Transfered']:
            if VERBOSE:
                print("Logical Issue, more pups born than embryos transferred")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#G0 deletion event detected'] > row['#Founder Pups Born']:
            if VERBOSE:
                print("Logical Issue, more pups with mutation than pups born")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#Founders Assayed'] > row['#Founder Pups Born']:
            if VERBOSE:
                print("Logical Issue, more Founders Assayed than Founder Pups Born")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#Founders Selected For Breeding'] > row['#G0 deletion event detected']:
            if VERBOSE:
                print("Logical Issue, more Founders Selected for Breeding than #G0 with deletion event detected")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#Founders Selected For Breeding'] == 0 and str(row['Genotype Confirmed']) == 't':
            if VERBOSE:
                print("Logical Issue, no Founders Selected for Breeding but Genotype is Confirmed ")
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1
        elif row['#G0 deletion event detected'] == 0 and row['Allele Subtype'] in deletionSubtypes:
            if VERBOSE:
                print("Logical Issue, no #G0 with deletion, but Allele Subtype is " + str(row['Allele Subtype']))
                print(str(row['Mi Attempt URL']))
            fd.drop(index, inplace=True)
            logicalIssues += 1

        """
        elif str(row['Genotype Confirmed']) == 't' and pd.isnull(row['Mutant Fasta Sequence']):
            if str(row['Production Centre']) not in qc_issues:
                qc_issues[str(row['Production Centre'])] = []
            qc_issues[str(row['Production Centre'])].append(str(row['Mi Attempt URL']))
            if VERBOSE:
                print("Logical Issue, Genotype is Confirmed but no QC data provided")
            fd.drop(index, inplace=True)
            logicalIssues += 1
        """

    print(f'\nRemoved {logicalIssues} records with logical issues:')
    print('More Embryos Transferred than Injected')
    print('More Founder Pups Born than Embryos Transferred')
    print('More #G0 pups than #Pups Born')
    print('More Founders Assayed than Founder Pups Born')
    print('Allele subtype is a deletion subtype but no #G0 with deletion reported')
    print('No Founders Selected for Breeding but Genotype is Confirmed')
    print('More Founders Selected for Breeding than #G0 with deletion event detected')

    """
    print("QC ISSUES:")
    num_qc_issues = 0
    for centre, missing_qc_attempts in qc_issues.items():
        print(centre + ":")
        print(str(len(missing_qc_attempts))+" ATTEMPT(S) MISSING ALLELE QC DATA")
        for attempt in missing_qc_attempts:
            print(attempt)
        num_qc_issues += len(missing_qc_attempts)

    print(str(num_qc_issues)+" IDENTIFIED")
    """

    '''CALCULATING RATES AND RATIOS'''
    # interested in the percentage of embryos injected that survived
    fd['Ratio Embryos Survived to Transfer'] = (fd['#Embryos Transfered'] / fd['#Embryos Injected'])
    # reporting the proportion of #G0 with the desired mutation that were selected for breeding
    fd['Ratio of #G0 with Mutation Selected for Breeding'] = (fd['#Founders Selected For Breeding'] / fd['#G0 deletion event detected'])
    # interested in the # Pups Born per #Embryos Transferred
    fd['Birth Rate'] = (fd['#Founder Pups Born'] / fd['#Embryos Transfered'])
    # calculate the Founder Rates
    fd['Founder Rate (per Embryos Transferred)'] = fd['#G0 deletion event detected'] / fd['#Embryos Transfered']
    # set all the GLT results based on Status Name. 'Founder obtained' experiments are explicitly not set
    conditions = [
        (fd['Status Name'] == 'Genotype confirmed'),
        (fd['Status Name'] == 'Micro-injection aborted'),
        (fd['Status Name'] == 'Founder obtained')
    ]
    results = ['t', 'f', '']
    fd['GLT'] = np.select(conditions, results)

    # fill in the rates that were 0/0 (so were not set) with 0s
    fillBlanks = ['Ratio Embryos Survived to Transfer',
                  'Ratio of #G0 with Mutation Selected for Breeding',
                  'Birth Rate',
                  'Founder Rate (per Embryos Transferred)'
                  ]
    fd[fillBlanks] = fd[fillBlanks].fillna(value=0)

    # mark repeated genes after all the filtering is done
    markRepeated(fd)

    return fd


def calculateCutSize(row, guideColumns):
    guideCuts = []

    # get all the guides for the row
    for guideType, column in enumerate(guideColumns):
        guide = str(row[column])
        if guide != 'nan':
            if guide.endswith('GG') and guide.startswith('CC'):
                # this is ambiguous, randomly choose whether to be pos or neg
                strand = random.choice(['+', '-'])
            elif guide.endswith('GG'):
                strand = '+'
            elif guide.startswith('CC'):
                strand = '-'
            else:
                # unsupported PAM, return False to indicate that the row should be removed
                return False
            # know direction, use guide's numbered start column to calculate cut
            startCol = 'Start Co-od'
            if guideType == 0:
                if strand == '+':
                    guideCut = int(row[startCol]) + 17
                else:
                    guideCut = int(row[startCol]) + 6
            else:
                # determine the column name
                startCol = startCol + '.' + str(guideType)
                if strand == '+':
                    guideCut = int(row[startCol]) + 17
                else:
                    guideCut = int(row[startCol]) + 6

            guideCuts.append(guideCut)

    guideCuts = sorted(guideCuts)
    return guideCuts


def addViability(viabilityReport, filteredData):
    # get the columns we're interested in
    viabilityReport = viabilityReport[['marker_symbol', 'categories']]

    '''
    First we want to determine the counts of the viability calls on a 
    gene by gene basis. i.e. if one gene had multiple colonies with
    different viablity calls they should both be recorded
    '''
    viabilityDict = {}
    for index, row in viabilityReport.iterrows():
        gene = str(row['marker_symbol'])
        viaCat = str(row['categories'])
        # either add or increment
        if gene not in viabilityDict:
            viabilityDict[gene] = {}
            viabilityDict[gene][viaCat] = 1
        else:
            if viaCat in viabilityDict[gene]:
                viabilityDict[gene][viaCat] += 1
            else:
                viabilityDict[gene][viaCat] = 1

    '''
    Now we want to add the counts for the viability calls and the consensus
    (if there are conflicts) to the filtered data. Tracking the values for the
    counts and consensus in separate lists which are added to the DataFrame
    after they're calculated for every row 
    '''
    import operator
    import collections
    viabCounts = []
    viabConsen = []
    noConsensusGenes = []  # track the genes with no consensus so they're not reported twice
    for index, row in filteredData.iterrows():
        gene = row['Gene Marker Symbol']
        if str(gene) in viabilityDict:
            # ensure that the gene has a consensus viability call
            value_to_key = collections.defaultdict(list)
            for k, v in viabilityDict[gene].items():
                # add all categories with non-zero counts to list
                if v > 0:
                    value_to_key[v].append(k)
                    # report categories that have the same # of calls
                    if len(value_to_key[v]) > 1:
                        if gene not in noConsensusGenes:
                            print(gene + " has no viability consensus: " + str(v) + " calls each for " + str(value_to_key[v]))
                            print("Will use '" + str(value_to_key[v][0]) + "' as the consensus")
                            noConsensusGenes.append(gene)
            # add the viability categories for the gene
            viabCounts.append(str(viabilityDict[gene]))
            # get the category with most number of counts
            maxC = max(viabilityDict[gene].items(), key=operator.itemgetter(1))[0]
            # report just Viable, Subviable, and Lethal (remove Homozygous and Hemizygous qualifiers)
            maxC = maxC.split(" - ")[1]
            # modification 03/27/2019: split based on Lethal v Non-Lethal (i.e. merge Viable and Subviable categories)
            if maxC in ['Viable', 'Subviable']:
                maxC = 'Non-Lethal'
            else:
                maxC = 'Lethal'

            viabConsen.append(maxC)
        else:
            # need a blank cell for genes with no viability calls
            viabCounts.append('')
            viabConsen.append('')

    # append the columns to the filtered data
    filteredData['Viability Counts'] = viabCounts
    filteredData['Viability Consensus'] = viabConsen


def joinData(imitsData):
    """
    Using the Mi Attempt URL, join the Mouse Production information with the F1 Colony information.
    Note that the merged count will be higher than Mouse Production alone since any Attempts that produce
    multiple colonies will have their rows duplicated in the merge. After filtering based on Allele type
    has been done, these duplicates will be removed. It's performed after the Allele type filtering since 
    if a single attempt produces colonies with different allele types, we only want the GLT status from
    the Deletion experiments (Not Indel, HDR, or HR)
    """
    # separate the Mouse Production and Colony sheets
    #attempts = pd.read_excel(imitsData, 'Mouse Production', header=5)  # column headings in fifth row
    # column headings were fixed March 2020, no need to specify header
    attempts = pd.read_excel(imitsData, 'Mouse Production')  # column headings in fifth row
    #colonies = pd.read_excel(imitsData, 'F1 Colonies', header=3)  # column headings in 3rd row
    colonies = pd.read_excel(imitsData, 'F1 Colonies')  # column headings in 3rd row

    # join the two sheets on their MI Attempt URL etc. also include other fields that are
    # reported in both sheets so they aren't duplicated (have ensured that these fields all match)
    joinedData = pd.merge(attempts, colonies, on=['Mi Attempt URL',
                                                  'Mi Attempt External Ref',
                                                  'Gene Marker Symbol',
                                                  'Gene MGI Accession ID',
                                                  'Consortium',
                                                  'Production Centre'], how='left')

    # remove fields from the merged dataset we have no interest in analyzing and don't need for filtering
    remove = ['IMPC Mutant Strain Donor', 'MF External Ref', 'Comments', 'Mi Attempt External Ref', 'F1 Colony Name',
              'Genotyping Comment', 'Uploaded Trace file', 'MGI Allele Symbol Superscript', 'MGI Allele Accession ID',
              'Centres Allele Description', 'Released from Genotyping (WTSI)']
    joinedData = joinedData.drop(columns=remove)

    if VERBOSE:
        # print out the shape/attributes of the data
        numAttempts = len(attempts.index)
        print("Number of Mouse Production Columns: " + str(len(attempts.columns)))
        print("Number of Attempts: " + str(numAttempts))
        numColonies = len(colonies.index)
        print("Number of F1 Colony Columns: " + str(len(colonies.columns)))
        print("Number of Colonies: " + str(numColonies))
        numMerged = len(joinedData.index)
        print("Number of columns retained in merged dataset: " + str(len(joinedData.columns)))
        print("Number of records when merged: " + str(numMerged))

    return joinedData


def removeFalseRepeats(repeatedAttempts, repeats, successAttempts):
    '''
    If a non-GLT (either MI Aborted or Founders Obtained experiment) that DID have #G0 with deletion > 0 has an
    Mi Date within 6 weeks of a GLT experiment, we can't conclude that it didn't go GLT due to failure since it may have
    been aborted simply because the other experiment already produced a colony.

    We also want to filter out those genes that have gone GLT at a different centre. We use the previously imported
    list of all GLT genes to do this.
    '''

    from dateutil.relativedelta import relativedelta

    toRemove = []

    falseRemoved = repeatedAttempts.copy(deep=True)
    for repeat in repeats:
        repeatSubset = falseRemoved.loc[falseRemoved['GeneSymbol_ProductionCentre'] == repeat]

        repeatSubset.sort_values('datetime')

        # if the first experiment in the group is successful, the group shouldn't be considered in the repeat analysis
        if repeatSubset.iloc[0]['GLT'] == 't':
            if VERBOSE:
                print("Removing group, First Attempt was Successful: " + str(repeat))
            toRemove.append(repeat)
            continue  # group already marked for removal, don't need to check specific experiments

        # check all the non-GLT experiments within the subset with Founders against the GLT experiments
        falseSubset = repeatSubset.loc[(repeatSubset['#G0 deletion event detected'] > 0) & (repeatSubset['GLT'] != 't')]
        gltSubset = repeatSubset.loc[repeatSubset['GLT'] == 't']
        for index, row in falseSubset.iterrows():
            # check if the 'false failures' were completed within w weeks of a GLT experiment from the 'repeat subset'
            exp_date = pd.to_datetime(row['Mi Date'])
            w = 6
            before = exp_date - relativedelta(weeks=+w)
            after = exp_date + relativedelta(weeks=+w)

            # get the subset of experiments performed within 6 weeks of a GLT experiment
            withinDate = gltSubset.loc[(pd.to_datetime(gltSubset['Mi Date']) <= after) & (pd.to_datetime(gltSubset['Mi Date']) >= before)]

            # check the length of subset of successes within date, if > 0 then remove the record
            if len(withinDate.index) > 0:
                if VERBOSE:
                    print("Removing non-GLT experiment " + str(row['Mi Attempt URL'] + " from " + str(row['GeneSymbol_ProductionCentre']) + " group as it was completed within 6 weeks of a GLT experiment"))
                falseRemoved.drop(index, inplace=True)

        # now, if there are no successes in the entire repeat group, ensure that the gene has not been successfully
        # knocked out in any experiment in the iMits dataset
        if len(gltSubset.index) == 0:
            # get the mgi id and centre from the first element of the repeat subset
            mgi_id = repeatSubset.iloc[0]['Gene MGI Accession ID']
            for attempt in successAttempts:
                if str(mgi_id) == str(attempt[0]):
                    # this means that all of the repeats should be removed as a successful attempt exists
                    toRemove.append(repeat)

    # now apply the dropping of the filtered experiments
    falseRemoved = falseRemoved[~falseRemoved.GeneSymbol_ProductionCentre.isin(toRemove)]

    return falseRemoved


def addCorrections(cleanData, correctionData):
    '''
    After sending some experiments we had questions about, BCM sent back some corrections to the counts for
    #Founders Selected for Breeding. Adding this information before exporting
    '''
    # add a blank column to store why the GLT failed
    cleanData['Reason GLT Failed'] = ''
    # for every row in the correction data, find the attempt by url in the cleanData and then update the columns
    for index, row in correctionData.iterrows():
        attemptURL = str(row['Mi Attempt URL'])
        attemptRow = cleanData.loc[cleanData['Mi Attempt URL'] == attemptURL]
        # ensure only one row is matched before updating it
        if attemptRow.shape[0] == 1:
            # only update the founders selected count if it's different
            if attemptRow.iloc[0]['#Founders Selected For Breeding'] != row['#Founders Selected For Breeding']:
                cleanData.loc[cleanData['Mi Attempt URL'] == attemptURL, '#Founders Selected For Breeding'] = row['#Founders Selected For Breeding']
        else:
            print("Error: correction row not found in clean data set, will not update")

    return


def addReasonGLTFailed(cleanData, gltData):
    '''
    Asked several centres about the reason GLT failed, integrating their responses here
    '''
    # add a blank column to store why the GLT failed
    # for every row in the data from centres, if it has a value in the Reason GLT Failed column, find the attempt by url
    for index, row in gltData.iterrows():
        if not pd.isnull(row['Reason GLT Failed']):
            attemptURL = str(row['Mi Attempt URL'])
            attemptRow = cleanData.loc[cleanData['Mi Attempt URL'] == attemptURL]
            # ensure only one row is matched before updating it
            if attemptRow.shape[0] == 1:
                cleanData.loc[cleanData['Mi Attempt URL'] == attemptURL, 'Reason GLT Failed'] = str(row['Reason GLT Failed'])
            else:
                print("Error: GLT Failure note row not found in clean data set, will not update")

    return


def repeatAnalysis(cleanData, successAttempts):
    '''
    Compile info about repeats and the variables changed between them. Outputs as a dataframe
    An experiment is considered a repeat if the gene symbol and production centre match
    '''

    # make a deep copy of the clean data to modify & output
    repeatedAttempts = cleanData.copy(deep=True)
    repeatedAttempts.sort_values('datetime')

    # remove the founder obtained experiments as these shouldn't be considered in the repeat analysis
    repeatedAttempts = repeatedAttempts[repeatedAttempts['Status Name'] != 'Founder obtained']

    # record the attempts that are included in the repeat analysis
    experiments_repeated = []

    # add a column that's the concatenation of the Gene Symbol and Production Centre
    repeatedAttempts['GeneSymbol_ProductionCentre'] = repeatedAttempts['Gene Marker Symbol'] + "-" + repeatedAttempts['Production Centre']

    prodCentre_gene = repeatedAttempts["GeneSymbol_ProductionCentre"]  # column to identify repeats at a centre
    repeatedAttempts = repeatedAttempts[prodCentre_gene.isin(prodCentre_gene[prodCentre_gene.duplicated()])]

    # variable to track which centre/gene combinations were in the repeated attempt list
    repeats = repeatedAttempts["GeneSymbol_ProductionCentre"].tolist()
    repeats = list(set(repeats))

    # remove the 'false' failures and sets where the first attempt is successful from the repeated attempts
    falseRemoved = removeFalseRepeats(repeatedAttempts, repeats, successAttempts)

    # refresh these variables after the 'false repeats' are removed
    prodCentre_gene = falseRemoved["GeneSymbol_ProductionCentre"]  # column to identify repeats at a centre
    falseRemoved = falseRemoved[prodCentre_gene.isin(prodCentre_gene[prodCentre_gene.duplicated()])]

    repeats = falseRemoved["GeneSymbol_ProductionCentre"].tolist()
    repeatDict = {repeat: 0 for repeat in repeats}  # to track which ones to do comparison for

    # sort by gene symbol/production centre, Mi date
    falseRemoved = falseRemoved.sort_values(by=['GeneSymbol_ProductionCentre', 'Mi Date'])
    falseRemoved = falseRemoved.reset_index(drop=False)

    # want to determine the differences in controlled variables
    # using sorted cuts to quickly see if the guides are different (Max Cut could potentially not report if
    # different guides had the same cut size by chance)
    controlledVariables = ['Delivery Method',
                           'mRNA Concentration',
                           'Protein Concentration',
                           'Cas9 Type',
                           'Sorted Cut-Sites',
                           'Num Guides']

    # 2d 'differences' matrix tracks which variables were changed between one experiment and its 'predecessor'
    differences = []
    for index, row in falseRemoved.iterrows():
        experiments_repeated.append(row['Mi Attempt URL'])
        row_differences = []
        centreGene = str(row['GeneSymbol_ProductionCentre'])
        instance = int(repeatDict[centreGene])
        if instance != 0:
            for variable in controlledVariables:
                # compare each variable to previous attempts same variables
                curr_val = str(falseRemoved.loc[index, variable])
                prev_val = str(falseRemoved.loc[index - 1, variable])
                if curr_val == 'nan' and prev_val == 'nan':
                    # nan's aren't ==, so need to explicitly check
                    row_differences.append('')
                elif curr_val != prev_val:
                    row_differences.append(variable)
                else:
                    # add blanks for spreadsheet
                    row_differences.append('')
        else:
            # blank row of differences for spreadsheet
            row_differences = [''] * len(controlledVariables)

        differences.append(row_differences)
        repeatDict[centreGene] = instance + 1

    controlledVariables = [u'Δ {0}'.format(element) for element in controlledVariables]  # cols need unique names
    falseRemoved[controlledVariables] = pd.DataFrame(differences, index=falseRemoved.index)
    falseRemoved = falseRemoved.drop(columns=['RepeatedGene'])  # they're all repeated so this is an unnecessary field
    # also drop the residual index column
    falseRemoved = falseRemoved.drop(columns=['index'])

    cleanData['ExperimentRepeated'] = ''
    cleanData.loc[cleanData['Mi Attempt URL'].isin(experiments_repeated), ['ExperimentRepeated']] = 't'

    return falseRemoved


def removeRepeatedGenes(derivedData, successfulGenes):
    """
    Want to keep only one experiment per gene, selected using the following criteria:
        - If only one was successful, keep that and remove the rest
        - If multiple were successful, keep the first success
        - If none were successful, keep the most recent failure
        - If none were failures (i.e. all 'Founder obtained'), keep the most recent

    However, if none were successful in this filtered dataset, this can falsely give the impression that no successful
    knockout exists for this gene. Using the successfulGenes list from the pre-filtered dataset, genes with any
    successful attempt in the initial iMits dataset are marked 't' in the column 'GLT attempt exits for gene'
    """
    repeatedGenesRemoved = derivedData.copy(deep=True)

    # get the 'groups' of attempts with the same gene
    geneColumn = repeatedGenesRemoved['Gene Marker Symbol']

    repeatedGenes = repeatedGenesRemoved[geneColumn.isin(geneColumn[geneColumn.duplicated()])]['Gene Marker Symbol'].tolist()
    repeatedGenes = list(set(repeatedGenes))

    for gene in repeatedGenes:

        # get the subset of attempts with that gene
        repeatedGeneSubset = repeatedGenesRemoved.loc[repeatedGenesRemoved['Gene Marker Symbol'] == str(gene)]
        repeatedGeneSubset.sort_values('datetime')  # sort by date
        status_counts = repeatedGeneSubset['Status Name'].value_counts()
        success_count = 0
        failure_count = 0
        progress_count = 0
        # assigning these to variables via loop to avoid missing key errors
        for val, cnt in status_counts.iteritems():
            if val == 'Genotype confirmed':
                success_count = int(cnt)
            elif val == 'Micro-injection aborted':
                failure_count = int(cnt)
            elif val == 'Founder obtained':
                progress_count = int(cnt)

        if success_count == 1:
            # remove all the failures
            for index, row in repeatedGeneSubset.iterrows():
                if row['Status Name'] != 'Genotype confirmed':
                    if VERBOSE:
                        print("Removing " + str(row['Mi Attempt URL']) + " as there is a successful attempt with the same gene")
                    repeatedGenesRemoved.drop(index, inplace=True)
        elif success_count > 1:
            # dataframe already sorted by date, set flag when first success seen and delete all others (and failures)
            seenSuccess = False
            for index, row in repeatedGeneSubset.iterrows():
                if row['Status Name'] == 'Genotype confirmed':
                    if not seenSuccess:
                        # don't remove but set the flag so the next one will be
                        seenSuccess = True
                    else:
                        if VERBOSE:
                            print("Removing " + str(row['Mi Attempt URL']) + " as there is an earlier successful attempt")
                        repeatedGenesRemoved.drop(index, inplace=True)
                else:
                    # also remove everything that's not a success
                    if VERBOSE:
                        print("Removing " + str(row['Mi Attempt URL']) + " as there are multiple other successful attempts")
                    repeatedGenesRemoved.drop(index, inplace=True)
        else:
            # if there are failures, keep the most recent
            if failure_count > 0:
                # sorted by date, remove all except latest failure (use failure count)
                failure_counter = 0
                for index, row in repeatedGeneSubset.iterrows():
                    # remove all the founder obtained
                    if row['Status Name'] == 'Founder obtained':
                        repeatedGenesRemoved.drop(index, inplace=True)
                    else:
                        failure_counter += 1  # it's a failure, increment
                        if failure_counter != failure_count:
                            # not the most recent, remove it
                            repeatedGenesRemoved.drop(index, inplace=True)
                        else:
                            # this is the most recent of the group, keep it
                            pass
            else:
                # ok there are no failures, keep the most recent founder obtained
                progress_counter = 0
                for index, row in repeatedGeneSubset.iterrows():
                    progress_counter += 1
                    if progress_counter != progress_count:
                        # not the most recent, remove it
                        repeatedGenesRemoved.drop(index, inplace=True)
                    else:
                        # this is the most recent of the group, keep it
                        pass

    # finally, mark all those for which a successful attempt exists for the gene
    repeatedGenesRemoved['SuccessfulAttemptExists'] = ''
    repeatedGenesRemoved.loc[repeatedGenesRemoved['Gene MGI Accession ID'].isin(successfulGenes), ['SuccessfulAttemptExists']] = 't'

    return repeatedGenesRemoved


def addGeneAnalysis(derivedData, geneInfo):
    '''
    Add in the annotation information related to the gene for an attempt
    Length	GCcontent	CpGsites	PercentageCpG
    '''

    gene_lengths = []
    GC_contents = []
    CpG_sites = []
    CpG_percents = []
    missing_count = 0
    if VERBOSE:
        print("Missing gene information for: ")

    for index, row in derivedData.iterrows():
        mgi_id = str(row['Gene MGI Accession ID'])
        geneRow = geneInfo.loc[geneInfo['MGI_ID'] == mgi_id]
        if geneRow.shape[0] == 1:
            gene_lengths.append(geneRow.iloc[0]['Length'])
            GC_contents.append(geneRow.iloc[0]['GCcontent'])
            CpG_sites.append(geneRow.iloc[0]['CpGsites'])
            CpG_percents.append(geneRow.iloc[0]['PercentageCpG'])
        elif geneRow.shape[0] > 1:
            gene_lengths.append('')
            GC_contents.append('')
            CpG_sites.append('')
            CpG_percents.append('')
            print("Error: duplicate entry for " + mgi_id)
        else:
            missing_count += 1
            gene_lengths.append('')
            GC_contents.append('')
            CpG_sites.append('')
            CpG_percents.append('')
            # option to print the symbols and mgi ids for genes we're missing
            if VERBOSE:
                print(str(row['Gene Marker Symbol']) + "\t" + str(mgi_id))

    if missing_count:
        print("Missing gene information for " + str(missing_count) + " genes")

    derivedData['Length'] = gene_lengths
    derivedData['GCcontent'] = GC_contents
    derivedData['CpGsites'] = CpG_sites
    derivedData['PercentageCpG'] = CpG_percents

    return


def addHumanOrthologs(derivedData, humanOrthologs):
    '''
    Using the accepted spreadsheet of mouse-human orthologs, add the human ortholog Symbol and ENSID to the dataset.
    All matching done on MGI ID
    '''
    no_ortholog_count = 0
    human_ortholog_symbols = []
    human_ortholog_ensids = []
    multiple_orthologs = []
    for index, row in derivedData.iterrows():
        mgi_id = str(row['Gene MGI Accession ID'])
        mouse_symbol = str(row['Gene Marker Symbol'])
        # find all the human genes orthologous to the mouse gene for this row
        orthologSubset = humanOrthologs.loc[humanOrthologs['mgi_id'] == mgi_id]
        if orthologSubset.shape[0] == 1:
            # exactly one human ortholog
            human_ortholog_symbols.append(str(orthologSubset.iloc[0]['symbol']))
            human_ortholog_ensids.append(str(orthologSubset.iloc[0]['ensembl_gene_id']))
        elif orthologSubset.shape[0] == 0:
            # print("NO HUMAN ORTHOLOGS FOR : " + mouse_symbol)
            no_ortholog_count += 1
            human_ortholog_symbols.append('')
            human_ortholog_ensids.append('')
        else:
            if VERBOSE: print("MULTIPLE HUMAN ORTHOLOGS FOR : " + mouse_symbol)
            multiple_orthologs.append(mouse_symbol)
            multiple_symbols = []
            multiple_ensids = []
            # see if any match on exact same symbol name (case-insensitive)
            foundMatch = False
            for subindex, subrow in orthologSubset.iterrows():
                if VERBOSE: print(str(subrow['symbol']))
                human_symbol = str(subrow['symbol'])
                human_ensid = str(subrow['ensembl_gene_id'])
                multiple_symbols.append(human_symbol)
                multiple_ensids.append(human_ensid)
                if human_symbol.lower() == mouse_symbol.lower():
                    #if VERBOSE: print("USING HUMAN SYMBOL WITH EXACT NAME MATCH: " + human_symbol)
                    human_ortholog_symbols.append(human_symbol)
                    human_ortholog_ensids.append(human_ensid)
                    foundMatch = True
            if not foundMatch:
                print("MULTIPLE HUMAN ORTHOLOGS FOR: " + str(mouse_symbol) + " (none with exact name match)")
                multiple_orthologs.append(mouse_symbol)
                print(", ".join(multiple_symbols))
                #print("Using first in list (" + str(multiple_symbols[0]) + ")")
                human_ortholog_symbols.append(str(multiple_symbols[0]))
                human_ortholog_ensids.append(str(multiple_ensids[0]))

    print("GET THIS: " + str(list(set(multiple_orthologs))))
    print("No human orthologs found for " + str(no_ortholog_count) + " mouse genes")
    derivedData['Human Ortholog'] = human_ortholog_symbols
    derivedData['Human ENSID'] = human_ortholog_ensids

    return


def addOrthologInfo(derivedData, humanOrthologs, gnomad):
    '''
    First map the mouse genes to their human orthologs. Record the relationship (i.e. one-to-one etc.) as well.
    Then, using the gnomad dataset, add in the pLI score and o/e score.
    '''

    human_ortholog_symbols = []
    human_ortholog_ensids = []
    human_ortholog_mapping = []
    pLI_scores = []
    oe_scores = []
    for index, row in derivedData.iterrows():
        mgi_id = str(row['Gene MGI Accession ID'])
        mouse_symbol = str(row['Gene Marker Symbol'])

        # find all the human genes orthologous to the mouse gene for this row
        orthologSubset = humanOrthologs.loc[humanOrthologs['mgi_id'] == mgi_id]
        if orthologSubset.shape[0] == 1 and str(orthologSubset.iloc[0]['human-to-mouse']) != 'one-to-many':

            # only one ortholog, add its pLI and o/e score from the gnomad dataset
            human_ensid = str(orthologSubset.iloc[0]['ensembl_gene_id'])  # used to get the pLI
            ortholog_row = gnomad.loc[gnomad['gene_id'] == human_ensid]
            if ortholog_row.shape[0] == 1:  # should always be one or 0
                pLI_scores.append(ortholog_row.iloc[0]['pLI'])
                oe_scores.append(ortholog_row.iloc[0]['oe_lof'])
            elif ortholog_row.shape[0] == 0:
                pLI_scores.append('')
                oe_scores.append('')
            else:
                print("ERROR: multiple results in gnomad dataset for gene with ensid: " + human_ensid)
            human_ortholog_symbols.append(str(orthologSubset.iloc[0]['symbol']))
            human_ortholog_ensids.append(human_ensid)
            human_ortholog_mapping.append(str(orthologSubset.iloc[0]['concatenate']))
        else:
            human_ortholog_symbols.append('NO ONE-TO-ONE ORTHOLOG')
            human_ortholog_ensids.append('')
            human_ortholog_mapping.append('')
            pLI_scores.append('')
            oe_scores.append('')

    derivedData['Human Ortholog'] = human_ortholog_symbols
    derivedData['Human ENSID'] = human_ortholog_ensids
    derivedData['Ortholog Relationship (human-to-mouse_mouse-to-human)'] = human_ortholog_mapping
    derivedData['pLI of Orthologs'] = pLI_scores
    derivedData['oe of Orthologs'] = oe_scores

    return


def addScores(derivedData, gnomadLof):
    '''
    Add the pLI scores from the gnomad dataset
    '''
    pLI_column = []
    for index, row in derivedData.iterrows():
        human_ensid = str(row['Human ENSID'])
        match_subset = gnomadLof.loc[gnomadLof['gene_id'] == human_ensid]
        if match_subset.shape[0] == 1:
            pLI = match_subset.iloc[0]['pLI']
            #pLI = round(pLI, 3)
            pLI_column.append(pLI)
        elif match_subset.shape[0] == 0:
            pLI_column.append('')
        else:
            print("ERROR: multiple results for gene with ensid: " + human_ensid)

    derivedData['pLI'] = pLI_column


def main():
    global VERBOSE
    VERBOSE = True
    engine='openpyxl'

    datePrefix = datetime.now().strftime("%m-%d-%Y")

    # using new data
    outputFile = datePrefix + "_IMPC_Cas9_2020-10-09"

    repeatedAttemptFile = datePrefix + "_Repeats-IMPC_Cas9_2020-10-09"
    # precautionary check, only written if duplicate attempts exist
    duplicateAttemptsFile = datePrefix + "_DuplicateAttemptsRemaining.xlsx"

    # read in raw data downloaded from IMITS (new as of May 29th, 2019) (replaced January 2020)
    imitsData = pd.ExcelFile(os.path.join(os.getcwd(), 'iMTS_AllCas9_20201009.xlsx'), engine='openpyxl')
    joinedData = joinData(imitsData)

    # read in data on all successful experiments from iMits
    gltExcel = pd.ExcelFile('input/20190610_mi_attempts_list.xlsx', engine='openpyxl')
    gltData = pd.read_excel(gltExcel, "All GLT")
    successGenes, successAttempts = getSuccesses(gltData)

    # add a datetime column to sort on
    joinedData['datetime'] = pd.to_datetime(joinedData['Mi Date'], format="%Y-%m-%d")
    joinedData = joinedData.sort_values('datetime')

    miURLS = joinedData["Mi Attempt URL"]
    duplicateAttempts = joinedData[miURLS.isin(miURLS[miURLS.duplicated()])]
    if VERBOSE:
        print("Initial # of records with duplicate attempts: " + str(len(duplicateAttempts.index)))

    print(f'{len(joinedData.index)} experiments initially ({len(duplicateAttempts.index)} are duplicates from merge)')

    # import the viability report which is used to add a viability column
    viabilityReport = pd.read_csv('input/all_genes_viability.csv')

    # add in corrections to data from BCM
    informationBCMExcel = pd.ExcelFile("input/NoFounder_Reasons/20190327_ExperimentsWithFoundersNoGLT_BCM.xlsx", engine='openpyxl')
    informationBCM = pd.read_excel(informationBCMExcel, '1+ G0 no GLT')
    addCorrections(joinedData, informationBCM)

    # filter out data we're not interested in or are missing preliminary data
    filteredData = filterData(joinedData, viabilityReport)

    # remove the non protein-coding genes
    proteinCodingGenes = pd.read_csv('input/protein-coding_MGI_20190604.csv')
    filteredData = removeNonProteinCoding(filteredData, proteinCodingGenes)

    # remove duplicates (complete duplicates this is, exact same attempt and colony
    filteredData_duplicatesRemoved = filteredData.loc[filteredData.astype(str).drop_duplicates().index]

    # import the report on essentiality
    essentialityFile = pd.ExcelFile('input/HumanOrtholog_Info/HumanEssentiality_MouseOrthologs.xlsx', engine='openpyxl')
    essentialityReport = pd.read_excel(essentialityFile, 'Sheet1')

    # add derived and calculated columns for analysis
    derivedData = addCalculatedColumns(filteredData_duplicatesRemoved, viabilityReport, essentialityReport)

    # add gene-related information
    geneInfoFile = pd.ExcelFile('input/AnnotationInfo.xlsx', engine='openpyxl')
    geneInfo = pd.read_excel(geneInfoFile, 'Sheet1')
    addGeneAnalysis(derivedData, geneInfo)

    # add human orthologs and ensembl gene ids for the human genes
    humanOrthologsExcel = pd.ExcelFile("input/HumanOrtholog_Info/Human_mouse_orthologues_for_IMPC.xlsx", engine='openpyxl')
    humanOrthologs = pd.read_excel(humanOrthologsExcel, "All_protCoding_score>=5")

    # add gnomad dataset for pLI and oe scores
    gnomadExcel = pd.ExcelFile("input/HumanOrtholog_Info/gnomad_v2_1_1_lof_metrics_by_gene.xlsx", engine='openpyxl')
    gnomad = pd.read_excel(gnomadExcel, 'Sheet1')

    # add the human orthologs and pLI/oe scores for the orthologs from the gnomad dataset
    # for cases where there are multiple orthologs, take the average of the scores
    addOrthologInfo(derivedData, humanOrthologs, gnomad)

    # outputting experiments with founder pups born but no GLT for questions prior to removal of columns
    foundersNoGLT = derivedData.loc[(derivedData['#G0 deletion event detected'] > 0) & (derivedData['GLT'] == 'f')]
    foundersNoGLT = foundersNoGLT.sort_values(by=['Production Centre', 'Mi Date'])
    foundersNoGLT.to_excel(os.path.join("output", "ExperimentsWithFoundersNoGLT.xlsx"), index=False)

    noFounderDir = "input/NoFounder_Reasons"
    derivedData['Reason GLT Failed'] = ''  # add column to track reason
    for GLT_Failure_Excel in os.listdir(noFounderDir):
        gltInfo = pd.read_excel(pd.ExcelFile(os.path.join(noFounderDir, GLT_Failure_Excel)), '1+ G0 no GLT', engine='openpyxl')
        addReasonGLTFailed(derivedData, gltInfo)

    # drop the columns that may have been used for filtering but aren't relevant to analysis
    fullDrop = ['Consortium',
                'IS Active',
                'Report Micro Injection Progress To Public',
                'Experimental (exclude from production reports)',
                'Plasmid Generated Guides',
                'Truncated_guide',
                'Individually set gRNA concentrations',
                'Truncated_guide.1',
                'Individually set gRNA concentrations.1',
                'Truncated_guide.2',
                'Individually set gRNA concentrations.2',
                'Truncated_guide.3',
                'Individually set gRNA concentrations.3',
                'mRNA Nuclease',
                'Protein Nuclease',
                'Report F1 Colony To Public',
                'Voltage',
                '#Pulses',
                'Embryo Transfer Day',
                '#Embryos Survived to 2 cell stage',
                'Assay Carried Out',
                '#G0 with detected mutation',
                '#G0 NHEJ event detected',
                '#G0 HR event detected',
                '#G0 HDR event detected',
                '#G0 all donor insertions detected',
                '#G0 subset of donors inserted detected',
                'Genotype Confirmed',
                'Background Strain',
                'Viability Counts',
                # 'datetime'  # this was used to sort the attempts
                ]
    cleanData = derivedData.drop(columns=fullDrop)
    print(f'\nRemoved {len(fullDrop)} columns irrelevant to analysis\n')
    if VERBOSE:
        for column in fullDrop:
            print(column)

    # export the attempts that produced two colonies after filtering
    miURLS = cleanData["Mi Attempt URL"]
    duplicateAttempts = cleanData[miURLS.isin(miURLS[miURLS.duplicated()])]
    print("Records with duplicate attempts after filtering: " + str(len(duplicateAttempts.index)))
    # only export if duplicates exist
    if len(duplicateAttempts.index) > 0:
        duplicateAttempts.to_excel("output/" + duplicateAttemptsFile, index=False)
        # now remove them, they will be identical records
        # excepting potentially the Allele Type/Subtype which aren't used in analysis anyways
        print("(They will be removed)")
        cleanData = cleanData.drop_duplicates(subset='Mi Attempt URL')
        # N.B. this is done before the repeated genes are removed since these attempts will be incorrectly marked
        # as repeated otherwise

    # the experiments that are repeats are exported in an excel which highlights what was changed
    repeatedAttempts = repeatAnalysis(cleanData, successAttempts)  # cleanData also has attempts that are repeat experiments marked

    repeatedAttempts = repeatedAttempts.drop(['datetime'], axis=1)  # remove the datetime column
    repeatedAttempts.to_excel(os.path.join("output", repeatedAttemptFile + ".xlsx"), index=False)
    repeatedAttempts.to_csv(os.path.join("output", repeatedAttemptFile + ".csv"), index=False)

    # perform filtering to only keep one attempt per gene, also mark those which have a filtered success
    repeatedGenesRemoved = removeRepeatedGenes(cleanData, successGenes)

    # export clean data with the repeats removed
    repeatedGenesRemoved = repeatedGenesRemoved.drop(['datetime'], axis=1)  # remove the datetime column
    repeatedGenesRemoved.to_excel(os.path.join("output", outputFile + ".xlsx"), index=False)
    repeatedGenesRemoved.to_csv(os.path.join("output", outputFile + ".csv"), index=False)
    print(f'\n{len(repeatedGenesRemoved.index)} experiments remaining after all cleaning and filtering')

    # also export the attempts that produced two colonies after filtering
    miURLS = cleanData["Mi Attempt URL"]
    duplicateAttempts = cleanData[miURLS.isin(miURLS[miURLS.duplicated()])]
    print("Records with duplicate attempts after filtering: " + str(len(duplicateAttempts.index)))
    # only export if duplicates exist
    if len(duplicateAttempts.index) > 0:
        duplicateAttempts.to_excel("output/" + duplicateAttemptsFile, index=False)


if __name__ == "__main__":
    main()
