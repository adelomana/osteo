'''

This script checks if regulons have consistent expression.

The steps are the following:

iterate over regulons.
plot row expression.
plot the log2 expression.
scale the data. Consider z-scores and real scaling within a defined range.
run PCA. Check how much variance is explained by PC1. This is the key to accept regulon tightness.
use scores as values for heatmap.
optionally, check the correlation of TF and PC1 to reveal if the TF is transcriptionally or posttranscriptionally regulated, e.g., phosphorylation.

'''

import pandas, numpy

import scipy, scipy.stats
import sklearn, sklearn.decomposition

import matplotlib, matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20, 'font.family':'FreeSans', 'xtick.labelsize':20, 'ytick.labelsize':20, 'figure.figsize':(8, 4)})

import pyensembl
annotation = pyensembl.EnsemblRelease(86) # better matching than version 100


def metadata_reader():

    metadata = {}
    donors = []
    days = []

    with open(metadata_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split(',')
            sample = v[0]
            day = v[-2]
            donor = v[-3]
            metadata[sample] = (donor, day)

            if donor not in donors:
                donors.append(donor)
            if day not in days:
                days.append(day)

    return metadata, donors, days
# 1. get all TFs that regulate DETs

# read all

# for each TF, do we have a hypergeometric enrichment?


###
### 0. user-defined functions
###

regulon_input_file = '/home/adrian/projects/osteo/results/grn/regulons.txt'
expression_data_input_file = '/home/adrian/projects/osteo/results/tpm/DESeq2_TPM_values.tsv'
metadata_file = '/home/adrian/projects/osteo/metadata/MSC_Sample list.csv'
DETs_input_file = '/home/adrian/projects/osteo/results/DETs_filtered/entire_DET_set.txt'
sub_regulon_file = '/home/adrian/projects/osteo/results/grn/subregulons.txt'

donors = ['6', '16', '17']
days = ['Day_3', 'Day_6', 'Day_9', 'Day_16', 'Day_28']


###
### 1. read data
###

##
## 1.1. read regulon memberships
##

df = pandas.read_csv(regulon_input_file, sep='\t', index_col=0)

##
## 1.2. read expression data
##

metadata, donors, days = metadata_reader()

expression = {}
ensembl_IDs = []

for sample in metadata:
    expression[sample] = {}

with open(expression_data_input_file, 'r') as f:

    header = f.readline()
    v = header.split('\t')
    sample_IDs = v[1:]
    sample_IDs[-1] = sample_IDs[-1].replace('\n', '')

    for line in f:
        v = line.split('\t')

        gene_name = v[0]
        ensembl_IDs.append(gene_name)
        v = [float(element) for element in v[1:]]

        for i in range(len(v)):
            value = v[i]
            sampleID = sample_IDs[i]
            expression[sampleID][gene_name] = value

###
### Let's work with expression trajectories for only the genes we care: DETs and TFs
###

# retrieve DETs for expression
DETs = []
with open(DETs_input_file, 'r') as f:
    for line in f:
        element = line.replace('\n', '')
        DETs.append(element)
print('DETs found: {}'.format(len(DETs)))
print(len(list(set(DETs))))

# add TFs, in case we work with their expression later
TFs = df.index.to_list()
print(len(TFs), TFs[:10])
TF_ensembl_ids = [annotation.gene_ids_of_gene_name(TF)[0] for TF in TFs]
print(len(TF_ensembl_ids), TF_ensembl_ids)

working_genes = list(set(DETs + TF_ensembl_ids))

print(len(working_genes))


# we will work from this point with gene symbols only, converting ENSEMBL to gene symbols
manual_annotation = {
    'ENSG00000274619':'CFD',
    'ENSG00000278843':'MMP28',
    'ENSG00000227746':'C4A',
    'ENSG00000282854':'ASAP3',
    'ENSG00000283009':'ZNF436',
    'ENSG00000283106':'ZNF436-AS1',
    'ENSG00000263238':'CTSO',
    'ENSG00000275482':'EPHB6',
    'ENSG00000262862':'KRTAP2-3',
    'ENSG00000277909':'SPC25',
    'ENSG00000281166':'NLRP10',
    'ENSG00000280610':'KIF15',
    'ENSG00000280682':'HYOU1',
    'ENSG00000280641':'CDH4',
    'ENSG00000273707':'CDKN1C',
    'ENSG00000233450':'KIFC1',
    'ENSG00000262634':'SKA1',
    'ENSG00000275125':'PRXL2B',
    'ENSG00000281165':'SLC2A6',
    'ENSG00000276657':'PYCR3',
    'ENSG00000282147':'FAM20C'
}

# convert ensembl to gene symbol
for working_gene in working_genes:

    # convert ensembl to gene symbol
    try:
        name = annotation.gene_name_of_gene_id(working_gene)
    except:
        print('Ensembl ID {} not found in annotation.'.format(working_gene))
        if working_gene in manual_annotation:
            name = manual_annotation[working_gene]
            print('\t present in manual annotation, all OK.')
        else:
            print('WE HAVE A PROBLEM')
        print()

# obtain a TPM trajectory (discrete for very low) across time over replicates and patients
expression_trajectories = {}

for working_gene in working_genes:

    # convert ensembl to gene symbol
    try:
        name = annotation.gene_name_of_gene_id(working_gene)
    except:
        print('Ensembl ID {} not found in annotation.'.format(working_gene))
        if working_gene in manual_annotation:
            name = manual_annotation[working_gene]
            print('\t present in manual annotation, all OK.')
        else:
            print('WE HAVE A PROBLEM')
    ############################

    ensembl_trajectory = []
    for day in days:
        expression_across_donors = []
        for donor in donors:

            # retrieve sample labels
            sample_labels = []
            for sample in sample_IDs:
                if metadata[sample] == (donor, day):
                    sample_labels.append(sample)

            # retrieve expression
            try:
                exp = [expression[label][working_gene] for label in sample_labels]
            except:
                exp = [0 for label in sample_labels]
                print('ENSEMBL ID {} not found in expression data set. Setting expression to zero.'.format(working_gene))

            average_expression = numpy.mean(exp)
            if numpy.isnan(average_expression) == False:
                expression_across_donors.append(average_expression)

        # compute median across patients
        median_expression = numpy.median(expression_across_donors)

        # round values, useful for very low values
        round_value = numpy.around(median_expression) + 1

        # add round value to trajectory
        ensembl_trajectory.append(round_value)

    # add trajectory to dictionary
    expression_trajectories[name] = ensembl_trajectory

    print(name, ensembl_trajectory)

###
### 2. analysis
###

time_trajectory = [int(day.split('_')[-1]) for day in days]
print(time_trajectory)
theColors = ['blue', 'green', 'gold', 'orange', 'red']
print(theColors)

for TF, row in df.iterrows():

    number_of_targets = int(row['TF targets DEGs / all TF targets'].split('/')[0])
    print(TF, number_of_targets)

    # plot raw expression of all genes
    targets = row['Observed TF target DEGs'].split(', ')
    print(targets)

    regulon_expression = []
    for target in targets:
        regulon_expression.append(expression_trajectories[target])
    #print(regulon_expression)

    log2_regulon_expression = numpy.log2(regulon_expression)
    #print('log2 reg expression')
    #print(log2_regulon_expression)

    # zscores
    #print('zscores')
    zscores_regulon = scipy.stats.zscore(log2_regulon_expression, axis=1)
    #print(zscores_regulon)
    for gene_trajectory in zscores_regulon:
        matplotlib.pyplot.plot(time_trajectory, gene_trajectory, 'o-', alpha=0.5)
    matplotlib.pyplot.xlabel('Time (day)')
    matplotlib.pyplot.ylabel('zscore')
    matplotlib.pyplot.xticks(time_trajectory)
    matplotlib.pyplot.title('original')
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('original.pdf')

    ###
    ### exhaustive PCA
    ###

    working_data = numpy.transpose(zscores_regulon)
    working_targets = targets
    print('original data shape:', working_data.shape)

    subregulons = []
    consistency_threshold = 0.85
    exhausted = False
    next_round_data = []
    next_round_targets = []
    while exhausted == False:

        iteration = 0
        PC1_explained_variance = 0



        # remove one by one until consistency threshold
        new_round = False
        while new_round == False:
            #print('iteration {}; shape: {}'.format(iteration, working_data.shape))
            # find lowest loading gene and remove it
            pca = sklearn.decomposition.PCA(n_components=3)
            projection = pca.fit_transform(working_data)
            PC1_explained_variance = pca.explained_variance_ratio_[0]
            #print('\t PC1 variance explained: {}'.format(PC1_explained_variance))
            PC1_loadings = pca.components_[0]
            lowest_loading_index = numpy.argmin(PC1_loadings)
            lowest_loading_gene = working_targets[lowest_loading_index]

            next_round_data.append(list(working_data[:, lowest_loading_index]))
            next_round_targets.append(lowest_loading_gene)
            working_data = numpy.delete(working_data, lowest_loading_index, 1)
            working_targets.remove(lowest_loading_gene)

            #print(working_data.shape)
            #print(next_round_data.shape)

            '''
            for gene_trajectory in numpy.transpose(working_data):
                matplotlib.pyplot.plot(time_trajectory, gene_trajectory, 'o-', alpha=0.5)
            matplotlib.pyplot.xlabel('Time (day)')
            matplotlib.pyplot.ylabel('zscore')
            matplotlib.pyplot.xticks(time_trajectory)
            matplotlib.pyplot.grid(alpha=0.5, ls=':')
            matplotlib.pyplot.ylim([-2, 2])
            matplotlib.pyplot.title(PC1_explained_variance)
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.show()
            '''

            iteration = iteration + 1

            if (PC1_explained_variance >= consistency_threshold) | (len(working_targets) <= 4):
                new_round = True

                if PC1_explained_variance >= consistency_threshold:
                    subregulons.append(working_targets)
                    print('found subregulon: {} {}'.format(working_targets, len(working_targets)))
                else:
                    print('WARNING: These genes are not tight, discarding.'.format(working_genes))

                working_data = numpy.transpose(numpy.array(next_round_data))
                working_targets = next_round_targets
                next_round_data = []
                next_round_targets = []



        #
        if len(working_targets) <= 5:
            exhausted = True


    print('################################## FINAL')
    working_data = numpy.transpose(zscores_regulon)
    working_targets = targets = row['Observed TF target DEGs'].split(', ')
    original = len(working_targets)
    new = sum([len(sub) for sub in subregulons])
    lost = original - new
    print('original size: {}; new size {}; lost: {}'.format(original, new, lost))


    for i in range(len(subregulons)):

        sub=subregulons[i]
        print(sub)
        indices = [working_targets.index(element) for element in sub]
        plotting = working_data[:, indices]

        for gene_trajectory in numpy.transpose(plotting):
                matplotlib.pyplot.plot(time_trajectory, gene_trajectory, 'o-', alpha=0.5)
        matplotlib.pyplot.xlabel('Time (day)')
        matplotlib.pyplot.ylabel('zscore')
        matplotlib.pyplot.xticks(time_trajectory)
        matplotlib.pyplot.grid(alpha=0.5, ls=':')
        matplotlib.pyplot.ylim([-2, 2])
        matplotlib.pyplot.title(len(sub))
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig('sub.{}{}.pdf'.format(TF, i))

        label = '{} | {}'.format(TF, str(i+1))
        average = numpy.mean(plotting, axis=1)
        print(label, average)

        #heatmap[label] = average
    # write results
