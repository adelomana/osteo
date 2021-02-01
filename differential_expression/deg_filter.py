import sys

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
            metadata[sample] = {'donor':donor, 'day':day}

            if donor not in donors:
                donors.append(donor)
            if day not in days:
                days.append(day)

    return metadata, donors, days

def read_DEGs(DEGs, donor, day, trend):

    '''
    Returns a dictionary of DEGs.
    '''

    working_file = DESeq2_folder + donor + '_' + day + '_' + '_' + trend + '.tsv'
    print(working_file)
    sys.exit()

    with open(working_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')

            ensembl = v[0]
            gene_name = v[1]
            biotype = v[2]
            description = v[3]
            basemean = float(v[4])
            logFC = float(v[5])
            pvalue = float(v[6])
            adjusted = float(v[7])

            info = (ensembl, gene_name, biotype, description, basemean, logFC, pvalue, adjusted)

            DEGs[experiment][concentration][time][trend].append(info)

    print('{} {} {} {} \t DEGs: {}'.format(experiment, concentration, time, trend, len(DEGs[experiment][concentration][time][trend])))

    return DEGs

###
### MAIN
###

#
# 0. use-defined variables
#

# paths
project_dir = '/home/adrian/projects/osteo/'
DESeq2_folder = project_dir + 'results/DETs_DESeq2/'
metadata_file = project_dir +'metadata/MSC_Sample list.csv'
expression_file = project_dir + 'tpm/DESeq2_TPM_values.tsv'
filtered_folder = project_dir + 'results/DETs_filtered/'

donors = ['6', '16', '17']
days = ['Day_3', 'Day_6', 'Day_9', 'Day_16', 'Day_28']
trend_tags = ['up', 'down']

#
# 1. read data
#

# 1.1. define the DEGs across experimental design
print('define DEGs')
DEGs = {}

for donor in donors:
    DEGs[donor] = {}
    for day in days[1:]:
        DEGs[donor][day] = {}
        for trend in trend_tags:
            DEGs[donor][day][trend] = []
            DEGs = read_DEGs(DEGs, donor, day, trend)

# 1.1. read metadata
print('read metadata')
metadata, donors, days = metadata_reader()

# 1.2. read expression data
print('read expression')

expression = {}
for sample in metadata:
    expression[sample] = {}

with open(expression_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')

        gene_name = v[0]
        v = [float(element) for element in v[1:]]

        for i in range(len(v)):
            value = v[i]
            sampleID = sample_IDs[i]
            expression[sampleID][gene_name] = value
