import sys, os, numpy

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

def read_DEGs(DEGs, donor, day, trend):

    '''
    Returns a dictionary of DEGs.
    '''

    working_file = DESeq2_folder + 'donor_' + donor + '_' + day + '_' + trend + '.tsv'

    if os.path.exists(working_file) == False:
        print('\t skipping donor {} {} {} because it does not exist'.format(donor, day, trend))
        del DEGs[donor][day][trend]
    else:
        all_ensembls = []
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

                if ensembl not in all_ensembls:
                    DEGs[donor][day][trend].append(info)
                    all_ensembls.append(ensembl)
                else:
                    print('\t\t ensembl {} found multiple times, skipping'.format(ensembl))

        print('\t donor {} {} {} \t DEGs: {}'.format(donor, day, trend, len(DEGs[donor][day][trend])))

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
expression_file = project_dir + 'results/tpm/DESeq2_TPM_values.tsv'
filtered_folder = project_dir + 'results/DETs_filtered/'

donors = ['6', '16', '17']
days = ['Day_3', 'Day_6', 'Day_9', 'Day_16', 'Day_28']
trend_tags = ['up', 'down']

expression_threshold = 2
discrete_fc_threshold = 1
noise_threshold = 1/2

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

# 1.2. read metadata
print('read metadata')
metadata, donors, days = metadata_reader()

# 1.3. read expression data
print('read expression')

expression = {}
for sample in metadata:
    expression[sample] = {}

with open(expression_file, 'r') as f:

    header = f.readline()
    v = header.split('\t')
    sample_IDs = v[1:]
    sample_IDs[-1] = sample_IDs[-1].replace('\n', '')

    for line in f:
        v = line.split('\t')

        gene_name = v[0]
        v = [float(element) for element in v[1:]]

        for i in range(len(v)):
            value = v[i]
            sampleID = sample_IDs[i]
            expression[sampleID][gene_name] = value

#
# 2. filter genes
#
print('filter options:')
print('\t expression threshold: {}'.format(expression_threshold))
print('\t discrete FC threshold: {}'.format(discrete_fc_threshold))
print('\t noise threshold: {}'.format(noise_threshold))

print('filter DEGs')
union = {}
for donor in DEGs:
    for day in DEGs[donor]:
        if day is not 'Day_3':
            if day not in union:
                union[day] = {}
            for trend in DEGs[donor][day]:
                if trend not in union[day]:
                    union[day][trend] = []

                # define sample and reference sample labels
                sample_labels = []; reference_labels = []
                for sample in sample_IDs:
                    if metadata[sample] == (donor, day):
                        sample_labels.append(sample)

                    if metadata[sample] == (donor, 'Day_3'):
                        reference_labels.append(sample)

                print('INFO: donor {}, day {}, trend {}'.format(donor, day, trend))
                print('INFO: sample labels {}'.format(sample_labels))
                print('INFO: reference labels {}'.format(reference_labels))
                print()

                # filters
                container = []
                before = len(DEGs[donor][day][trend])
                for case in DEGs[donor][day][trend]:
                    including = True
                    ensembl = case[0]

                    # gather TPM expression
                    ref = [expression[label][ensembl] for label in reference_labels]
                    sam = [expression[label][ensembl] for label in sample_labels]

                    # filter 1: identify low-expressed genes
                    r = numpy.median(ref); s = numpy.median(sam)
                    top = numpy.max([r, s])

                    # filter 2: identify fold-changes using discrete values
                    ###
                    ###            [round(x, epsilon)/epsilon ] + 1
                    ###  FC = abs  -------------------------------- > 1
                    ###            [round(y, epsilon)/epsilon ] + 1
                    ###
                    ###
                    ###  epsilon = 1
                    num = numpy.around(s) + 1
                    den = numpy.around(r) + 1
                    fc = num/den
                    abs_log2FC = numpy.abs(numpy.log2(fc))

                    # filter 3: identify noisy genes
                    ref_int = numpy.around(ref) + 1
                    sam_int = numpy.around(sam) + 1
                    sem_ref = numpy.std(ref_int) / numpy.sqrt(len(ref_int))
                    rsem_ref = sem_ref / numpy.mean(ref_int)
                    sem_sam = numpy.std(sam_int) / numpy.sqrt(len(sam_int))
                    rsem_sam = sem_sam / numpy.mean(sam_int)
                    noise = numpy.max([rsem_ref, rsem_sam])

                    # selection
                    if abs_log2FC < discrete_fc_threshold:
                        including = False
                        info = '\t WARNING: small change gene discarded. Expression changes from {:.3f} ({}) to {:.3f} ({}), resulting in abs_log2FC {:.3f}. {}, {}'.format(r, den, s, num, abs_log2FC, case[1], case[3])
                        print(info)

                    if (including == True) and (top < expression_threshold):
                        including = False
                        info = '\t WARNING: low-expression gene discarded. Expression changes from {:.3f} to {:.3f}. {}, {}'.format(r, s, case[1], case[3])
                        print(info)

                    if (including == True) and (noise > noise_threshold):
                        including = False
                        info = '\t WARNING: noisy gene. Ref: {}, RSEM {:.3f}; Sam: {}, RSEM {:.3f}. {}, {}'.format(ref, rsem_ref, sam, rsem_sam, case[1], case[3])
                        print(info)

                    if including == True:
                        content = list(case)
                        content.append(r); content.append(s); content.append(abs_log2FC)
                        container.append(content)
                        union[day][trend].append(case[0])

                # info about low-expression cases
                after = len(container)
                print('{} {} {} | DEGs final reduction from {} to {}\n'.format(donor, day, trend, before, after))

                # store reduced version
                DEGs[donor][day][trend] = container


### 2.4. learn genes responding in at least two patients
intersection = {}
for day in union:
    intersection[day] = {}
    for trend in union[day]:
        intersection[day][trend] = []
        unique = list(set(union[day][trend]))

        for element in unique:
            if union[day][trend].count(element) > 1:
                intersection[day][trend].append(element)
        print('{}; {}; all genes found:{}, unique list: {}; found in two patients: {}'.format(day, trend, len(union[day][trend]), len(unique), len(intersection[day][trend])))

# 3. store
## 3.1. store filtered DEGs per donor and condition
for donor in DEGs:
    for day in DEGs[donor]:

        for trend in DEGs[donor][day]:

            storage = filtered_folder + '{}_{}_{}_filtered.tsv'.format(donor, day, trend)
            with open(storage, 'w') as f:
                f.write('{}_{}_{}\n'.format(donor, day, trend))
                f.write('ENSEMBL\tGene name\tBiotype\tDescription\tBase mean\tlog2FC\tP value\tAdjusted P-value\tReference expression (TPM)\tSample expression (TPM)\tDiscrete abs(log2FC)\tIn at least two patients\n')
                for content in DEGs[donor][day][trend]:

                    if content[0] in intersection[day][trend]:
                        content.append('yes')
                    else:
                        content.append('no')

                    line = ''
                    for element in content:
                        if isinstance(element, str) == False:
                            if numpy.abs(element) < 0.05 and element != 0.:
                                sub = '{:.4E}'.format(element)
                            else:
                                sub = '{:.4f}'.format(element)
                        else:
                            sub = element
                        line=line+sub+'\t'
                    line=line+'\n'
                    f.write(line)
