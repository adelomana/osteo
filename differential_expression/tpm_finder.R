#
# 1. load libraries
#
library(DESeq2)                
library(tximport)  

# genome annotation
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

metadata_file = '/home/adrian/projects/osteo/metadata/MSC_Sample list.csv'
results_dir = '/home/adrian/projects/osteo/results/tpm/'
kallisto_dir = "/home/adrian/projects/osteo/data/kallisto_run"

#
# 2. read metadata
#
metadata = read.table(metadata_file, header = TRUE, sep = ",")
metadata

#
# 3. build annotation reference
# 
k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
tx2gene = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")

#
# 4. read files
#
files = file.path(kallisto_dir, metadata$FILE.NAME, 'abundance.h5')
print(files)
txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)


# 5. find abundance
tpm = txi$abundance
colnames(tpm) = metadata$FILE.NAME
dim(tpm)

# 6. convert annotation into genes
ensembl_ids = rownames(tpm)
annotation = select(EnsDb.Hsapiens.v86, ensembl_ids, c("GENEBIOTYPE", "GENENAME"), "GENEID")

# 7. store
store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

store = paste(results_dir, 'annotation.tsv', sep='')
write.table(annotation, file=store, quote=FALSE, sep='\t', col.names=NA)