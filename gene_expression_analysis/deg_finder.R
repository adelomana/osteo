###
### 0. installation
###
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::version()
# 
# BiocManager::install('openssl')
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("readr")
# BiocManager::install("rlist")
# BiocManager::install('EnsDb.Hsapiens.v86')
# BiocManager::install('org.Hs.eg.db')
# BiocManager::install('rhdf5')

###
### 1. load libraries
###
library(DESeq2)                
library(tximport)         
library(readr)
library(rlist) # required to append to lists
library(rhdf5)

# genome annotation
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# performance
library(BiocParallel)
library(tictoc)

###
### 2. user defined variables
###
tic()
register(MulticoreParam(20))
setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/osteo/data/kallisto_run"
metadata_file = '/home/adrian/projects/osteo/metadata/MSC_Sample list.csv'
results_dir = '/home/adrian/projects/osteo/results/DETs_DESeq2/'

###
### 4. build annotation reference
### 
k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
tx2gene = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")

###
### 3. read metadata
###
metadata = read.table(metadata_file, header = TRUE, sep = ",")
metadata

###
### 6. build hypotheses
### 
donors = c(6, 16, 17)
days = c('Day_3', 'Day_6', 'Day_9', 'Day_16', 'Day_28')

for (donor in donors) {
  for (day in days[2:length(days)]){
    
    tag=paste('donor', donor, day, sep='_')
    print(tag)
    
    reference_samples = metadata[metadata$TIME == 'Day_3' & metadata$DONOR == donor, ]$FILE.NAME
    print(reference_samples)
    testing_samples = metadata[metadata$TIME == day & metadata$DONOR == donor, ]$FILE.NAME
    print(testing_samples)
    all_samples = c(reference_samples, testing_samples)
    print(all_samples)
    if (length(all_samples) < 5) {
      message('bypassing because of lack of data...')
      next
    } 
    hypothesis_metadata = subset(metadata, FILE.NAME %in% all_samples, select = c("FILE.NAME", "TIME"))
    
    print(hypothesis_metadata)

    ### 6.1. read files
    files = file.path(kallisto_dir, all_samples, 'abundance.h5')
    print(files)
    txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)
    
    ### 6.2. build object
    dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~TIME)
    dds$TIME = relevel(dds$TIME, ref = "Day_3")
    
    ### 6.3. filtering
    threshold = 10
    keep = rowSums(counts(dds)) >= threshold  
    print(paste('dimensions before filtering', dim(dds)[1], dim(dds)[2]))
    dds = dds[keep, ]
    print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))
    
    ### 6.4. analysis
    print('analysis')
    dds = DESeq(dds, parallel=TRUE)
    
    ### 6.5. filter, annotate, format and store
    print('filter')
    res = results(dds, lfcThreshold=1, parallel=TRUE)
    filt1 = res[which(res$pvalue < 0.05), ]
    filt2 = filt1[which(filt1$padj < 0.1), ]
    print(paste('DEGs found', dim(filt2)[1], sep=' '))
    
    print('annotate...')
    df = as.data.frame(filt2)
    df['common'] = rownames(df)
    selected = rownames(df)
    info = select(EnsDb.Hsapiens.v86, selected, c("GENEBIOTYPE", "GENENAME"), "GENEID")
    info['common'] = info$GENEID
    
    descriptions = tryCatch({
      descriptions = select(org.Hs.eg.db, keys=selected, columns=c("GENENAME"), keytype="ENSEMBL")
    }, error = function(e) {
      print('Warning: no description found for ENSEMBL IDs')
      descriptions = data.frame('ENSEMBL'=selected, 'GENENAME'=rep('Not found', each=length(selected)))
    })
    #print(descriptions)
    
    names(descriptions)[names(descriptions) == "GENENAME"] <- "DESCRIPTION" # arrow is needed here!
    descriptions['common'] = descriptions$ENSEMBL
    dh = merge(df, info, by='common')
    di = merge(dh, descriptions, by='common')
    
    print('format')
    formatted = di[ , c(8,10,9,12,2,3,6,7)]
    up = formatted[formatted$log2FoldChange > 0, ]
    down = formatted[formatted$log2FoldChange < 0, ]
    sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
    sorted_down = down[order(down$log2FoldChange), ]
    
    print('store')
    store = paste(results_dir, tag, '_up', '.tsv', sep='')
    write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
    store = paste(results_dir, tag, '_down', '.tsv', sep='')
    write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)
    
    print('---')
    
  }
  print('-----------------------------------------')
} 
toc()

