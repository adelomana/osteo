# install from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dorothea")

# install the development version from GitHub
# install.packages("devtools")
devtools::install_github("saezlab/dorothea")

library(dorothea)
df = dorothea_hs
write.table(df, "/home/adrian/projects/osteo/data/dorothea/database.txt", sep="\t", quote=FALSE)


viper_regulons = df2regulon(dorothea_hs)
