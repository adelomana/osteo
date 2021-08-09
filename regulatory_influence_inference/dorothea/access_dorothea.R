# install from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dorothea")
# 
# install the development version from GitHub
#install.packages("devtools")
#devtools::install_github("saezlab/dorothea")

library(dorothea)

# 0. user-defined variables
database_output_dir = '/Users/adrian/gd15/hi/research/osteo/data/dorothea/'
regulons_output_dir = '/Users/adrian/gd15/hi/research/osteo/data/dorothea/regulons/'

# 1. write all interactions
df = dorothea_hs
output_file = paste(database_output_dir, 'database.txt', sep='')
write.table(df, output_file, sep="\t", quote=FALSE)

### 2. write database in the form of regulons
viper_regulons = df2regulon(dorothea_hs)

for (tf in names(viper_regulons)) {
  print(tf)
  targets = names(viper_regulons[[tf]]$tfmode)
  output_file = paste(regulons_output_dir, tf, '.txt', sep='')
  write(targets, output_file)
}
