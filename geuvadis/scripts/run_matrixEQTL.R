# main matrixEQTL script

# adapted from the FastGxC code

# libraries
# chooseCRANmirror(78)
# install.packages("MatrixEQTL")


# command line arguments
args=commandArgs(trailingOnly=TRUE)
pop=args[1]
N_PCs=args[2]
var_method=args[3]
exp_PCs=args[4]
covariate_type=args[5]
# expression_file_name = args[2]
# gene_location_file_name = args[3]
# SNP_file_name = args[3]
# snps_location_file_name = args[4]
# covariate_file_name = args[4]
# output_file_name_cis = args[5]
# SNP_file_name = args[6]

library(MatrixEQTL);
library(data.table);
library(glue);

# file names
expression_file_name = glue('../expressionData/{pop}_{var_method}_{exp_PCs}_expression_df_T.tsv'); # corrected expression
gene_location_file_name = glue('../expressionData/{pop}_gene_loc.tsv'); # gene positions
SNP_file_name = glue('../genotypeData/{pop}_T_snps.txt'); # SNP by sample, genotype file
snps_location_file_name = glue('../genotypeData/{pop}_snps_loc.txt'); # SNP location file

G_eigenvec_file_name = glue('../genotypeData/{pop}.eigenvec.T')
covariate_file_name = glue('../genotypeData/{pop}.{covariate_type}_covariates.T')
output_file_name_cis = glue('../results/{pop}_cis_eQTLs_{N_PCs}gPCs_{var_method}_{exp_PCs}ePCs_{covariate_type}'); #output for cis eQTLs

output_file_name_tra = tempfile();

# MatrixEQTL arguments
useModel = modelLINEAR; 

pvOutputThreshold_cis = .05;   # store ALL cis eQTLs less than .1 
pvOutputThreshold_tra = 0;   # don't find trans eQTLs

cisDist = 1e6; # distance for which to test SNPs and gene expression change


expression_mat = as.matrix(data.frame(fread(input = expression_file_name, sep = "\t", header = T),row.names = 1, check.names = F));
G_eigenvec_mat = as.matrix(data.frame(fread(input = G_eigenvec_file_name, sep = "\t", header = T),row.names = 1, check.names = F))[1:N_PCs,];
covariates_extra_mat = as.matrix(data.frame(fread(input = covariate_file_name, sep = "\t", header = T),row.names = 1, check.names = F));

covariates_mat = rbind(G_eigenvec_mat,covariates_extra_mat)


print('Expression and covariate matrices are fine')

genepos = read.table(file = gene_location_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE);
snpspos = read.table(file = snps_location_file_name, header = FALSE, stringsAsFactors = FALSE);

print('gene pos and snppos are fine')

# gene expression data
gene = SlicedData$new();
gene$CreateFromMatrix(expression_mat);

print('gene exp from SlicedData is fine')


# covariates! 
covariates = SlicedData$new();
covariates$CreateFromMatrix(covariates_mat);

print('covariate SlicedData is fine')


# SNPs
## Genotype data with snp position
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

print('SNP sliced data is fine')

# meh = Matrix_eQTL_main(
#   snps = snps,
#   gene = gene,
#   cvrt = covariates,
#   output_file_name  = output_file_name_tra,
#   pvOutputThreshold  = pvOutputThreshold_tra,
#   useModel = useModel,
# #   errorCovariance = errorCovariance,
#   verbose = TRUE,
#   output_file_name.cis = output_file_name_cis,
#   pvOutputThreshold.cis = pvOutputThreshold_cis,
#   snpspos = snpspos,
#   genepos = genepos,
#   cisDist = cisDist,
#   pvalue.hist = 100,
#   min.pv.by.genesnp = FALSE,
#   noFDRsaveMemory = FALSE);



# png(file=glue('./{pop}_cis_eQTLS_{N_PCs}PCs_pval_dist.png'))
# plot(meh)
# dev.off();


meq = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = covariates,
  output_file_name  = output_file_name_tra,
  pvOutputThreshold  = pvOutputThreshold_tra,
  useModel = useModel,
#   errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


png(file=glue('../results/{pop}_cis_eQTLs_{N_PCs}gPCs_{var_method}_{exp_PCs}ePCs_{covariate_type}_qqplot.png'))
plot(meq)
dev.off();
