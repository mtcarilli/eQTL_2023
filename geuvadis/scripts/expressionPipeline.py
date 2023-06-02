# expression normalization pipeline for Geuvadis data

# import
import sys, os
import argparse

# numbers
import numpy as np

# PCA
from sklearn.decomposition import PCA
import scipy.stats as stats

# text processing
import pandas as pd
from collections.abc import Iterable

# process zipped gtf file
import functools
import gzip


# anndata
import anndata


# pipeline specifications:
parser = argparse.ArgumentParser()


parser.add_argument('--name', type=str)
parser.add_argument('--population')
parser.add_argument('--var_method', default = 'log1p', type=str)
parser.add_argument('--num_PCs', default = 10, type=int)

args = parser.parse_args()

TPM = None   # make all counts sum to this? 
var_method = args.var_method
num_PCs = args.num_PCs
population = args.population
name = args.name


# if not isinstance(populations, Iterable):
#     populations = [populations]
    

out_path = f'../expressionData/{name}_{var_method}_{num_PCs}_expression_df_T.tsv'
gene_path = f'../expressionData/{name}_gene_loc.tsv'
samp_path = f'../expressionData/{name}_sample_ids.tsv'
SNP_path = f'../genotypeData/'

# TO DO: ADD FILE THAT REMOVES SEX AND MHC GENES
expression_path = '../expressionData/expression_noSex_noMHC.adata'

# functions
def quantile_normalize(df):
    """Quantile normalizes a pandas DataFrame.
      From PCCA github
      https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe
      Args:
        df: A pandas DataFrame
      Returns:
        The dataframe where all columns have the same distribution.
      """
    rank_mean = df.stack().groupby(
    df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

def inverse_normal_transform(array):
    ''' Transform each row of expression (sample expression) to the inverse normal distribution. 
    USES BLOM METHOD: https://www.statsdirect.com/help/data_preparation/transform_normal_scores.htm 
    '''
    
    # rank each array 
    new_array = np.ones_like(array,dtype = float)
    for i,row in enumerate(array):
        for j,val in enumerate(row):
            val_quantile =  (sum(row<val) + .375)/(len(row) + 1.4)
            inv_val = stats.norm.ppf(val_quantile)
            new_array[i,j] = inv_val
            
    return(new_array)

# useful function for unpacking arrays
def nd(arr):
    return np.asarray(arr).reshape(-1)




# populations_string = ''
# for p in populations:
#     populations_string += p 
# out_path = args.out_path if args.out_path not None else f'./{populations_string}_expression_df_T.tsv'

# load in data
adata = anndata.read(expression_path)

# and metadata
sample_err_pop_overlap = pd.read_csv(SNP_path + 'err_id_sample_pop_map_unique.csv')
err_overlap = sample_err_pop_overlap['ERR ID']
sample_overlap = sample_err_pop_overlap['Sample']
pop_overlap = sample_err_pop_overlap['Population'];


# unique samples -- one ERR per genotyped sample, can change this, for now saves first ERR corresponding to sample ID
samps_unique = []
err_to_samps_unique = []
pop_to_samps_unique = []

for i,err in enumerate(err_overlap):
    samp = sample_overlap[i]
    pop = pop_overlap[i]
    if samp in samps_unique:
        continue
    else:
        samps_unique.append(samp)
        err_to_samps_unique.append(err)
        pop_to_samps_unique.append(pop)
        
# subset     
print(population)
err_subset = []
samps_subset = []
err_subset += [err for i,err in enumerate(err_to_samps_unique) if population in pop_to_samps_unique[i]]
samps_subset += [samp for i,samp in enumerate(samps_unique) if population in pop_to_samps_unique[i]]
adata = adata[err_subset]


# save the unique samples
samps_save = ['0_' + samp for samp in samps_subset]
s = len(samps_subset)
print(f'THE LENGTH OF SAMP SUBSET IS {s}')
samp_save_df = pd.DataFrame({'fam' : [0]*len(samps_subset), 'samp' : samps_subset})
samp_save_df.to_csv(samp_path,sep='\t',header=None,index=None)

# remove all transcripts with abundance < 0.1 TPM
# SHOULD BE DONE EVERY TIME
avg_tpm  = nd((adata.X).sum(axis=0))/len(adata)
adata = adata[:,avg_tpm > 0.1]

# normalize? -- each sample sum to N TPM FLAG, unnecessary i think
if TPM != None:
    adata.X = TPM*adata.X/(adata.X.sum(axis=1))

# variance stabilize

if var_method == 'log1p':
    adata.X = np.log1p(np.asarray(adata.X))

elif var_method == 'quantile':
    df = pd.DataFrame(adata.X.todense())
    df_qnorm = quantile_normalize(df)
    adata.X = df_qnorm
    
elif var_method == 'INT':
    data = np.asarray(adata.X.todense())
    adata.X = inverse_normal_transform(data)
    
elif var_method == 'None':
    adata.X = adata.X

# mean center to mean 0 
adata.X = adata.X.todense() - adata.X.mean(axis=0)
expression = np.asarray(adata.X.todense())

# regress out the top num_pca PCA components

pca_reduce = PCA(n_components = num_PCs)

# fit 
pca_reduce.fit(expression)
x_reduced = pca_reduce.transform(expression)

# QUESTION -- will the returned values be centered at 0? because PCA 0-centers the data (why, i have forgotten)
x_regress = pca_reduce.inverse_transform(x_reduced)

# remove residuals: 
expression_corrected = expression - x_regress

# set observed counts to corrected expression. Now, this matrix should be ready for eQTL analysis!
adata.X = expression_corrected


expression_df_T = pd.DataFrame(expression_corrected.T,
                               index = [g.split('.')[0] for g in adata.var['gene_name']],
                                        columns = samps_save)

expression_df_T.to_csv(out_path,sep='\t')



# GENES to perform eQTL analysis on: only those retained in the matrix

def build_gene_location_dict(gtf_file_name, genes_to_use):
    genes = []
    chromosomes = []
    s1_list = []
    s2_list = []
    
    if gtf_file_name[-3:] == '.gz':
        open_fn = functools.partial(gzip.open, mode='rt')
    else:
        open_fn = open
    with open_fn(gtf_file_name) as gtf_file:
        for line in gtf_file:
            if line[0] == '#':  # skip comments
                continue
            l = line.strip().split()
            feat_type = l[2]
            gene_name = l[9][1:-2]
            if feat_type == 'gene' and gene_name in genes_to_use:
                strand = l[6]
                chrom = str(l[0])
                if strand == '+' and chrom not in ['X', 'Y', 'M']:
                    s1 = l[3]
                    s2 = l[4]
                    genes.append(gene_name)
                    chromosomes.append(chrom)
                    s1_list.append(s1)
                    s2_list.append(s2)
                elif chrom not in ['X', 'Y', 'M']:
                    s1 = l[3]
                    s2 = l[4]
                    genes.append(gene_name)
                    chromosomes.append(chrom)
                    s1_list.append(s1)
                    s2_list.append(s2)
                    
    return pd.DataFrame({'geneid' : genes, 'chr' : chromosomes, 's1' : s1_list, 's2' : s2_list})

keep = [g.split('.')[0] for g in adata.var.gene_id.tolist()]
gene_data_frame = build_gene_location_dict('../kallisto_alignment/Homo_sapiens.GRCh38.104.gtf.gz',keep)


gene_data_frame = gene_data_frame.set_index('geneid')
gene_data_frame.to_csv(gene_path,sep='\t')


