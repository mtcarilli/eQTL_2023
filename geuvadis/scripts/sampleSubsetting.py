# sampleSubsetting
# must have all necessary metadata 


import pandas as pd
import numpy as np

SNP_path = '../genotypeData/'
exp_path = '../expressionData/'

# load in metadata, sample IDs, and ERR IDs (with aligned mRNA)
samples_full = pd.read_table(SNP_path + '/igsr_Geuvadis.tsv') # metadata, links sample IDs to ERR IDs
samples_lab_metadata = pd.read_table(SNP_path + 'sample_metadata.tsv',header=None)
samples_sex_metadata = pd.read_table(SNP_path + 'sample_sex.tsv')
err_aligned = pd.read_table(exp_path + 'sample_table.txt',header=None)[0].tolist() # ERR IDs

samples_genotype = pd.read_table(data_path + '/samples.txt',header=None)[0].to_list() # samples with available genotypes
samples_genotype = [s[2:] for s in samples_genotype]


err_overlap = []
url_overlap = []
sample_id_overlap = []
pop_overlap = []


# match ERR to sample IDs
for i,err in enumerate(err_aligned):
    for j,url in enumerate(samples_full['url'].tolist()):
        
        if err in err_overlap:
            continue
        elif ((err in url) and (samples_full['Sample'][j] in samples_genotype) and (samples_full['Analysis group'][j] == 'mRNA')):
                sample = samples_full['Sample'][j]
                sample_id_overlap.append(sample)
                url_overlap.append(url)
                err_overlap.append(err)
                pop_overlap.append(samples_full['Population'][j])   

# data frame of mapping
df = pd.DataFrame({'ERR ID' : err_overlap, 'Sample' : sample_id_overlap, 'Population' : pop_overlap})
# save the mapping dataframe
df.to_csv(SNP_path + 'err_id_sample_pop_map.csv',header=None,index=None)


# pick one ERR ID per genotyped sample ID
# can change this, for now saves first ERR corresponding to sample ID
samps_unique = []
err_to_samps_unique = []
pop_unique = []
lab_unique = []
sex_unique = []

for i,err in enumerate(err_overlap):
    samp = sample_id_overlap[i]
    pop = pop_overlap[i]
    if samp in samps_unique:
        continue
    elif samp == 'NA07000':
        samps_unique.append(samp)
        err_to_samps_unique.append(err)
        pop_unique.append(pop)
        lab_unique.append(1.0)
        sex = 'female'
        sex_unique.append(sex)
    else:
        samps_unique.append(samp)
        err_to_samps_unique.append(err)
        pop_unique.append(pop)
        lab_unique.append(lab)
        sex = samples_sex_metadata[samples_sex_metadata['Sample name'] == samp]['Sex'].values[0]
        sex_unique.append(sex)

# make and save df
unique_df = pd.DataFrame({'ERR ID' : err_to_samps_unique, 'Sample' : samps_unique, 'Population' : pop_unique})
unique_df.to_csv(SNP_path + '/err_id_sample_pop_map_unique.csv')

# sort 
sorted_samples_full_meta = samples_full.sort_values('Sample')
ind = [(samp in samps_unique) for samp in samples_full['Sample']]
sorted_samples_subset_meta = samples_full[ind]

# save
sorted_samples_subset_meta.to_csv(SNP_path + 'sorted_sample_subset_metadata.csv')
