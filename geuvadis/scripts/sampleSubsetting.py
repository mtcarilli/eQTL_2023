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





# also store covariates 
samps_to_store = ['0_' + s for s in samps_unique]


# convert sexes to 0 and 1 
sexes = [1 if s == 'female' else 0 for s in sex_unique]

# one hot encode lab
labs = np.zeros((7,len(lab_unique)))
for i in range(1,8):
    labs[i-1,:] = [1 if l == i else 0 for l in lab_unique]
    
# make and save df
covariate_df = pd.DataFrame({'Lab 1' : labs[0,:], 
                             'Lab 2' : labs[1,:], 
                             'Lab 3' : labs[2,:], 
                             'Lab 4' : labs[3,:], 
                             'Lab 5' : labs[4,:], 
                             'Lab 6' : labs[5,:], 
                             'Lab 7' : labs[6,:], 
                             'Sex' : sexes})

covariate_df_T = covariate_df.T
covariate_df_T.columns = samps_to_store

covariate_df_T.to_csv('./geuvadis.covariates.T',header=True,sep = '\t',index=True)


# separate out covariates by population
populations = ['Fin', 'Tos', 'Yor', 'Ut', 'Bri']

for pop in populations:
    
    index = [True if pop in population else False for population in pop_unique]
    
    pop_labs_unique = np.unique(np.array(lab_unique)[index])
    
    if 7.0 in pop_labs_unique:
        pop_labs_unique = np.setdiff1d(pop_labs_unique, [7.0])
    

    
    cov_dict_ = {f'Lab {i}' : labs[int(i-1),index] for i in pop_labs_unique}
    cov_dict_['Sex'] = np.array(sexes)[index]
    
    
    # make and save df
    covariate_df_ = pd.DataFrame(cov_dict_)

    covariate_df_T_ = covariate_df_.T
    covariate_df_T_.columns = np.array(samps_to_store)[index]
    
    covariate_df_T_.to_csv(f'./{pop}.lab_sex_covariates.T',header=True,sep = '\t',index=True)
    
    # make and save LAB df
    lab_dict_ = {f'Lab {i}' : labs[int(i-1),index] for i in pop_labs_unique}
    lab_df_ = pd.DataFrame(lab_dict_)

    lab_df_T_ = lab_df_.T
    lab_df_T_.columns = np.array(samps_to_store)[index]
    
    lab_df_T_.to_csv(f'./{pop}.lab_covariates.T',header=True,sep = '\t',index=True)
    
    # make and save SEX df
    sex_df_ = pd.DataFrame({'Sex' : np.array(sexes)[index]})

    sex_df_T_ = sex_df_.T
    sex_df_T_.columns = np.array(samps_to_store)[index]
    
    sex_df_T_.to_csv(f'./{pop}.sex_covariates.T',header=True,sep = '\t',index=True)
