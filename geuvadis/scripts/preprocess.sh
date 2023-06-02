# process genotype and gene expression data for each population before running MatrixEQTL

pops_to_analyze=('Tos' 'Fin' 'Yor' 'Ut')
exp_path='../expressionData'
SNP_path='../genotypeData'
script_path='./'
num_ePCs=10
num_gPCs=10
var_method='quantile'
COV='lab_sex'

for POP in "${pops_to_analyze[@]}"
do
    # expression processing
    python3 $script_path/expressionPipeline.py --name $POP --population $POP --var_method $var_method --num_PCs num_ePCs

    # genotype processing
    # subset matrix
    $SNP_path/plink --vcf $SNP_path/plink_subset.vcf --keep $exp_path/${POP}_sample_ids.tsv --maf 0.05 --not-chr XY --make-bed -out $SNP_path/${POP}_plink


#     # PCA of genotype
    $SNP_path/plink --bfile $SNP_path/${POP}_plink --chr 1-22 --snps-only --pca --maf 0.05 --make-bed --out $SNP_path/${POP}

    # transpose and add alleles
    $SNP_path/plink --bfile $SNP_path/${POP} --recode A-transpose --out $SNP_path/${POP}_T
    
    # remove all but variants 
    # 430 is just a large enough number to access all individuals
    cut -f2,7-430 $SNP_path/${POP}_T.traw > $SNP_path/${POP}_T_snps.txt

    # get the SNPS location file in order
    cut -f1,4 $SNP_path/${POP}.bim > $SNP_path/${POP}_chr_pos.txt
    cut -f2 $SNP_path/${POP}.bim > $SNP_path/${POP}_snp.txt
    paste $SNP_path/${POP}_snp.txt $SNP_path/${POP}_chr_pos.txt > $SNP_path/${POP}_snps_loc.txt

#     # reformat eigenvecs
#     python3 $SNP_path/reformatEigenVecs.py --name $POP
    
    # save covariates
    python3 $script_path/reformatCovariates.py --name $POP --covariates $COV
done
