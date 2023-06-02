pops_to_analyze=("Tos" "Fin" "Ut" "Yor" "Gbr")

for POP in "${pops_to_analyze[@]}"
do
    
    COV_TYPE="lab_sex"
    # vary number of exp PCs with 10 genotype PCs 
    exp_PCs=10
    var_method='quantile'
    Rscript run_matrixEQTL.R $POP 10 $var_method 1 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 10 $var_method 5 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 10 $var_method 10 $COV_TYPE
    
    # vary normalization
    Rscript run_matrixEQTL.R $POP 10 "log1p" 10 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 10 "None" 10 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 10 "quantile" 10 $COV_TYPE
    
    # vary number of genotype PCs
    Rscript run_matrixEQTL.R $POP 1 $var_method 10 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 5 $var_method 10 $COV_TYPE
    Rscript run_matrixEQTL.R $POP 10 $var_method 10 $COV_TYPE
    
    # vary covariate type 
    Rscript run_matrixEQTL.R $POP 10 $var_method 10 "lab"
    Rscript run_matrixEQTL.R $POP 10 $var_method 10 "sex"
    
done
