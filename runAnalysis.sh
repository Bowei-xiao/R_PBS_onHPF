# !/bin/bash
# Syntax: runAnalysis.sh data1 data2 phenotype ChrNumber
# This is a dependent file for gwas_association_binary.pbs

f1=$1   # data from GW1
f2=$2   # data from GW2
# Tor siblings were added to JHU cohort so that 
# 1. we could run logistic using Toronto samples, reducing calculation time
# 2. More importantly, Toronto sample may not converge since there are so few clusters
# To sum up, now only run JHU using GEE but run the other 3 using logistics
# Actual code was in R


Rscript Meta_log.R ${f1} ${f2} 'fam_fr.fam'
Rscript Meta_log.R ${f1} ${f2} 'fam_tor.fam'
Rscript Meta_log.R ${f1} ${f2} 'fam_unc.fam'
Rscript Meta_gee.R ${f1} ${f2} 'fam_jhu.fam'
echo 'Analysis is Done' 
