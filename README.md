# R_PBS_onHPF
PBS script and related R script used on HPF(High-Peformancing computing Facility). Mainly for GWAS data, designed to parallel on 23 chromosome.
The files in this repo are:

* gwas_association_binary.pbs: FIRST STEP to run on HPF,  unzip data and formatting if necessary, passed down the data to 
runAnalysis.sh to start the analysis and collect the results and bind back to proper format for future use. The output file name format would be like `chr1.block1_fr_out.txt'
   - runAnalysis.sh: intermediate script that pass down data from gwas_association_binary.pbs to a suitable  Meta_log/gee.R
   - Meta_log/gee.R: R script that do the corresponding analysis (logistic/GEE). This runs the usual association phenotype~genotype using an additive model. There are conditional analysis code avaiable with minor changes of this main code

* runMetaRE2.pbs: SECOND STEP to run on HPF, take in the output file from the FIRST STEP, and using MetaSoft to run the Meta-Analysis on 4 cohorts with Han & Eskin random effect. This would be the final result in the format as `chr1.block1_metaOut.txt'
   - runFileForMetasoft.R: Dependence file for the SECOND STEP, it takes in the output files from 4 cohorts and re-format to the      required format Metasoft needs
   
* runfixFormatMeta.pbs: OPTIONAL STEP to run, takes in the final result and makes it more R-friendly format so that it would be easier to do later work like Manhattan Plot, Forest Plot etc. Note that this file only keeps the Han & Eskin Random Effect result, so need to refer to the previous one for the complete results. The output file name format would be like `formatFixed_chr11.block1_metaOut.txt'

