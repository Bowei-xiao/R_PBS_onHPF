# R_PBS_onHPF
PBS script and related R script used on HPF(High-Peformancing computing Facility). Mainly for GWAS data, designed to parallel on 23 chromosome.
The files in this repo are:

--gwas_association_binary.pbs: Main file to run on HPF,  unzip data and formatting if necessary, passed down the data to 
runAnalysis.sh to start the analysis and collect the results and bind back to proper format for future use.

--runAnalysis.sh: intermediate script that pass down data from gwas_association_binary.pbs to a suitable  Meta_log/gee.R

--Meta_log/gee.R: R script that do the corresponding analysis (logistic/GEE)

