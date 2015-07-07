#....
#......
#........
# ======== La funcion principal ===========
logResult = function(dosage,pheno,m, snpname){
  # remove MAF that are < 0.02 
  if (sum(dosage)/(2*m) < 0.02){
    return(c(SNP=snpname, Beta=NA, SE=NA, P_value=NA))
  }else{
    asso= cbind(pheno, dosage)
    # Removing the individuals with missing phenotype status
    asso$MI[which(asso$MI == '-9')] = NA 
    asso = na.omit(asso)
    
    # logistic regression
    # covariates including platform and PC correcting for ethnicity
    res = glm(as.formula(paste0('MI~dosage+as.factor(PLATFORM)+'
                                   ,c(paste0('PC',1:6,collapse = '+'))))
                          ,family='binomial', data = asso)
    
    if (is.na(res$coefficients['dosage'])){
      return(c(SNP=snpname, estimates = NA, SE = NA, P_value = NA))
    }else{
    return(c(SNP=snpname, Beta = as.character(coef(summary(res))["dosage",1]),
             SE = as.character(coef(summary(res))["dosage",2]),
             P_value=as.character(coef(summary(res))["dosage",4])))
  }}
}
# ============ El Fin =================================

library('doMC')
library('foreach')
library('data.table')
registerDoMC(30) # This better matched with PBS file setup  

args=(commandArgs(TRUE))
genoFile1 = as.character(args[1]) # file name for gw1, should be like chr1.block1.dosage
genoFile2 = as.character(args[2]) # file name for gw2, should be like chr1.block1.dosage
phenoFile = as.character(args[3]) # file name for phenotype file

pheno = read.table(phenoFile, header=T, stringsAsFactors = F)
geno1 = fread(genoFile1,data.table=F)
geno2 = fread(genoFile2,data.table=F)

m1 = dim(geno1)[2]-5 # N of indvs
m2 = dim(geno2)[2]-5
n = dim(geno1)[1] # N of SNPS in this block

if (FALSE){ # This is the alternative code using bigmemory package
library('bigmemory')
geno1 = read.big.matrix(genoFile1, type='double', header=F,sep='\t')
geno2 = read.big.matrix(genoFile2, type='double', header=F,sep='\t')
# Since gwas1 and gwas2 have the exact SNPs in the same order
# only need to record once
col_class= c('NULL','character',rep('NULL',dim(geno1)[2]-2))
# The SNP names
name = read.table(genoFile1, header=F, colClasses=col_class,sep='\t')
}
# SNP names
name = geno1$V2
print('read-in Sucessfully')

# Run logistic regression for association analysis 
rstList = foreach(i=1:n, .combine=rbind) %dopar% {
logResult(dosage=as.numeric(c(geno1[i,-c(1:5)],geno2[i,-c(1:5)])),m=m1+m2,pheno=pheno, snpname=name[i])
}
print('p-val success')

# Save output
out_name = paste0(strsplit(path1,split='.dosage')[[1]][1],'_',strsplit(path3,split='[[:punct:]]')[[1]][2],'_out.txt')
# filename would be like `chr1.block1_fr.out.txt'

out_dir = paste0('./results/',strsplit(out_dir,split='/')[[1]][4])
print('******************The file is at')
print(out_dir)

write.table(rstList,out_dir,row.names = F,col.names=F,quote = F)
# Should have saved column names, but in some cases that a big block was separated
# and analysis was done separately, having columns would results in the column names
# showed up in the middle of the file in the re-binded file.

# COULD SOLVE using bash command, will be implemented later
print('file written')


