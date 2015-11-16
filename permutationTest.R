#setwd('/Users/bowei/Desktop/SLC9A3/0918geneTest/')
library('doMC')
library('foreach')
registerDoMC(30)

args=(commandArgs(TRUE))
theRepeatCount=as.numeric(args[1])
g1_file = as.character(args[2])
g2_file = as.character(args[3])
phenoAll = read.table('fam_all.fam',header=T,stringsAsFactors = F)
ge1= read.table(g1_file, header=F, sep='\t', stringsAsFactors=F)
ge2= read.table(g2_file, header=F, sep='\t', stringsAsFactors=F)

snpList = as.character(ge1$V2)
dose =  do.call(cbind, lapply(FUN=function(t) c(as.matrix(ge1[t,-c(1:5)]),as.matrix(ge2[t,-c(1:5)])), 1:nrow(ge1))) #Snps per column
# First Check, how many pairs/triplets we have
colnames(dose) = c(paste0('dose',1:nrow(ge1)))
phenoAll = cbind(phenoAll,dose)

# 4458 singlton, 750 pairs, and 26 >=2 (25 tri- and 1 quad-)
print('Finish loading data')

getUnrelated = function(pheno){
set.seed(1020)
p1=NULL;p2=NULL;pg2=NULL;
allFID = unique(pheno$FID)
for (i in 1:length(allFID)){
  p = pheno[grep(paste0(allFID[i],'$'),pheno$FID),]
  if (nrow(p) == 1){
    p1 = rbind(p1,p)
  } else if (nrow(p) ==2){
    p2 = rbind(p2,p)
  } else {
    pg2 = rbind(pg2,p)
  }}

# Method 1: Throw away pairs and triplets
# Pick the individual with MI if possible
pickMIind = function(inx,pair,sibNum,jump=9999999){
  if (sibNum >2 & !(is.na(k))){jump=k}
  MIs = pair$MI[sibNum*inx+c(1:sibNum)+(inx>(jump-1))]
  if (sum(MIs)==0){
    return(sibNum*inx+sample(sibNum,1)+(inx>(jump-1)))
  } else {
    return((sibNum*inx+which(MIs == 1))[1]+(inx>(jump-1)))
  }
}
k = (grep('JH0088',pg2$FID)[1]-1)/3+1
#set.seed(1020)
remainp2 = p2[mapply(FUN=function(t) pickMIind(t,p2,2), 0:(nrow(p2)/2-1)),]
# JH0088 is a quad
if (is.na(k)){
remainp3 = pg2[mapply(FUN=function(t) pickMIind(t,pg2,3), 0:(nrow(pg2)/3-1)),]  
remainList = rbind(p1, remainp2,remainp3)
 }else {
remainp3 = pg2[c(mapply(FUN=function(t) pickMIind(t,pg2,3),0:(k-2)) 
                 ,mapply(FUN=function(t) pickMIind(t,pg2,3),k:((nrow(pg2)-1)/3-1))),] #triplets without JH0088
remainList = rbind(p1, remainp2,remainp3,pg2[(k-1)*3+sample(4,1),])
}
return(remainList)
}
print('finish subsetting Individuals')

unrelatPermu = function(dset,snpList){
# Permute within all possible platform-cohort conbination
  dset$subgroup = paste0(dset$COHORT,'_',dset$PLATFORM)
  dset$newMI = 0
  for (sub in names(table(dset$subgroup))){
    indexSub = which(dset$subgroup == sub)
    dset$newMI[indexSub] = dset$MI[sample(indexSub)] 
  }
  
z2 = NULL;z2Ori=NULL; 
for (i in 1:length(snpList)){
if (length(levels(as.factor(dset$PLATFORM)))==1){
    res = glm(as.formula(paste0('newMI~dose',i,'+as.factor(COHORT)+'
                                ,c(paste0('PC',1:6,collapse = '+'))))
              ,family='binomial', data = dset)
   res2 = glm(as.formula(paste0('MI~dose',i,'+as.factor(COHORT)+'
                                ,c(paste0('PC',1:6,collapse = '+'))))
              ,family='binomial', data = dset)
}else{
   res = glm(as.formula(paste0('newMI~dose',i,'+as.factor(PLATFORM)+as.factor(COHORT)+',c(paste0('PC',1:6,collapse = '+')))),family='binomial', data = dset)
res2 = glm(as.formula(paste0('MI~dose',i,'+as.factor(PLATFORM)+as.factor(COHORT)+',c(paste0('PC',1:6,collapse = '+')))),family='binomial', data = dset)
}
  z2 = c(z2, as.numeric(coef(summary(res))[paste0("dose",i),3])^2)
  z2Ori = c(z2Ori, as.numeric(coef(summary(res2))[paste0("dose",i),3])^2)
}

return(list(z2,z2Ori))
}
if (FALSE){
remainList1 = getUnrelated(pheno = phenoAll[phenoAll$MI !=-9,])
chiListAll = foreach(i=1:10000, .combine=rbind) %dopar% {
   unrelatPermu(remainList1,snpList)[[1]]
}
chiListOri = unrelatPermu(remainList1,snpList)[[2]]
write.table(chiListAll, paste0('EmpChi_T2allApical',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)
write.table(chiListOri, paste0('EmpChi_T2allApicalOri',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)

remainList2 = getUnrelated(pheno = phenoAll[phenoAll$MI !=-9 & phenoAll$PLATFORM == 'GWAS1',])
chiListG1 = foreach(i=1:10000, .combine=rbind) %dopar% {
   unrelatPermu(remainList2,snpList)[[1]]
}
chiListOriG1 = unrelatPermu(remainList2,snpList)[[2]]
write.table(chiListG1, paste0('EmpChi_T2G1',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)
write.table(chiListOriG1, paste0('EmpChi_T2G1Ori',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)
}

remainList3 = getUnrelated(pheno = phenoAll[phenoAll$MI !=-9 & phenoAll$PLATFORM != 'GWAS1' & phenoAll$COHORT != 'FR',])
chiListG2 = foreach(i=1:10000, .combine=rbind) %dopar% {
   unrelatPermu(remainList3,snpList)[[1]]
}
chiListOriG2 = unrelatPermu(remainList3,snpList)[[2]]
write.table(chiListG2, paste0('EmpChi_T2G2NA',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)
write.table(chiListOriG2, paste0('EmpChi_T2G2NAOri',theRepeatCount,'.txt'),col.names=F,row.names=F,quote=F)
