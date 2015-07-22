# Syntax: fileForMetasoft.R ${theFile} 
# The reads-in variable `file' is the file name as chrxx.blockxx
# This function reads in the 4 separate cohort created before
# And then created the file in metasoft format:
# rs1 b1 se1 b2 se2 b3 se3 b4 se4
# rs2 b1 se1 b2 se2 b3 se3 b4 se4....
# Write the result file as ${theFile}_forMetasoft.txt in the same directory as this R Script
args=(commandArgs(TRUE))
file = as.character(args[1])
print(file)
file_list = list.files(pattern = paste0('^',file,'_+'))
print(file_list)
final_list = read.table(file_list[1],header=F,stringsAsFactors = F)[1:3]
for (i in (2:length(file_list))){
  cohort_info = read.table(file_list[i],header=F,stringsAsFactors = F)
  final_list = cbind(final_list, cohort_info[,2:3])
}
write.table(final_list,paste0(file,'_forMetasoft.txt'),quote = F,col.names = F,row.names = F)
