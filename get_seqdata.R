getTCGAmatrix2 <- function(files.loc){
i=files.loc[1]
#print(i)

seqdata <- read.csv2(i, sep = "\t", skip = 1, header = T)[-c(1:4),c(1,4)]
colnames(seqdata)[2] <- str_split_fixed (i, "/",4)[,4]

for (i in files.loc[2: length(files.loc)]) {
seqdata_temp <- read.csv2(i, sep = "\t", skip = 1, header = T)[-c(1:4),c(1,4)]
if (  all (all(seqdata[,1] == seqdata_temp[,1])) ) {
        seqdata[, str_split_fixed (i, "/",4)[,4]] <- seqdata_temp[,2]
 }
 }
return(seqdata)
}

TCGA <- dir( pattern = "TCGA-")
library(stringr)
for ( i in TCGA ) {
#  print(i)
print(paste0(i, "/", i, "__seqdata.rds"))
  if ( file.exists(paste0(i, "/", i, "__seqdata.rds")  )) {
  print("file exists!")
  }
    else
  {

  folders <- dir(paste0(i, "/", str_replace(i, "TCGA-", "") ))
  folders <- folders[!grepl("Normal", folders)]
  loc <- paste0(i, "/", str_replace(i, "TCGA-", ""), "/",  folders, "/")
 # print(loc)
  files.loc <- c()
  for (j in loc)
   # print(dir(j))
  files.loc <-c( files.loc, 
                 paste0(j,  dir(j) )
                 )

#files.loc <- files.loc[1:10]
seqdata <- getTCGAmatrix2(files.loc)
print(dim(seqdata))

saveRDS(seqdata, file= paste0(i, "/", i, "__seqdata.rds"))
}

}
