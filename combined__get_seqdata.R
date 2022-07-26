library(stringr)



getTCGAmatrix3 <- function(files.loc){
i=files.loc[1]
#print(i)

seqdata <- read.csv2(i, sep = "\t", skip = 1, header = T)[-c(1:4),c(1,4)]
colnames(seqdata)[2] <- str_split_fixed (i, "/",3)[,3]
for (i in files.loc[1 : length(files.loc)]) {
#print(i)
seqdata_temp <- read.csv2(i, sep = "\t", skip = 1, header = T)[-c(1:4),c(1,4)]
     if (all (all(seqdata[,1] == seqdata_temp[,1])) ) { 
        seqdata[, str_split_fixed (i, "/",3)[,3]] <- seqdata_temp[,2]
    }
 }
return(seqdata)
}


 loc <- dir("NONTCGA")
 files.loc <- c()
 for (j in loc) {
   print(j)
  files.loc <-c( files.loc, 
                 paste0("NONTCGA/", j,"/",  dir(paste0("NONTCGA/",j)) )
                 )
}
  length(files.loc)
files.loc[1]
files.loc[9193]
files.data <- data.frame(str_split_fixed(files.loc, "/", 3))
files.data$loc <- files.loc
 
 
 
sample <- read.csv2(paste0( dir( pattern = "gdc_sample_sheet") ), sep = "\t")

for ( i in unique(sample$Project.ID)) {
files.loc.n <- files.data[files.data$X3 %in% sample[sample$Project.ID %in% i,  ]$File.Name,]$loc

seqdata <- getTCGAmatrix3(files.loc.n)
print(dim(seqdata))
saveRDS(seqdata, file= paste0("NON_TCGA", "/'", i, "'__seqdata.rds"))
}
