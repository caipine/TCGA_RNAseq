library(dplyr)
library(DESeq2)

TCGA <- dir(pattern = "TCGA-")
  for (i in TCGA) {  
  
  print(paste0(i, "/", i, "__deseq2.RDA"))
  if ( file.exists(paste0(i, "/", i, "__deseq2.RDA")  )) {
  print("file exists!")
  }
    else
  {
 
  
#i = TCGA[1]
seqdata <- readRDS(paste0(i, "/", i, "__seqdata.rds"))
sample <- read.csv2(paste0(i, "/", dir(i, pattern = "gdc_sample_sheet") ), sep = "\t")
cli <- read.csv2(paste0(i, "/", 'clinical.tsv'), sep = "\t")
cli$Case.ID <- cli$case_submitter_id
ff <- merge( sample, cli, by = "Case.ID" )
cli2 <- unique(ff[,c("Case.ID", "days_to_last_follow_up", "vital_status"  ,"days_to_death", "File.Name")])
cli2$Case.ID <-  paste0(cli2$Case.ID,"_" , 1:nrow(cli2))
rownames(cli2) <- cli2$File.Name


rownames(seqdata) <- seqdata[,1]
#seqdata_x <- seqdata_x[,-1]
cli2x <- cli2[colnames(seqdata)[-1],]
colnames(seqdata)[-1] <-  cli2x$Case.ID
#colnames(seqdata_x) ==  cli2x$Case.ID
cli2x <- cli2x[!is.na(colnames(seqdata)[-1]),]
rownames(cli2x) <- cli2x$Case.ID
sampleinfo <- cli2
     seqdata <- seqdata [, -1]    
                  
keep <- rowSums(seqdata) > 5
countdata <- seqdata[keep,]
countdata <- as.matrix(countdata)


sampleinfo <- sampleinfo[!sampleinfo$vital_status %in% "Not Reported",]
dim(countdata)
countdata <- countdata[, colnames( countdata)%in%  sampleinfo$Case.ID]
dim(countdata)
sampleinfo <- sampleinfo[sampleinfo$Case.ID %in% colnames( countdata), ]
countdata <- countdata[, sampleinfo$Case.ID]

sampleinfo$Condition <- factor(sampleinfo$vital_status, levels = c("Alive", "Dead") )
sampleinfo$librarySizes <-  colSums(countdata)

#rlogcounts <- vst(countdata) 
all(colnames(countdata) == sampleinfo$Case.ID)
rownames(sampleinfo) <- sampleinfo$Case.ID

#ggplot(data=sampleinfo, aes(x=File.Name, y=librarySizes, color=Condition )) +
#  geom_bar(stat="identity")+
#theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

design <- as.formula(~ Condition )
#design 
modelMatrix <- model.matrix(design, data = sampleinfo)
#modelMatrix

# create the DESeqDataSet object
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design)
#nrow(ddsObj.raw)
keep <- rowSums(counts(ddsObj.raw)) > 1
ddsObj.raw <- ddsObj.raw[keep,]

#relevel  if nesessery
#ddsObj.raw$Condition
ddsObj.raw$Condition <- relevel(ddsObj.raw$Condition, ref = "Alive")

ddsObj <- DESeq(ddsObj.raw)
#resultsNames(ddsObj)
#rld <- rlog(ddsObj, blind=FALSE)

vst<- vst(ddsObj, blind=FALSE)

save(vst, sampleinfo, ddsObj, file =  paste0(i, "/", i, "__deseq2.RDA") )
  }
  
  }
  
