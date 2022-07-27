library(dplyr)
library(DESeq2)


sample <- read.csv2(paste0( dir( pattern = "gdc_sample_sheet") ), sep = "\t")
cli <- read.csv2(paste0( 'clinical.tsv'), sep = "\t")
cli$Condition <- cli$vital_status
cli[!cli$Condition %in% "Dead",]$Condition <- "Alive"

cli$Case.ID <- cli$case_submitter_id
ff <- merge( sample, cli, by = "Case.ID" )
#cli2$Case.ID <-  paste0(cli2$Case.ID,"_" , 1:nrow(cli2))
cli2 <- unique(ff[,c("Case.ID", "days_to_last_follow_up", "vital_status"  ,"days_to_death", "File.Name", "Project.ID" , "Condition" )])

cli2$Case.ID2 <- cli2$Case.ID 
cli2$Case.ID <- paste0("X_",cli2$Case.ID)
cli2$File.Name <- paste0("X_",cli2$File.Name)
#rownames(cli2) <- cli2$File.Name
print(unique(cli2$Project.ID))

n <-  1
for ( i in unique(sample$Project.ID)) {
print (n)
if (file.exists(  paste0("'", i, "'__deseq2.RDA") ) ) {
print( paste0("'", i, "'__deseq2.RDA  exists") )
 } else {
 #(file.exists( paste0("'", i, "'__seqdata.rds"))      & !file.exists( paste0("'", i, "'__deseq2.RDA")) ) {

#i =  unique(sample$Project.ID)[3]
print(paste0("Working DEseq2 on ", i) )

seqdata <- readRDS( paste0("'", i, "'__seqdata.rds"))

rownames(seqdata) <- seqdata[,1]
seqdata <- seqdata[, -1]
colnames(seqdata) <- paste0("X_",colnames(seqdata))

cli2x <- cli2[ cli2$File.Name %in%  colnames(seqdata),]
cli2x[cli2x$days_to_death %in% "'--",]$days_to_death <- cli2x[cli2x$days_to_death %in% "'--",]$days_to_last_follow_up
dim(cli2x)

if ( nrow(cli2x[!cli2x$days_to_death %in% "'--",]) > 0 ) {
cli2x <- cli2x[!cli2x$days_to_death %in% "'--",]


dim(cli2x)
dim(cli2x[ cli2x$File.Name %in% cli2x$File.Name[duplicated(cli2x$File.Name)],])

seqdata <- seqdata[,cli2x$File.Name]

sampleinfo <- cli2x

rownames(sampleinfo) <- sampleinfo$File.Name
sampleinfo <- sampleinfo[colnames(seqdata),]

keep <- rowSums(seqdata) > 5
countdata <- seqdata[keep,]
countdata <- as.matrix(countdata)


#sampleinfo <- sampleinfo[!sampleinfo$vital_status %in% "Not Reported",]
dim(countdata)
countdata <- countdata[, colnames( countdata)%in%  sampleinfo$File.Name]
rownames(sampleinfo) <- sampleinfo$File.Name
dim(countdata)
sampleinfo <- sampleinfo[colnames( countdata), ]
countdata <- countdata[, sampleinfo$File.Name]
dim(sampleinfo)
all(colnames(countdata) == sampleinfo$File.Name)

if (length(unique(sampleinfo$Condition)) == 1) {
sampleinfo$Condition[1:4] <- c("Alive", "Dead", "Alive", "Dead")
}

sampleinfo$Condition <- factor(sampleinfo$Condition, levels = c("Alive", "Dead") )
sampleinfo$librarySizes <-  colSums(countdata)


#rlogcounts <- vst(countdata) 
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
print("starting DESeq")
ddsObj <- DESeq(ddsObj.raw)
print("Done")

#resultsNames(ddsObj)
#rld <- rlog(ddsObj, blind=FALSE)

vst<- vst(ddsObj, blind=FALSE)

save(vst, sampleinfo, ddsObj, file =  paste0("'", i, "'__deseq2.RDA") )
save(vst, sampleinfo, file =  paste0("'", i, "'__deseq2-2.RDA") )
   
} 
else { print(" Survival data not available") }
}
n =n +1
}

