gene_sur3 <- function(genelist ) {

 # genelist <-  c("ACLY" , "ACACA", "FASN") # ALCY_TCA_PURINR 
  #load( paste0("'",5, "'__deseq2-2.RDA") )
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

geneidnames <-  read.csv2("../TCGA-LIHC/LIHC/Primary_Tumor/00be7d10-fd2d-4483-a712-72f24846610f.rna_seq.augmented_star_gene_counts.tsv", sep = "\t", skip = 1, header = T)[-c(1:4),c(1,2)]
#geneidnames[geneidnames$gene_name %in% genelist, ]
geneids <- geneidnames[geneidnames$gene_name %in% genelist, ]$gene_id
geneids <- geneids[geneids %in% rownames(vst@assays@data@listData[[1]])]

data_subset_norm <- data.frame(vst@assays@data@listData[[1]][geneids,]  )
rownames(data_subset_norm)
 data_subset_norm2 <-   t(apply(as.matrix(data_subset_norm ), 1, cal_z_score))
 a <- data.frame(colSums(data_subset_norm2))

 
 if (length(genelist) == 1) {
   print("here")
 # data_subset_norm <-  data.frame(matrix(vst@assays@data@listData[[1]][geneids,], nrow =1   ))
#(  data_subset_norm ) <- colnames(vst@assays@data@listData[[1]])
a <- data.frame(vst@assays@data@listData[[1]][geneids,])
}

colnames(a) <- "A3"
a$name<- rownames(a)
#a$name <- str_replace(a$name, "-", ".")
#a$name <- str_replace(a$name, "-", ".")
#a$name <- str_replace(a$name, "-", ".")


sinfo <- sampleinfo[!( sampleinfo$days_to_last_follow_up %in% "'--" &  sampleinfo$days_to_death  %in% "'--"),]
sinfo[sinfo$days_to_death %in% "'--",]$days_to_death <- sinfo[sinfo$days_to_death %in% "'--",]$days_to_last_follow_up

#sinfo$days_to_death <- sinfo$days_to_last_follow_up ###############################################################################################?
sinfo[sinfo$vital_status %in% "Alive",]$days_to_death <- sinfo[sinfo$vital_status %in% "Alive",]$days_to_last_follow_up
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")
rownames(sinfo) <- str_replace(rownames(sinfo), "-", ".")

rownames(sinfo)[1]
a$name[1]
rownames(sinfo)[1] ==a$name[1]
s <- sinfo[rownames(sinfo) %in% a$name,]
s


s$name <- rownames(s)

s$sur <- as.numeric(s$days_to_death)
s$alive <- 1
if (nrow(s[s$vital_status %in% "Alive",]) >0 )
s[s$vital_status %in% "Alive",]$alive <- 0  # dead = 1

#s$sur2 <- as.numeric(s$sur2)
rownames(s) == rownames(a)
#a[rownames(s) ,]


d1 <- merge(a, s, by= "name" )
colnames(d1)

d1 <- d1[order(d1[,"A3"]),]

return(d1)
}
