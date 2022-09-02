
```{r}
read_unmatch_maf <- function(sample, loc) {

b1 <- read.table(paste0(loc,"tmp/wanted_", sample,"_chr1.unmatched_filtered.vcf"), head=T)
b  <- read.table(paste0(loc,"tmp/wanted_", sample, "_24chr.unmatched_filtered.vcf"), head=F )
b[,1] <- data.frame(str_split_fixed(b[,1],":",2))[,"X2"]
colnames(b) <- colnames(b1)
colnames(b)[10] <- "tumor"

ann <- read.csv(paste0(loc, sample, "/" , sample, ".annovar.hg38_multianno.csv"))


colnames(b)[c(1,2,4,5)] <- c("Chr", "Start", "Ref", "Alt")
#colnames(ann)[1:5] <- c("CHROM", "POS", "End", "REF", "ALT")
#ann1 <-   ann[ (ann$Func.refGene %in% "exonic"  & ann$SIFT_pred %in% "D" & ann$MutationTaster_pred %in% "D" ) | (ann$Func.refGene %in% "exonic"  & !ann$cosmic70 %in% ".")
#a <- merge(b, ann, by= c("CHROM", "POS", "REF", "ALT"))


a <- merge(b, ann, by= c("Chr", "Start", "Ref", "Alt"))
c <- a[grepl("exonic", a$Func.refGene),]
c <- c[!grepl("ncRNA", c$Func.refGene),]
c <- c[!c$ExonicFunc.refGene %in% c( "synonymous SNV","unknown"  ),]

write.table( c, file = paste0(loc, sample, "/" , sample, ".annovar.for_maftool.txt"), sep="\t",  row.names = F)

d <- annovarToMaf(annovar = paste0(loc, sample, "/" , sample, ".annovar.for_maftool.txt"), Center = 'CSI-NUS', refBuild = 'hg38', table = 'refGene')
return(d)
}
```
```{r}
PDXplusMCL.list <- read.table("PDXplusMCL.list")$V1

loc <- "../called_TUMORONLY_MCL_nops_ann/"
sample <- i <- PDXplusMCL.list[1]
```
read_unmatch_maf
```{r}
PDXplusMCL.list <- read.table("PDXplusMCL.list")$V1

  i <- PDXplusMCL.list[1]
  d <- read_unmatch_maf(i,"../called_TUMORONLY_MCL_nops_ann/")
  d_all <- d

for ( i in PDXplusMCL.list[2: length(PDXplusMCL.list)]) {   # 
  print(i)
   d <- read_unmatch_maf(i,"../called_TUMORONLY_MCL_nops_ann/")
  d_all <-rbind(d_all, d) 
}

  
  unique(d_all$Tumor_Sample_Barcode)
  saveRDS(d_all, file = "2022.0830.maf.mclpdx.rds")
```


```{r}
list <- read.table("samples.STM.list")$V1
list <- list[!grepl("V", list)]
list <- list[!grepl("GCART", list)]
list <- list[!grepl("CD19", list)]
list <- list[!grepl("MPN", list)]
list <- list[!grepl("678180N", list)]



 i <- list[1]
  d <- read_unmatch_maf(i,"../called_TUMORONLYstm_nops_ann/")
  d_all <- d


for ( i in list[2: length(list) ]) {
   d <- read_unmatch_maf(i,"../called_TUMORONLYstm_nops_ann/")
  d_all <-rbind(d_all, d) 
}

 unique(d_all$ExonicFunc.refGene)
  unique(d_all$Sample)
 saveRDS(d_all, file= "2022.0830.maf.STM.rds")
```




```{r}
rm(list = ls())
```


```{r}
library(stringr)
d_all_malpdx <- readRDS("2022.0830.maf.mclpdx.rds")
d_all_stm <- readRDS("2022.0830.maf.STM.rds")
d_all <- rbind(d_all_malpdx, d_all_stm)
rm(d_all_malpdx)
rm(d_all_stm)
t1 <-  str_split_fixed(d_all[,"tumor"],":",5)
saveRDS(t1,  file= "2022.0830.maf.PDXMCLSTM_VAF.rds")
```


```{r}
d_all_malpdx <- readRDS("2022.0830.maf.mclpdx.rds")
d_all_stm <- readRDS("2022.0830.maf.STM.rds")
d_all <- rbind(d_all_malpdx, d_all_stm)
rm(d_all_malpdx)
rm(d_all_stm)
t1 <-  readRDS("2022.0830.maf.PDXMCLSTM_VAF.rds")
d_all$VAF <-  t1[,3]
d_all$VAF <-  as.data.frame(sapply(d_all$VAF, as.numeric))[,1]
d_all$DEPTH <-  t1[,4]
d_all$DEPTH <-  as.data.frame(sapply(d_all$DEPTH, as.numeric))[,1]
unique(d_all$ExonicFunc.refGene)

saveRDS(d_all,  file= "2022.0831.maf.PDXMCLSTM.rds")
```
```{r}
rm(list = ls())
```

# from here


```{r}
d_all <- readRDS( "2022.0831.maf.PDXMCLSTM.rds")
load("../sc145_DNA/2022.0824.data251.rdata") # ann5, data251, tmp4

tmp4[,c("WES", "WES4", "RNA")]
rm(data251, ann5)
for ( i  in unique(d_all$Tumor_Sample_Barcode)[unique(d_all$Tumor_Sample_Barcode) %in% tmp4$WES ]) {
  d_all[d_all$Tumor_Sample_Barcode %in% i,]$Tumor_Sample_Barcode <- tmp4[tmp4$WES %in% i, ]$WES4
}
unique(d_all$Tumor_Sample_Barcode)
```


```{r}
d_all[d_all$ExonicFunc.refGene %in% c("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift insertion", "nonframeshift substitution"),]$Variant_Classification <- "MISSENSE" 
d_all[(d_all$SIFT_pred %in% "D" | d_all$SIFT4G_pred %in% "D")  & d_all$ExonicFunc.refGene %in% c("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift insertion", "nonframeshift substitution"),]$Variant_Classification <- "DelMISSENSE"
d_all[d_all$ExonicFunc.refGene %in% c("stoploss", "startloss", "stopgain","frameshift insertion"),]$Variant_Classification <-  "TRUNC"

unique(d_all$Variant_Classification)

d_all <- d_all[d_all$Variant_Classification %in% c("DelMISSENSE", "TRUNC"),]
unique(d_all$Variant_Classification)
```

`
```{r}
ninfo902 <- read_xlsx("\\\\d1prprsh3ccifs/scratch/lym_myl_rsch/qcai1/MCLdata/PDX_WES_TACC2/RMD/info902.merged.edited2.xlsx")
ninfo902$WES
ninfo902
ninfo902[, c("RNA", "WES",
             "BTKi.status", 
             "STM.status",
             "Date.collected", 
             "BTKi.start.date", 
             "BTKi.off.date2",  
             "DateDeathorLastVisit", 
             "sample.on.IBN.time",
             "PT.on.IBN.time",
             "IBN.time.after.collect",
             "alive") ]

ninfo902.wes <- ninfo902[ninfo902$WES %in% unique(d_all$Tumor_Sample_Barcode),]
ninfo902.wes$overall.sur <- ninfo902.wes $DateDeathorLastVisit  - ninfo902.wes$Date.collected

cli <- ninfo902.wes[, c("RNA",
             "WES",
             "sample.on.IBN.time",
             "PT.on.IBN.time",
             "IBN.time.after.collect",
             "overall.sur", 
             "alive",
             "BTKi.status" ), ]
cli$Tumor_Sample_Barcode  <- cli$WES
```

```{r}
mat1 <- d_all[ !grepl("GMCL", d_all$Tumor_Sample_Barcode) &
                !grepl("VJC", d_all$Tumor_Sample_Barcode) &
                !grepl("CD19", d_all$Tumor_Sample_Barcode) &
                 !grepl("0N", d_all$Tumor_Sample_Barcode) &
                 !grepl("MPN", d_all$Tumor_Sample_Barcode) &
                 !grepl("1135027V1", d_all$Tumor_Sample_Barcode) &
                !grepl("1138372V2", d_all$Tumor_Sample_Barcode) &
                !grepl("1135647V7", d_all$Tumor_Sample_Barcode) &
                !grepl("1144000V25", d_all$Tumor_Sample_Barcode) &
                !grepl("GCART", d_all$Tumor_Sample_Barcode) 
               & !grepl("PDX50", d_all$Tumor_Sample_Barcode) 
              & !grepl("PDX51", d_all$Tumor_Sample_Barcode) 
              & !grepl("PDX52", d_all$Tumor_Sample_Barcode) 
              & !grepl("PDX53", d_all$Tumor_Sample_Barcode) 
              & !grepl("PDX68", d_all$Tumor_Sample_Barcode) 
               & !grepl("PDX71", d_all$Tumor_Sample_Barcode) 
               & !grepl("PDX11", d_all$Tumor_Sample_Barcode)
              & !grepl("PDX12", d_all$Tumor_Sample_Barcode)
              & !grepl("PDX1883", d_all$Tumor_Sample_Barcode)
              & !grepl("PDX1884", d_all$Tumor_Sample_Barcode)
              & !grepl("MCL57", d_all$Tumor_Sample_Barcode)
              & !grepl("MCL58", d_all$Tumor_Sample_Barcode) #56, 57,58,59
              & !grepl("MCL78", d_all$Tumor_Sample_Barcode) #40,78
              ,]


dim(mat1[mat1$VAF > 0.2 & mat1$DEPTH >50,])
mat2 <- mat1[mat1$VAF > 0.2& mat1$DEPTH > 50,]
dim(mat2)
```


```{r}
unique(mat2[ grepl("__A0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)
unique(mat2[ grepl("__B0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)
unique(mat2[ grepl("__C0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)
unique(mat2[ grepl("__D0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)
unique(mat2[ grepl("__E0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)
unique(mat2[ grepl("__F0", mat2$Tumor_Sample_Barcode),]$Tumor_Sample_Barcode)

```

```{r,  fig.width=20, fig.height=25}
mat2_PDX_A2 <- mat2[ grepl("__A", mat2$Tumor_Sample_Barcode),]
mat2_PDX_A <- mat2[ grepl("__A", mat2$Tumor_Sample_Barcode) & !grepl("__A0", mat2$Tumor_Sample_Barcode),]
mat2_PDX_B <- mat2[ grepl("__B", mat2$Tumor_Sample_Barcode) & !grepl("__B0", mat2$Tumor_Sample_Barcode),]
mat2_PDX_C <- mat2[ grepl("__C", mat2$Tumor_Sample_Barcode) & !grepl("__C0", mat2$Tumor_Sample_Barcode),]
mat2_PDX_D <- mat2[ grepl("__D", mat2$Tumor_Sample_Barcode) & !grepl("__D0", mat2$Tumor_Sample_Barcode),]
mat2_PDX_E <- mat2[ grepl("__E", mat2$Tumor_Sample_Barcode) & !grepl("__E0", mat2$Tumor_Sample_Barcode),]
mat2_PDX_F <- mat2[ grepl("__F", mat2$Tumor_Sample_Barcode) & !grepl("__F0", mat2$Tumor_Sample_Barcode),]

PDX_names = unique(c(mat2_PDX_A$Tumor_Sample_Barcode,
               mat2_PDX_B$Tumor_Sample_Barcode,
               mat2_PDX_C$Tumor_Sample_Barcode,
               mat2_PDX_D$Tumor_Sample_Barcode,
               mat2_PDX_E$Tumor_Sample_Barcode,
               mat2_PDX_F$Tumor_Sample_Barcode
           ))
 
mat2_RS <- mat2[ !mat2$Tumor_Sample_Barcode %in%  PDX_names]
mat2_R <-  mat2_RS[ mat2_RS$Tumor_Sample_Barcode %in% (cli[cli$BTKi.status %in% "R",]$WES),]
mat2_S <-  mat2_RS[ mat2_RS$Tumor_Sample_Barcode %in% (cli[cli$BTKi.status %in% "S",]$WES),]
unique(mat2_R$Tumor_Sample_Barcode)
unique(mat2_S$Tumor_Sample_Barcode)



cli2 <- data.frame (Tumor_Sample_Barcode = mat2$Tumor_Sample_Barcode)
cli2$Group <- NA
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_A$Tumor_Sample_Barcode),]$Group <- "PDX_A"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_B$Tumor_Sample_Barcode),]$Group <- "PDX_B"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_C$Tumor_Sample_Barcode),]$Group <- "PDX_C"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_D$Tumor_Sample_Barcode),]$Group <- "PDX_D"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_E$Tumor_Sample_Barcode),]$Group <- "PDX_E"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_PDX_F$Tumor_Sample_Barcode),]$Group <- "PDX_F"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_R$Tumor_Sample_Barcode),]$Group <- "R"
cli2[cli2$Tumor_Sample_Barcode %in% unique(mat2_S$Tumor_Sample_Barcode),]$Group <- "S"
cli2[ cli2$Tumor_Sample_Barcode %in%  unique(cli2[is.na(cli2$Group),]$Tumor_Sample_Barcode),]$Group <- "NA"
cli2$BTKi.status <- NA
cli2[cli2$Tumor_Sample_Barcode %in% (cli[cli$BTKi.status %in% "R", ]$WES),]$BTKi.status <- "R"
cli2[cli2$Tumor_Sample_Barcode %in% (cli[cli$BTKi.status %in% "S", ]$WES),]$BTKi.status <- "S"

cli2$Group2 <- cli2$Group
cli2[grepl("PDX_",cli2$Group ) & !is.na(cli2$BTKi.status),]$Group2 <- cli2[grepl("PDX_",cli2$Group ) & !is.na(cli2$BTKi.status),]$BTKi.status
cli2 <- unique(cli2)
#write.table( cli2, file = "2022.0830.cli2.txt", sep="\t",  row.names = F)
```


```{r,  fig.width=20, fig.height=25}
mat2_R2 <- mat2_R[! mat2_R$Tumor_Sample_Barcode %in% c("PT1106365__A0", "1106365B4__A0",
                                                         "MCL63__B0",  "MCL64__B0" ,
                                                         "MCL67__C0",
                                                         "MCL92__D0" ,
                                                         "MCL93__E0"
                                                         ), ]
mat2_S2 <- mat2_S[! mat2_S$Tumor_Sample_Barcode %in% c("PT1106365__A0", "1106365B4__A0",
                                                         "MCL63__B0",  "MCL64__B0" ,
                                                         "MCL67__C0",
                                                         "MCL92__D0" ,
                                                         "MCL93__E0"
                                                         ), ]
unique(mat2_R$Tumor_Sample_Barcode)
unique(mat2_S$Tumor_Sample_Barcode)


unique(mat2_R2$Tumor_Sample_Barcode)
unique(mat2_S2$Tumor_Sample_Barcode)
mat2_RS2 <- rbind(mat2_R2, mat2_S2)


cli2A <- cli2
cli2A$Group3 <- "others"
cli2A[cli2A$Tumor_Sample_Barcode %in% unique(mat2_PDX_A$Tumor_Sample_Barcode)[grepl("AS", unique(mat2_PDX_A$Tumor_Sample_Barcode) )] ,]$Group3 <- "A_S"
cli2A[cli2A$Tumor_Sample_Barcode %in% unique(mat2_PDX_A$Tumor_Sample_Barcode)[grepl("AV", unique(mat2_PDX_A$Tumor_Sample_Barcode) )] ,]$Group3 <- "A_V"
cli2A[cli2A$Tumor_Sample_Barcode %in% unique(mat2_PDX_A$Tumor_Sample_Barcode)[grepl("G1", unique(mat2_PDX_A$Tumor_Sample_Barcode) )] ,]$Group3 <- "A_G1"
cli2A[grepl("__A0", cli2A$Tumor_Sample_Barcode)  ,]$Group3 <- "A_PT"

```

```{r}
cnv2 <- readRDS("../RMD/2022.0705.cnv.rds")  # from 2022.0705.short.VCF.CNV.clustering.RMD

colnames(cnv2)[1:2] <- c("Sample_name", "Gene")
cnv3 <- cnv2[abs(cnv2$SegmentMean) >= 1,]
cnv3$CN <- NA
cnv3[cnv3$SegmentMean <= -1, ]$CN <- "Del"
cnv3[cnv3$SegmentMean <= -2, ]$CN <- "DeepDel"

cnv3[cnv3$SegmentMean >= 1, ]$CN <- "ShallowAmp"
cnv3[cnv3$SegmentMean >= 2, ]$CN <- "Amp"
cnv3 <- cnv3[,c(2,1,5)]

for ( i  in unique(cnv3$Sample_name)[unique(cnv3$Sample_name) %in% tmp4$WES ]) {
  cnv3[cnv3$Sample_name %in% i,]$Sample_name <- tmp4[tmp4$WES %in% i, ]$WES4
}
unique(cnv3$Sample_name)
rm(cnv2)
```

```{r, fig.width=20, fig.height=10}
cli3 <- cli
cli3$status <- 1
cli3[cli3$alive %in% "alive" ,]$status  <- 0
data_tmp <- mat2_RS2
maf.RS2_sur <- read.maf(maf =data_tmp ,                  ,
                     clinicalData =  cli3,
                     vc_nonSyn = c( "DelMISSENSE", "MISSENSE" ,   "TRUNC"),
                     cnTable = cnv3[cnv3$Sample_name %in% unique(data_tmp$Tumor_Sample_Barcode) ,],
                     verbose = FALSE)
#fab.ce.RS2 = clinicalEnrichment(maf = maf.RS2, clinicalFeature = 'BTKi.status')
```


```{r, fig.width=5, fig.height=5}
mafSurvival(maf = maf.RS2_sur, genes = 'TP53', time = 'overall.sur', Status = "status",isTCGA = F)
```
```{r}
prog_geneset = survGroup(maf = maf.RS2_sur, top = 20, geneSetSize = 3, time = 'overall.sur', Status = "status", verbose = FALSE)
prog_geneset[P_value < 0.05]$Gene_combination[!grepl("MUC", prog_geneset[P_value < 0.05]$Gene_combination)]
```
```{r}
mafSurvGroup(maf =maf.RS2_sur, geneSet = c("LENG9", "TARM1", "HLA-DRB5"),time = 'overall.sur', Status = "status",)
```

```{r, fig.width=20, fig.height=10}
data_tmp <- mat2_RS2
maf.RS2 <- read.maf(maf =data_tmp ,                  ,
                     clinicalData =  cli2A,
                     vc_nonSyn = c( "DelMISSENSE", "MISSENSE" ,   "TRUNC"),
                     cnTable = cnv3[cnv3$Sample_name %in% unique(data_tmp$Tumor_Sample_Barcode) ,],
                     verbose = FALSE)
fab.ce.RS2 = clinicalEnrichment(maf = maf.RS2, clinicalFeature = 'BTKi.status')
```

```{r}
fab.ce.RS2$groupwise_comparision[Hugo_Symbol %in% c("ATM", "TP53"),]
```


```{r, fig.width=20, fig.height=10}
t2 <- fab.ce.RS2$groupwise_comparision# [n_mutated_group2 %in% c("0 of 31", "0 of 33"),]
t2R <- t2[Group1 %in% "R",]
t2R_1 <- str_split_fixed(t2R$n_mutated_group1, " ", 3)
t2R_2 <- str_split_fixed(t2R$n_mutated_group2, " ", 3)

t2R$R <- as.data.frame(sapply(t2R_1[,1], as.numeric))[,1]/ as.data.frame(sapply(t2R_1[,3], as.numeric))[,1]
t2R$S <- as.data.frame(sapply(t2R_2[,1], as.numeric))[,1]/ as.data.frame(sapply(t2R_2[,3], as.numeric))[,1]


fab.ce.RS2$groupwise_comparision[p_value < 0.05]
unique(fab.ce.RS2$groupwise_comparision[p_value < 0.2]$Hugo_Symbol)
 plotEnrichmentResults(enrich_res = fab.ce.RS2, pVal = 0.2, geneFontSize = 0.5, annoFontSize = 0.6)

#write.table(t2R[p_value < 02 & R > S]$Hugo_Symbol,              file="2022.0902.RS.txt", quote = F, row.names = F)
```


```{r}
genelist <- t2R[p_value < 0.5 & R > S]$Hugo_Symbol
pathways.all <- readRDS("C:/Users/qcai1/OneDrive - Inside MD Anderson/007analysis/RNAseq/007B1/results/925.pathways.all.rds")

#names(pathways.all)
a <- c()
for (i in names(pathways.all) ){
 a <- c(a, length(pathways.all[[i]][pathways.all[[i]] %in%  genelist]))
}

length(a)
a1 <- data.frame(S = 1: length(a), N = a)

 gl<- names(pathways.all)[a1[a1$N > 5,]$S]
 gl[!grepl("TARGET", gl) & grepl("REAC", gl)]
  gl[grepl("KERA", gl)]
 
 for (i in gl[!grepl("TARGET", gl) & grepl("REAC", gl)]){
   print(i)
      print(intersect(pathways.all[[i]],  genelist))
 }
```



```{r}
alter_fun <- list(background = function(x, y, w, h){
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = "#F0F0F0", col = NA))
    },
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                  gp = gpar(fill = "red", col = NA))
    },
    DEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                  gp = gpar(fill = "#021AFF", col = NA))
    },
    LOH = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"),
                  gp = gpar(col = "purple"))
    },
    

    GAIN = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                  gp = gpar(fill = "#FB9FB5", col = NA))
    },
    TRUNC = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33,
                  gp = gpar(fill = "black", col = NA))
    },
    MISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33,
                  gp = gpar(fill = "#34A02C", col = NA))
    },
    DelMISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "grey", col = NA))     # "#89419D"
        grid.points(x, y, pch = 16, size= unit(0.3, "char") )
    }
   )


col <- c("DEL" = "#021AFF", 
         "LOH" = "purple", 
         "AMP" = "red",
          "MISSENSE" = "#34A02C",
         "TRUNC" = "black",
         "DelMISSENSE" = "grey")
#"GAIN" = "#FB9FB5",

 # gl <- genelist
```


```{r}

mn <-  mat2[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification" )]
colnames(mn) <- c("sample", "Gene.refGene", "input")
mmn <- rbind(mn[, c("sample", "Gene.refGene", "input")],cnv[, c("sample", "Gene.refGene", "input")])
mmn <- mmn[mmn$Gene.refGene %in% genelist &
           mmn$sample %in%  unique(mat2$Tumor_Sample_Barcode),]

mat.unique <- unique(mmn[ ,c("sample", "Gene.refGene")])
colnames(mat.unique) <- c("tumor_name","gene")
sample <- unique(mmn$sample)
gene <-   unique(mmn$Gene.refGene)
mat.u <- matrix(NA, ncol = length(sample), nrow = length(gene))
mmn <- data.frame(mmn)

for (i in 1:length(gene)){
    for (j in 1:length(sample)){
        data2 <- mmn[mmn$Gene.refGene %in% gene[i] & mmn$sample %in% sample[j], ]
        input <- unlist(data2$input)
        if (length(input) > 0){
            INPUT <- paste(input[1], ";", sep = "")
            if (length(input) > 1){
                for (k in 2:length(input)){
                    INPUT <- paste(INPUT, paste(input[k], ";", sep = ""), sep = "")
                }
            }
                 #  print(INPUT)
        mat.u[i, j] <- INPUT
        }

     } 
}
 

rownames(mat.u) <- gene
colnames(mat.u) <- sample
mat.u[is.na(mat.u)] <- ""
dim(mat.u)
mat.2<- mat.u
dim(mat.2)
```

```{r,  fig.width=10, fig.height=8}
#onco(mat.2 )

mat.2_PDX_A <- mat.2[,unique(mat2_PDX_A$Tumor_Sample_Barcode)]
dim(mat.2_PDX_A)
mat.2_PDX_B <- mat.2[,unique(mat2_PDX_B$Tumor_Sample_Barcode)]
mat.2_PDX_C <- mat.2[,unique(mat2_PDX_C$Tumor_Sample_Barcode)]
mat.2_PDX_D <- mat.2[,unique(mat2_PDX_D$Tumor_Sample_Barcode)]
mat.2_PDX_E <- mat.2[,unique(mat2_PDX_E$Tumor_Sample_Barcode)]
mat.2_PDX_F <- mat.2[,unique(mat2_PDX_F$Tumor_Sample_Barcode)]

mat.2_R <- mat.2[,colnames(mat.2) %in% unique(mat2_R$Tumor_Sample_Barcode)]
mat.2_S <- mat.2[,colnames(mat.2) %in% unique(mat2_S$Tumor_Sample_Barcode)]

mat.2_R2 <- mat.2[,colnames(mat.2) %in%unique(mat2_R2$Tumor_Sample_Barcode)]
mat.2_S2 <- mat.2[,colnames(mat.2) %in%unique(mat2_S2$Tumor_Sample_Barcode)]
```

```{r,  fig.width=20, fig.height=25}
onco2 <- function(mat){
return (
  oncoPrint(mat,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_names_gp = gpar(fontsize = 9, fontface = "italic"),
          column_names_gp = gpar(fontsize =6),
          pct_gp = gpar(fontsize = 6),
         # labels_gp = gpar(fontsize = 12),
          # row_order = row.order,
          # column_order = sample.order$name,
          show_column_names = TRUE,
          show_row_names = TRUE,
         
          #column_title = "Landscape of MCL mutation",
          remove_empty_columns = FALSE,
          heatmap_legend_param = list(title = "Alterations",
                                      at = c("AMP",
                                             "DEL",
                                             "LOH" ,
                                             "TRUNC",
                                             
                                             "MISSENSE",
                                             "DelMISSENSE" ),
                                      labels = c(  "Amplification", 
                                        "Deletion",
                                                  "LOH",
                                                  "Truncating",
                                                   "Missense",
                                                 "Deleterious \n missense") ) )
)
}

onco <- function(mat){
return (
  oncoPrint(mat,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_names_gp = gpar(fontsize = 9, fontface = "italic"),
          column_names_gp = gpar(fontsize =6),
          pct_gp = gpar(fontsize = 6),
         # labels_gp = gpar(fontsize = 12),
          # row_order = row.order,
          # column_order = sample.order$name,
          show_column_names = TRUE,
          show_row_names = F,
         
          #column_title = "Landscape of MCL mutation",
          remove_empty_columns = FALSE,
          heatmap_legend_param = list(title = "Alterations",
                                      at = c("AMP",
                                             "DEL",
                                             "LOH" ,
                                             "TRUNC",
                                             
                                             "MISSENSE",
                                             "DelMISSENSE" ),
                                      labels = c(  "Amplification", 
                                        "Deletion",
                                                  "LOH",
                                                  "Truncating",
                                                   "Missense",
                                                 "Deleterious \n missense") ) )
)
}




```



```{r,  fig.width=20, fig.height=20}
gl2 <-   genelist
  onco(mat.2_R[gl2,] ) + 
  onco2(mat.2_S[gl2,] ) +
  onco(mat.2_PDX_A[gl2,] ) +
  onco(mat.2_PDX_B[gl2,] ) +
  onco2(mat.2_PDX_C[gl2,] ) +
  onco(mat.2_PDX_D[gl2,] ) +
  onco(mat.2_PDX_E[gl2,] ) +
  onco(mat.2_PDX_F[gl2,] ) 
```


```{r,  fig.width=20, fig.height=4}
gl2 <-  intersect(pathways.all[["REACTOME_DEVELOPMENTAL_BIOLOGY"]],  genelist)
onco(mat.2_R[gl2,] ) + onco(mat.2_S[gl2,] ) +
  onco(mat.2_PDX_A[gl2,] ) +
  onco(mat.2_PDX_B[gl2,] ) +
  onco(mat.2_PDX_C[gl2,] ) +
  onco(mat.2_PDX_D[gl2,] ) +
  onco(mat.2_PDX_E[gl2,] ) +
  onco2(mat.2_PDX_F[gl2,] ) 
```
```{r}
d_all[d_all$Hugo_Symbol %in% "TUBA3C",]
```

```{r,  fig.width=20, fig.height=4}
gl2 <-  intersect(pathways.all[["REACTOME_CELL_CYCLE"]],  genelist)
onco(mat.2_R[gl2,] ) + onco(mat.2_S[gl2,] ) +
  onco(mat.2_PDX_A[gl2,] ) +
  onco(mat.2_PDX_B[gl2,] ) +
  onco(mat.2_PDX_C[gl2,] ) +
  #onco(mat.2_PDX_D[gl2,] ) +
  onco(mat.2_PDX_E[gl2,] ) +
  onco2(mat.2_PDX_F[gl2,] ) 
```


```{r,  fig.width=20, fig.height=20}
gl2 <- genelist
onco2(mat.2_R[gl2,] ) + onco(mat.2_S[gl2,] ) +
  onco(mat.2_PDX_A[gl2,] ) +
  onco2(mat.2_PDX_B[gl2,] ) +
  onco(mat.2_PDX_C[gl2,] ) +
  onco(mat.2_PDX_D[gl2,] ) +
  onco(mat.2_PDX_E[gl2,] ) +
  onco(mat.2_PDX_F[gl2,] ) 
```


```{r,  fig.width=20, fig.height=7}
gl2 <-  intersect(c(pathways.all[["REACTOME_ADAPTIVE_IMMUNE_SYSTEM"]],
                            pathways.all[["REACTOME_INNATE_IMMUNE_SYSTEM"]]),
                            genelist)
onco(mat.2_R[gl2,] ) + onco(mat.2_S[gl2,] ) +
  onco(mat.2_PDX_A[gl2,] ) +
  onco(mat.2_PDX_B[gl2,] ) +
  onco(mat.2_PDX_C[gl2,] ) +
  onco(mat.2_PDX_D[gl2,] ) +
  onco(mat.2_PDX_E[gl2,] ) +
  onco2(mat.2_PDX_F[gl2,] ) 
```

