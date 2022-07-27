
library("survival")
library("survminer")
library(stringr)

source("https://raw.githubusercontent.com/caipine/TCGA_RNAseq/main/gene_sur3.fucntion.R")
pv <- c()
pvn <- c()
n <- 1
sample <- read.csv2(paste0( dir( pattern = "gdc_sample_sheet") ), sep = "\t")

for (i in unique(sample$Project.ID) )   {
  print(n)
  if (file.exists(paste0("'", i, "'__deseq2-2.RDA") )) {
#i <- unique(sample$Project.ID)[3]
  print(paste0("'", i, "'__deseq2-2.RDA exists"))
  load( paste0("'", i, "'__deseq2-2.RDA") )
d1 <- gene_sur3(c("ADSL", "ATIC"  ))
d1$Group <- c(sample("Low", round(nrow(d1)/2) , replace = T),
              sample("High", nrow(d1) - round(nrow(d1)/2) , replace = T))

fit <- survfit(Surv(sur, alive) ~ Group, data = d1)
pv <- c(pv, surv_pvalue(fit)$pval)
pvn <- c(pvn, i)

print(
  ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE, title =paste0( "ACLY_PUR_in__", i),
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(base_size = 15))
)
  }
n = n+1
}
pv2 <- data.frame(TCGA = pvn, H50_L50 = pv)
