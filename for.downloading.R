
TCGA <- dir(pattern = "TCGA-")
  for (i in TCGA) {
    #i = TCGA[1 ]
    print(i)
gdc_manifest <-  dir(i, pattern = "gdc_manifest")

m1 <- c( "import os",
         "import sys",
         "sys.path.append('../') ", 
         "from tcga_downloader import * ",
         "ids=get_ids('XXX.txt')", 
         "payload=prepare_payload(ids,data_type='Gene Expression Quantification')",
         "metadata=get_metadata(payload)",
         "download_data(metadata,sep='\t', outdir='XXX')"
  )
m1[5] <- paste0("ids=get_ids('", gdc_manifest, "')" )
m1[8] <- paste0("download_data(metadata,sep='\\t', outdir='", str_replace(i, "TCGA-", "" ), "')" )
write.table(matrix(m1), file = paste0(i,"/", "downloading.", str_replace(i, "TCGA-", "" ),".py"), row.names =F, quote = F,col.names = F )
  }
