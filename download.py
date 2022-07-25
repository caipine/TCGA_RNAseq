import os
import sys        
sys.path.append('../') 
from tcga_downloader import *
ids=get_ids('gdc_manifest_20220725_200519.txt')
payload=prepare_payload(ids,data_type='Gene Expression Quantification')
metadata=get_metadata(payload)
download_data(metadata,sep='\t', outdir='SKCM')

