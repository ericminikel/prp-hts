# Eric Minikel
# CureFFI.org
# 2013-11-17
# Analysis of high throughput screening data from two screens

require(rcdk)
options(stringsAsFactors=FALSE)
setwd('c:/sci/042prphts/data/')

# 1. PrP-FEHTA screening results:
# Screen summary: http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=720596

# primary screen data: # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=720596&q=expdata_csvsave
fehta_primary = read.csv('prpfehta/AID_720596_data.csv') 


# 2. PrP 5'UTR inhibitor screening results:
# Screen summary: http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=488894&loc=ea_ras

# primary screen data: # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=488862&q=expdata_csvsave
utr_primary = read.csv('5utr/AID_488862_data.csv') 
