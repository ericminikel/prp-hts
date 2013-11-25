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
dim(fehta_primary) # 370276 7

# 2. PrP 5'UTR inhibitor screening results:
# Screen summary: http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=488894&loc=ea_ras

utr_primary = read.csv('5utr/AID_488862_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=488862&q=expdata_csvsave
utr_cherrypick_dose = read.csv('5utr/AID_504539_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=504539&q=expdata_csvsave
utr_cherrypick_sing = read.csv('5utr/AID_504592_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=504592&q=expdata_csvsave
utr_drypowder_dose = read.csv('5utr/AID_504932_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=504932&q=expdata_csvsave
utr_control = read.csv('5utr/AID_540255_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=540255&q=expdata_csvsave
utr_counter_gfp = read.csv('5utr/AID_540343_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=540343&q=expdata_csvsave
utr_counter_app = read.csv('5utr/AID_588451_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=588451&q=expdata_csvsave
utr_counter_luc = read.csv('5utr/AID_588507_data.csv')  # http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=588507&q=expdata_csvsave

dim(utr_primary) # 335011 11

#####
# figure out some basics of how the data are organized

# are they the same compounds in both screens?
# using PUBCHEM_SID
sum(fehta_primary$PUBCHEM_SID %in% utr_primary$PUBCHEM_SID)   # 323084
sum(!(fehta_primary$PUBCHEM_SID %in% utr_primary$PUBCHEM_SID)) # 47192
sum(!(utr_primary$PUBCHEM_SID %in% fehta_primary$PUBCHEM_SID)) # 11927

# using PUBCHEM_CID
sum(fehta_primary$PUBCHEM_CID %in% utr_primary$PUBCHEM_CID)   # 328985
sum(!(fehta_primary$PUBCHEM_CID %in% utr_primary$PUBCHEM_CID)) # 41291
sum(!(utr_primary$PUBCHEM_CID %in% fehta_primary$PUBCHEM_CID)) #  5942

# what is the difference between the SID and CID
sum(is.na(fehta_primary$PUBCHEM_SID)) #0
sum(is.na(fehta_primary$PUBCHEM_CID)) #2
dim(fehta_primary) # 370726 7
length(unique(fehta_primary$PUBCHEM_SID)) #370726 - this is a unique identifier
length(unique(fehta_primary$PUBCHEM_CID)) #369940 - not quite unique

# same for the utr screen?
dim(utr_primary) # 335011 11
length(unique(utr_primary$PUBCHEM_SID)) #335011 - this is a unique identifier
length(unique(utr_primary$PUBCHEM_CID)) #334622 - not quite unique


# how is the 5'UTR screen data stored?
head(utr_primary)
names(utr_primary)

# is the activity score the readout?
hist(utr_primary$PUBCHEM_ACTIVITY_SCORE)
# does it differ between compounds dubbed "Active" vs. "Inactive"?
t.test(utr_primary$PUBCHEM_ACTIVITY_SCORE ~ utr_primary$PUBCHEM_ACTIVITY_OUTCOME)
# yes

# what is the cutoff?
hist(utr_primary$PUBCHEM_ACTIVITY_SCORE[utr_primary$PUBCHEM_ACTIVITY_OUTCOME=='Active'])
# apparently a score of > 50 gets you labeled "Active"
max(utr_primary$PUBCHEM_ACTIVITY_SCORE[utr_primary$PUBCHEM_ACTIVITY_OUTCOME=='Inactive'])
# and all of <= 50 are labeled "Inactive"

# what are these other columns?
head(unique(utr_primary$PUBCHEM_ASSAYDATA_COMMENT))
# all NA

plot(utr_primary$REPLICATE_A_ACTIVITY_SCORE, utr_primary$REPLICATE_B_ACTIVITY_SCORE)
# largely reproducible

hist(utr_primary$REPRODUCIBILITY_COSINE_TRANSFORM)
# ranges 0 to 1, heavier at the 1 end

# is the reproducibility score just capturing how well the 4 replicates match each other?
# plot variance of 4 replicates vs. reproducibility score
plot(apply(utr_primary[,8:11],1,var,na.rm=TRUE),utr_primary$REPRODUCIBILITY_COSINE_TRANSFORM)
# no, that's not it at all.

# how about the fehta screen?
head(fehta_primary)
# is this column the readout?
hist(fehta_primary$Inhibition.at.13.8.uM) # yes, perhaps
length(unique(fehta_primary$PUBCHEM_ACTIVITY_OUTCOME)) # 2 - Inactive and Active
t.test(fehta_primary$Inhibition.at.13.8.uM ~ fehta_primary$PUBCHEM_ACTIVITY_OUTCOME)
# yes, big difference

plot(fehta_primary$Inhibition.at.13.8.uM, fehta_primary$PUBCHEM_ACTIVITY_SCORE)
# activity score ranges 0 - 100.

max(fehta_primary$Inhibition.at.13.8.uM) # 135.29
# so definitely not a percent reduction in signal

# are hits in one assay more likely to be hits in the other assay than by chance?
shared = merge(fehta_primary[,c("PUBCHEM_SID","PUBCHEM_ACTIVITY_OUTCOME")],
               utr_primary[,c("PUBCHEM_SID","PUBCHEM_ACTIVITY_OUTCOME")],
               by="PUBCHEM_SID")
colnames(shared) = c("PUBCHEM_SID", "FEHTA", "UTR")
table(shared[,2:3])
#            UTR
# FEHTA      Active Inactive
# Active       70     1263
# Inactive    994   320757
fisher.test(table(shared[,2:3]),alternative="two.sided")
# p < 2e-16, odds ratio 17
# so yes, hits from one assay are very strongly enriched in the other

# how were the compounds chosen for followup assays

followup_sids = unique(c(utr_cherrypick_dose$PUBCHEM_SID,
                  utr_cherrypick_sing$PUBCHEM_SID,
                  utr_drypowder_dose$PUBCHEM_SID,
                  utr_control$PUBCHEM_SID,
                  utr_counter_gfp$PUBCHEM_SID,
                  utr_counter_app$PUBCHEM_SID,
                  utr_counter_luc$PUBCHEM_SID))
length(followup_sids) # 1144
mean(followup_sids %in% utr_primary$PUBCHEM_SID) # 98.07% were in the primary
# wonder why the other 2% were counterscreened despite not being in primary.

table(utr_primary$PUBCHEM_ACTIVITY_OUTCOME)
# Active 1169, Inactive 333842
# so they did followup work on almost all of the actives.

table(utr_primary$PUBCHEM_ACTIVITY_OUTCOME[utr_primary$PUBCHEM_SID %in% followup_sids])
# Active 1121, Inactive 1
# hmm, so of the 1144 substances from the followup screens,
# 1121 were actives, 1 was inactive and 22 were not screened in the primary.

# fyi - all of the other 22 are in the control screen
sum(!(utr_control$PUBCHEM_SID %in% utr_primary$PUBCHEM_SID))
# so if I can get the structures for the primary + the control, I've got all.

# finding the structure data for these
utr_control_smiles = read.table('5utr/utr-control-smiles-386951438944468372.txt',comment.char='',quote='',sep='\t')
utr_primary_smiles = read.table('5utr/utr-primary-smiles-3818840491230252424.txt',comment.char='',quote='',sep='\t')
fehta_primary_smiles = read.table('prpfehta/fehta-primary-smiles-4280565902671935407.txt',comment.char='',quote='',sep='\t')

all_smiles = unique(rbind(utr_control_smiles,utr_primary_smiles,fehta_primary_smiles))
dim(all_smiles) # 382225 2
colnames(all_smiles) = c("PUBCHEM_SID","SMILES")

