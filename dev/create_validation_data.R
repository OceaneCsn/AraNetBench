######## dap seq ###########


# data from http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php

library(DBI )
library( RSQLite )
library( dplyr )
library( dbplyr )
con <- dbConnect(RSQLite::SQLite(), "D:/These/DAPSeqData/databaseDapSeq.db")
src_dbi(con)

query <- dbSendQuery(con, "SELECT * FROM binding WHERE amplified = 0")
dbBind(query, list())
res <- data.frame(dbFetch(query))
res <- res[1:2]
colnames(res) <- c("from", "to")
dap_seq <- res
dap_seq$type <- "DAPSeq"

####### chip seq ##########

# data from https://connectf.org/query, with all_in_planta_bound


chip_seq <- read.csv("D:/These/ValidationRegulation/connecTF_chipSeq.csv")
chip_seq <- chip_seq[c("gene_id", "TARGET")]
colnames(chip_seq) <- c("from", "to")
chip_seq$type <- "CHIPSeq"


######### TARGET ###################

# downloaded from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09522-1/MediaObjects/41467_2019_9522_MOESM5_ESM.xlsx

target <- read.csv("D:/These/ValidationRegulation/target.csv", 
                   skip = 2, h = T, sep = ';')
target <- target[c("TF.GeneID", "Target")]
colnames(target) <- c("from", "to")
target$type <- "TARGET"


######### ATRM (litterature) ###################

# downloaded from http://atrm.gao-lab.org/for_download/Regulations_in_ATRM.xlsx

atrm <- read.csv("D:/These/ValidationRegulation/ATRM.csv", 
                   h = T, sep = ',')
atrm <- atrm[c("TF.ID", "Target.ID")]
colnames(atrm) <- c("from", "to")
atrm$type <- "Litterature"



######### ATREGNET ###################

# downloaded from https://agris-knowledgebase.org/downloads.html



# 1. The TF binds directly to the regulatory region of the target gene (as shown by electromobility shift assay (EMSA), yeast one-hybrid analysis, and/or chromatin immunoprecipitation (ChIP)).
# OR
# 2. The TF directly regulates the target gene, based on use of transgenic plants expressing an inducible TF-GR (glucocorticoid receptor) fusion protein. Fusion of TFs to the GR domain has been demonstrated to cause the fusion protein to be retained in the cytoplasm in the absence of steroid hormone. Upon application of the synthetic hormone dexamethasome (DEX), the protein moves into the nucleus and activates expression. Direct targets are those that are either induced or repressed by a given TF-GR upon DEX induction even in the presence of the translation inhibitor cycloheximide (CYC).
# AND
# 3. In vivo evidence of regulation: Expression of the target gene is affected by either loss of function mutations in the TF or ectopic expression of that TF in the plant.
# 
# atregnet <- read.csv("D:/These/ValidationRegulation/AtRegNet.csv")
# 
# atrm <- atrm[c("TF.ID", "Target.ID")]
# colnames(atrm) <- c("from", "to")
# atrm$type <- "Litterature"



validated_edges <- rbind.data.frame(chip_seq, target, atrm, dap_seq)

usethis::use_data(validated_edges, version = 3, overwrite = T)


#save(validated_edges, file = "D:/These/ValidationRegulation/validated_edges_all_data_sources.rdata")

table(validated_edges$type)





