library(openxlsx)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(pROC)
library(xgboost)
library(plotly)
library(stringi)
library(dummies)
library(dummy)

dataChr<- read.xlsx('sources_r/clinvar_conflicting.xlsx', sheet = 1, startRow = 1, colNames = TRUE)


dataChr<- read.csv('sources_r/clinvar_conflicting.csv', sep=",", header = TRUE)



dataChr[, "CHROM"] <- as.character(dataChr[, "CHROM"])

for (ln in nrow(dataChr):1){
  if (grepl("MT", dataChr[ln,"CHROM"]) == TRUE) {
    dataChr <- dataChr[-ln,]
  }
}

for (ln in nrow(dataChr):40000){
  if (grepl("X", dataChr[ln,"CHROM"]) == TRUE) {
    dataChr[ln,"CHROM"] <- "23"
  }
}

dataChr[, "CHROM"] <- as.factor(dataChr[, "CHROM"])

str(dataChr)

#check missing values
sum(sapply(dataChr, function(y) sum(is.na(y) | is.null(y)))) != 0

sum(sapply(dataChr, function(y) sum(is.na(y) | is.null(y))))

# missing value by column
apply(dataChr, 2, function(x) sum(is.na(x) | is.null(x)))


## CHROM            POS                REF                    ALT             AF_ESP           AF_EXAC         AF_TGP 
## 0                0                  0                      0               0                0               0 
##
## CLNDISDB         CLNDISDBINCL       CLNDN                  CLNDNINCL       CLNHGVS          CLNSIGINCL      CLNVC 
## 0                65112              0                      65112           0                65112           0 
##
## CLNVI            MC                 ORIGIN                 SSR             CLASS            Allele          Consequence 
## 37529            6969               6123                   65084           0                0               0 
##
## IMPACT           SYMBOL             Feature_type           Feature         BIOTYPE          EXON            INTRON 
## 0                16                 14                     14              16               8893            56385 
##
## cDNA_position    CDS_position       Protein_position       Amino_acids     Codons           DISTANCE        STRAND 
## 8884             9955               9955                   10004           10004            65080           14 
##
## BAM_EDIT         SIFT               PolyPhen               MOTIF_NAME      MOTIF_POS       HIGH_INF_POS     MOTIF_SCORE_CHANGE 
## 33219            40352              40392                  65186           65186           65186            65186 
##
## LoFtool          CADD_PHRED         CADD_RAW               BLOSUM62 
## 4213             1092               1092                   39595

unique(dataChr[, "CHROM"])
unique(dataChr[, "REF"]) # 866 lvl (faut à mon avis drop les lignes qui ont plus de 1 lettre)
unique(dataChr[, "ALT"]) # 458 lvl (pareil que pour REF)
unique(dataChr[, "Allele"]) # 374 lvl
unique(dataChr[, "Feature_type"]) # 2 lvl 
unique(dataChr[, "CLNDISDB"]) #9234 Levels
unique(dataChr[, "CLNDISDBINCL"]) #49 Levels => 65112 missing val
unique(dataChr[, "CLNDN"]) #9260 Levels
unique(dataChr[, "CLNDNINCL"]) #55 levels => 65112 missing val
unique(dataChr[, "CLNHGVS"]) #65188 Levels => à retrier en enlevant POS, REF, ALT de la feature (NC_000001.10:g.40770175C>T => après le "g.") 
unique(dataChr[, "CLNSIGINCL"])# 69 Levels => 65112 missing val
unique(dataChr[, "CLNVC"])# 7 Levels => à séparer en 7 colonnes
unique(dataChr[, "CLNVI"]) # 26290 Levels
unique(dataChr[, "BIOTYPE"]) # 2 lvl 
unique(dataChr[, "MC"]) # 90 Levels => faire un multi découpage car beaucoup de lvl viennent du fait de plusieurs ref en 1 ligne 
#(ex SO:0001583 et SO:0001819 présents sur la même ligne => rajoute 1 lvl au lieu de compter chacun une fois)
# => sans doute à découper sur une dizaine de colonnes
unique(dataChr[, "ORIGIN"]) #27 Levels
unique(dataChr[, "SSR"]) #2 Levels => 65084 missing val
unique(dataChr[, "Consequence"]) #48 Levels (pareil que MC, c'est du à des "doublons" de cas sur une même ligne)
unique(dataChr[, "IMPACT"]) #4 lvl
unique(dataChr[, "SYMBOL"]) #2329 LVL
unique(dataChr[, "STRAND"]) #


##### CHROM #####
prop.table(table(dataChr$CHROM,dataChr$CLASS), 1)
table(dataChr$CHROM) 
chisq.test(dataChr$CHROM,dataChr$CLASS) # seems ok to conserve


##### SSR #####
prop.table(table(dataChr$SSR,dataChr$CLASS), 1)
table(dataChr$SSR) 
chisq.test(dataChr$SSR,dataChr$CLASS) # 0.61 => pas retenue


##### CLNVC #####

prop.table(table(dataChr[!(dataChr$CLNVC %in% c("Microsatellite", "Inversion", "Insertion")), "CLNVC"],dataChr[!(dataChr$CLNVC %in% c("Microsatellite", "Inversion", "Insertion")),"CLASS"]), 1) # intéressant aussi
table(dataChr$CLNVC)
chisq.test(dataChr[!(dataChr$CLNVC %in% c("Microsatellite", "Inversion", "Insertion")),"CLNVC"],dataChr[!(dataChr$CLNVC %in% c("Microsatellite", "Inversion", "Insertion")),"CLASS"])
# après analyse => pas retenu



##### ORIGIN ######

prop.table(table(dataChr$ORIGIN,dataChr$CLASS), 1)
sum(table(dataChr$ORIGIN))
chisq.test(dataChr$ORIGIN,dataChr$CLASS) # ORIGIN discriminant 

# après analyse la catégorie 1 comporte plus de 95% des valeurs (si on ne prend pas en compte les missings values)


##### BIOTYPE #####
for(ln in nrow(dataChr):1) {
  if (sapply(dataChr[ln,"BIOTYPE"], function(y) is.na(y) | is.null(y)) == TRUE) {
    dataChr <- dataChr[-ln, ]
  }
}

prop.table(table(dataChr$BIOTYPE, dataChr$CLASS), 1)
table(dataChr$BIOTYPE) 
chisq.test(dataChr$BIOTYPE,dataChr$CLASS) # bof


###### Allele #####

dataChr$Allele <- as.character(dataChr$Allele)
unique(dataChr[, "Allele"]) # 

for (ln in nrow(dataChr):1){
  if (nchar(dataChr[ln, "Allele"], type = "chars")>1){
    dataChr <- dataChr[-ln, ]
  }
}



unique(dataChr[, "Allele"])

dataChr$Allele <- as.factor(dataChr$Allele)

temp_df <- dataChr[,c("Allele", "CLASS")]
temp_df2 <- dataChr[,c("Allele", "CLASS")]
temp_df <- dummy.data.frame(temp_df)
temp_df <- temp_df[, -5]

dataChr <- cbind(dataChr, temp_df)


prop.table(table(dataChr$AlleleA,dataChr$CLASS), 1)
table(dataChr$AlleleA)
chisq.test(dataChr$AlleleA, dataChr$CLASS)

prop.table(table(dataChr$AlleleC,dataChr$CLASS), 1)
table(dataChr$AlleleC)
chisq.test(dataChr$AlleleC, dataChr$CLASS)

prop.table(table(dataChr$AlleleT,dataChr$CLASS), 1)
table(dataChr$AlleleT)
chisq.test(dataChr$AlleleT, dataChr$CLASS)

prop.table(table(dataChr$AlleleG,dataChr$CLASS), 1)
table(dataChr$AlleleG)
chisq.test(dataChr$AlleleG, dataChr$CLASS)

###### feature type ######
prop.table(table(dataChr$Feature_type, dataChr$CLASS), 1)
table(dataChr$Feature_type) 
chisq.test(dataChr$Feature_type,dataChr$CLASS) # pas utilisable on dirait



###### Feature ######
unique(dataChr[, "Feature"])  

dataChr[, "New_Feature"] <- ""

for(ln in 1:nrow(dataChr)){
  dataChr[ln,"New_Feature"] <- stri_split_fixed(dataChr[ln,"Feature"], "_") [[1]][1]
}


dataChr[!(dataChr$Feature %in% c("_")), "New_Feature"]

unique(dataChr[, "New_Feature"])  
prop.table(table(dataChr$New_Feature, dataChr$CLASS), 1)
table(dataChr$New_Feature) 
chisq.test(dataChr$New_Feature,dataChr$CLASS) 
# pas convainquant pour le moment


######## STRAND ########

unique(dataChr[, "STRAND"])  

for (ln in 1:nrow(dataChr)){
  if ((dataChr[ln,"STRAND"]) == -1) {
    dataChr[ln, "STRAND"] <- 0
  }
}

prop.table(table(dataChr$STRAND, dataChr$CLASS), 1)
table(dataChr$STRAND) 
chisq.test(dataChr$STRAND,dataChr$CLASS) # ça à l'air OK


######## traitement de "Consequence" ######

unique(dataChr[, "Consequence"])
#découpage de la colonne "Consequence"
val1 <- "missense_variant"
dataChr[, "missense_variant"] <- 0
val2 <- "synonymous_variant"
dataChr[, "synonymous_variant"] <- 0
val3 <- "splice_region_variant"
dataChr[, "splice_region_variant"] <- 0
val4 <- "intron_variant"
dataChr[, "intron_variant"] <- 0
val5 <- "3_prime_UTR_variant"
dataChr[, "3_prime_UTR_variant"] <- 0
val6 <- "frameshift_variant"
dataChr[, "frameshift_variant"] <- 0
val7 <- "inframe_insertion"
dataChr[, "inframe_insertion"] <- 0
val8 <- "inframe_deletion"
dataChr[, "inframe_deletion"] <- 0
val9 <- "stop_lost"
dataChr[, "stop_lost"] <- 0
val10 <- "stop_gained"
dataChr[, "stop_gained"] <- 0
val11 <- "5_prime_UTR_variant"
dataChr[, "5_prime_UTR_variant"] <- 0
val12 <- "splice_acceptor_variant"
dataChr[, "splice_acceptor_variant"] <- 0
val13 <- "coding_sequence_variant"
dataChr[, "coding_sequence_variant"] <- 0
val14 <- "start_lost"
dataChr[, "start_lost"] <- 0
val15 <- "downstream_gene_variant"
dataChr[, "downstream_gene_variant"] <- 0
val16 <- "splice_donor_variant"
dataChr[, "splice_donor_variant"] <- 0
val17 <- "stop_retained_variant"
dataChr[, "stop_retained_variant"] <- 0
val18 <- "non_coding_transcript_variant"
dataChr[, "non_coding_transcript_variant"] <- 0
val19 <- "start_retained_variant"
dataChr[, "start_retained_variant"] <- 0
val20 <- "protein_altering_variant"
dataChr[, "protein_altering_variant"] <- 0
val21 <- "TF_binding_site_variant"
dataChr[, "TF_binding_site_variant"] <- 0
val22 <- "intergenic_variant"
dataChr[, "intergenic_variant"] <- 0



#val1:22
for (ln in 1:nrow(dataChr)){
  if (grepl(val1, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val1] <- 1
  }
  if (grepl(val2, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val2] <- 1
  }
  if (grepl(val3, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val3] <- 1
  }
  if (grepl(val4, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val4] <- 1
  }
  if (grepl(val5, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val5] <- 1
  }
  if (grepl(val6, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val6] <- 1
  }
  if (grepl(val7, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val7] <- 1
  }
  if (grepl(val8, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val8] <- 1
  }
  if (grepl(val9, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val9] <- 1
  }
  if (grepl(val10, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val10] <- 1
  }
  if (grepl(val11, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val11] <- 1
  }
  if (grepl(val12, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val12] <- 1
  }
  if (grepl(val13, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val13] <- 1
  }
  if (grepl(val14, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val14] <- 1
  }
  if (grepl(val15, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val15] <- 1
  }
  if (grepl(val16, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val16] <- 1
  }
  if (grepl(val17, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val17] <- 1
  }
  if (grepl(val18, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val18] <- 1
  }
  if (grepl(val19, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val19] <- 1
  }
  if (grepl(val20, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val20] <- 1
  }
  if (grepl(val21, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val21] <- 1
  }
  if (grepl(val22, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val22] <- 1
  }
}

for (ln in 1:nrow(dataChr)){
  if (grepl(val4, dataChr[ln,"Consequence"]) == TRUE) {
    dataChr[ln, val4] <- 1
    number_of_1 <<- number_of_1 + 1 
  }
}
number_of_1

for (valu in 47:68){
  for (ln in 1:nrow(dataChr)){
    print(unique(dataChr[,valu]))
  }
}


######## traitement de MC #########
unique(dataChr[, "MC"]) 
not_in_consequence_MC <- " "

val23 <- "nonsense"
dataChr[, "nonsense"] <- 0
val24 <- "2KB_upstream_variant"
dataChr[, "2KB_upstream_variant"] <- 0

for (ln in 1:nrow(dataChr)){
  if (grepl(val23, dataChr[ln,"MC"])) {
    dataChr[ln, val23] <- 1
    }
  if (grepl(val24, dataChr[ln,"MC"])) {
    dataChr[ln, val24] <- 1
  }
 
}
#variables suivantes ajoutées
#SO:0001587|nonsense              => nonsense
#SO:0001636|2KB_upstream_variant  => 2KB_upstream_variant

######## EXON + INTRON ######

dataChr$EXON <- as.character(dataChr$EXON)
dataChr$INTRON <- as.character(dataChr$INTRON)
test2 <- ""
for (ln in 1:nrow(dataChr)){
  if ((dataChr[ln,"INTRON"] != "") && (dataChr[5,"EXON"] != "")) {
    test2 <<- ln
  }
}

dataChr[,"EXONINTRON"] <- ""

for (ln in 1:nrow(dataChr)){
  if(dataChr[ln,"EXONINTRON"] == ""){
    dataChr[ln,"EXONINTRON"] <- dataChr[ln,"INTRON"]
  }
}

for (ln in 1:nrow(dataChr)){
  if(dataChr[ln,"EXONINTRON"] == ""){
    dataChr[ln,"EXONINTRON"] <- dataChr[ln,"EXON"]
  }
}

nbr <- 0
for (ln in 1:nrow(dataChr)){
  if(dataChr[ln,"EXONINTRON"] == ""){
    nbr <<- nbr+1
  }
}
nbr

for(ln in 1:nrow(dataChr)){
  dataChr[ln,"EXONINTRON"] <- as.numeric(stri_split_fixed(dataChr[ln,"EXONINTRON"], "/") [[1]][1]) / as.numeric(stri_split_fixed(dataChr[ln,"EXONINTRON"], "/") [[1]][2])
}

dataChr$EXONINTRON <- as.numeric(dataChr$EXONINTRON)

chisq.test(dataChr$EXONINTRON,dataChr$CLASS) # à garder => après analyse boxplot c'est pas discriminant


######## NEW col Impact ########

unique(dataChr[, "IMPACT"]) #4 lvl

prop.table(table(dataChr$IMPACT,dataChr$CLASS), 1) # => très intéressant
table(dataChr$IMPACT)
chisq.test(dataChr$IMPACT,dataChr$CLASS) # rejet de l'hypothèse nulle => dépendance entre les colonnes donc à garder.

var_NewImp <- "HIGH"
dataChr[, "New_IMPACT"] <- 0

for (ln in 1:nrow(dataChr)){
  if (grepl(var_NewImp, dataChr[ln,"IMPACT"]) == TRUE) {
    dataChr[ln, "New_IMPACT"] <- 1
  }
}


unique(dataChr[, "New_IMPACT"])
prop.table(table(dataChr$New_IMPACT,dataChr$CLASS), 1) # => très intéressant
table(dataChr$New_IMPACT)
chisq.test(dataChr$New_IMPACT,dataChr$CLASS) # => on garde (a remplacé "IMPACT")


###### REF et ALT #########


dataChr$ALT <- as.character(dataChr$ALT)
dataChr$REF <- as.character(dataChr$REF)
unique(dataChr[, "ALT"]) # 

for (ln in nrow(dataChr):1){
  if (nchar(dataChr[ln, "ALT"], type = "chars")>1){
    dataChr <- dataChr[-ln, ]
  }
}

for (ln in nrow(dataChr):1){
  if (nchar(dataChr[ln, "REF"], type = "chars")>1){
    dataChr <- dataChr[-ln, ]
  }
}

dataChr$ALT <- as.factor(dataChr$ALT)
dataChr$REF <- as.factor(dataChr$REF)

prop.table(table(dataChr[,"ALT"],dataChr$CLASS), 1)
table(dataChr[,"ALT"]) 
chisq.test(dataChr[,"ALT"],dataChr$CLASS) # OK

prop.table(table(dataChr[,"REF"],dataChr$CLASS), 1)
table(dataChr[,"REF"]) 
chisq.test(dataChr[,"REF"],dataChr$CLASS) # OK


temp_df <- dataChr[,c("REF", "CLASS")]
temp_df2 <- dataChr[,c("REF", "CLASS")]
temp_df <- dummy.data.frame(temp_df)
temp_df <- temp_df[, -5]

dataChr <- cbind(dataChr, temp_df)


prop.table(table(dataChr$REFA,dataChr$CLASS), 1)
table(dataChr$REFA)
chisq.test(dataChr$REFA, dataChr$CLASS)

prop.table(table(dataChr$REFC,dataChr$CLASS), 1)
table(dataChr$REFC)
chisq.test(dataChr$REFC, dataChr$CLASS)

prop.table(table(dataChr$REFT,dataChr$CLASS), 1)
table(dataChr$REFT)
chisq.test(dataChr$REFT, dataChr$CLASS)

prop.table(table(dataChr$REFG,dataChr$CLASS), 1)
table(dataChr$REFG)
chisq.test(dataChr$REFG, dataChr$CLASS)


########## new col CLNDISDB #########

unique(dataChr[, "MedGen"])

dataChr[, "MedGen"] <- ""

for(ln in 1:nrow(dataChr)){
  dataChr[ln,"MedGen"] <- stri_split_fixed(dataChr[ln,"CLNDISDB"], ",") [[1]][1]
}
for(ln in 1:nrow(dataChr)){
  dataChr[ln,"MedGen"] <- stri_split_fixed(dataChr[ln,"MedGen"], "|") [[1]][1]
}
for(ln in 1:nrow(dataChr)){
  dataChr[ln,"MedGen"] <- stri_split_fixed(dataChr[ln,"MedGen"], ":") [[1]][1]
}

var1 <- "[.]"

for(ln in 1:nrow(dataChr)){
  if(grepl(var1, dataChr[ln,"MedGen"])){
    dataChr[ln,"MedGen"] <- "EFO"
  }
}

var2 <- "EFO"

for(ln in nrow(dataChr):1){
  if(grepl(var2, dataChr[ln,"MedGen"])){
    dataChr <- dataChr[-ln, ]
  }
}


temp_df <- dataChr[,c("MedGen", "CLASS")]
temp_df2 <- dataChr[,c("MedGen", "CLASS")]
temp_df <- dummy.data.frame(temp_df)
temp_df <- temp_df[, -5]

dataChr <- cbind(dataChr, temp_df)

unique(dataChr[, "MedGen"])
prop.table(table(dataChr$MedGen,dataChr$CLASS), 1)
table(dataChr$MedGen)
chisq.test(dataChr$MedGen, dataChr$CLASS)


prop.table(table(dataChr$MedGenGene,dataChr$CLASS), 1)
table(dataChr$MedGenGene)
chisq.test(dataChr$MedGenGene, dataChr$CLASS) #pas intéressant

prop.table(table(dataChr$MedGenHuman_Phenotype_Ontology,dataChr$CLASS), 1)
table(dataChr$MedGenHuman_Phenotype_Ontology)
chisq.test(dataChr$MedGenHuman_Phenotype_Ontology, dataChr$CLASS)# très intéressant

prop.table(table(dataChr$MedGenMedGen,dataChr$CLASS), 1)
table(dataChr$MedGenMedGen)
chisq.test(dataChr$MedGenMedGen, dataChr$CLASS) # oui aussi

prop.table(table(dataChr$MedGenMeSH,dataChr$CLASS), 1)
table(dataChr$MedGenMeSH)
chisq.test(dataChr$MedGenMeSH, dataChr$CLASS) # bof



  
######## new col CLNDN #########

dataChr[, "New_CLNDN"] <- ""

for(ln in 1:nrow(dataChr)){
  dataChr[ln,"New_CLNDN"] <- stri_split_fixed(dataChr[ln,"CLNDN"], ";") [[1]][1]
}
for(ln in 1:nrow(dataChr)){
  dataChr[ln,"New_CLNDN"] <- stri_split_fixed(dataChr[ln,"New_CLNDN"], "|") [[1]][1]
}
for(ln in 1:nrow(dataChr)){
  dataChr[ln,"New_CLNDN"] <- stri_split_fixed(dataChr[ln,"New_CLNDN"], ",") [[1]][1]
}
for(ln in 1:nrow(dataChr)){
  dataChr[ln,"New_CLNDN"] <- stri_split_fixed(dataChr[ln,"New_CLNDN"], "_") [[1]][1]
}

unique(dataChr[, "New_CLNDN"])
prop.table(table(dataChr[, "New_CLNDN"],dataChr[, "CLASS"]), 1)
table(dataChr$New_CLNDN)
chisq.test(dataChr$New_CLNDN,dataChr$CLASS)


######## SYMBOL ######

temp_ds <- dataChr
unique(dataChr[, "SYMBOL"]) #2329 LVL
prop.table(table(dataChr[, "SYMBOL"],dataChr[, "CLASS"]), 1)
table(dataChr$SYMBOL)
chisq.test(dataChr$SYMBOL,dataChr$CLASS)

toto <- 0
nbr <-0
for(name in 1:2293){
  if (table(dataChr$SYMBOL)[[name]]>50){
    if (table(dataChr$SYMBOL)[[name]] > toto){
      toto <<- table(dataChr$SYMBOL)[[name]]
    }
    nbr <<- (nbr + 1)
  }
}
toto
nbr

# 212 variables de + de 50 valeurs (sur 2329) le max d'une variable est à 2601
# intéressant de la garder du coup ?

########### études colonnes AF_... et tests #########
# faire étude discriminance sur les features numériques (AF_XX et LoF, CADD)
# 
# tot <- !is.na(dataChr$CLNDISDBINCL)
# dataChr$CLNDISDBINCL[tot]
# 
# 
# is.numeric(typeof(dataChr[,1]))
# length(dataChr)
# 
# num_col = c()
# for (colu in 1:length(dataChr)){
#   if (typeof(dataChr[,colu]) == "double" | typeof(dataChr[,colu]) == "numeric") {
#     print(paste(colu))
#   }
# }
# 
# dataChr[,c(5)] <- as.numeric(dataChr[,c(5)])
# dataChr[,c(6)] <- as.numeric(dataChr[,c(6)])
# dataChr[,c(7)] <- as.numeric(dataChr[,c(7)])
# dataChr[,c(43)] <- as.numeric(dataChr[,c(43)])
# dataChr[,c(44)] <- as.numeric(dataChr[,c(44)])
# dataChr[,c(45)] <- as.numeric(dataChr[,c(45)])
# dataChr[,c(19)] <- as.integer(dataChr[,c(19)])
# 
# data_to_plot <- dataChr[,c(5, 6, 7, 43, 44, 45, 19)]
# colnames(data_to_plot)[7] <- "class"
# data_to_plot[,c(7)] <- as.integer(data_to_plot[,c(7)])
# 
# 
# 
# str(data_to_plot)
# 
# apply(data_to_plot, 2, function(x) sum(is.na(x) || is.null(x)))
# 
# is.na(data_to_plot[5,2])
# 
# # essayer de changer les distances par un row_number ou dans le style
# 
# for (col in 1:length(data_to_plot)){
#   for (ln in nrow(data_to_plot):1){
#     if (is.na(data_to_plot[ln,col]) || is.null(data_to_plot[ln,col]) == TRUE) {
#       data_to_plot <- data_to_plot[-ln,]
#     }
#   }
# }
# 
# 
# ggplot(data = melt(data_to_plot,id.vars = "class", variable.name = "field")) + 
#   geom_boxplot(aes(x = as.factor(class), y=value, color = as.factor(class))) + 
#   facet_wrap(~field,scale = "free") +
#   scale_color_manual(values= c("green", "red")) +
#   theme(legend.position = "none") 




####### plot ######

data_to_plot <- data.table(dataChr[, c("REF","ALT", "CLASS")])
data_to_plot[, count:= 1] 
data_to_plot %<>%
  .[, CLASS := factor(CLASS, levels = c(0,1), labels = c("sane", "ill"))] %>%
  dcast(REF + ALT ~ CLASS,fun.aggregate = sum, value.var = "count")

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1, 2, NA),
                        c(3, 3, 4))
)


p1 <- ggplotGrob(ggplot(data = data_to_plot) +
  geom_bar(aes(x = ALT, y = sane, fill = REF), stat = "identity") +
  facet_wrap(. ~ REF, nrow = 4, strip.position = "right"))

p2 <- ggplotGrob(ggplot(data = data_to_plot) +
  geom_bar(aes(x = ALT, y = ill, fill = REF), stat = "identity") +
  facet_wrap(. ~ REF, nrow = 4, strip.position = "right"))

p <- rbind(p1, p2, size = "first")
p$widths <- unit.pmax(p1$widths, p2$widths)
grid.newpage()
grid.draw(p)
  

data_to_plot <- data.table(dataChr[, c("MedGen","CLASS")])
data_to_plot[, count:= 1] 
data_to_plot %<>% 
  dcast(CLASS ~ MedGen,fun.aggregate = sum, value.var = "count") %>%
  melt(id.vars = "CLASS", variable.name = "Genotype") %>%
  .[, CLASS := as.factor(CLASS)]

ggplot(data = data_to_plot) +
  geom_bar(aes(x = CLASS, y = value, fill = missense_variant), stat = "identity", position = "stack")
ggplot(data = data_to_plot) +
  geom_boxplot(aes(x = CLASS, y = value, colour = CLASS))

data_to_plot <- data.table(dataChr[, c("missense_variant","CLASS")])
data_to_plot[, count:= 1] 
data_to_plot %<>% 
  dcast(CLASS ~ missense_variant,fun.aggregate = sum, value.var = "count") %>%
  melt(id.vars = "CLASS", variable.name = "missense_variant") %>%
  .[, CLASS := as.factor(CLASS)]

ggplot(data = data_to_plot) +
  geom_boxplot(aes(x = CLASS, y = value, colour = CLASS))



data_to_plot <- data.table(dataChr[, c("New_CLNDN","CLASS")])
data_to_plot[, count:= 1] 
data_to_plot %<>% 
  dcast(CLASS ~ New_CLNDN,fun.aggregate = sum, value.var = "count") %>%
  melt(id.vars = "CLASS", variable.name = "New_CLNDN") %>%
  .[, CLASS := as.factor(CLASS)]

ggplot(data = data_to_plot) +
  geom_bar(aes(x = CLASS, y = value, fill = New_CLNDN), stat = "identity", position = "stack")


ggplot(dataChr[1:10000,], aes(x="", y=CLNVC, fill=REF))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)

test <- data.table(dataChr) %>%
  .[, .(REF, CLNVC)] %>%
  dcast(REF ~ CLNVC) %>%
  data.table %>%
  .[, REF:= as.character(REF)]

sapply(test, function(x) if(class(x)!="factor") sum(x))

plot_ly(test[single_nucleotide_variant!=0, ]) %>%
  add_bars(x=~REF, y=~single_nucleotide_variant, marker = list(color = rainbow(4)))

data.table(dataChr) %>%
  .[, .(CLASS, CLNVC)] %>%
  dcast(CLASS ~ CLNVC) %>%
  plot_ly() %>%
  add_bars(x=~CLASS, y=~single_nucleotide_variant, marker = list(color = rainbow(2)))


  plot_ly() %>%
  add_pie(values = ~Duplication, labels = ~REF)
  
  
# plot du ration exon et intron en fonction de la cible
data_to_plot <- dataChr[,c("EXONINTRON", "CLASS")]
  
ggplot(data = melt(data_to_plot,id.vars = "CLASS", variable.name = "field")) + 
   geom_boxplot(aes(x = as.factor(CLASS), y=value, color = as.factor(CLASS))) +
   facet_wrap(~field,scale = "free") +
   scale_color_manual(values= c("green", "red")) +
   theme(legend.position = "none")

#visiblement c'est pas discriminant


data_to_plot <- data.table(dataChr[, c("ORIGIN","CLASS")])
data_to_plot[, count:= 1] 
data_to_plot %<>% 
  dcast(CLASS ~ ORIGIN,fun.aggregate = sum, value.var = "count") %>%
  melt(id.vars = "CLASS", variable.name = "ORIGIN") %>%
  .[, CLASS := as.factor(CLASS)]

ggplot(data = data_to_plot) +
  geom_bar(aes(x = CLASS, y = value, fill = ORIGIN), stat = "identity", position = "stack")



# chrom

data_to_conserve$CHROM <- as.integer(data_to_conserve$CHROM) 
data_to_conserve$CLASS <- as.factor(data_to_conserve$CLASS)

ggplot(data = melt(data_to_conserve[, c("CHROM", "CLASS")],id.vars = "CLASS", variable.name = "field")) + 
  geom_boxplot(aes(x = CLASS, y=value, color = CLASS)) +
  facet_wrap(~field,scale = "free") +
  scale_color_manual(values= c("green", "red")) +
  theme(legend.position = "none")

table(data_to_conserve$CHROM, data_to_conserve$CLASS)
prop.table(table(dataChr[, "CHROM"],dataChr[, "CLASS"]), 1)

###### data to be conserved #######
prop.table(table(dataChr[,val23],dataChr[,"CLASS"]), 1)


str(dataChr)
prop.table(table(dataChr[,70],dataChr$CLASS), 1)
table(dataChr[,70]) 
chisq.test(dataChr[,70],dataChr$CLASS) #  inframe deletion à 0.036 => ça passe mais on va pas la garder pour le moment
# pareil pour start_loss qui est à 0.045
# downstream_gene_variant à 0.007 donc ça passe mais mis de côté pour le moment
table(dataChr$CLNVC)
# data_to_conserve <- dataChr[ , c("CLASS", "missense_variant", "splice_region_variant", "intron_variant", "3_prime_UTR_variant", 
#                                  "frameshift_variant", "stop_gained", "5_prime_UTR_variant", "splice_acceptor_variant", "splice_donor_variant",
#                                  "intergenic_variant", "CHROM", "CLNVC", "New_IMPACT", "nonsense", "ALT", "REF", "ORIGIN", "MedGen",
#                                  "STRAND", "Allele")] # , "ORIGIN" mis de coté pour le moment
# 
unique(dataChr[, "ORIGIN"])


for(ln in 1:nrow(dataChr)){
  if(sapply(dataChr[ln,"ORIGIN"], function(y) is.na(y) | is.null(y)) == TRUE){
    dataChr[ln,"ORIGIN"] <- mean(dataChr$ORIGIN)
  }
}

table(dataChr$ORIGIN)

dataChr$Allele <- as.character(dataChr$Allele)
dataChr$Allele <- factor(dataChr$Allele)
data_to_conserve <- dataChr[ , c("CLASS", "missense_variant", "splice_region_variant", "intron_variant", "3_prime_UTR_variant", 
                                 "stop_gained", "5_prime_UTR_variant", "splice_acceptor_variant", "splice_donor_variant", 
                                 "nonsense", "AlleleA", "AlleleT", "AlleleC", "AlleleG", "STRAND", "New_IMPACT", "REFA", "REFC",
                                 "REFT", "REFG", "MedGenHuman_Phenotype_Ontology", "MedGenMedGen")]


str(data_to_conserve)

sapply(data_to_conserve, function(y) sum(is.null(y))) 


for(ln in 1:nrow(data_to_conserve)){
  if(sapply(data_to_conserve[ln,"ORIGIN"], function(y) is.na(y) | is.null(y)) == TRUE){
    dataChr[ln,"ORIGIN"] <- mean(dataChr$ORIGIN)
  }
}


data_to_conserve$CLASS <- as.factor(data_to_conserve$CLASS)
# data_to_conserve$missense_variant <- as.factor(data_to_conserve$missense_variant)
# data_to_conserve$splice_region_variant <- as.factor(data_to_conserve$splice_region_variant)
# data_to_conserve$intron_variant <- as.factor(data_to_conserve$intron_variant)
# data_to_conserve[,"3_prime_UTR_variant"] <- as.factor(data_to_conserve[,"3_prime_UTR_variant"])
# data_to_conserve$stop_gained <- as.factor(data_to_conserve$stop_gained)
# data_to_conserve[,"5_prime_UTR_variant"] <- as.factor(data_to_conserve[,"5_prime_UTR_variant"])
# data_to_conserve$splice_acceptor_variant <- as.factor(data_to_conserve$splice_acceptor_variant)
# data_to_conserve$splice_donor_variant <- as.factor(data_to_conserve$splice_donor_variant)
# data_to_conserve$nonsense <- as.factor(data_to_conserve$nonsense)
# data_to_conserve$STRAND <- as.factor(data_to_conserve$STRAND)
# data_to_conserve$New_IMPACT <- as.factor(data_to_conserve$New_IMPACT)
# #data_to_conserve$REF <- as.factor(data_to_conserve$REF)
# #data_to_conserve$ALT <- as.factor(data_to_conserve$ALT)
# data_to_conserve$MedGen <- as.factor(data_to_conserve$MedGen)


str(data_to_conserve)
colnames(data_to_conserve)[5] <- "prime_UTR_Variant_3"
colnames(data_to_conserve)[7] <- "prime_UTR_Variant_5"


lapply(c("missense_variant", "splice_region_variant", "intron_variant", "prime_UTR_Variant_3",
         "prime_UTR_Variant_5", "splice_acceptor_variant", "splice_donor_variant", "nonsense"), 
       function(col) {
         column <- as.integer(data.table(data_to_conserve)[, get(col)])
         sum(column)
       })

sapply(data_to_conserve, function(col) sum(col == 1))

####### mise en place du train et test #######
set.seed(525)

eval <- sample(1:nrow(data_to_conserve), 0.35*nrow(data_to_conserve))
test <- data_to_conserve[eval,]
train <- data_to_conserve[-eval,]




temp_dataset <- data_to_conserve[data_to_conserve$CLASS == 1, ][1:8000,] 

eval <- sample(1:nrow(data_to_conserve[data_to_conserve$CLASS == 0, ]), nrow(temp_dataset) + 120)
temp_dataset <- rbind(temp_dataset, data_to_conserve[data_to_conserve$CLASS == 0, ][eval,])
temp_true_dataset <- data_to_conserve[data_to_conserve$CLASS == 0, ][-eval,]
temp_true_dataset <- rbind(temp_true_dataset, data_to_conserve[data_to_conserve$CLASS == 1, ][8001: nrow(data_to_conserve[data_to_conserve$CLASS == 1, ]),]) 


eval <- sample(1:nrow(temp_dataset), 0.35*nrow(temp_dataset))
test <- temp_dataset[eval,]
train <- temp_dataset[-eval,]



########## Random Forest ##########


RF_classifier <- randomForest(CLASS ~ ., data = train, do.trace=20, ntree=100, mtry=9)

predict1<-predict(RF_classifier, type='response', newdata=train)
confusionMatrix(data = predict1, reference = train$CLASS, positive = "0")

predict2<-predict(RF_classifier, type='response', newdata=test)
confusionMatrix(data = predict2, reference = test$CLASS, positive = "0")


mod.rf.tri = randomForest(CLASS~., data = temp_dataset, do.trace=20, ntree=100, mtry=9)


predictR1<-predict(mod.rf.tri, type='response', newdata=temp_dataset)
confusionMatrix(data = predictR1, reference = temp_dataset$CLASS, positive = "0")


predictR2<-predict(mod.rf.tri, type='response', newdata=temp_true_dataset)
confusionMatrix(data = predictR2, reference = temp_true_dataset$CLASS, positive = "0")

#confusionMatrix(data = RF_predicted, reference = temp_true_dataset$class, positive = "0")

rf_roc_curve <- roc(temp_true_dataset$CLASS, as.integer(as.vector(predictR2)))
plot(rf_roc_curve)


########## Logisitque Regression ##########

lr_classifier_t <- glm(CLASS ~ ., family = binomial(), data = train)
lr_predictor_t <- predict(lr_classifier_t, train, type = "response")
confusionMatrix(data = factor(ifelse(lr_predictor_t>0.5, 1,0), levels = c(0, 1)), reference = as.factor(train$CLASS), positive = "0")

lr_predictor_test <- predict(lr_classifier_t, test, type = "response")
confusionMatrix(data = factor(ifelse(lr_predictor_test>0.5, 1,0), levels = c(0, 1)), reference = as.factor(test$CLASS), positive = "0")


lr_classifier <- glm(CLASS ~ ., family = binomial(), data = temp_dataset)
lr_predictor <- predict(lr_classifier, temp_dataset, type = "response")
confusionMatrix(data = factor(ifelse(lr_predictor>0.5, 1,0), levels = c(0, 1)), reference = as.factor(temp_dataset$CLASS), positive = "0")



lr_predicted <- predict(lr_classifier, temp_true_dataset, type = "response")
confusionMatrix(data = factor(ifelse(lr_predicted>0.5, 1,0), levels = c(0, 1)), reference = as.factor(temp_true_dataset$CLASS), positive = "0")

lr_roc_curve <- roc(temp_true_dataset$CLASS  , as.integer(as.vector(factor(ifelse(lr_predicted>0.5, 1,0), levels = c(0, 1)))))
plot(lr_roc_curve)



######### Data Transformation #######

norm_data <- setDT(data_to_conserve)
normalizer <- bestNormalize(norm_data[[2]])

lapply(colnames(norm_data), function(col){
  if(col != "CLASS"){
    norm_data[[col]] <<- predict(normalizer, norm_data[, get(col)])
  }
})

par(mfrow = c(3,7))
colors <- rainbow(22)
i <- 0
for(col in colnames(norm_data)){
  i <- i+1
  if(col !="CLASS"){
    plot(density(norm_data[[col]]), main = col, type = "o", col = colors[i])
  }
}

data_to_plot <- norm_data %>%
  melt(value.name = "value", variable.name = "field", id.vars = "CLASS")

ggplot(data = data_to_plot) +
  geom_boxplot(aes(x = CLASS, y = value, colour = CLASS)) +
  facet_wrap(~field,scale = "free") +
  scale_color_manual(values= c("green", "red")) +
  theme(legend.position = "none")

#######

set.seed(5) 
  
  
eval <- sample(1:nrow(norm_data), 0.35*nrow(norm_data))
test <- norm_data[eval,]
train <- norm_data[-eval,]

sapply(train[, .(CLASS)], function(col) sum(col == 1)/(sum(col ==1)+sum(col ==0)))

#Add weights to the dataset
train[, weight := rep(1, nrow(train))]
train[CLASS == 1, weight := 2]
weights = train[, weight]
train$weight <- NULL


lr_classifier_t <- glm(CLASS ~ ., family = binomial("logit"), weights = weights, data = train)
summary(lr_classifier_t)


CLASS <- as.numeric(as.vector(train[, CLASS]))
plot(lr_classifier_t, which = 1)
plot(predict(lr_classifier_t),residuals(lr_classifier_t),col=c("blue","red")[1+CLASS])
abline(h=0,lty=2,col="grey")

plot(train[, intron_variant],residuals(lr_classifier_t),col=c("blue","red")[1+CLASS])
lines(lowess(train[, intron_variant],residuals(lr_classifier_t)),col="black",lwd=2)
lines(lowess(train[, intron_variant][CLASS==0],residuals(lr_classifier_t)[CLASS==0]),col="blue")
lines(lowess(train[, intron_variant][CLASS==1],residuals(lr_classifier_t)[CLASS==1]),col="red")
abline(h=0,lty=2,col="grey")

lr_predictor_t <- predict(lr_classifier_t, train, type = "response")
confusionMatrix(data = factor(ifelse(lr_predictor_t>0.5, 1,0), levels = c(0, 1)), reference = as.factor(train$CLASS), positive = "0")


lr_predictor_t <- predict(lr_classifier_t, test, type = "response")
confusionMatrix(data = factor(ifelse(lr_predictor_t>0.5, 1,0), levels = c(0, 1)), reference = as.factor(test$CLASS), positive = "0")


RF_classifier <- randomForest(CLASS ~ ., data = train, do.trace=20, ntree=100, mtry=9)

predict1<-predict(RF_classifier, type='response', newdata=train)
confusionMatrix(data = predict1, reference = train$CLASS, positive = "0")

predict2<-predict(RF_classifier, type='response', newdata=test)
confusionMatrix(data = predict2, reference = test$CLASS, positive = "0")
