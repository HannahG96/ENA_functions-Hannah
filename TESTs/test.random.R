setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")

#Create input data for testing the function:
library("data.table")
A1<-fread(file = "Network data/A1.csv", data.table = FALSE)
A1<-data.matrix(A1)
A2<-fread(file = "Network data/A2.csv", data.table = FALSE)
A2<-data.matrix(A2)
D1<-fread(file = "Network data/D1.csv", data.table = FALSE)
D1<-data.matrix(D1)
D2<-fread(file = "Network data/D2.csv", data.table = FALSE)
D2<-data.matrix(D2)
B1<-fread(file = "Network data/B1.csv", data.table = FALSE)
B1<-data.matrix(B1)
B2<-fread(file = "Network data/B2.csv", data.table = FALSE)
B2<-data.matrix(B2)
E1<-fread(file = "Network data/E1.csv", data.table = FALSE)
E1<-data.matrix(E1)
E2<-fread(file = "Network data/E2.csv", data.table = FALSE)
E2<-data.matrix(E2)
C1<-fread(file = "Network data/C1.csv", data.table = FALSE)
C1<-data.matrix(C1)
C2<-fread(file = "Network data/C2.csv", data.table = FALSE)
C2<-data.matrix(C2)
F1<-fread(file = "Network data/F1.csv", data.table = FALSE)
F1<-data.matrix(F1)
F2<-fread(file = "Network data/F2.csv", data.table = FALSE)
F2<-data.matrix(F2)
list.FW<-list(A1,A2,B1,B2,C1,C2,D1,D2,E1,E2,F1,F2)
names(list.FW)<-c("A1","A2","B1","B2","C1","C2","D1","D2","E1","E2","F1","F2")
#-->convert networks in extended transfer matrices:
for(i in 1:length(list.FW)){
  web<-list.FW[[i]]
  Z_fw<-web[,"Z"]
  E_fw<-c(web[,"E"],0,0,0)
  R_fw<-c(web[,"R"],0,0,0)
  zeros<-rep(0,length(Z_fw))
  web<-web[c(1:length(Z_fw)),c(1:length(Z_fw))]
  web<-rbind(web,zeros,zeros,Z_fw)
  web<-cbind(web,E_fw,R_fw,rep(0,length(E_fw)))
  rownames(web)<-colnames(web)<-c("1","2","3","4","5","6","7","8","9","10","11","12",
                                  "13","14","15","16","EXPORT","RESP","IMPORT")
  list.FW[[i]]<-web}

#This is additionally needed for function:
#balance
#constrain
#runs
#design<-c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", "3Hw", "3Hw", "1Hw", "1Hw")
#exclude.extremes
#nl<-3
#attributes
#-->Require certain informations about nods:
#-->Compartment numbers of Primary producers need to be indicated in a vector
#PP <- c(1:5)
## Compartment number of sediment detritus needs to be indicated
#SDet<-14
## Compartment number of water detritus needs to be indicated
#WDet<-15
## Compartment number of dissolved organic carbon needs to be indicated
#DOC<-16
###############################################################################

#Full randomization=only constrain link number
#without extreme values:
TEST1<-random(list.FW,
nl=3,
design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                               "3Hw", "3Hw", "1Hw", "1Hw"),
exclude.extremes=TRUE,
runs=999,
balance="unbalanced",
constrain="link number",
attributes=FALSE,
PP=c(1:5),
SDet=14,
WDet=15,
DOC=16,
trophic.analysis=FALSE)
#with extreme values:
TEST2<-random(list.FW,
              nl=3,
              design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                       "3Hw", "3Hw", "1Hw", "1Hw"),
              exclude.extremes=FALSE,
              runs=999,
              balance="unbalanced",
              constrain="link number",
              attributes=FALSE,
              PP=c(1:5),
              SDet=14,
              WDet=15,
              DOC=16,
              trophic.analysis=FALSE)
#-->RECOMMENDED to use with attributes=FALSE & balance="unbalanced"
#   as setting these TRUE leads to messy values, bugs and errors
   

#Randomization with constrained topology
#without extreme values:
TEST3<-random(list.FW,
              nl=3,
              design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                       "3Hw", "3Hw", "1Hw", "1Hw"),
              exclude.extremes=TRUE,
              runs=999,
              balance="unbalanced",
              constrain="topology",
              attributes=FALSE,
              PP=c(1:5),
              SDet=14,
              WDet=15,
              DOC=16,
              trophic.analysis=FALSE)
#with extreme values:
TEST4<-random(list.FW,
              nl=3,
              design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                       "3Hw", "3Hw", "1Hw", "1Hw"),
              exclude.extremes=FALSE,
              runs=999,
              balance="unbalanced",
              constrain="topology",
              attributes=FALSE,
              PP=c(1:5),
              SDet=14,
              WDet=15,
              DOC=16,
              trophic.analysis=FALSE)
#-->RECOMMENDED to use with attributes=FALSE & balance="unbalanced"
#   as setting these TRUE leads to messy values, bugs and errors

#Randomization with constrained link strengths
#without extreme values:
TEST5<-random(list.FW,
              nl=3,
              design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                       "3Hw", "3Hw", "1Hw", "1Hw"),
              exclude.extremes=TRUE,
              runs=999,
              balance="avg2",
              constrain="link strength",
              attributes=TRUE,
              PP=c(1:5),
              SDet=14,
              WDet=15,
              DOC=16,
              trophic.analysis=FALSE)
#with extreme values:
TEST6<-random(list.FW,
              nl=3,
              design=c("C", "C", "3Hw", "3Hw", "1Hw", "1Hw", "C", "C", 
                       "3Hw", "3Hw", "1Hw", "1Hw"),
              exclude.extremes=FALSE,
              runs=99,
              balance="avg2",
              constrain="link strength",
              attributes=TRUE,
              PP=c(1:5),
              SDet=14,
              WDet=15,
              DOC=16,
              trophic.analysis=FALSE)
#PROBLEMS: FUNCTION DOES NOT SHOW PLOTS (-->in step by step running it works!)
#-->Server problem?