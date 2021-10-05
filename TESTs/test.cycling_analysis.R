setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")
source("OLD_ENAfuncs.R")
##TEST CYCLING ANALYSIS FUNCTION

#function to identify negative values in acyclic interaction matrix:
negresults <- function(myfoodweb){
  
  neg.Res<-as.data.frame(matrix(rep(NA), nrow=0, ncol=3,dimnames=list(NULL, c("Out", "In", "Flow"))))
  neg.res<-which(myfoodweb[[7]]<0)
  if(length(neg.res)==0){neg.Res<-NULL}
  else{
  for(i in 1:length(neg.res)){
    row<-neg.res[i]%%nrow(myfoodweb[[7]])
    if(row==0)row<-nrow(myfoodweb[[7]])
    Out<-rownames(myfoodweb[[7]])[row]
    In<-colnames(myfoodweb[[7]])[ceiling(neg.res[i]/nrow(myfoodweb[[7]]))]
    Flow<-myfoodweb[[7]][neg.res[i]]
    neg.Res<-rbind(neg.Res,c(Out, In, Flow))
    colnames(neg.Res)<-c("Out", "In", "Flow")} }
  
  return(neg.Res)
}

#TEST 1: Cone Spring
names_cs <- c("Plants", "Bacteria", "Detritivores", "Carnivores", "Detritus")
com_cs <- length(names_cs)
liv_cs <- 4
Z_cs <- c(11184, 0, 0, 0, 635)
T_cs <- matrix(c(0, 0, 0, 0, 8881,
                 0, 0, 75, 0, 1600,
                 0, 0, 0, 370, 200,
                 0, 0, 0, 0, 167,
                 0, 5205, 2309, 0, 0), nrow = com_cs, byrow = TRUE)
E_cs <- c(300, 255, 0, 0, 860)
R_cs <- c(2003, 3275, 1814, 203, 3109)
names(Z_cs) <- names(E_cs) <- names(R_cs) <- names_cs
colnames(T_cs) <- rownames(T_cs)  <- names_cs

ConeSpring1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs)#default multiple weak arc selection
negresults(ConeSpring1) #-->no negative results in acyclic matrix
ConeSpring2<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, select="weighted")#weighted multiple weak arc selection
negresults(ConeSpring2) #-->no negative results in acyclic matrix

##COMPARE RESULTS
ConeSpring3<-cyc.analysis(Z_cs,T_cs,E_cs,R_cs)
#Compare cycle number:
ConeSpring1[[1]]
ConeSpring2[[1]]
ConeSpring3[[1]]
#   -->SAME
#Compare cyclelength distribution:
ConeSpring1[[2]]
ConeSpring2[[2]]
ConeSpring3[[2]]
#   -->SAME
#Compare weak arcs:
ConeSpring1[[4]]
ConeSpring2[[4]]
ConeSpring3[[4]]
#   -->SAME
#Compare acyclic matrix:
ConeSpring1[[7]] - ConeSpring3[[8]] 
ConeSpring2[[7]] - ConeSpring3[[8]] 
#   -->SAME
#Compare residual matrix:
ConeSpring1[[8]] - ConeSpring3[[9]] 
ConeSpring2[[8]] - ConeSpring3[[9]] 
#   -->SAME
#Compare cycle distribution:
ConeSpring1[[5]]
ConeSpring2[[5]]
ConeSpring3[[6]]
#   -->SAME
#Compare normalized cycle distribution:
ConeSpring1[[6]]
ConeSpring2[[6]]
ConeSpring3[[7]]
#   -->SAME

##SEARCH FOR CYCLES IN THE ACYCLIC NETWORK
T_cs<-ConeSpring1[[7]]
ConeSpring_acyc1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
ConeSpring_acyc1[[1]]
T_cs<-ConeSpring2[[7]]
ConeSpring_acyc2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
ConeSpring_acyc2[[1]]
T_cs<-ConeSpring3[[8]]
ConeSpring_acyc3<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
ConeSpring_acyc3[[1]]
#   -->OKE

###########################
#TEST2: Crystal River Creek

source("Network data/cryscon.R")
T_cs<-T
R_cs<-R
E_cs<-E
Z_cs<-Z

CrystalRiverCreek1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)#default multiple weak arc selection
negresults(CrystalRiverCreek1) #negative results of acyclic matrix range from 10^-16 to 10^-19
CrystalRiverCreek2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs, select="weighted")#weighted multiple weak arc selection
negresults(CrystalRiverCreek2) #negative results of acyclic matrix range from 10^-16 to 10^-19

##COMPARE RESULTS
CrystalRiverCreek3<-cyc.analysis(Z_cs,T_cs,E_cs,R_cs)
#Compare cycle number:
CrystalRiverCreek1[[1]]
CrystalRiverCreek2[[1]]
CrystalRiverCreek3[[1]]
#   -->SAME
#Compare cyclelength distribution:
CrystalRiverCreek1[[2]]
CrystalRiverCreek2[[2]]
CrystalRiverCreek3[[2]]
#   -->SAME
#Compare weak arcs:
nrow(CrystalRiverCreek1[[4]])
CrystalRiverCreek1[[4]][, "Flow.removed"]
CrystalRiverCreek1[[4]][, "Flow"]
nrow(CrystalRiverCreek2[[4]])
CrystalRiverCreek2[[4]][, "Flow.removed"]
CrystalRiverCreek2[[4]][, "Flow"]
length(CrystalRiverCreek3[[4]])
CrystalRiverCreek3[[4]]
#   -->SAME
#Compare acyclic matrix:
CrystalRiverCreek1[[7]] - CrystalRiverCreek3[[8]] 
CrystalRiverCreek2[[7]] - CrystalRiverCreek3[[8]] 
#   -->largest differences are in the range of 10^-3 (To be worried?...)
#Compare residual matrix:
CrystalRiverCreek1[[8]] - CrystalRiverCreek3[[9]] 
CrystalRiverCreek2[[8]] - CrystalRiverCreek3[[9]] 
#   -->largest differences are in the range of 10^-3 (To be worried?...)
#Compare cycle distribution:
CrystalRiverCreek1[[5]]
CrystalRiverCreek2[[5]]
CrystalRiverCreek3[[6]]
#   -->SAME
#Compare normalized cycle distribution:
CrystalRiverCreek1[[6]]
CrystalRiverCreek2[[6]]
CrystalRiverCreek3[[7]]
#   -->SAME

##SEARCH FOR CYCLES IN THE ACYCLIC NETWORK
T_cs<-CrystalRiverCreek1[[7]]
CrystalRiverCreek_acyc1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
CrystalRiverCreek_acyc1[[1]]
T_cs<-CrystalRiverCreek2[[7]]
CrystalRiverCreek_acyc2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
CrystalRiverCreek_acyc2[[1]]
T_cs<-CrystalRiverCreek3[[8]]
CrystalRiverCreek_acyc3<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
CrystalRiverCreek_acyc3[[1]]
#   -->OKE


#TEST3: Chesa
chesa <- as.matrix(read.table(file = "Network data/chesa.txt", header = TRUE, sep = "\t"))
rownames(chesa) <- colnames(chesa)
T_cs<-chesa[c(2:37),c(2:37)]
R_cs<-chesa[c(2:37), 39]
E_cs<-chesa[c(2:37), 38]
Z_cs<-chesa[c(2:37), 1]
names(R_cs)<-names(E_cs)<-names(Z_cs)<-colnames(T_cs)

Chesa1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)#default multiple weak arc selection
negresults(Chesa1)#negative results of acyclic matrix range from 10^-13 to 10^-16
Chesa2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs, select="weighted")#weighted multiple weak arc selection
negresults(Chesa2)#negative results of acyclic matrix range from 10^-13 to 10^-16

##COMPARE RESULTS
Chesa3<-cyc.analysis(Z_cs,T_cs,E_cs,R_cs)
#Compare cycle number:
Chesa1[[1]]
Chesa2[[1]]
Chesa3[[1]]
#   -->SAME
#Compare cyclelength distribution:
Chesa1[[2]]
Chesa2[[2]]
Chesa3[[2]]
#   -->SAME
#Compare weak arcs:
nrow(Chesa1[[4]])
Chesa1[[4]][, "Flow.removed"]
Chesa1[[4]][, "Flow"]
nrow(Chesa2[[4]])
Chesa2[[4]][, "Flow.removed"]
Chesa2[[4]][, "Flow"]
length(Chesa3[[4]])
Chesa3[[4]]
sum(Chesa1[[4]][, "Flow"]) - sum(Chesa3[[4]])
sum(Chesa2[[4]][, "Flow"]) - sum(Chesa3[[4]])
#   -->SAME
#Compare acyclic matrix:
Chesa_acyc1<-Chesa1[[7]] - Chesa3[[8]]
      #Largest differences:
      #difference zooplankton-->sea nettle=0.392
      #difference ctenophore-->sea nettle=-0.392
      #difference zooplankton-->ctenophore=-6.92
      #difference ctenophore-->ctenophore=6.52
      #difference ciliates-->suspended POC=6.92
      #difference suspended POC-->zooplankton=6.92
Chesa_acyc2<-Chesa2[[7]] - Chesa3[[8]]
      #Largest differences:
      #difference suspended POC-->zooplankton=6.92
#   -->SOME LARGE DIFFERENCES
#Compare residual matrix:
Chesa_res1<-Chesa1[[8]] - Chesa3[[9]] 
Chesa_res2<-Chesa2[[8]] - Chesa3[[9]] 
#   -->SAME LARGE DIFFERENCES AS FOR ACYCLIC MATRIX (of course, just sign of differences is reversed)
#Compare cycle distribution:
Chesa1[[5]]
Chesa2[[5]]
Chesa3[[6]]
#   -->Similar values BUT NOT the same!
#Compare normalized cycle distribution:
Chesa1[[6]]
Chesa2[[6]]
Chesa3[[7]]
#   -->Similar values BUT NOT the same!

##SEARCH FOR CYCLES IN THE ACYCLIC NETWORK
T_cs<-Chesa1[[7]]
Chesa_acyc1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
Chesa_acyc1[[1]]
T_cs<-Chesa2[[7]]
Chesa_acyc2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
Chesa_acyc2[[1]]
T_cs<-Chesa3[[8]]
Chesa_acyc3<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)
Chesa_acyc3[[1]]
#   -->OKE


#TEST4: St Marks
stmarks <- as.matrix(read.table(file = "Network data/st_marks.txt", header = TRUE, sep = "\t"))
rownames(stmarks) <- colnames(stmarks)
T_cs<-stmarks[c(2:52),c(2:52)]
R_cs<-stmarks[c(2:52), 54]
E_cs<-stmarks[c(2:52), 53]
Z_cs<-stmarks[c(2:52), 1]
names(R_cs)<-names(E_cs)<-names(Z_cs)<-colnames(T_cs)

Stmarks1<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs)#default multiple weak arc selection
negresults(Stmarks1)#negative results of acyclic matrix range from 10^-16 to 10^-19
Stmarks2<-cycling.analysis(Z_cs, T_cs, E_cs, R_cs, select="weighted")#weighted multiple weak arc selection
negresults(Stmarks2)#negative results of acyclic matrix range from 10^-16 to 10^-19

##COMPARE RESULTS
Stmarks3<-cyc.analysis(Z_cs,T_cs,E_cs,R_cs)
#Compare cycle number:
Stmarks1[[1]]
Stmarks2[[1]]
Stmarks3[[1]]
#   -->SAME
#Compare cyclelength distribution:
Stmarks1[[2]]
Stmarks2[[2]]
Stmarks3[[2]]
#   -->SAME
#Compare weak arcs:
nrow(Stmarks1[[4]])
Stmarks1[[4]][, "Flow.removed"]
Stmarks1[[4]][, "Flow"]
nrow(Stmarks2[[4]])
Stmarks2[[4]][, "Flow.removed"]
Stmarks2[[4]][, "Flow"]
length(Stmarks3[[4]])
Stmarks3[[4]]
sum(Stmarks1[[4]][, "Flow"]) - sum(Stmarks3[[4]])
sum(Stmarks2[[4]][, "Flow"]) - sum(Stmarks3[[4]])
#   -->SAME
#Compare acyclic matrix:
Stmarks_acyc1<-Stmarks1[[7]] - Stmarks3[[8]]
      #Largest differences in the range of 10^-2
Stmarks_acyc2<-Stmarks2[[7]] - Stmarks3[[8]]
      #Largest differences in the range of 10^-2 (Seems there are LESS differences...)
#Compare residual matrix:
Stmarks_res1<-Stmarks1[[8]] - Stmarks3[[9]] 
Stmarks_res2<-Stmarks2[[8]] - Stmarks3[[9]] 
#   -->SAME LARGE DIFFERENCES AS FOR ACYCLIC MATRIX (of course, just sign of differences is reversed)
#Compare cycle distribution:
Stmarks1[[5]]
Stmarks2[[5]]
Stmarks3[[6]]
#   -->Similar values (BUT NOT the same)
#Compare normalized cycle distribution:
Stmarks1[[6]]
Stmarks2[[6]]
Stmarks3[[7]]
#   -->Similar values (BUT NOT exactly the same)

##SEARCH FOR CYCLES IN THE ACYCLIC NETWORK
T_cs<-Stmarks1[[7]]
Stmarks_acyc1<-test.function(Z_cs, T_cs, E_cs, R_cs)
Stmarks_acyc1[[1]]
T_cs<-Stmarks2[[7]]
Stmarks_acyc2<-test.function(Z_cs, T_cs, E_cs, R_cs)
Stmarks_acyc2[[1]]
T_cs<-Stmarks3[[8]]
Stmarks_acyc3<-test.function(Z_cs, T_cs, E_cs, R_cs)
Stmarks_acyc3[[1]]
#   -->OKE

##################
#Test function on a very large foodweb
#-->number of detected cycles is limited to nmax=10 000 (default)
Cypwet<-SCOR.convert("Network data/cypwet.dat")
Z_cs<-Cypwet[[3]]
E_cs<-Cypwet[[4]]
R_cs<-Cypwet[[5]]
T_cs<-Cypwet[[6]]

Cypwet_cs<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs)#takes ~4-5 min 
Cypwet_cs[[1]]
negresults(Cypwet_cs)#44 items in the range of10^-19 to 10^-22

Cypwet_cs.w<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, select="weighted")#somewhat longer 
Cypwet_cs.w[[1]]
negresults(Cypwet_cs.w)#44 items in the range of10^-19 to 10^-22
