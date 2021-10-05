setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")
source("OLD_ENAfuncs.R")
###TEST NETWORK BALANCING FUNCTIONS###

##TEST 1: Cypwet network=unbalanced network

Cypwet<-SCOR.convert("Network data/cypwet.dat") 
Z_cs<-Cypwet[[4]]#Import
E_cs<-Cypwet[[5]]#Export
R_cs<-Cypwet[[6]]#Respiration
T_cs<-Cypwet[[7]]#intercompartmental exchange

###INPUT-BASED APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="inp")
myINP<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
INP.correct<-IN.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myINP-(apply(INP.correct,1,sum)-apply(INP.correct,2,sum)))

###OUTPUT-BASED APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="out")
myOUT<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
OUT.correct<-OUT.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myOUT-(apply(OUT.correct,1,sum)-apply(OUT.correct,2,sum)))

###AVG APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="avg")
myAVG<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
AVG.correct<-AVG.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myAVG-(apply(AVG.correct,1,sum)-apply(AVG.correct,2,sum)))

###IO APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="io")
myIO<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
IO.correct<-IO.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myIO-(apply(IO.correct,1,sum)-apply(IO.correct,2,sum)))

###OI APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="oi")
myOI<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
OI.correct<-OI.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myOI-(apply(OI.correct,1,sum)-apply(OI.correct,2,sum)))

###AVG2 APPROACH
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="avg2")
myAVG2<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
AVG2.correct<-AVG2.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
unname(myAVG2-(apply(AVG2.correct,1,sum)-apply(AVG2.correct,2,sum)))

##TEST2: ConeSpring=balanced network
names_cs <- c("Plants", "Bacteria", "Detritivores", "Carnivores", "Detritus")
Z_cs <- c(11184, 0, 0, 0, 635)
T_cs <- matrix(c(0, 0, 0, 0, 8881,
                 0, 0, 75, 0, 1600,
                 0, 0, 0, 370, 200,
                 0, 0, 0, 0, 167,
                 0, 5205, 2309, 0, 0), nrow = 5, byrow = TRUE)
E_cs <- c(300, 255, 0, 0, 860)
R_cs <- c(2003, 3275, 1814, 203, 3109)
names(Z_cs) <- names(E_cs) <- names(R_cs) <- names_cs
colnames(T_cs) <- rownames(T_cs)  <- names_cs
a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="avg2")


##TEST3: Test function on 48 networks (SCOR format) and compare results

#Import networks:
source("test.convert_SCOR.R")

#Create a FOR loop to test network balancing functions on each network and to
#compare results with the ENA function:

for(i in 1:length(foodwebs.corrected)){
  
  print(paste("Network", i, ":", sep=""))
  
  Z_cs<-foodwebs.corrected[[i]][[4]]#Import
  E_cs<-foodwebs.corrected[[i]][[5]]#Export
  R_cs<-foodwebs.corrected[[i]][[6]]#Respiration
  T_cs<-foodwebs.corrected[[i]][[7]]#intercompartmental exchange
  
  ##INPUT-BASED
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="inp")
  if(length(a)!=0){
  myINP<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  INP.correct<-IN.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  INP.comp<-unname(myINP-(apply(INP.correct,1,sum)-apply(INP.correct,2,sum)))
      if(length(which(INP.comp==0))==length(INP.comp)){print("INP-OKE")
               }else{print("INP-PROBLEM!!!!!!!")} }
  
  ##OUTPUT-BASED
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="out")
  if(length(a)!=0){
  myOUT<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  OUT.correct<-OUT.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  OUT.comp<-unname(myOUT-(apply(OUT.correct,1,sum)-apply(OUT.correct,2,sum)))
      if(length(which(OUT.comp==0))==length(OUT.comp)){print("OUT-OKE")
              }else{print("OUT-PROBLEM!!!!!!!")} }
  
  ##AVG
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="avg")
  if(length(a)!=0){
  myAVG<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  AVG.correct<-AVG.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  AVG.comp<-unname(myAVG-(apply(AVG.correct,1,sum)-apply(AVG.correct,2,sum)))
     if(length(which(AVG.comp==0))==length(AVG.comp)){print("AVG-OKE")
             }else{print("AVG-PROBLEM!!!!!!!")} }
  
  ##IO
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="io")
  if(length(a)!=0){
  myIO<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  IO.correct<-IO.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  IO.comp<-unname(myIO-(apply(IO.correct,1,sum)-apply(IO.correct,2,sum)))
    if(length(which(IO.comp==0))==length(IO.comp)){print("IO-OKE")
            }else{print("IO-PROBLEM!!!!!!!")} }
  
  ##OI
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="oi")
  if(length(a)!=0){
  myOI<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  OI.correct<-OI.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  OI.comp<-unname(myOI-(apply(OI.correct,1,sum)-apply(OI.correct,2,sum)))
    if(length(which(OI.comp==0))==length(OI.comp)){print("OI-OKE")
           }else{print("OI-PROBLEM!!!!!!!")} }
  
  ##AVG2
  a<-network.balance(Z_cs, E_cs, R_cs, T_cs, method="avg2")
  if(length(a)!=0){
  myAVG2<-apply(a[[2]],1,sum)-apply(a[[2]],2,sum)
  AVG2.correct<-AVG2.balance(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs)
  AVG2.comp<-unname(myAVG2-(apply(AVG2.correct,1,sum)-apply(AVG2.correct,2,sum)))
    if(length(which(AVG2.comp==0))==length(AVG2.comp)){print("AVG2-OKE")
          }else{print("AVG2-PROBLEM!!!!!!!")} }

}#END OF FOR LOOP
