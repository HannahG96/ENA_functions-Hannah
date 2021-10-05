setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")
##########################################################################
############TEST NETWORK INDICES FUNCTIONS###############
#+COMPARE OUTPUTS OF partial.matrix=TRUE FALSE FUNCTIONS (and optimized indices)

##Import network:
Cypwet<-SCOR.convert("Network data/cypwet.dat") 
Z_cs<-Cypwet[[4]]#Import
E_cs<-Cypwet[[5]]#Export
R_cs<-Cypwet[[6]]#Respiration
T_cs<-Cypwet[[7]]#intercompartmental exchange
biomass<-Cypwet[[3]]#stock sizes (only necessary for biomass-incl. ascendency)
nl<-Cypwet[[2]][1]-Cypwet[[2]][2]#number of non-living compartments (only necessary for foodweb connectance)

##Shannon Diversity of flows
#partial.matrix=TRUE
shannon.flow(Z_cs, E_cs, R_cs, T_cs, type="whole")[[1]]
shannon.flow(T_cs=T_cs, type="intercompartmental")[[1]]
shannon.flow(T_cs=T_cs, nl=nl, type="foodweb")[[1]]
#partial.matrix=FALSE
shannon.flow(Z_cs, E_cs, R_cs, T_cs, type="whole",partial.matrix=FALSE)
shannon.flow(T_cs=T_cs, type="intercompartmental",partial.matrix=FALSE)
shannon.flow(T_cs=T_cs, nl=nl, type="foodweb",partial.matrix=FALSE)

##Average Mutual Information
#partial.matrix=TRUE
Ami(Z_cs, E_cs, R_cs, T_cs)[[1]]
#partial.matrix=FALSE
Ami(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)

##Development capacity
#partial.matrix=TRUE
dC(Z_cs, E_cs, R_cs, T_cs)[[1]]
#partial.matrix=FALSE
dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)

##Ascendency
#partial.matrix=TRUE
Asc(Z_cs, E_cs, R_cs, T_cs)[[1]]
#partial.matrix=FALSE
Asc(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)

##Overhead
#partial.matrix=TRUE
Overhead(Z_cs, E_cs, R_cs, T_cs)[[1]]
#partial.matrix=FALSE
Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)

##Internal Development Capacity
#internal.dC(Z_cs, E_cs, R_cs, T_cs)#slow
Internal.dC(Z_cs, E_cs, R_cs, T_cs)#fast

##Internal Ascendency
#internal.Asc(Z_cs, E_cs, R_cs, T_cs)#slow
Internal.Asc(Z_cs, E_cs, R_cs, T_cs)#fast

##Biomass-inclusive Ascendency
#partial.matrix=TRUE
Asc.biomass(Z_cs, E_cs, R_cs, T_cs, biomass)[[1]]
#partial.matrix=FALSE
Asc.biomass(Z_cs, E_cs, R_cs, T_cs, biomass, partial.matrix=FALSE)

###########################################################################
#COMPARE OUTPUTS OF OPTIMIZED AND ENA INDICES
############################################################################
#TEST 1: test functions on Cypwet food web

##Import network:
Cypwet<-SCOR.convert("Network data/cypwet.dat") 
Z_cs<-Cypwet[[4]]#Import
E_cs<-Cypwet[[5]]#Export
R_cs<-Cypwet[[6]]#Respiration
T_cs<-Cypwet[[7]]#intercompartmental exchange
biomass<-Cypwet[[3]]#stock sizes (only necessary for biomass-incl. ascendency)
nl<-Cypwet[[2]][1]-Cypwet[[2]][2]#number of non-living compartments (only necessary for foodweb connectance)

##Test indices functions:

###1.Total System Throughput (TST)
Tst(Z_cs, E_cs, R_cs, T_cs)#my function
TST(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)#correct

###1b.Shannon큦 Diversity of flows (H, derived from Shannon큦 information measure)
#shannon.flow(Z_cs, E_cs, R_cs, T_cs, type="whole")[[1]]
#shannon.flow(T_cs=T_cs, type="intercompartmental")[[1]]
#shannon.flow(T_cs=T_cs, nl=nl, type="foodweb")[[1]]

###2.Average Mutual Information (AMI)
Ami(Z_cs, E_cs, R_cs, T_cs)[[1]] #my function
AMI(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###3.Maximal Average Mutual Information
Ami.max(Z_cs, E_cs, R_cs, T_cs) #my function
AMI.max(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###4.Development Capacity (DC)
dC(Z_cs, E_cs, R_cs, T_cs)[[1]] #my function
DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###5.Internal development capacity (DCi)
Internal.dC(Z_cs, E_cs, R_cs, T_cs) #my function
internal.DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###6.Network Ascendency (A)
Asc(Z_cs, E_cs, R_cs, T_cs)[[1]] #my function
ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###7.Internal Ascendency (Ai)
Internal.Asc(Z_cs, E_cs, R_cs, T_cs) #my function
internal.ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###8.Overhead (OV)
Overhead(Z_cs, E_cs, R_cs, T_cs)[[1]] #my function
OVERHEAD(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct
REDUNDANCY(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.imports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.exports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.resp(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)

###9.Biomass-inclusive Ascendency (Ab)
Asc.biomass(Z_cs, E_cs, R_cs, T_cs, biomass)[[1]] #my function
ASC.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, bi=biomass) #correct

###10.Overall connectance (Ulanowicz & Wolff, 1991)
#-->A BIT SLOW!!
connectance(Z_cs, E_cs, R_cs, T_cs, type="whole") #my function
O.C(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###11.Intercompartmental Connectance
#-->A BIT SLOW!!
connectance(T_cs=T_cs, type="intercompartmental") #my function
I.C(inter=T_cs) #correct

###12.Foodweb Connectance
#-->A BIT SLOW!!
connectance(T_cs=T_cs,nl=nl,type="foodweb") #my function
FW.C(inter=T_cs, nl=nl) #correct

###13.Weighted Connectance (Zorach and Ulanowicz, 2003)
#-->A BIT SLOW!!
eff.C(Z_cs, E_cs, R_cs, T_cs) #my function
effconn.m(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###14.Link Density (Bersier et al., 2002)
LDq(Z_cs, E_cs, R_cs, T_cs) #my function
#-->A BIT SLOW!!
LD.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###15.AMI contribution per nod
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")#whole system
Ami.relative(T_cs=T_cs, type="intercompartmental")#only exchanges between compartments
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)#only exchanges between living compartments

#check outputs:
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)["AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)["AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]



#########################################################################
##TEST 2: Test function (fast version) on 48 networks (SCOR format) and compare results 
#+test speed!

#Import networks:
source("test.convert_SCOR.R")

#Create a FOR loop to test the indices functions on each network and to
#compare results with the ENA function:
#Numbering-->(i1)TST; (i2)AMI; (i3)max.AMI; (i4)C; (i5)internal C; (i6)A;
#            (i7)internal A; (i8)OV; (i9)redundancy; (i10)import OV;
#           (i11)export OV; (i12)respiration OV; (i13)biomass-incl. A;
#           (i14)overall Connectance; (i15)intercompartmental Connectance;
#           (i16)foodweb Connectance; (i17)Weighted Connectance;
#           (i18)Link density; (i19)ami contribution per nod of inflows (whole system);
#           (i20)ami contribution per nod of inflows (only living compartments);
#           (i21)ami contribution per nod of outflows (whole system);
#           (i22)ami contribution per nod of outflows (only living compartments)

for(i in 1:length(foodwebs.corrected)){
  
  print(paste("Network", i, ":", sep=""))
  
  #Import network i
  Z_cs<-foodwebs.corrected[[i]][[4]]#Import
  E_cs<-foodwebs.corrected[[i]][[5]]#Export
  R_cs<-foodwebs.corrected[[i]][[6]]#Respiration
  T_cs<-foodwebs.corrected[[i]][[7]]#intercompartmental exchange
  biomass<-foodwebs.corrected[[i]][[3]]#stock sizes 
  nl<-foodwebs.corrected[[i]][[2]][1]-foodwebs.corrected[[i]][[2]][2]#number of non-living compartments 
  
  #Compute indices and substract it from ENA-calculated indices
  #If value is NOT zero, print (PROBLEM)
  i1<-TST(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Tst(Z_cs, E_cs, R_cs, T_cs)
  i2<-AMI(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Ami(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)
  i3<-AMI.max(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Ami.max(Z_cs, E_cs, R_cs, T_cs) 
  i4<-DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)
  i5<-internal.DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Internal.dC(Z_cs, E_cs, R_cs, T_cs) #my function
  i6<-ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Asc(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)
  i6<-as.vector(round(i6, digits=7))
  if(length(which(i6==0))!=2)i6<-1 else i6<-0
  i7<-internal.ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)[1,]-Internal.Asc(Z_cs, E_cs, R_cs, T_cs)
  i7<-as.vector(round(i7, digits=7))
  if(length(which(i7==0))!=2)i7<-1 else i7<-0
  i8<-OVERHEAD(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)[1,]
  i8<-as.vector(round(i8, digits=7))
  if(length(which(i8==0))!=2)i8<-1 else i8<-0
  i9<-REDUNDANCY(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)[2,]
  i9<-as.vector(round(i9, digits=7))
  if(length(which(i9==0))!=2)i9<-1 else i9<-0  
  i10<-OVERHEAD.imports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)[3,]
  i10<-as.vector(round(i10, digits=7))
  if(length(which(i10==0))!=2)i10<-1 else i10<-0
  i11<-OVERHEAD.exports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)[4,]
  i11<-as.vector(round(i11, digits=7))
  if(length(which(i11==0))!=2)i11<-1 else i11<-0
  i12<-OVERHEAD.resp(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)[5,]
  i12<-as.vector(round(i12, digits=7))
  if(length(which(i12==0))!=2)i12<-1 else i12<-0
  i13<-ASC.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, bi=biomass)-Asc.biomass(Z_cs, E_cs, R_cs, T_cs, biomass, partial.matrix=FALSE)
  i14<-O.C(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-connectance(Z_cs, E_cs, R_cs, T_cs, type="whole") 
  i15<-I.C(inter=T_cs)-connectance(T_cs=T_cs, type="intercompartmental")
  if(nrow(T_cs)-nl<=1){print("i16: only 1 or less living compartment, return NA");i16<-NA}#include exception (which is not accounted for in ENA function but only in mine)
  else{i16<-FW.C(inter=T_cs, nl=nl)-connectance(T_cs=T_cs,nl=nl, type="foodweb")} 
  i17<-effconn.m(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-eff.C(Z_cs, E_cs, R_cs, T_cs)
  i18<-LD.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)-LDq(Z_cs, E_cs, R_cs, T_cs)
  if(nrow(T_cs)-nl<=1){
    print("i19-22: return NA")
    i19<-i20<-i21<-i22<-NA}#include exception (which is not accounted for in ENA function but only in mine)
  else{
    i19<-Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
    i19<-as.vector(round(i19,digits=7))
    if(length(which(i19!=0))>0)i19<-1 else i19<-0
    i20<-Ami.relative(T_cs=T_cs, nl=nl, type="foodweb")[,"AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]
    i20<-as.vector(round(i20,digits=7))
    if(length(which(i20!=0))>0)i20<-1 else i20<-0
    i21<-Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
    i21<-as.vector(round(i21,digits=7))
    if(length(which(i21!=0))>0)i21<-1 else i21<-0
    i22<-Ami.relative(T_cs=T_cs, nl=nl, type="foodweb")[,"AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]
    i22<-as.vector(round(i22,digits=7))
    if(length(which(i22!=0))>0)i22<-1 else i22<-0}
  
  RES<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22)
  RES<-unname(round(RES, digits=7))
  
  probs<-which(RES!=0)
  if(length(probs)==0)print("GOOD")else print(paste("PROBLEM in indice(s)",probs,"!!!",sep=" "))
  
}#END OF FOR LOOP

#Separate Check of Shannon큦 Diversity of Flows Index and relative ami of intercompartmental transfers:
#-->no function in ENA-script available to compare it
for(i in 1:length(foodwebs.corrected)){
  
  print(paste("Network", i, ":", sep=""))
  
  #Import network i
  Z_cs<-foodwebs.corrected[[i]][[4]]#Import
  E_cs<-foodwebs.corrected[[i]][[5]]#Export
  R_cs<-foodwebs.corrected[[i]][[6]]#Respiration
  T_cs<-foodwebs.corrected[[i]][[7]]#intercompartmental exchange
  biomass<-foodwebs.corrected[[i]][[3]]#stock sizes 
  nl<-foodwebs.corrected[[i]][[2]][1]-foodwebs.corrected[[i]][[2]][2]#number of non-living compartments 
  
  #Compute H:
  shannon1<-shannon.flow(Z_cs, E_cs, R_cs, T_cs, type="whole", partial.matrix=FALSE)
  shannon2<-shannon.flow(T_cs=T_cs, type="intercompartmental", partial.matrix=FALSE)
  shannon3<-shannon.flow(T_cs=T_cs, nl=nl, type="foodweb", partial.matrix=FALSE)
  
  Ami.relative(T_cs=T_cs, type="intercompartmental")
  
  print("OK")
}#END OF FOR LOOP


#########################################################################################
##TEST 3: Exceptional case=no intercompartmental exchanges
Z_cs<-Cypwet[[4]]#Import
E_cs<-Cypwet[[5]]#Export
R_cs<-Cypwet[[6]]#Respiration
T_cs<-matrix(0, nrow=length(Z_cs), ncol=length(Z_cs), 
             dimnames=list(names(Z_cs), names(Z_cs)))#Intercompartmental exchanges
biomass<-Cypwet[[3]]#stock sizes (only necessary for biomass-incl. ascendency)
nl<-Cypwet[[2]][1]-Cypwet[[2]][2]#number of non-living compartments (only necessary for foodweb connectance)


###1.Total System Throughput (TST)
Tst(Z_cs, E_cs, R_cs, T_cs)#my function
TST(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)#correct

###1b.Shannon큦 Diversity of Flows (H)
shannon.flow(Z_cs, E_cs, R_cs, T_cs, type="whole", partial.matrix=FALSE)
shannon.flow(T_cs=T_cs, type="intercompartmental", partial.matrix=FALSE)
shannon.flow(T_cs=T_cs, nl=nl, type="foodweb", partial.matrix=FALSE)

###2.Average Mutual Information (AMI)
Ami(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE) #my function
AMI(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###3.Maximal Average Mutual Information
Ami.max(Z_cs, E_cs, R_cs, T_cs) #my function
AMI.max(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###4.Development Capacity (DC)
dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)#my function
DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###5.Internal development capacity (DCi)
Internal.dC(Z_cs, E_cs, R_cs, T_cs) #my function 
internal.DC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###6.Network Ascendency (A)
Asc(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE) #my function
ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###7.Internal Ascendency (Ai)
Internal.Asc(Z_cs, E_cs, R_cs, T_cs) #my function 
internal.ASC(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###8.Overhead (OV)
Overhead(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE) #my function
OVERHEAD(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct
REDUNDANCY(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.imports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.exports(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)
OVERHEAD.resp(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs)

###9.Biomass-inclusive Ascendency (Ab)
Asc.biomass(Z_cs, E_cs, R_cs, T_cs, biomass, partial.matrix=FALSE) #my function 
ASC.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, bi=biomass) #correct

###10.Overall connectance (Ulanowicz & Wolff, 1991)
connectance(Z_cs, E_cs, R_cs, T_cs, type="whole") #my function 
O.C(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###11.Intercompartmental Connectance
connectance(T_cs=T_cs, type="intercompartmental") #my function
I.C(inter=T_cs) #correct

###12.Foodweb Connectance
connectance(T_cs=T_cs, nl=nl, type="foodweb") #my function
FW.C(inter=T_cs, nl=nl) #correct

###13.Weighted Connectance (Zorach and Ulanowicz, 2003)
eff.C(Z_cs, E_cs, R_cs, T_cs) #my function 
effconn.m(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###14.Link Density (Bersier et al., 2002)
LDq(Z_cs, E_cs, R_cs, T_cs) #my function 
LD.B(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs) #correct

###15.AMI contribution per nod
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")#whole system
Ami.relative(T_cs=T_cs, type="intercompartmental")#only exchanges between compartments
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)#only exchanges between living compartments

#check outputs:
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)["AMI.inflows"]-AMI.in(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]
Ami.relative(Z_cs, E_cs, R_cs, T_cs, type="whole")[,"AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[1]]
Ami.relative(T_cs=T_cs, type="foodweb", nl=nl)["AMI.outflows"]-AMI.out(inp=Z_cs, outt=E_cs, diss=R_cs, inter=T_cs, nl=nl)[[2]]



