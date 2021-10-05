#################################RANDOMIZATION ALGORITHM#################################################

#PROBLEMS:
#(1)if there are negative indice values log(i)=NA, these values are excluded for
#   calculation of quantiles
#(2)Randomization (especially if only link number is constrained) can lead to 
#  problems in balancing procedure (det(X)=0):it is recommended to keep the 
#  networks unbalanced if only link number is set as constrain.
#(3)FCI requires calculation of leontief matrix. Can lead to errors if det(X)=0

#Aim of the function:

#(1)Input is a series of networks originating from an experimental setting:
#-->each network defined as extended transfer matrix
#-->same species community but under hold different treatments (with replicates)
#   which replicates belong to which treatment must be indicated in function as
#   a vector (design=c(); "C" as control is necessary)
#   Experiment should be balanced (?????)
#-->matrix is expected to be filled by columns (byrow=FALSE)
#-->number of non-living compartments (nl) must be additionally specified

###########OPTION: exclude.extremes=TRUE/FALSE
###Exclusion of extreme values
#(1b)for each treatment identify extremes of every spp-spp transfer
#-->min and max flow value
#+pooled networks: identify extremes+calculate mean and sd
#(1c)For each spp-spp interaction we discard the 2 flows with the most extreme 
# values resulting in a set of (nb experimental foodwebs - 2) foodwebs to sample from.
###ELSE: use all foodweb flows to create random networks
#-->results in the total nb of experimental foodwebs to sample from.

#(2)For n runs:
#we simulate foodwebs by randomly shuffling intercompartmental 
#flows among treatments
#-->exclude.extremes=TRUE: (nb experimental foodwebs - 2) random foodwebs are created
#-->exclude.extremes=FALSE: (nb experimental foodwebs) random foodwebs are created

#(i)random topology and link strengths but same number of non-null links+
#   at least 1 non-null import + respiration non-null for all compartments; 
#(ii)constant topology, random link strengths;
#(iii)constant topology, random link strengths constrained by realized spp-spp exchanges
#-->number of runs must be specified (runs=)
#-->method must be specified in function as string expression
#   constrain=c("link number", "topology", "link strength")
#-->Use of Rbase function: sample()

#(3)Next (i)a balancing procedure is applied on the simulated networks or (ii)networks are 
#   kept unbalanced
#-->balancing method must be specified in function as string expression
#   balance=c("unbalanced","inp","out","avg","io","oi","avg2")

#(4)We calculate indices for simulated and true networks.
#-->To do so a treatment is randomly attributed to each random foodweb
#   Thereby we obtain a distribution of values for each indice: 
#-->enables to calculate log ratios, to evaluate significance of true log-ratios
#   and to attribute confidence intervals to each true log ratio
##################OPTION: attributes=TRUE/FALSE
#if set TRUE 18 attribute indices are computed in addition to the 
#information and connectance indices
#this requires knowledge of the primary producer compartment numbers (PP),
#the water & sediment detritus compartment numbers (WDet, SDet)
#-->PROBLEM: computation of grazing chain efficiency, Lindeman spine efficiency
#            and system omnivory index is VERY SLOW! 
#  SOLUTION: only calculate them if asked for-->trophic.analysis=TRUE/FALSE)
#-->PROBLEM: negative APL values lead to NA values of log ratios (This is problematic
#            for calculating quantiles)
#   SOLUTION: instead of APL defined by Wulff et al. (1989) -> APL = [T - Z]/Z we 
#             use the definition of Finn (1976) -> APL = TST/Z

#(5)We visualize results in log-ratio plots
# -->DEPENDENCY ON ggplot2-package

###################################################################################

###################################################################################
random<-function(list.FW,nl,design,exclude.extremes=TRUE,runs,balance="unbalanced",
                 constrain="link number",attributes=TRUE,PP,SDet=NA,WDet=NA,DOC=NA,
                 trophic.analysis=TRUE){

library("ggplot2")

###CREATE FOR FUNCTION RELEVANT OBJECTS#######################################

#number of compartments of the experimental foodwebs:
comps<-nrow(list.FW[[1]])-3 

#number of living compartments of the experimental foodwebs:
living<-comps-nl

#compartment names:
names.comps<-rownames(list.FW[[1]])[c(1:comps)]

#topology of experimental foodwebs:
topology<-which(list.FW[[1]]!=0)

#Calculated information/connectance indices:
#indices<-c("TST","Development Capacity (DC)","Ascendency (A)", "A/DC",
#"Average Mutual Information","Overhead on Imports","Overhead on Exports",
#"Dissipative Overhead","Redundancy","Total Overhead (OV)","OV/TST",
#"Internal Capacity","Internal Ascendency","Overall Connectance",
#"Intercompartmental Connectance","Food web Connectance") #names of all computed indices

if(attributes==FALSE){
indices<-c("TST", "DC", "H", "A", "A(%)", "AMI",									
           "OI", "OE", "OD", "R", "O", "OI(%)", "OE(%)", "OD(%)", "R(%)", 
           "O(%)", "Hc", "int_Hc","IC", "IA", "IR","OC", "ICC", "FWC") } #24

#Calculated attribute indices:
#c("Finn Cycling Index","Effective amount of carbon recycled","Herbivory",
#"Detrivory","Circulation within detrital pool","Detrivory/Herbivory","Export
#from sediment detritus","Export from water detritus","Export from water dissolved 
#organic carbon","Respiration:NPP","Grazing chain trophic efficiency","Lindeman spine trophic 
#efficiency", "Average Path Length","Herbivory/NPP","System Omnivory Index","Total
#gross Primary Productivity","Total export from primary producers","Total 
#respiration from living compartments")
if(attributes==TRUE){
  indices<- c("TST", "DC", "H", "A", "A(%)", "AMI",									
              "OI", "OE", "OD", "R", "O", "OI(%)", "OE(%)", "OD(%)", "R(%)", 
              "O(%)", "Hc", "int_Hc","IC", "IA", "IR","OC", "ICC", "FWC",
              "FCI","C_recycled","Herb","Det","Det_circulation","Det/Herb","E_SDet",
              "E_WDet","E_DOC","R/NPP","GrazingChain","LindemanSpine","APL","Herb/NPP",
              "SOI","GPP","E_PP","R_living")} #42


#########################################################################################
###COMPUTE TRUE INDICE VALUES OF EXPERIMENTAL FOOD WEBS

#matrix to store true indices values:
IND.TRUE<-matrix(NA, nrow=length(indices), ncol=length(list.FW),
                 dimnames=list(indices,names(list.FW)))

for(i in 1:length(list.FW)){
  
  #Check if matrix i is filled by columns:
  #-->if filled by rows print an error message 
  #-->TO DO!!?
  
  #Obtain in-, output, respiration vectors+matrix of intercompartmental transfers:
  Z_f1<-list.FW[[i]]["IMPORT",c(1:comps)]
  E_f1<-list.FW[[i]][c(1:comps),"EXPORT"]
  R_f1<-list.FW[[i]][c(1:comps),"RESP"]
  T_f1<-list.FW[[i]][c(1:comps),c(1:comps)]  
  
  #Balance experimental networks (i) or not (ii)
  if(balance!="unbalanced"){
    list.FW[[i]]<-network.balance(Z_cs=Z_f1,E_cs=E_f1,R_cs=R_f1,T_cs=T_f1,
                                  method=balance)[[2]]
    #Recalculate in-, output, respiration vectors+matrix of intercompartmental transfers:
    Z_f1<-list.FW[[i]]["IMPORT",c(1:comps)]
    E_f1<-list.FW[[i]][c(1:comps),"EXPORT"]
    R_f1<-list.FW[[i]][c(1:comps),"RESP"]
    T_f1<-list.FW[[i]][c(1:comps),c(1:comps)]  
  if(length(Z_f1)==0){
    print("ERROR in Network balancing procedure: det(X)=0. Abort function. Use `balance=unbalanced´ to avoid this error.")
    invokeRestart("abort")}}
  
  #Calculate true indices values of each foodweb and treatment:
  #-->must be specified in function as vector with string expressions(???)
  #-->food web connectance requires to indicate the number of non-living compartments (nl)
  
  ## --- (1) total system throughput
  i1 <- Tst(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
  ##
  ## --- (2) development capacity
  i2 <- dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
  ##
  ## --- (3) uncertainty (Shannon's index of diversity, H = DC/TST)
  i3 <- i2/i1
  ##
  #Ascendency
  output <- Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
  ## --- (4) ascendency (absolute value)
  i4<-output[1]
  ##
  ## --- (5) ascendency (ratio)
  i5 <- output[2]
  ##
  ## --- (6) average mutual information
  i6 <- Ami(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
  ##
  ##
  #Overhead
  output<-Overhead(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE, show.components=TRUE)
  ##
  ## --- (7) overhead on imports (absolute value)
  i7 <- output[3,1]
  ##
  ## --- (8) overhead on exports (absolute value)
  i8 <- output[4,1]
  ##
  ## --- (9) overhead on dissipations (absolute value)
  i9 <- output[5,1]
  ##
  ## --- (10) redundancy (absolute value)
  i10 <- output[2,1]
  ##
  ## --- (11) total overhead (absolute value)
  i11 <- output[1,1]
  ##
  ## --- (12) overhead on imports (relative value)
  i12 <- output[3,2]
  ##
  ## --- (13) overhead on exports (relative value)
  i13 <- output[4,2]
  ##
  ## --- (14) overhead on dissipations (relative value)
  i14 <- output[5,2]
  ##
  ## --- (15) redundancy (relative value)
  i15 <- output[2,2]
  ##
  ## --- (16) total overhead (relative value)
  i16 <- output[1,2]
  ##
  ## --- (17) residual diversity (Hc = O/TST)
  i17 <- i11/i1
  ##
  ## --- (18) internal diversity (int_Hc = R/TST)
  i18 <- i10/i1
  ##
  ##
  ##
  ##
  ## --- (19) internal capacity
  i19 <- Internal.dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
  ##
  ## --- (20) internal ascendency (relative value)
  i20 <- Internal.Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)[1]
  ##
  ## --- (21) internal redundancy (relative value)
  i21 <- i10/i20
  ##
  ##
  ##
  ##
  ## --- (22) overall connectance
  i22 <- connectance(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, type="whole")
  ##
  ## --- (23) intercompartmental connectance
  i23 <- connectance(T_cs=T_f1, type="intercompartmental")
  ##
  ## --- (24) food web connectance
  i24 <- connectance(T_cs=T_f1,nl=nl, type="foodweb")
  
if(attributes==FALSE){  
  inds<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,
          i21,i22,i23,i24)}
  
if(attributes==TRUE){
  
  ## --- (25) Finn cycling index (FCI, Finn 1980)
  att1 <- FCI(inp=Z_f1, inter=T_f1, outt=E_f1, diss=R_f1)[1]
  ##
  ## --- (26) effective amount of carbon recycled
  att2 <- att1 * i1
  ##
  ##
  ##
  ##
  ## --- (27) herbivory
  att3 <- sum(T_f1[PP,c(1:living)])
  
  ##
  ## --- (28) detritivory
  att4 <- sum(T_f1[c((living+1):comps),c(1:living)])
  ##
  ## --- (29) circulation within detrital pool
  att5 <- sum(T_f1[c((living+1):comps),c((living+1):comps)])
  ##
  ## --- (30) detritivory:herbivory
  att6 <- att4/att3
  ##
  ##
  ##
  ##if indicated:
  if(is.na(SDet)==FALSE){
    ## --- (31) export from sediment detritus (SD) = compartment SD (14)
    att7 <- E_f1[SDet]
  }else{att7<-NA}
  ##
  if(is.na(WDet)==FALSE){
    ## --- (32) export from water detritus (WD) = compartment WD (15)
    att8 <- E_f1[WDet]
  }else{att8<-NA}
  ##
  if(is.na(DOC)==FALSE){
    ## --- (33) export from water dissolved organic carbon (DOC) = compartment DOC (16)
    att9 <- E_f1[DOC]
  }else{att9<-NA}
  ##
  ## --- (34) respiration:NPP
  att10 <- sum(R_f1[c(1:living)])/(sum(Z_f1[PP])-sum(R_f1[PP]))#-->CAN LEAD TO NEG VALUES
  #att10<-NA
  ##
  #att11&att12 require TROPHIC ANALYSIS:
  if(trophic.analysis==TRUE){
  output<-TP.CTA(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl)#VERY SLOW!
  
  ## --- (35) grazing chain trophic efficiency (TL1 ---> TL2)
  s35_1 <- output[[5]][1]
  s35_2 <- output[[5]][2]
  att11 <- s35_2/s35_1
  #   
  ##
  ## --- (36) Lindeman spine trophic efficiency (TL1 ---> TL2)
  att12 <- output[[11]][1]}
  else{att11<-att12<-NA}
  ##
  ## --- (37) average path length (APL = [T - Z]/Z) - Wulff et al. (1989) at page 47
  #att13 <- (sum(T_f1) - sum(Z_f1[PP]))/sum(Z_f1[PP])-->CAN LEAD TO NEGATIVE VALUES
  #Alternative definition by Finn (1976)-->TST/Z
  att13<-i1/sum(Z_f1)
  ##
  ## --- (38) herbivory:NPP (herb:NPP)
  att14 <- att3/(sum(Z_f1[PP])-sum(R_f1[PP]))#-->CAN LEAD TO NEG VALUES
  #att14<-NA
  ##
  ## --- (39) system omnivory index (SOI), only for living compartments
  if(trophic.analysis==TRUE){
  att15<-SOI(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl) }#VERY SLOW!
  else{att15<-NA}
  ##
  ##
  ##
  ##
  ## --- (40) total gross primary productivity (GPP)
  att16 <- sum(Z_f1[PP])
  ##
  ## --- (41) total export from primary producers
  att17 <- sum(E_f1[PP])
  ##
  ## --- (42) total respiration from living compartments
  att18 <- sum(R_f1[c(1:living)])
  
  inds<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,
          i21,i22,i23,i24,
          att1,att2,att3,att4,att5,att6,att7,att8,att9,att10,att11,att12,att13,att14,
          att15,att16,att17,att18) }
  
  IND.TRUE[,i]<-inds
}#END OF FOR LOOP

inds.exclude<-which(is.na(inds)==TRUE)#exclude undefined indices (if there are)
if(length(inds.exclude)>0){inds<-inds[-inds.exclude]
IND.TRUE<-IND.TRUE[-inds.exclude,]
indices<-indices[-inds.exclude]}

#####################################################################################
###CALCULATE LOG-RATIOS FOR EACH TREATMENT
#based on mean values of replicates
#based on equations of Hedges et al. (1999)

cols.C<-which(design=="C")#indices columns of IND.TRUE assigned to control
reps.C<-length(cols.C)#replicate number of control
mean.C<-apply(IND.TRUE[,cols.C],1,mean,na.rm=TRUE)#mean indice values of control
sd.C<-apply(IND.TRUE[,cols.C],1,sd,na.rm=TRUE)#standard deviations of mean indices of control

treatment<-unique(design[which(design!="C")])#identify treatment levels

LOGR.TRUE<-matrix(NA,nrow=length(indices), ncol=(5*length(treatment)))#matrix to store log-ratio results
LOGR.names<-vector()
l<-1
for(i in 1:length(treatment)){
  
  cols.tr<-which(design==treatment[i])
  reps.tr<-length(cols.tr)
  mean.tr<-apply(IND.TRUE[,cols.tr],1,mean,na.rm=TRUE)
  sd.tr<-apply(IND.TRUE[,cols.tr],1,sd,na.rm=TRUE)
  
  #-->log-ratio: log(mean i_treatment)-log(mean i_control)
  LOGR.TRUE[,l]<-log(mean.tr)-log(mean.C)
  LRR.name<-paste("LRR_",treatment[i],sep="")
  
  #-->Variance: sd(mean i_treatment)^2/(Nb replicate treatment*(mean i_treatment)^2)+
  #      sd(mean i_control)^2/(Nb replicate control*(mean i_control)^2)   
  #(assumes normal distribution of mean i_control,mean i_treatment)
  LOGR.TRUE[,l+1]<-(sd(mean.tr,na.rm=TRUE)^2)/(reps.tr*(mean.tr)^2)+(sd(mean.C,na.rm=TRUE)^2)/(reps.C*(mean.C)^2)
  V.name<-paste("Variance_",treatment[i],sep="")  
  
  #-->Standard deviation: SD=sqrt(V) 
  LOGR.TRUE[,l+2]<-sqrt(LOGR.TRUE[,l+1])
  SD.name<-paste("SD_",treatment[i],sep="")
  
  #-->Standard error: SE=SD/sqrt(Number of replicates) 
  #(number of replicates should be the same for control and treatment???)
  mean.reps<-(reps.tr+reps.C)/2
  LOGR.TRUE[,l+3]<-LOGR.TRUE[,l+2]/sqrt(mean.reps)
  SE.name<-paste("SE_",treatment[i],sep="")
  
  #-->95% Confidence Interval: SE*1.96
  LOGR.TRUE[,l+4]<-LOGR.TRUE[,l+3]*1.96
  CI95.name<-paste("CI95_",treatment[i],sep="")
  
  LOGR.names<-c(LOGR.names,LRR.name,V.name,SD.name,SE.name,CI95.name)
  l<-l+5
}#END OF FOR LOOP

colnames(LOGR.TRUE)<-LOGR.names
rownames(LOGR.TRUE)<-indices

###COMPUTE RANDOM NETWORKS AND INDICES#####################################

#Re-organize total of food web flows in 1 matrix:
all.FW.links<-matrix(NA, nrow=length(topology), ncol=length(list.FW),
                     dimnames=list(topology, names(list.FW)))
for(i in 1:length(topology)){
  for(j in 1:length(list.FW)){
    all.FW.links[i,j]<-list.FW[[j]][topology[i]]}}

#Build matrix of total of flows to sample from:

###OPTION: exclude.extremes=TRUE
#-->for each spp-spp interaction the most extreme values are discarded

if(exclude.extremes==TRUE){
  sample.FW.links<-matrix(NA, nrow=length(topology), ncol=length(list.FW)-2,
                          dimnames=list(topology, NULL))
  for(i in 1:length(topology)){
    interaction<-all.FW.links[i,]
    max.flow<-which(interaction==max(interaction))
    if(length(max.flow)>1)max.flow<-max.flow[1]#if there are >1 maximal flows, discard first one
    min.flow<-which(interaction==min(interaction))
    if(length(min.flow)>1)min.flow<-min.flow[1]#if there are >1 minimal flows, discard first one
    IA.sample<-unname(interaction[-c(min.flow,max.flow)])
    sample.FW.links[i,]<-IA.sample } 
  
  ###ELSE: exclude.extremes=FALSE
  #-->use of all foodweb flows for sampling
}else{
  sample.FW.links<-all.FW.links}

#Special case: randomized topology
#define all topology positions that are not used for randomization
#-->RESPIRATION,EXPORT-rows and IMPORT-column 
#+constrain randomization of topology: 
#-->all respiration flows must be non-zero
if(constrain=="link number"){
  
  topo.outs<-list.FW[[1]]
  topo.outs[,c("RESP","IMPORT")]<-NA
  topo.outs[c("RESP","EXPORT"),]<-NA
  topo.outs["IMPORT",c((comps+1):ncol(topo.outs))]<-NA
  topo.outs<-which(is.na(topo.outs==TRUE))
  
  topo.resp<-list.FW[[1]]
  topo.resp[c(1:comps),"RESP"]<-"resp"
  topo.resp<-which(topo.resp=="resp")#these matrix positions are filled with 100% probability with non-zero values
}

#Create for simulation necessary objects:
list_INDS.RAND<-vector(mode="list", length=runs)#list to store calculated indices values of each simulation

LOGR.RAND<-matrix(NA,nrow=runs,ncol=length(indices), dimnames=list(NULL,indices))
list_LOGR.RAND<-vector(mode="list", length=length(treatment))#list to store log-ratios of each simulation:1 matrix per treatment
for(i in 1:length(treatment))list_LOGR.RAND[[i]]<-LOGR.RAND
names(list_LOGR.RAND)<-treatment

for(k in 1:runs){
  
  #1.Create set of (Nb of experimental networks-2) random matrices:
  if(exclude.extremes==TRUE){
    n<-length(list.FW)-2
  }else{n<-length(list.FW)}
  
  list.FW.rand<-vector(mode="list", length=n)#list to temporally store set of random networks
  INDS.RAND<-matrix(NA, nrow=length(indices), ncol=n,
                    dimnames=list(indices,NULL))#matrix to store indices values of set ofrandom networks
  
  for(p in 1:n){
    FW.rand<-matrix(0,ncol=comps+3, nrow=comps+3, 
                    dimnames=list(c(names.comps,"EXPORT","RESP","IMPORT")
                                  ,c(names.comps,"EXPORT","RESP","IMPORT")))
    
    #Fill matrix with randomized link strengths:
    if(constrain=="link number"){
      
      #constrain randomization of topology:
      #-->at least 1 import flow must be non-zero
      topo.imp<-FW.rand
      topo.imp["IMPORT",c(1:comps)]<-"import"
      topo.imp<-which(topo.imp=="import")
      imp.sample<-sample(topo.imp,1)#randomly select one import flow that is guaranteed non-zero
      topology.sample<-c(1:((comps+3)^2))[-c(topo.outs,imp.sample)]
      
      #sample links:
      all.links<-as.vector(unname(sample.FW.links))#extract all flows to sample from as vector
      links.rand<-sample(all.links, length(topology), replace=FALSE)
      
      #sample topology:
      #-->must include all respiration flows+at least 1 import flow
      topo.sample<-sample(topology.sample, size=(length(topology)-comps-1), replace=FALSE)
      topo.sample<-c(topo.sample, topo.resp, imp.sample)
      topology.rand<-sample(topo.sample, replace=FALSE)#randomly permutate topology vector
      FW.rand[topology.rand]<-links.rand}
    
    if(constrain=="topology"){
      #sample links:
      all.links<-as.vector(unname(sample.FW.links))#extract all flows to sample from as vector
      links.rand<-sample(all.links, length(topology), replace=FALSE)
      FW.rand[topology]<-links.rand }
    
    if(constrain=="link strength"){
      #sample links constrained by true transfers realized between species across treaments:
      for(i in 1:length(topology)){
        link.sample<-unname(sample.FW.links[i,])
        link.rand<-sample(link.sample,1,replace=FALSE)
        FW.rand[topology[i]]<-link.rand } }
    
    #Obtain in-, output, respiration vectors+matrix of intercompartmental transfers:
    Z_f1<-FW.rand["IMPORT",c(1:comps)]
    E_f1<-FW.rand[c(1:comps),"EXPORT"]
    R_f1<-FW.rand[c(1:comps),"RESP"]
    T_f1<-FW.rand[c(1:comps),c(1:comps)]  
    
    #2.Balance network (or not)
    if(balance!="unbalanced"){
      FW.rand<-network.balance(Z_cs=Z_f1,E_cs=E_f1,R_cs=R_f1,T_cs=T_f1,
                               method=balance)[[2]]
      #Recalculate in-, output, respiration vectors+matrix of intercompartmental transfers:
      Z_f1<-FW.rand["IMPORT",c(1:comps)]
      E_f1<-FW.rand[c(1:comps),"EXPORT"]
      R_f1<-FW.rand[c(1:comps),"RESP"]
      T_f1<-FW.rand[c(1:comps),c(1:comps)]  
      if(length(Z_f1)==0){
        print("ERROR in Network balancing procedure: det(X)=0. Abort function. Use `balance=unbalanced´ to avoid this error.")
        invokeRestart("abort")}}
    
    #3.Compute indices and store them in matrix
    
    ## --- (1) total system throughput
    i1.rand <- Tst(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
    ##
    ## --- (2) development capacity
    i2.rand <- dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
    ##
    ## --- (3) uncertainty (Shannon's index of diversity, H = DC/TST)
    i3.rand <- i2/i1
    ##
    #Ascendency
    output <- Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
    ## --- (4) ascendency (absolute value)
    i4.rand<-output[1]
    ##
    ## --- (5) ascendency (ratio)
    i5.rand <- output[2]
    ##
    ## --- (6) average mutual information
    i6.rand <- Ami(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
    ##
    ##
    #Overhead
    output<-Overhead(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE, show.components=TRUE)
    ##
    ## --- (7) overhead on imports (absolute value)
    i7.rand <- output[3,1]
    ##
    ## --- (8) overhead on exports (absolute value)
    i8.rand <- output[4,1]
    ##
    ## --- (9) overhead on dissipations (absolute value)
    i9.rand <- output[5,1]
    ##
    ## --- (10) redundancy (absolute value)
    i10.rand <- output[2,1]
    ##
    ## --- (11) total overhead (absolute value)
    i11.rand <- output[1,1]
    ##
    ## --- (12) overhead on imports (relative value)
    i12.rand <- output[3,2]
    ##
    ## --- (13) overhead on exports (relative value)
    i13.rand <- output[4,2]
    ##
    ## --- (14) overhead on dissipations (relative value)
    i14.rand <- output[5,2]
    ##
    ## --- (15) redundancy (relative value)
    i15.rand <- output[2,2]
    ##
    ## --- (16) total overhead (relative value)
    i16.rand <- output[1,2]
    ##
    ## --- (17) residual diversity (Hc = O/TST)
    i17.rand <- i11/i1
    ##
    ## --- (18) internal diversity (int_Hc = R/TST)
    i18.rand <- i10/i1
    ##
    ##
    ##
    ##
    ## --- (19) internal capacity
    i19.rand <- Internal.dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
    ##
    ## --- (20) internal ascendency (relative value)
    i20.rand <- Internal.Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)[1]
    ##
    ## --- (21) internal redundancy (relative value)
    i21.rand <- i10/i20
    ##
    ##
    ##
    ##
    ## --- (22) overall connectance
    i22.rand <- connectance(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, type="whole")
    ##
    ## --- (23) intercompartmental connectance
    i23.rand <- connectance(T_cs=T_f1, type="intercompartmental")
    ##
    ## --- (24) food web connectance
    i24.rand <- connectance(T_cs=T_f1,nl=nl, type="foodweb")

    if(attributes==FALSE){    
    INDS.RAND[,p]<-c(i1.rand,i2.rand,i3.rand,i4.rand,i5.rand,i6.rand,i7.rand,i8.rand,
                     i9.rand,i10.rand,i11.rand,i12.rand,i13.rand,i14.rand,i15.rand,
                     i16.rand,i17.rand,i18.rand,i19.rand,i20.rand,
                     i21.rand,i22.rand,i23.rand,i24.rand)}
    
    if(attributes==TRUE){  
    ## --- (25) Finn cycling index (FCI, Finn 1980)
    att1.rand <- FCI(inp=Z_f1, inter=T_f1, outt=E_f1, diss=R_f1)[1]
    ##
    ## --- (26) effective amount of carbon recycled
    att2.rand <- att1.rand * i1
    ##
    ##
    ##
    ##
    ## --- (27) herbivory
    att3.rand <- sum(T_f1[PP,c(1:living)])
    
    ##
    ## --- (28) detritivory
    att4.rand <- sum(T_f1[c((living+1):comps),c(1:living)])
    ##
    ## --- (29) circulation within detrital pool
    att5.rand <- sum(T_f1[c((living+1):comps),c((living+1):comps)])
    ##
    ## --- (30) detritivory:herbivory
    att6.rand <- att4.rand/att3.rand
    ##
    ##
    ##
    ##if indicated:
    if(is.na(SDet)==FALSE){
      ## --- (31) export from sediment detritus (SD) = compartment SD (14)
      att7.rand <- E_f1[SDet]
    }else{att7.rand<-NA}
    ##
    if(is.na(WDet)==FALSE){
      ## --- (32) export from water detritus (WD) = compartment WD (15)
      att8.rand <- E_f1[WDet]
    }else{att8.rand<-NA}
    ##
    if(is.na(DOC)==FALSE){
      ## --- (33) export from water dissolved organic carbon (DOC) = compartment DOC (16)
      att9.rand <- E_f1[DOC]
    }else{att9.rand<-NA}
    ##
    ## --- (34) respiration:NPP
    att10.rand <- sum(R_f1[c(1:living)])/(sum(Z_f1[PP])-sum(R_f1[PP]))#-->CAN LEAD TO NEG VALUES
    #att10.rand<-NA
    ##
    #att11&att12 require TROPHIC ANALYSIS:
    if(trophic.analysis==TRUE){
      output<-TP.CTA(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl)#VERY SLOW!
      
      ## --- (35) grazing chain trophic efficiency (TL1 ---> TL2)
      s35_1 <- output[[5]][1]
      s35_2 <- output[[5]][2]
      att11.rand <- s35_2/s35_1
      #   
      ##
      ## --- (36) Lindeman spine trophic efficiency (TL1 ---> TL2)
      att12.rand <- output[[11]][1]}
    else{att11.rand<-att12.rand<-NA}
    ##
    ## --- (37) average path length (APL = [T - Z]/Z) - Wulff et al. (1989) at page 47
    #att13.rand <- (sum(T_f1) - sum(Z_f1[PP]))/sum(Z_f1[PP])-->CAN LEAD TO NEGATIVE VALUES
    #Alternative definition by Finn (1976)-->TST/Z
    att13.rand<-i1/sum(Z_f1)
    ##
    ## --- (38) herbivory:NPP (herb:NPP)
    att14.rand <- att3.rand/(sum(Z_f1[PP])-sum(R_f1[PP]))#-->CAN LEAD TO NEG VALUE
    #att14.rand<-NA
    ##
    ## --- (39) system omnivory index (SOI), only for living compartments
    if(trophic.analysis==TRUE){
      att15.rand<-SOI(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl) }#VERY SLOW!
    else{att15.rand<-NA}
    ##
    ##
    ##
    ##
    ## --- (40) total gross primary productivity (GPP)
    att16.rand <- sum(Z_f1[PP])
    ##
    ## --- (41) total export from primary producers
    att17.rand <- sum(E_f1[PP])
    ##
    ## --- (42) total respiration from living compartments
    att18.rand <- sum(R_f1[c(1:living)])
    
    inds.rand<-c(i1.rand,i2.rand,i3.rand,i4.rand,i5.rand,i6.rand,i7.rand,i8.rand,
            i9.rand,i10.rand,i11.rand,i12.rand,i13.rand,i14.rand,i15.rand,
            i16.rand,i17.rand,i18.rand,i19.rand,i20.rand,
            i21.rand,i22.rand,i23.rand,i24.rand,
            att1.rand,att2.rand,att3.rand,att4.rand,att5.rand,att6.rand,att7.rand,
            att8.rand,att9.rand,att10.rand,att11.rand,att12.rand,att13.rand,
            att14.rand,att15.rand,att16.rand,att17.rand,att18.rand)
    if(length(inds.exclude)>0)inds.rand<-inds.rand[-inds.exclude]
    INDS.RAND[,p]<-inds.rand}
  }#END OF SET OF SET OF RANDOM NETWORKS FOR-LOOP
  
  list_INDS.RAND[[k]]<-INDS.RAND #store indices results of this simulation in list
  
  #4.Calculate log ratios and store them in matrix:
  
  #Randomly attribute a treatment to each random network of the set:
  
  ###OPTION: exclude.extremes=TRUE
  #-->Since 2 networks were excluded we randomly 2x1replicate/treatment
  #-->condition: excluded replicates are of different treatments
  if(exclude.extremes==TRUE){
    out1<-sample(design,1)
    which1<-which(design==out1)
    reps1<-length(which1)
    out2<-sample(design[-which1],1)
    which2<-which(design==out2)
    reps2<-length(which2)
    design.sample<-c(design[-c(which1,which2)],rep(out1,(reps1-1)),rep(out2,(reps2-1)))
    design.rand<-sample(design.sample,replace=FALSE)
  
  ###ELSE: exclude.extremes=FALSE
  }else{
    design.rand<-sample(design,replace=FALSE)  } #randomly permutate treatment vector
  
  #Compute log-ratios for each indice and store results in list
  cols.C.rand<-which(design.rand=="C")#indices columns of IND.RAND assigned to control
  mean.C.rand<-apply(INDS.RAND[,cols.C.rand],1,mean,na.rm=TRUE)#mean indice values of control
  
  for(i in 1:length(treatment)){
    cols.tr.rand<-which(design.rand==treatment[i])
    mean.tr.rand<-apply(INDS.RAND[,cols.tr.rand],1,mean,na.rm=TRUE)
    list_LOGR.RAND[[i]][k,]<-log(mean.tr.rand)-log(mean.C.rand) }
  
}#END OF SIMULATION-FOR LOOP

###CALCULATE QUANTILES FOR UPPER/LOWER LIMIT OF LOG-RATIO DISTRIBUTIONS###############
#10%-threshold, 5%-threshold, 2.5%-threshold, 1%-threshold
#0% and 100% quantiles
#AND
###CALCULATE P-VALUES##############################################################
#-->Identify which indice log-ratio are significantly lower or larger across treatments

LOGR.QUANTILES<-matrix(NA, nrow=length(indices), ncol=length(treatment)*5*2,
                       dimnames=list(indices,NULL))#empty matrix to store quantile results
names.quantiles<-c("LRR0%","LRR1%","LRR2.5%","LRR5%","LRR10%","LRR90%","LRR95%","LRR97.5%","LRR99%","LRR100%")
colnames.LOGR.QUANTILES<-vector()#column names are assigned during computation of quantiles

INDS.TRUE.p_val<-matrix(NA,nrow=length(indices), ncol=length(treatment)*2,
                        dimnames=list(indices,NULL))#matrix to store p-value results
names.p_val<-c("P_Left.tailed.test","P_Right.tailed.test")
colnames.INDS.TRUE.p_val<-vector()

l<-1
k<-1
g<-1

for(i in 1:length(treatment)){
  
  #Pool random log-ratios and true log-ratios:
  all_LOGR<-rbind(list_LOGR.RAND[[i]],LOGR.TRUE[,l])
  
  #1.Calculate quantiles and store results:

  ##(for 1-tailed testing)
  #lower limit (zero-quantile):
  Q_zero<-apply(all_LOGR,2,function(x)quantile(x,probs=0,na.rm=TRUE))
  #upper limit (100%-quantile):
  Q_100<-apply(all_LOGR,2,function(x)quantile(x,probs=1,na.rm=TRUE))
  
  ##10%:
  upper_10<-apply(all_LOGR,2,function(x)quantile(x,probs=0.9,na.rm=TRUE))
  lower_10<-apply(all_LOGR,2,function(x)quantile(x,probs=0.1,na.rm=TRUE))
  ##5%:
  upper_5<-apply(all_LOGR,2,function(x)quantile(x,probs=0.95,na.rm=TRUE))
  lower_5<-apply(all_LOGR,2,function(x)quantile(x,probs=0.05,na.rm=TRUE))
  ##2.5%:
  upper_2.5<-apply(all_LOGR,2,function(x)quantile(x,probs=0.975,na.rm=TRUE))
  lower_2.5<-apply(all_LOGR,2,function(x)quantile(x,probs=0.025,na.rm=TRUE))
  ##1%:
  upper_1<-apply(all_LOGR,2,function(x)quantile(x,probs=0.99,na.rm=TRUE))
  lower_1<-apply(all_LOGR,2,function(x)quantile(x,probs=0.01,na.rm=TRUE))
  
  
  LOGR.QUANTILES[,k]<-unname(Q_zero)
  LOGR.QUANTILES[,k+1]<-unname(lower_1)
  LOGR.QUANTILES[,k+2]<-unname(lower_2.5)
  LOGR.QUANTILES[,k+3]<-unname(lower_5)
  LOGR.QUANTILES[,k+4]<-unname(lower_10)
  LOGR.QUANTILES[,k+5]<-unname(upper_10)
  LOGR.QUANTILES[,k+6]<-unname(upper_5)
  LOGR.QUANTILES[,k+7]<-unname(upper_2.5)
  LOGR.QUANTILES[,k+8]<-unname(upper_1)
  LOGR.QUANTILES[,k+9]<-unname(Q_100)
  
  colnames<-paste(names.quantiles,treatment[i],sep="_")
  colnames.LOGR.QUANTILES<-c(colnames.LOGR.QUANTILES,colnames)
  
  #2.calculate p-values of each indice:
  
  #1-TAILED TEST:
  p_inf<-p_sup<-rep(NA,length(indices))
  
  for(z in 1:length(indices)){
    true.logr<-LOGR.TRUE[z,l]
    p_inf[z]<-length(which(all_LOGR[,z]<=true.logr))/nrow(all_LOGR)
    p_sup[z]<-length(which(all_LOGR[,z]>=true.logr))/nrow(all_LOGR)  }
  
  INDS.TRUE.p_val[,g]<-p_inf
  INDS.TRUE.p_val[,g+1]<-p_sup
  names.p<-paste(names.p_val,treatment[i],sep="_")
  colnames.INDS.TRUE.p_val<-c(colnames.INDS.TRUE.p_val,names.p)
  
  k<-k+10
  l<-l+5
  g<-g+2
  
  #2-TAILED TEST?
  
}#END OF FOR LOOP
colnames(LOGR.QUANTILES)<-colnames.LOGR.QUANTILES
colnames(INDS.TRUE.p_val)<-colnames.INDS.TRUE.p_val

###PRODUCE BARPLOT FIGURE####################################################################

#1.Reorganize log-ratio data
BARPLOT.data<-as.data.frame(matrix(NA,nrow=length(LOGR.RAND)*2,ncol=6,
                                   dimnames=list(NULL,c("treatment","indice","LOGR","LOGR.mean","Q_2.5","Q_97.5"))))
indice.barplot<-treatment.barplot<-LOGR.barplot<-means.barplot<-lowerQ.barplot<-upperQ.barplot<-vector()
l<-3
for(i in 1:length(treatment)){
  means<-apply(list_LOGR.RAND[[i]],2,mean,na.rm=TRUE)
  for(z in 1:length(indices)){
    LOGR.barplot<-c(LOGR.barplot,list_LOGR.RAND[[i]][,z])
    means.barplot<-c(means.barplot,rep(means[z],runs))
    indice.barplot<-c(indice.barplot,rep(indices[z],runs))
    treatment.barplot<-c(treatment.barplot,rep(treatment[i],runs))
    lowerQ.barplot<-c(lowerQ.barplot,rep(LOGR.QUANTILES[z,l],runs))
    upperQ.barplot<-c(upperQ.barplot,rep(LOGR.QUANTILES[z,l+5],runs))}
  l<-l+10}
BARPLOT.data[,"treatment"]<-treatment.barplot
BARPLOT.data[,"indice"]<-indice.barplot
BARPLOT.data[,"LOGR"]<-LOGR.barplot
BARPLOT.data[,"LOGR.mean"]<-means.barplot
BARPLOT.data[which(is.infinite(BARPLOT.data[,"LOGR.mean"])==TRUE),"LOGR.mean"]<-NA
BARPLOT.data[,"Q_2.5"]<-lowerQ.barplot
BARPLOT.data[,"Q_97.5"]<-upperQ.barplot
BARPLOT.data[which(is.infinite(BARPLOT.data[,"Q_2.5"])==TRUE),"Q_2.5"]<-NA
BARPLOT.data[which(is.infinite(BARPLOT.data[,"Q_97.5"])==TRUE),"Q_97.5"]<-NA
neg<-BARPLOT.data
neg[which(BARPLOT.data[,"LOGR.mean"]>0),"LOGR.mean"]<-0
neg[which(BARPLOT.data[,"LOGR.mean"]>0),"Q_2.5"]<-0
neg[which(BARPLOT.data[,"LOGR.mean"]>0),"Q_97.5"]<-0
pos<-BARPLOT.data
pos[which(BARPLOT.data[,"LOGR.mean"]<0),"LOGR.mean"]<-0
pos[which(BARPLOT.data[,"LOGR.mean"]<0),"Q_2.5"]<-0
pos[which(BARPLOT.data[,"LOGR.mean"]<0),"Q_97.5"]<-0

if(attributes==FALSE){
#2.Produce ggplot-figure
BARPLOT.data$indice<-factor(BARPLOT.data$indice,levels=indices)
neg$indice<-factor(neg$indice,levels=indices)
pos$indice<-factor(pos$indice,levels=indices)
plot<-
ggplot(BARPLOT.data, aes(x = LOGR.mean, y = indice,xmin=Q_2.5, xmax=Q_97.5, fill = treatment)) + 
  theme_bw()+
  geom_bar(data = neg,position="dodge", stat = "identity") +
  geom_bar(data = pos,position="dodge", stat = "identity") + 
  geom_errorbarh(color="gray48",position =position_dodge(.9), height=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  xlab("LRR")+
  ylab("Information/Connectance Indices")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        legend.title = element_blank(),
        axis.title.x = element_text(face="bold", size=10, color="grey28"),
        axis.title.y = element_text(face="bold", size=10, color="grey48")) }

if(attributes==TRUE){
  #2a.Subset data in information/connectance indices & attribute indices
  neg.split<-split(neg,neg[,"indice"])
  pos.split<-split(pos,pos[,"indice"])
  ICi<-indices[c(1:24)]#information and connectance indices
  Ai<-indices[c(25:length(indices))]
  neg.ICi<-neg.Ai<-pos.ICi<-pos.Ai<-data.frame()
  for(i in 1:length(ICi)){
    index.neg.ICi<-neg.split[[which(names(neg.split)==ICi[i])]]
    index.pos.ICi<-pos.split[[which(names(pos.split)==ICi[i])]]
    if(nrow(neg.ICi)==0){neg.ICi<-index.neg.ICi}
    else{neg.ICi<-rbind(neg.ICi,index.neg.ICi)}
    if(nrow(pos.ICi)==0){pos.ICi<-index.pos.ICi}
    else{pos.ICi<-rbind(pos.ICi,index.pos.ICi)} }
  BARPLOT.data.ICi<-rbind(neg.ICi,pos.ICi)
  for(i in 1:length(Ai)){
    index.neg.Ai<-neg.split[[which(names(neg.split)==Ai[i])]]
    index.pos.Ai<-pos.split[[which(names(pos.split)==Ai[i])]]
    if(nrow(neg.Ai)==0){neg.Ai<-index.neg.Ai}
    else{neg.Ai<-rbind(neg.Ai,index.neg.Ai)}
    if(nrow(pos.Ai)==0){pos.Ai<-index.pos.Ai}
    else{pos.Ai<-rbind(pos.Ai,index.pos.Ai)} }
  BARPLOT.data.Ai<-rbind(neg.Ai,pos.Ai)
  
  #2b.Produce 2 ggplot-figures
  #-->information/connectance indices
  #-->attribute indices
  #Plot1:
  neg.ICi$indice<-factor(neg.ICi$indice,levels=ICi)
  pos.ICi$indice<-factor(pos.ICi$indice,levels=ICi)
  BARPLOT.data.ICi$indice<-factor(BARPLOT.data.ICi$indice,levels=ICi)
  plot1<-
  ggplot(BARPLOT.data.ICi, aes(x = LOGR.mean, y = indice,xmin=Q_2.5, xmax=Q_97.5, fill = treatment)) + 
    theme_bw()+
    geom_bar(data = neg.ICi,position="dodge", stat = "identity") +
    geom_bar(data = pos.ICi,position="dodge", stat = "identity") + 
    geom_errorbarh(color="gray48",position =position_dodge(.9), height=0.25)+
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
    xlab("LRR")+
    ylab("Information/Connectance Indices")+
    theme(panel.grid.major.x = element_blank() ,
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="grey" ),
          legend.title = element_blank(),
          axis.title.x = element_text(face="bold", size=10, color="grey28"),
          axis.title.y = element_text(face="bold", size=10, color="grey48"))
  #Plot2:
  neg.Ai$indice<-factor(neg.Ai$indice,levels=Ai)
  pos.Ai$indice<-factor(pos.Ai$indice,levels=Ai)
  BARPLOT.data.Ai$indice<-factor(BARPLOT.data.Ai$indice,levels=Ai)
  plot2<-
  ggplot(BARPLOT.data.Ai, aes(x = LOGR.mean, y = indice,xmin=Q_2.5, xmax=Q_97.5, fill = treatment)) + 
    theme_bw()+
    geom_bar(data = neg.Ai,position="dodge", stat = "identity") +
    geom_bar(data = pos.Ai,position="dodge", stat = "identity") + 
    geom_errorbarh(color="gray48",position =position_dodge(.9), height=0.25)+
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
    xlab("LRR")+
    ylab("Attribute Indices")+
    theme(panel.grid.major.x = element_blank() ,
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="grey" ),
          legend.title = element_blank(),
          axis.title.x = element_text(face="bold", size=10, color="grey28"),
          axis.title.y = element_text(face="bold", size=10, color="grey48")) }


###FUNCTION OUTPUT#################################################################
#-->True indices
#-->True log-ratios
#-->Random log-ratios (+ indices)
#-->Quantiles
#-->p-values

OUTPUT.FINAL<-list(IND.TRUE, LOGR.TRUE, list_INDS.RAND, list_LOGR.RAND,
                   LOGR.QUANTILES,INDS.TRUE.p_val)
names(OUTPUT.FINAL)<-c("Indices of true networks", "Log-ratios (LRR) of 
                       true networks","Indices of simulated networks", 
                       "Log-ratios of simulated networks", "Quantiles for upper 
                       and lower limit of LRR distributions", "Summary of 
                       One Tailed Significance Tests")
if(attributes==FALSE){
  OUTPUT.FINAL[["Figure"]]<-plot}
if(attributes==TRUE){
  OUTPUT.FINAL[["Figure 1"]]<-plot1
  OUTPUT.FINAL[["Figure 2"]]<-plot2}
if(balance!="unbalanced"){
  OUTPUT.FINAL[["Balanced true networks"]]<-list.FW }


#if(exclude.extremes==TRUE)-->SHOULD FUNCTION RETURN EDGELIST OF EXCLUDED VALUES?

return(OUTPUT.FINAL)
} #END OF FUNCTION
##################################################################################

