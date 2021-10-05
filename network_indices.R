###FUNCTIONS FOR NETWORK INDICES-->SPEEDED UP############################################
#Approach: we set the option partial.matrix=TRUE/FALSE
#If TRUE-->output returns index value and partial matrix (slow function)
#If FALSE-->output only returns index value (hopefully rapid function)

#-->computed based on matrix of intercompartmental exchanges, input vector,
#   export vector and respiration vector

########################################################################################
#INFORMATION INDICES####################################################################
########################################################################################

##########
##Total System Throughput
Tst<-function(Z_cs, E_cs, R_cs, T_cs){
  TST<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)
  return(TST)
}#END OF TST FUNCTION

##########
##Shannon´s Diversity of flows (H, derived from Shannon´s information measure)

shannon.flow<-function(Z_cs=NA, E_cs=NA, R_cs=NA, T_cs, nl, type, partial.matrix=TRUE, k=1){
  
  if(type=="whole"){
    #-->takes into account non-living nodes and exchanges with surroundings (i.e., import,
    #   export and respiration vectors)  
    T_matrix<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
    T_matrix<-cbind(T_matrix,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))#define matrix of transfers
    rownames(T_matrix)<-colnames(T_matrix)<-c(rownames(T_cs), "EXPORT", "RESP", "IMPORT")
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="intercompartmental"){
    #-->takes into account only exchanges between compartments (living and non-living)
    T_matrix<-T_cs #define matrix of transfers
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="foodweb"){
    #-->takes into account only exchanges between living compartments  
    #-->requires additional information about number of non-living compartments (nl)
    living<-nrow(T_cs)-nl
    T_matrix<-T_cs[c(1:living), c(1:living)] #define matrix of transfers
    comps<-living} 
  
  if(living==0 || 
     living==1 && T_matrix==0){ #NO LIVING COMPARTMENTS OR ONE LIVING COMPARTMENT & NO FLOWS
    list_h<-list(NA,NA)
    h.FINAL<-NA
    
  }else if(living==1 && T_matrix>0){ #ONE LIVING COMPARTMENT & FLOW(=CANNIBALISM)
    list_h<-list(0,NA) 
    h.FINAL<-0
  
  }else{
    
    #2.Compute H
    Tdotdot<-sum(T_matrix)#total system throughput
    whichflows<-which(T_matrix!=0)
    flows<-T_matrix[which(T_matrix!=0)]
    
    if(partial.matrix==TRUE){
      T_matrix.h<-matrix(0,nrow=nrow(T_matrix),ncol=ncol(T_matrix), 
                                               dimnames =list(rownames(T_matrix), colnames(T_matrix)))}
    if(partial.matrix==FALSE){h.FINAL<-0}
    
    if(length(flows)>0){
      
      for(i in 1:length(flows)){
        
        if(partial.matrix==TRUE){
        T_matrix.h[whichflows[i]]<-NA
        row<-which(is.na(T_matrix.h) == TRUE, arr.ind=TRUE)[1]
        col<-which(is.na(T_matrix.h) == TRUE, arr.ind=TRUE)[2]
        h<-(-1)*k*(T_matrix[row,col]/Tdotdot)*
          log2(T_matrix[row,col]/Tdotdot)
        T_matrix.h[row,col]<-h }
        
        if(partial.matrix==FALSE){
          h<-(-1)*k*(T_matrix[whichflows[i]]/Tdotdot)*
            log2(T_matrix[whichflows[i]]/Tdotdot) 
          h.FINAL<-h.FINAL+h }}}
    
    if(partial.matrix==TRUE)list_h<-list(sum(T_matrix.h), T_matrix.h) }#END OF ELSE CLAUSE
  
  if(partial.matrix==TRUE){
  names(list_h)<-c("Shannon Diversity of flows (H)","Partial H Matrix")
  return(list_h)}
  
  if(partial.matrix==FALSE){return(h.FINAL)}
}#END OF H FUNCTION


##########
##Average Mutual Information
Ami<-function(Z_cs, E_cs, R_cs, T_cs, partial.matrix=TRUE, k=1){
  
  #1.Produce extended transfer matrix T*:
  Tstar<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
  Tstar<-cbind(Tstar,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))
  
  #2.Compute AMI
  Tdotdot<-sum(Tstar)#total system throughput
  totIN<-apply(Tstar, 2,sum)
  totOUT<-apply(Tstar,1,sum)
  whichflows<-which(Tstar!=0)
  flows<-Tstar[which(Tstar!=0)]
  
  if(partial.matrix==TRUE){
  Tstar.ami<-matrix(0,nrow=nrow(Tstar),ncol=ncol(Tstar), 
                    dimnames =list(c(rownames(T_cs), "EXPORT", "RESP", "IMPORT"),
                                   c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")))}
  if(partial.matrix==FALSE){
    ami.FINAL<-0
    Tstar.OUT<-matrix(rep(totOUT, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                      byrow=FALSE)
    Tstar.IN<-matrix(rep(totIN, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                     byrow=TRUE) }
  
 for(i in 1:length(flows)){
    
 if(partial.matrix==TRUE){
    Tstar.ami[whichflows[i]]<-NA
    row<-which(is.na(Tstar.ami) == TRUE, arr.ind=TRUE)[1]
    col<-which(is.na(Tstar.ami) == TRUE, arr.ind=TRUE)[2]
    if(totOUT[row]==0){
      Tstar.ami[row, col]<-0}
    else if(totIN[col]==0){
      Tstar.ami[row, col]<-0}
    else{ami<-k*(Tstar[row,col]/Tdotdot)*
      log2((Tstar[row,col]*Tdotdot)/(totOUT[row]*totIN[col]))
    Tstar.ami[row,col]<-ami }}
    
  if(partial.matrix==FALSE){
    ami<-k*(Tstar[whichflows[i]]/Tdotdot)*
      log2((Tstar[whichflows[i]]*Tdotdot)/(Tstar.OUT[whichflows[i]]*Tstar.IN[whichflows[i]]))
    if(is.infinite(ami)==TRUE)ami<-0
    ami.FINAL<-ami.FINAL+ami} }
  
  if(partial.matrix==TRUE){
  list_ami<-list(sum(Tstar.ami), Tstar.ami)
  names(list_ami)<-c("Average Mutual Information (AMI)","Partial AMI Matrix")
  return(list_ami)}
  
  if(partial.matrix==FALSE){return(ami.FINAL)}
}#END OF AMI FUNCTION

################
##Maximal Average Mutual Information
#equation: AMI.max= DC/TST
Ami.max<-function(Z_cs, E_cs, R_cs, T_cs){
  Tdotdot<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)
  DC<-dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)
  ami.max<-DC/Tdotdot
  return(ami.max)
}#END OF Ami.max FUNCTION

################
##Development Capacity (DC)
dC<-function(Z_cs, E_cs, R_cs, T_cs, partial.matrix=TRUE){
  
  #1.Produce extended transfer matrix T*:
  Tstar<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
  Tstar<-cbind(Tstar,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))
  
  #2.Compute DC
  Tdotdot<-sum(Tstar)#total system throughput
  whichflows<-which(Tstar!=0)
  flows<-Tstar[which(Tstar!=0)]
  
  if(partial.matrix==TRUE){
  Tstar.DC<-matrix(0,nrow=nrow(Tstar),ncol=ncol(Tstar), 
                   dimnames =list(c(rownames(T_cs), "EXPORT", "RESP", "IMPORT"),
                                  c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")))}
  
  if(partial.matrix==FALSE){DC.FINAL<-0}
  
  for(i in 1:length(flows)){
    
  if(partial.matrix==TRUE){
    Tstar.DC[whichflows[i]]<-NA
    row<-which(is.na(Tstar.DC) == TRUE, arr.ind=TRUE)[1]
    col<-which(is.na(Tstar.DC) == TRUE, arr.ind=TRUE)[2]
    DC<-(-1)*Tstar[row,col] * log2(Tstar[row,col]/Tdotdot)
    Tstar.DC[row,col]<-DC }
    
  if(partial.matrix==FALSE){
    DC<-(-1)*Tstar[whichflows[i]] * log2(Tstar[whichflows[i]]/Tdotdot)
    DC.FINAL<-DC.FINAL+DC}}
  
  if(partial.matrix==TRUE){
  list_DC<-list(sum(Tstar.DC), Tstar.DC)
  names(list_DC)<-c("Development Capacity (DC)","Partial DC Matrix")
  return(list_DC)}
  
  if(partial.matrix==FALSE){return(DC.FINAL)}
}#END OF DC FUNCTION

############
##Network Ascendency (A)
Asc<-function(Z_cs, E_cs, R_cs, T_cs, partial.matrix=TRUE){
  
  #1.Produce extended transfer matrix T*:
  Tstar<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
  Tstar<-cbind(Tstar,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))
  
  #2.Compute A
  Tdotdot<-sum(Tstar)#total system throughput
  totIN<-apply(Tstar, 2,sum)
  totOUT<-apply(Tstar,1,sum)
  whichflows<-which(Tstar!=0)
  flows<-Tstar[which(Tstar!=0)]
  
  if(partial.matrix==TRUE){
  Tstar.a<-matrix(0,nrow=nrow(Tstar),ncol=ncol(Tstar), 
                  dimnames =list(c(rownames(T_cs), "EXPORT", "RESP", "IMPORT"),
                                 c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")))} 
  
  if(partial.matrix==FALSE){
    A.FINAL<-0
    Tstar.OUT<-matrix(rep(totOUT, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                      byrow=FALSE)
    Tstar.IN<-matrix(rep(totIN, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                     byrow=TRUE) }
  
  for(i in 1:length(flows)){
    
  if(partial.matrix==TRUE){
    Tstar.a[whichflows[i]]<-NA
    row<-which(is.na(Tstar.a) == TRUE, arr.ind=TRUE)[1]
    col<-which(is.na(Tstar.a) == TRUE, arr.ind=TRUE)[2]
    if(totOUT[row]==0){
      Tstar.a[row, col]<-0}
    else if(totIN[col]==0){
      Tstar.a[row, col]<-0}
    else{A<-Tstar[row,col]*
      log2((Tstar[row,col]*Tdotdot)/(totOUT[row]*totIN[col]))
    Tstar.a[row,col]<-A }}
    
  if(partial.matrix==FALSE){
    A<-Tstar[whichflows[i]]*
      log2((Tstar[whichflows[i]]*Tdotdot)/(Tstar.OUT[whichflows[i]]*Tstar.IN[whichflows[i]])) 
    if(is.infinite(A)==TRUE)A<-0  
    A.FINAL<-A.FINAL+A  }}
  
  if(partial.matrix==TRUE)A.FINAL<-sum(Tstar.a)
  
  DC<-dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)#calculate development capacity for scaling
  A.ratio<-A.FINAL/DC #ascendency scaled by development capacity
  output<-c(A.FINAL, A.ratio)
  names(output)<-c("value", "ratio")
  
  if(partial.matrix==TRUE){
  list_a<-list(output, Tstar.a)
  names(list_a)<-c("Network Ascendency (A)", "Partial A Matrix")
  return(list_a)}
  
  if(partial.matrix==FALSE){return(output)}
}#END OF A FUNCTION

##############
##Overhead (OV)
Overhead<-function(Z_cs, E_cs, R_cs, T_cs, partial.matrix=TRUE, show.components=TRUE){
  
  #1.Produce extended transfer matrix T*:
  Tstar<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
  Tstar<-cbind(Tstar,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))
  
  #2.Compute OV
  totIN<-apply(Tstar, 2,sum)
  totOUT<-apply(Tstar,1,sum)
  whichflows<-which(Tstar!=0)
  flows<-Tstar[which(Tstar!=0)]
  
  if(partial.matrix==FALSE && show.components==FALSE){
    OV.FINAL<-0
    Tstar.OUT<-matrix(rep(totOUT, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                      byrow=FALSE)
    Tstar.IN<-matrix(rep(totIN, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                     byrow=TRUE) }
  
  else{
  Tstar.ov<-matrix(0,nrow=nrow(Tstar),ncol=ncol(Tstar), 
                   dimnames =list(c(rownames(T_cs), "EXPORT", "RESP", "IMPORT"),
                                  c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")))} 
  

  
  for(i in 1:length(flows)){
    
  if(partial.matrix==FALSE && show.components==FALSE){
      ov<-(-1)*Tstar[whichflows[i]]*
        log2((Tstar[whichflows[i]])^2/(Tstar.OUT[whichflows[i]]*Tstar.IN[whichflows[i]]))
      if(is.infinite(ov)==TRUE)ov<-0
      OV.FINAL<-OV.FINAL+ov }
    
  else{
    Tstar.ov[whichflows[i]]<-NA
    row<-which(is.na(Tstar.ov) == TRUE, arr.ind=TRUE)[1]
    col<-which(is.na(Tstar.ov) == TRUE, arr.ind=TRUE)[2]
    if(totOUT[row]==0){
      Tstar.ov[row, col]<-0}
    else if(totIN[col]==0){
      Tstar.ov[row, col]<-0}
    else{ov<-(-1)*Tstar[row,col]*
      log2((Tstar[row,col])^2/(totOUT[row]*totIN[col]))
    Tstar.ov[row,col]<-ov} } }
  
  if(partial.matrix==TRUE ||
     partial.matrix==FALSE && show.components==TRUE){OV.FINAL<-sum(Tstar.ov)}
  DC<-dC(Z_cs, E_cs, R_cs, T_cs, partial.matrix=FALSE)#calculate development capacity for scaling
  OV.ratio<-OV.FINAL/DC #overhead scaled by development capacity
  output<-c(OV.FINAL,OV.ratio)
  names(output)<-c("value", "ratio")
  
  if(partial.matrix==TRUE){
  list_ov<-list(output, Tstar.ov)
  names(list_ov)<-c("Overhead (OV)", "Partial OV Matrix")}
  
  if(show.components==TRUE){
    
    #####################
    ##Flow Redundancy (OVr)
    Tstar.OVr<-Tstar.ov[c(1:nrow(T_cs)), c(1:ncol(T_cs))]
    OVr.value<-sum(Tstar.OVr)
    OVr.ratio<-OVr.value/DC #redundancy scaled by DC
    #####################    
    ##Overhead due to Import (OVi)
    Tstar.OVi<-Tstar.ov[nrow(Tstar.ov), c(1:ncol(Tstar.ov))]
    OVi.value<-sum(Tstar.OVi)
    OVi.ratio<-OVi.value/DC
    #####################
    ##Overhead due to Export (OVe)
    Tstar.OVe<-Tstar.ov[c(1:nrow(Tstar.ov)), ncol(Tstar.ov)-2]
    OVe.value<-sum(Tstar.OVe)
    OVe.ratio<-OVe.value/DC
    #####################
    ##Overhead due to Dissipation (OVd)
    Tstar.OVd<-Tstar.ov[c(1:nrow(Tstar.ov)), ncol(Tstar.ov)-1]
    OVd.value<-sum(Tstar.OVd)
    OVd.ratio<-OVd.value/DC
    
    output2<-matrix(c(OVr.value, OVr.ratio,
                      OVi.value, OVi.ratio,
                      OVe.value, OVe.ratio,
                      OVd.value, OVd.ratio), nrow=4, ncol=2, byrow=TRUE,
                    dimnames=list(NULL, c("value", "ratio")))
    output.final<-rbind(output, output2)
    rownames(output.final)<-c("Total Overhead","Redundancy", "Overhead of Imports",
                              "Overhead of Export", "Dissipative Overhead")
    }#END OF show.components=TRUE CLAUSE 
  
  if(partial.matrix==FALSE && show.components==FALSE){return(output)}
  if(partial.matrix==TRUE && show.components==TRUE){
    list_ov[[1]]<-output.final
    return(list_ov)}
  else{return(output.final)}
  
}#END OF OV FUNCTION

############
##Internal development capacity (DCi)
#internal.dC<-function(Z_cs, E_cs, R_cs, T_cs){
  
  #Obtain Partial DC Matrix:
  #Tstar.DC<-dC(Z_cs, E_cs, R_cs, T_cs)[[2]]
  
  #Compute DCi:
  #Tstar.DCi<-Tstar.DC[c(1:nrow(T_cs)), c(1:ncol(T_cs))]
  #DCi<-sum(Tstar.DCi)
  #return(DCi)
#}#END OF DCi FUNCTION

#OR (faster):
Internal.dC<-function(Z_cs, E_cs, R_cs, T_cs){
  
  Tdotdot<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)#total system throughput
  whichflows<-which(T_cs!=0)
  flows<-T_cs[which(T_cs!=0)]
  DCi.FINAL<-0
  
  if(length(flows>0)){
  for(i in 1:length(flows)){
    DCi<-(-1)*T_cs[whichflows[i]] * log2(T_cs[whichflows[i]]/Tdotdot)
    DCi.FINAL<-DCi.FINAL+DCi}}

  return(DCi.FINAL)
}#END OF DCi FUNCTION


############
##Internal Ascendency (Ai)
#internal.Asc<-function(Z_cs, E_cs, R_cs, T_cs){
  
  #Obtain Partial A Matrix:
  #Tstar.A<-Asc(Z_cs, E_cs, R_cs, T_cs)[[2]]
  
  #Compute Ai:
  #Tstar.Ai<-Tstar.A[c(1:nrow(T_cs)), c(1:ncol(T_cs))]
  #Ai.value<-sum(Tstar.Ai)
  #DCi<-internal.dC(Z_cs, E_cs, R_cs, T_cs)#calculate internal DC for scaling
  #Ai.ratio<-Ai.value/DCi #internal ascendency scaled by internal DC
  
  #output<-matrix(c(Ai.value, Ai.ratio), nrow=1, ncol=2, 
                 #dimnames=list("Internal Ascendency (Ai)", c("value","ratio")))
  
 # return(output)
#}#END OF Ai FUNCTION
#-->In ENA function of internal ascendency also redundancy is returned (??)
#-->redundancy is scaled by internal development capacity
  
#OR (faster):

Internal.Asc<-function(Z_cs, E_cs, R_cs, T_cs){
  
  #1.Produce extended transfer matrix T*:
  Tstar<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
  Tstar<-cbind(Tstar,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))
  
  #2.Compute Ai
  Tdotdot<-sum(Tstar)#total system throughput
  totIN<-apply(Tstar, 2,sum)
  totOUT<-apply(Tstar,1,sum)
  Tstar.OUT<-matrix(rep(totOUT, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                    byrow=FALSE)
  Tstar.OUT<-Tstar.OUT[c(1:nrow(T_cs)),c(1:ncol(T_cs))]
  Tstar.IN<-matrix(rep(totIN, ncol(Tstar)),nrow=nrow(Tstar),ncol=ncol(Tstar),
                   byrow=TRUE)
  Tstar.IN<-Tstar.IN[c(1:nrow(T_cs)),c(1:ncol(T_cs))]
  whichflows<-which(T_cs!=0)
  flows<-T_cs[which(T_cs!=0)]
  Ai.FINAL<-0
  
  if(length(flows)>0){
  for(i in 1:length(flows)){
    Ai<-T_cs[whichflows[i]]*
      log2((T_cs[whichflows[i]]*Tdotdot)/(Tstar.OUT[whichflows[i]]*Tstar.IN[whichflows[i]])) 
    if(is.infinite(Ai)==TRUE)Ai<-0
    Ai.FINAL<-Ai.FINAL+Ai  }}
  
  DCi<-Internal.dC(Z_cs, E_cs, R_cs, T_cs)#calculate internal DC for scaling
  Ai.ratio<-Ai.FINAL/DCi #internal ascendency scaled by internal DC
  
  output<-matrix(c(Ai.FINAL, Ai.ratio), nrow=1, ncol=2, 
                 dimnames=list("Internal Ascendency (Ai)", c("value","ratio")))
  return(output)
}#END OF Ai FUNCTION



#####################
##Biomass-inclusive Ascendency (Ab)
#-->only based on intercompartmental exchanges

Asc.biomass<-function(Z_cs, E_cs, R_cs, T_cs, biomass, partial.matrix=TRUE){
  
  Tdotdot<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)#total system throughput
  Bdot<-sum(biomass)#total biomass
  whichflows<-which(T_cs!=0)
  flows<-T_cs[which(T_cs!=0)]
  
  if(partial.matrix==TRUE){
  Tstar.ab<-matrix(0,nrow=nrow(T_cs),ncol=ncol(T_cs), 
                   dimnames =list(rownames(T_cs), colnames(T_cs)))}
  if(partial.matrix==FALSE){
    B.matrix.row<-matrix(rep(biomass, ncol(T_cs)),nrow=nrow(T_cs),ncol=ncol(T_cs),
                              byrow=FALSE)
    B.matrix.col<-matrix(rep(biomass, ncol(T_cs)),nrow=nrow(T_cs),ncol=ncol(T_cs),
                         byrow=TRUE)
    Ab.FINAL<-0}
  
  if(length(flows)>0){
    for(i in 1:length(flows)){
      
    if(partial.matrix==TRUE){
      Tstar.ab[whichflows[i]]<-NA
      row<-which(is.na(Tstar.ab) == TRUE, arr.ind=TRUE)[1]
      col<-which(is.na(Tstar.ab) == TRUE, arr.ind=TRUE)[2]
      if(biomass[row]==0){
        Tstar.ab[row, col]<-0}
      else if(biomass[col]==0){
        Tstar.ab[row, col]<-0}
      else{Ab<-T_cs[row,col]*
        log2((T_cs[row,col]*(Bdot^2))/(Tdotdot*biomass[row]*biomass[col]))
      Tstar.ab[row,col]<-Ab }}
      
    if(partial.matrix==FALSE){
      Ab<-T_cs[whichflows[i]]*
        log2((T_cs[whichflows[i]]*(Bdot^2))/(Tdotdot*B.matrix.row[whichflows[i]]*
                                         B.matrix.col[whichflows[i]]))
      Ab.FINAL<-Ab.FINAL+Ab }  }}
  
  if(partial.matrix==TRUE){
  list_ab<-list(sum(Tstar.ab), Tstar.ab)
  names(list_ab)<-c("Biomass-inclusive Ascendency (Ab)","Partial Ab Matrix")
  return(list_ab)}
  
  if(partial.matrix==FALSE){return(Ab.FINAL)}
}#END OF Ab FUNCTION

########################################################################################
#CONNECTANCE INDICES####################################################################
########################################################################################

##Overall/intercompartmental/foodweb connectance (=after Ulanowicz & Wolff, 1991)

connectance<-function(Z_cs=NA, E_cs=NA, R_cs=NA, T_cs, nl=NA, type, k=1){
  
  if(type=="whole"){
    #-->takes into account non-living nodes and exchanges with surroundings (i.e., import,
    #   export and respiration vectors)  
    T_matrix<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
    T_matrix<-cbind(T_matrix,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))#define matrix of transfers
    rownames(T_matrix)<-colnames(T_matrix)<-c(rownames(T_cs), "EXPORT", "RESP", "IMPORT")
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="intercompartmental"){
    #-->takes into account only exchanges between compartments (living and non-living)
    T_matrix<-T_cs #define matrix of transfers
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="foodweb"){
    #-->takes into account only exchanges between living compartments  
    #-->requires additional information about number of non-living compartments (nl)
    living<-nrow(T_cs)-nl
    T_matrix<-T_cs[c(1:living), c(1:living)] #define matrix of transfers
    comps<-living} 
  
  if(living<=1){
    Hc.final<-0
  }else{
    
    #1.Calculate conditional diversity Hc
    Tdotdot<-sum(T_matrix)
    totOUT<-apply(T_matrix,1,sum)
    totIN<-apply(T_matrix,2,sum)
    
    #for intercompartmental exchanges:
    #-->overall, intercompartmental, foodweb connectance
    Hc.matrix<-matrix(0, nrow=comps, ncol=comps)
    T_matrix.inter<-T_matrix[c(1:comps),c(1:comps)]
    
    whichflows<-which(T_matrix.inter!=0)
    flows<-T_matrix.inter[which(T_matrix.inter!=0)]
    
    if(length(flows)>0){
      for(i in 1:length(flows)){
        Hc.matrix[whichflows[i]]<-NA
        row<-which(is.na(Hc.matrix) == TRUE, arr.ind=TRUE)[1]
        col<-which(is.na(Hc.matrix) == TRUE, arr.ind=TRUE)[2]
        
        if(totIN[col]==0 && totOUT[row]==0 ||
           totIN[col]==0 || totOUT[row]==0){
          Hc<-0}
        else{
          Hc<-(-1)*k*(T_matrix.inter[row,col]/Tdotdot)*
            log2((T_matrix.inter[row,col]^2)/(totIN[col]*totOUT[row]))}
        
        Hc.matrix[whichflows[i]]<-Hc}}
    
    if(type=="intercompartmental" || type=="foodweb"){
      Hc.final<-sum(Hc.matrix)}
    
    if(type=="whole"){
      
      #Calculate Hc of import/export/dissipation:
      #-->only for overall connectance
      Himp.vec<-rep(0,length(Z_cs))
      Hexp.vec<-rep(0,length(Z_cs))
      Hresp.vec<-rep(0,length(Z_cs))
      for(i in 1:length(Z_cs)){
        if(totIN[i]==0 && totOUT[i]!=0){
          Hexp.vec[i]<-(-2)*k*(E_cs[i]/Tdotdot)*log2(E_cs[i]/totOUT[i])
          Hresp.vec[i]<-(-2)*k*(R_cs[i]/Tdotdot)*log2(R_cs[i]/totOUT[i])}
        else if(totIN[i]!=0 && totOUT[i]==0){
          Himp.vec[i]<-(-2)*k*(Z_cs[i]/Tdotdot)*log2(Z_cs[i]/totIN[i])}
        else if(totIN[i]!=0 && totOUT[i]!=0){
          Himp.vec[i]<-(-2)*k*(Z_cs[i]/Tdotdot)*log2(Z_cs[i]/totIN[i])
          Hexp.vec[i]<-(-2)*k*(E_cs[i]/Tdotdot)*log2(E_cs[i]/totOUT[i])
          Hresp.vec[i]<-(-2)*k*(R_cs[i]/Tdotdot)*log2(R_cs[i]/totOUT[i])} }
      Himp.vec[which(is.na(Himp.vec)==TRUE)]<-0
      Hexp.vec[which(is.na(Hexp.vec)==TRUE)]<-0
      Hresp.vec[which(is.na(Hresp.vec)==TRUE)]<-0
      
      Hc.final<-sum(Hc.matrix)+sum(Himp.vec)+sum(Hexp.vec)+sum(Hresp.vec)}
  }#END OF ELSE CLAUSE (=living>1)
  
  #2.Compute connectance index
  connectance<-2^(Hc.final/2)
  
  return(connectance)
}#END OF connectance FUNCTION


####################
##Weighted Connectance, Zorach and Ulanowicz (2003)
eff.C<-function(Z_cs, E_cs, R_cs, T_cs){
  
  Tdotdot<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)
  totOUT<-apply(T_cs,1,sum)+E_cs+R_cs
  totIN<-apply(T_cs,2,sum)+Z_cs
  effC.matrix<-matrix(1, nrow=nrow(T_cs), ncol=ncol(T_cs))
  
  whichflows<-which(T_cs!=0)
  flows<-T_cs[which(T_cs!=0)]
  
  if(length(flows)>0){
    for(i in 1:length(flows)){
      effC.matrix[whichflows[i]]<-NA
      row<-which(is.na(effC.matrix) == TRUE, arr.ind=TRUE)[1]
      col<-which(is.na(effC.matrix) == TRUE, arr.ind=TRUE)[2]
      
      if(totIN[col]==0 && totOUT[row]==0 ||
         totIN[col]==0 || totOUT[row]==0){effC.matrix[whichflows[i]]<-1}
      else{
        effC<-((T_cs[row,col]^2)/(totIN[col]*totOUT[row]))^
          ((-0.5)*(T_cs[row,col]/Tdotdot))
        effC.matrix[whichflows[i]]<-effC } }}
  
  effC.final<-prod(effC.matrix)
  return(effC.final)
}#END OF eff.C FUNCTION

####################
##Link Density, Bersier et al. (2002)
LDq<-function(Z_cs, E_cs, R_cs, T_cs){
  
  Tdotdot<-sum(T_cs)+sum(Z_cs)+sum(E_cs)+sum(R_cs)
  totOUT<-apply(T_cs,1,sum)+E_cs+R_cs
  totIN<-apply(T_cs,2,sum)+Z_cs
  
  #1.Calculate diversity of in- and outflows
  Hin.matrix<-matrix(0, ncol=ncol(T_cs), nrow=nrow(T_cs))
  Hout.matrix<-matrix(0, ncol=ncol(T_cs), nrow=nrow(T_cs))
  
  whichflows<-which(T_cs!=0)
  flows<-T_cs[which(T_cs!=0)]
  
  if(length(flows)>0){
    for(i in 1:length(flows)){
      Hin.matrix[whichflows[i]]<-NA
      row<-which(is.na(Hin.matrix) == TRUE, arr.ind=TRUE)[1]
      col<-which(is.na(Hin.matrix) == TRUE, arr.ind=TRUE)[2]
      
      if(totIN[col]==0 && totOUT[row]==0){
        Hin.matrix[whichflows[i]]<-0}
      else if(totIN[col]==0 && totOUT[row]!=0){
        Hin.matrix[whichflows[i]]<-0
        Hout<-(-1)*(T_cs[row,col]/totOUT[row])*log2(T_cs[row,col]/totOUT[row])
        Hout.matrix[whichflows[i]]<-Hout}
      else if(totIN[col]!=0 && totOUT[row]==0){
        Hin<-(-1)*(T_cs[row,col]/totIN[col])*log2(T_cs[row,col]/totIN[col])
        Hin.matrix[whichflows[i]]<-Hin }
      else{
        Hin<-(-1)*(T_cs[row,col]/totIN[col])*log2(T_cs[row,col]/totIN[col])
        Hin.matrix[whichflows[i]]<-Hin 
        Hout<-(-1)*(T_cs[row,col]/totOUT[row])*log2(T_cs[row,col]/totOUT[row])
        Hout.matrix[whichflows[i]]<-Hout } }}
  
  H.IN<-apply(Hin.matrix,2,sum)
  H.OUT<-apply(Hout.matrix,1,sum)
  
  #2.Calculate link density
  LD.in<-rep(0,ncol(T_cs))
  LD.out<-LD.in
  for(i in 1:length(LD.in)){
    LD.in[i]<-(totIN[i]/Tdotdot)*2^H.IN[i]
    LD.out[i]<-(totOUT[i]/Tdotdot)*2^H.OUT[i]}
  LDq<-0.5*(sum(LD.in)+sum(LD.out))
  return(LDq)
}#END OF LDq FUNCTION

################
##Average Mutual Information relative to flow contribution entering/exiting each nod
##(1)calculate based on incoming links
##(2)calculated based on outgoing links
Ami.relative<-function(Z_cs=NA, E_cs=NA, R_cs=NA, T_cs, type, nl=NA, k=1){
  
  if(type=="whole"){
    #-->takes into account non-living nodes and exchanges with surroundings (i.e., import,
    #   export and respiration vectors)  
    T_matrix<-rbind(T_cs, rep(0, ncol(T_cs)), rep(0, ncol(T_cs)), Z_cs)
    T_matrix<-cbind(T_matrix,c(E_cs,rep(0,3)), c(R_cs,rep(0,3)), rep(0,ncol(T_cs)+3))#define matrix of transfers
    rownames(T_matrix)<-colnames(T_matrix)<-c(rownames(T_cs), "EXPORT", "RESP", "IMPORT")
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="intercompartmental"){
    #-->takes into account only exchanges between compartments (living and non-living)
    T_matrix<-T_cs #define matrix of transfers
    comps<-nrow(T_cs)
    living<-2}
  
  if(type=="foodweb"){
    #-->takes into account only exchanges between living compartments  
    #-->requires additional information about number of non-living compartments (nl)
    living<-nrow(T_cs)-nl
    T_matrix<-T_cs[c(1:living), c(1:living)] #define matrix of transfers
    comps<-living} 
  
  if(living>1){
    
    #Compute relative AMI for in- and outflows of T_matrix:
    Tdotdot<-sum(T_matrix)#total system throughput of T_matrix
    totIN<-apply(T_matrix, 2,sum)#sum of ingoing flows per nod
    totOUT<-apply(T_matrix,1,sum)#sum of outgoing flows per nod
    whichflows<-which(T_matrix!=0)
    flows<-T_cs[which(T_matrix!=0)]
    T_matrix.amiIN<-matrix(0,nrow=nrow(T_matrix),ncol=ncol(T_matrix), 
                           dimnames =list(rownames(T_matrix),colnames(T_matrix))) 
    T_matrix.amiOUT<-T_matrix.amiIN
    
    if(length(flows)>0){    
      for(i in 1:length(flows)){
        T_matrix.amiIN[whichflows[i]]<-NA
        row<-which(is.na(T_matrix.amiIN) == TRUE, arr.ind=TRUE)[1]
        col<-which(is.na(T_matrix.amiIN) == TRUE, arr.ind=TRUE)[2]
        if(totOUT[row]==0){
          T_matrix.amiIN[row, col]<-0}
        else if(totIN[col]==0){
          T_matrix.amiIN[row, col]<-0}
        else{
          T_matrix.amiIN[row, col]<-k*(T_matrix[row,col]/totIN[col])*
            log2((T_matrix[row,col]*Tdotdot)/(totOUT[row]*totIN[col]))
          T_matrix.amiOUT[row, col]<-k*(T_matrix[row,col]/totOUT[row])*
            log2((T_matrix[row,col]*Tdotdot)/(totOUT[row]*totIN[col]))}}}
    
    amiIN<-unname(apply(T_matrix.amiIN,2,sum))
    amiOUT<-unname(apply(T_matrix.amiOUT,1,sum)) }
  
  #Function output:
  
  if(living==0){
    output<-data.frame(compartment.ID=NA, AMI.inflows=NA, 
                       AMI.outflows=NA)}
  if(living==1 && T_matrix==0){
    output<-data.frame(compartment.ID=rownames(T_cs)[1], AMI.inflows=NA, 
                       AMI.outflows=NA)}
  if(living==1 && T_matrix>0){ 
    output<-data.frame(compartment.ID=rownames(T_cs)[1], AMI.inflows=0, 
                       AMI.outflows=0)
  }else{
    output<-data.frame(compartment.ID=rownames(T_matrix)[c(1:comps)], 
                       AMI.inflows=amiIN[c(1:comps)], 
                       AMI.outflows=amiOUT[c(1:comps)])}
  
  return(output)
}#END OF RELATIVE AMI FUNCTION
