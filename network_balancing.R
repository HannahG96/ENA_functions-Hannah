###########BALANCING PROCEDURE###########################################################
#including special case: det(matrix)=0 -->no matrix inversion performed ; return NA

network.balance<-function(Z_cs, E_cs, R_cs, T_cs, method){

##1.Convert Z_cs, E_cs, R_cs, T_cs in extended transfer matrix:

  emptyrow<-rep(0, length(Z_cs))
  ZT_cs<-rbind(T_cs, emptyrow, emptyrow, Z_cs)
  exp<-c(E_cs, rep(0,3))
  resp<-c(R_cs, rep(0,3))
  emptycol<-rep(0, length(Z_cs)+3)
  Tstar<-cbind(ZT_cs, exp, resp, emptycol)
  rownames(Tstar)<-colnames(Tstar)<-c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")
  unbalanced.network<-Tstar#store T* in a second matrix-->needed for avg2 approach

##2.Check steady state of system
  
  INminusOUT<-unname(apply(Tstar,1,sum)-apply(Tstar,2,sum))[c(1:ncol(T_cs))]
  INminusOUT<-round(INminusOUT, digits=9)#values are rounded till 9th digit
  
####SYSTEM IS AT STEADY_STATE###########################
  
  if(length(which(INminusOUT==0))==length(INminusOUT)){ 
    print("Network is balanced. No balancing procedure was applied.")
    list_bal<-NULL

    
####SYSTEM IS NOT AT STEADY STATE#######################
  }else{ 
    
##3.Produce summary of unbalances:
    
    nods<-colnames(T_cs)
    nodIN.sum<-unname(apply(Tstar, 1, sum))[c(1:length(nods))]
    nodOUT.sum<-unname(apply(Tstar, 2, sum))[c(1:length(nods))]
    nod.DIF<-nodIN.sum - nodOUT.sum
    imbalSUM<-data.frame(compartment_ID=nods, IN.sum=nodIN.sum, OUT.sum=nodOUT.sum, 
                         Diff=nod.DIF)
    

    
########INPUT-BASED APPROACH###########
    
if(method=="inp" ||
   method=="avg" ||
   method=="io"  ||
   method=="avg2"){
      
      totnods<-ncol(Tstar)-3 #total number of nods of the network
      
#Step1: divide each coef tij with i[1,...,N] & j[1,...,N+3] by the i-th row sum
#       and store results in matrix F*
      Fstar<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                    dimnames=list(rownames(Tstar),colnames(Tstar)))
      OUTsum<-apply(Tstar[c(1:totnods),], 1, sum)
      for(t in 1:length(OUTsum)){
        if(OUTsum[t]!=0)Fstar[t,]<- Tstar[t,] / OUTsum[t] }
      
#Step2: transpose the NxN part of matrix, substract identity matrix and store results in
#       matrix R
      Imatrix<-diag(totnods)
      R<-t(Fstar[c(1:totnods), c(1:totnods)]) - Imatrix
      
#Step3: invert matrix R
     
      ddd <- det(R)
#Include exception of det(matrix)=0:
  if(ddd == 0){ Tstar_balIN<- NA
####      
  }else{
      R<-solve(R)

      
#Step4: multiply each coef rij by the corresponding j-th input of T* and change its sign
      for(j in 1:ncol(R)) R[,j]<- -R[,j] * Z_cs[j]
      
#Step5: sum i-th row to build vector U
      U<-apply(R,1,sum)
      
#Step6: multiply each fij by the corresponding ui to obtain a balanced form of T* matrix and
#       and store results in Tstar_balIN
      Tstar_balIN<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                          dimnames=list(rownames(Tstar),colnames(Tstar)))
      for(i in 1:totnods)Tstar_balIN[i,]<- Fstar[i,] * U[i] 
      for(i in (totnods+1):nrow(Tstar_balIN))Tstar_balIN[i,]<-Tstar[i,]
      
    }#END OF ELSE CLAUSE (Condition: det(matrix)!=0)
      
  if(method=="inp")T_star_bal<-Tstar_balIN
}#END OF INPUT-BASED APPROACH
    
    
    
########OUTPUT-BASED APPROACH###########
    
if(method=="out" ||
   method=="avg" ||
   method=="oi"  ||
   method=="avg2"){
    
  totnods<-ncol(Tstar)-3 #total number of nods of the network  
  
#Step1:transpose matrix T*
     Tstar<-t(Tstar)
      
      
#Step2: divide each coef tij with i[1,...,N] & j[1,...,N+3] by the i-th row sum
#       and store results in matrix F*
      Fstar<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                    dimnames=list(rownames(Tstar),colnames(Tstar)))
      OUTsum<-apply(Tstar[c(1:totnods),], 1, sum)
      for(t in 1:length(OUTsum)){
        if(OUTsum[t]!=0)Fstar[t,]<- Tstar[t,] / OUTsum[t] }
      
#Step3: transpose the NxN part of matrix, substract identity matrix and store results in
#       matrix R
      Imatrix<-diag(totnods)
      R<-t(Fstar[c(1:totnods), c(1:totnods)]) - Imatrix
      
#Step4: invert matrix R
      
      ddd <- det(R)
#Include exception of det(matrix)=0:
  if(ddd == 0){ Tstar_balOUT<- NA
####      
  }else{
        R<-solve(R)
        
        
#Step5: multiply each coef rij by the corresponding j-th output of T* and change its sign
        OUT <- E_cs + R_cs
        for(j in 1:ncol(R)) R[,j]<- -R[,j] * OUT[j]
        
#Step6: sum i-th row to build vector U
        U<-apply(R,1,sum)
        
#Step6: multiply each fij by the corresponding ui to obtain a balanced form of T* matrix and
#       and store results in Tstar_balIN
        Tstar_balout<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                            dimnames=list(rownames(Tstar),colnames(Tstar)))
        for(i in 1:totnods)Tstar_balout[i,]<- Fstar[i,] * U[i] 
        for(i in (totnods+1):nrow(Tstar_balout))Tstar_balout[i,]<-Tstar[i,]
        
#Step7: transpose the balanced matrix
        Tstar_balOUT<-t(Tstar_balout)
        
      }#END OF ELSE CLAUSE (Condition: det(matrix)!=0)
        
  if(method=="out")T_star_bal<-Tstar_balOUT
}#END OF OUTPUT-BASED APPROACH
    

########AVG APPROACH###########
    
if(method=="avg"){

#Include exception of det(matrix)=0:  
if(length(Tstar_balIN)==1){T_star_bal<-NA}
if(length(Tstar_balOUT)==1){T_star_bal<-NA
###  
}else{

#Produce balanced network based on average of input-based and output-based balanced network:
    T_star_bal<-0.5*(Tstar_balIN+Tstar_balOUT)}
}#END OF AVERAGE APPROACH
    

########IO APPROACH###########

if(method=="io" ||
   method=="avg2") {
  
#Include exception of det(matrix)=0:  
if(length(Tstar_balIN)==1){Tstar_balIO<-NA

###  
}else{
  
##Step1: produce extended transfer matrix to be balanced based on average of input-based 
#        balanced network and unbalanced network
  Tstar.io<-0.5*(Tstar_balIN+unbalanced.network)
  OUT <- Tstar.io[,"EXPORT"] + Tstar.io[,"RESP"]
  
##Step2: apply the output-based approach to balance the network
  
############
  totnods<-ncol(Tstar.io)-3 #total number of nods of the network  
  
  #Step1:transpose matrix T*
  Tstar<-t(Tstar.io)
  
  
  #Step2: divide each coef tij with i[1,...,N] & j[1,...,N+3] by the i-th row sum
  #       and store results in matrix F*
  Fstar<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                dimnames=list(rownames(Tstar),colnames(Tstar)))
  OUTsum<-apply(Tstar[c(1:totnods),], 1, sum)
  for(t in 1:length(OUTsum)){
    if(OUTsum[t]!=0)Fstar[t,]<- Tstar[t,] / OUTsum[t] }
  
  #Step3: transpose the NxN part of matrix, substract identity matrix and store results in
  #       matrix R
  Imatrix<-diag(totnods)
  R<-t(Fstar[c(1:totnods), c(1:totnods)]) - Imatrix
  
  #Step4: invert matrix R
  
  ddd <- det(R)
  #Include exception of det(matrix)=0:
  if(ddd == 0){ Tstar_balIO<- NA
  ####      
  }else{
    R<-solve(R)
    
    
  #Step5: multiply each coef rij by the corresponding j-th output of T* and change its sign

    for(j in 1:ncol(R)) R[,j]<- -R[,j] * OUT[j]
    
  #Step6: sum i-th row to build vector U
    U<-apply(R,1,sum)
    
  #Step7: multiply each fij by the corresponding ui to obtain a balanced form of T* matrix and
  #       and store results in Tstar_balIN
    Tstar_balio<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                         dimnames=list(rownames(Tstar),colnames(Tstar)))
    for(i in 1:totnods)Tstar_balio[i,]<- Fstar[i,] * U[i] 
    for(i in (totnods+1):nrow(Tstar_balio))Tstar_balio[i,]<-Tstar[i,]
    
  #Step8: transpose the balanced matrix
    Tstar_balIO<-t(Tstar_balio)

  }#END OF ELSE CLAUSE (Condition: det(matrix)!=0)
########  
  }#END OF OUTER ELSE CLAUSE (Condition: det(matrix)!=0)
  
if(method=="io")T_star_bal<-Tstar_balIO
}#END OF IO APPROACH
    

########OI APPROACH###########
    
if(method=="oi" ||
   method=="avg2") {
      
#Include exception of det(matrix)=0:  
if(length(Tstar_balOUT)==1){Tstar_balOI<-NA
      
###  
}else{
  
##Step1: produce extended transfer matrix to be balanced based on average of output-based 
#        balanced network and unbalanced network
Tstar.oi<-0.5*(Tstar_balOUT+unbalanced.network)
  IN <- Tstar.oi["IMPORT",] 
  
##Step2: apply the input-based approach to balance the network
  Tstar<-Tstar.oi
############
  totnods<-ncol(Tstar)-3 #total number of nods of the network
  
  #Step1: divide each coef tij with i[1,...,N] & j[1,...,N+3] by the i-th row sum
  #       and store results in matrix F*
  Fstar<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                dimnames=list(rownames(Tstar),colnames(Tstar)))
  OUTsum<-apply(Tstar[c(1:totnods),], 1, sum)
  for(t in 1:length(OUTsum)){
    if(OUTsum[t]!=0)Fstar[t,]<- Tstar[t,] / OUTsum[t] }
  
  #Step2: transpose the NxN part of matrix, substract identity matrix and store results in
  #       matrix R
  Imatrix<-diag(totnods)
  R<-t(Fstar[c(1:totnods), c(1:totnods)]) - Imatrix
  
  #Step3: invert matrix R
  
  ddd <- det(R)
  #Include exception of det(matrix)=0:
  if(ddd == 0){ Tstar_balOI<- NA
  ####      
  }else{
    R<-solve(R)
    
    
    #Step4: multiply each coef rij by the corresponding j-th input of T* and change its sign
    for(j in 1:ncol(R)) R[,j]<- -R[,j] * IN[j]
    
    #Step5: sum i-th row to build vector U
    U<-apply(R,1,sum)
    
    #Step6: multiply each fij by the corresponding ui to obtain a balanced form of T* matrix and
    #       and store results in Tstar_balIN
    Tstar_balOI<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                        dimnames=list(rownames(Tstar),colnames(Tstar)))
    for(i in 1:totnods)Tstar_balOI[i,]<- Fstar[i,] * U[i] 
    for(i in (totnods+1):nrow(Tstar_balOI))Tstar_balOI[i,]<-Tstar[i,]
    
  }#END OF INNER ELSE CLAUSE (Condition: det(matrix)!=0)
############
  
}#END OF OUTER ELSE CLAUSE (Condition: det(matrix)!=0)

if(method=="oi")T_star_bal<-Tstar_balOI
}#END OF OI APPROACH

    
########AVG2 APPROACH###########
    
if(method=="avg2"){
      
#Include exception of det(matrix)=0:  
if(length(Tstar_balIO)==1){T_star_bal<-NA}
if(length(Tstar_balOI)==1){T_star_bal<-NA
###  
}else{
        
#Produce balanced network based on average of input-based and output-based balanced network:
        T_star_bal<-0.5*(Tstar_balIO+Tstar_balOI)}
}#END OF AVERAGE2 APPROACH

  
########################################################################  
#Include exception of det(matrix)=0:
if(length(T_star_bal)==1){
  print("Error in network balancing: det(X)=0")#print error message
  list_bal<-NULL
  
}else{  
##FUNCTION OUTPUT:
#-->list of (1)summary of imbalances, (2)balanced Z_cs, (3)balanced E_cs,
#             (4)balanced R_cs, (5)balanced T_cs
  
  Z_cs_bal<-T_star_bal[nrow(T_star_bal), c(1:length(nods))] #import
  T_cs_bal<-T_star_bal[c(1:length(nods)), c(1:length(nods))] #matrix of intercompartmental exchanges
  E_cs_bal<-T_star_bal[c(1:length(nods)), length(nods)+1] #export
  R_cs_bal<-T_star_bal[c(1:length(nods)), length(nods)+2] #respiration
  
  list_bal<-list(imbalSUM, T_star_bal, Z_cs_bal, E_cs_bal, R_cs_bal, T_cs_bal)
  names(list_bal)<-c("Summary of imbalances", "Extended Transfer Matrix", "Import", "Export", 
                     "Respiration", "Intercompartmental exchanges")
}#END OF ELSE CLAUSE (Condition:det(matrix)!=0)
    
} #END OF OUTER ELSE CLAUSE=no steady-state condition

  return(list_bal)  
  
} #END OF FUNCTION