###########BALANCING PROCEDURE#########

#NOTE: include exceptions in case of length(Tstar)=1 (???)

network.balance<-function(Z_cs, E_cs, R_cs, T_cs, method){
  
  
##1.function to convert Z_cs, E_cs, R_cs, T_cs in extended transfer matrix:

extTM<-function(Z_cs, E_cs, R_cs, T_cs){
  emptyrow<-rep(0, length(Z_cs))
  ZT_cs<-rbind(T_cs, emptyrow, emptyrow, Z_cs)
  exp<-c(E_cs, rep(0,3))
  resp<-c(R_cs, rep(0,3))
  emptycol<-rep(0, length(Z_cs)+3)
  Tstar<-cbind(ZT_cs, exp, resp, emptycol)
  rownames(Tstar)<-colnames(Tstar)<-c(colnames(T_cs), "EXPORT", "RESP", "IMPORT")
  return(Tstar) }

##2.Produce the extended matrix of network (T*):
Tstar<-extTM(Z_cs, E_cs, R_cs, T_cs)



##3.Check steady state of system

INminusOUT<-unname(apply(Tstar,1,sum)-apply(Tstar,2,sum))[c(1:ncol(T_cs))]
INminusOUT<-round(INminusOUT, digits=9)#values are rounded till 9th digit

####SYSTEM IS AT STEADY_STATE

if(length(which(INminusOUT==0))==length(INminusOUT)){ 
  print("Network is balanced. No balancing procedure was applied.")
        list_bal<-NULL
        return(list_bal)
  
####SYSTEM IS NOT AT STEADY STATE
}else{ 

##4.Produce summary of unbalances:

nods<-colnames(T_cs)
nodIN.sum<-unname(apply(Tstar, 1, sum))[c(1:length(nods))]
nodOUT.sum<-unname(apply(Tstar, 2, sum))[c(1:length(nods))]
nod.DIF<-nodIN.sum - nodOUT.sum
imbalSUM<-data.frame(compartment_ID=nods, IN.sum=nodIN.sum, OUT.sum=nodOUT.sum, Diff=nod.DIF)


##5.Produce functions to apply a network balancing approach:

########INPUT-BASED APPROACH###########

INP.bal <- function(Tstar, flow=Z_cs){
  
  totnods<-ncol(Tstar)-3 #total number of nods of the network
  
#Step1: divide each coef tij with i[1,...,N] & j[1,...,N+3] by the i-th row sum
#       and store results in matrix F*
 Fstar<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
               dimnames=list(rownames(Tstar),colnames(Tstar)))
  OUTsum<-apply(Tstar[c(1:totnods),], 1, sum)
  for(t in 1:length(OUTsum)){
    if(OUTsum[t]!=0)Fstar[t,]<- Tstar[t,] / OUTsum[t] }
#NOTE: problem if rowsums[t]=0!!

#Step2: transpose the NxN part of matrix, substract identity matrix and store results in
#       matrix R
  Imatrix<-diag(totnods)
  R<-t(Fstar[c(1:totnods), c(1:totnods)]) - Imatrix
  
#Step3: invert matrix R
  R<-solve(R)
  #NOTE:  ddd <- det(R)
  #       if(ddd == 0)return(NA) ????????????????????
  
#Step4: multiply each coef rij by the corresponding j-th input of T* and change its sign
  for(j in 1:ncol(R)) R[,j]<- -R[,j] * flow[j]

#Step5: sum i-th row to build vector U
  U<-apply(R,1,sum)

#Step6: multiply each fij by the corresponding ui to obtain a balanced form of T* matrix and
#       and store results in Tstar_balIN
  Tstar_balIN<-matrix(0, nrow=totnods+3, ncol=totnods+3, 
                      dimnames=list(rownames(Tstar),colnames(Tstar)))
  for(i in 1:totnods)Tstar_balIN[i,]<- Fstar[i,] * U[i] 
  for(i in (totnods+1):nrow(Tstar_balIN))Tstar_balIN[i,]<-Tstar[i,]
  

return(Tstar_balIN)

} #END OF INPUT-BASED BALANCING FUNCTION


########OUTPUT-BASED APPROACH###########
#Same method as for output-based approach except:
#-->T* is transposed before Step1 
#-->Step4: multiply each coef rij by the corresponding j-th Output of T* and change its sign
#-->balanced T* is transposed after Step6

OUT.bal <- function(Tstar){
  
  OUT<-Tstar[,"EXPORT"] + Tstar[,"RESP"]
  Tstar.inv<-t(Tstar)
  Tstar_bal<-INP.bal(Tstar=Tstar.inv, flow=OUT)
  Tstar_balOUT<-t(Tstar_bal)
  
return(Tstar_balOUT)
  
}  #END OF OUTPUT-BASED BALANCING FUNCTION


########AVG APPROACH#########################

AVG.bal <- function(Tstar){
  
  Tstar_balIN<-INP.bal(Tstar)
  Tstar_balOUT<-OUT.bal(Tstar)
  Tstar_balAVG<-0.5*(Tstar_balIN+Tstar_balOUT)

  return(Tstar_balAVG)

} #END OF AVG-BASED BALANCING FUNCTION


########IO APPROACH##########################

IO.bal <- function(Tstar){
  
  Tstar_balIN<-INP.bal(Tstar)
  Tstar.io<-0.5*(Tstar_balIN+Tstar)
  Tstar_balIO<-OUT.bal(Tstar.io)
  
  return(Tstar_balIO)
  
} #END OF IO-BASED BALANCING FUNCTION


########OI APPROACH##########################

OI.bal <- function(Tstar){
  
  Tstar_balOUT<-OUT.bal(Tstar)
  Tstar.oi<-0.5*(Tstar_balOUT+Tstar)#CORRECT
  INP.oi<-Tstar.oi["IMPORT",]
  Tstar_balOI<-INP.bal(Tstar.oi, flow=INP.oi)
  
  return(Tstar_balOI)
  
} #END OF IO-BASED BALANCING FUNCTION


########AVG2 APPROACH##########################

AVG2.bal <- function(Tstar){
  
  Tstar_balIO<-IO.bal(Tstar)
  Tstar_balOI<-OI.bal(Tstar)
  Tstar_balAVG2<-0.5*(Tstar_balIO+Tstar_balOI)
  
  return(Tstar_balAVG2)
  
} #END OF AVG-BASED BALANCING FUNCTION


##6.Apply a balancing procedure depending on the chosen method:

if(method == "inp"){
  T_star_bal<-INP.bal(Tstar) }

if(method == "out"){
  T_star_bal<-OUT.bal(Tstar) }

if(method == "avg"){
  T_star_bal<-AVG.bal(Tstar) }
  
if(method == "io"){
  T_star_bal<-IO.bal(Tstar) }

if(method == "oi"){
  T_star_bal<-OI.bal(Tstar) }  

if(method == "avg2"){
  T_star_bal<-AVG2.bal(Tstar) }

##7.Function Output:
#-->list of (1)summary of imbalances, (2)balanced Z_cs, (3)balanced E_cs,
#             (4)balanced R_cs, (5)balanced T_cs

Z_cs_bal<-T_star_bal[nrow(T_star_bal), c(1:length(nods))] #import
T_cs_bal<-T_star_bal[c(1:length(nods)), c(1:length(nods))] #matrix of intercompartmental exchanges
E_cs_bal<-T_star_bal[c(1:length(nods)), length(nods)+1] #export
R_cs_bal<-T_star_bal[c(1:length(nods)), length(nods)+2] #respiration

list_bal<-list(imbalSUM, T_star_bal, Z_cs_bal, E_cs_bal, R_cs_bal, T_cs_bal)
names(list_bal)<-c("Summary of imbalances", "Extended Transfer Matrix", "Import", "Export", 
                   "Respiration", "Intercompartmental exchanges")


return(list_bal)


} #END OF ELSE CLAUSE
} #END OF FUNCTION



