setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/SCOR format")

###CONVERSION OF SCOR-FORMAT NETWORK INTO:
#####-->information of living/total number of compartments
#####-->vector of compartment biomass
#####-->vector of import flows (Z_cs)
#####-->vector of export flows (E_cs)
#####-->vector of energy dissipation flows=respiration (R_cs)
#####-->matrix of intercompartamental exchanges (T_cs) 

SCOR.convert<-function(file, headline=TRUE){

#headline=TRUE/FALSE
if(headline==FALSE){header<-0 #no headline in SCOR file=info is read from 1st line on
}else{header<-1} #Default: there is a headline=info is read from line 2 on

#1.Info of non-living/living compartments:
comps<-scan(file, skip=header, nlines=1)
nb_comps<-comps[which(comps==max(comps))]#Total number of compartments
names(comps)<-c("total", "living")

#2.vector of compartment names
lines_skip<-header+1
names_comps<-scan(file, skip=lines_skip, what="character", nlines=nb_comps, sep="\n")

#3.Import material flows
lines_skip<-lines_skip+nb_comps
materialflows<-scan(file, skip=lines_skip, what="character", sep="\n")

#Extract compartment biomass
biomass<-rep(0, nb_comps)
names(biomass)<-names_comps
repeat{
  biom<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  biom<-biom[!is.na(biom)]
  if(biom[1]<0){
    materialflows<-materialflows[-1]  
    break
  }else{
    biomass[biom[1]]<-biom[2]
    materialflows<-materialflows[-1]} }

#Extract import flows
Z_cs<-rep(0, nb_comps)
names(Z_cs)<-names_comps
repeat{
  imp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  imp_flow<-imp_flow[!is.na(imp_flow)]
  if(imp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break
  }else{
    Z_cs[imp_flow[1]]<-imp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract export flows
E_cs<-rep(0, nb_comps)
names(E_cs)<-names_comps
repeat{
  exp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  exp_flow<-exp_flow[!is.na(exp_flow)]
  if(exp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break
  }else{
    E_cs[exp_flow[1]]<-exp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract respiration flows
R_cs<-rep(0, nb_comps)
names(R_cs)<-names_comps
repeat{
  resp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  resp_flow<-resp_flow[!is.na(resp_flow)]
  if(resp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break
  }else{
    R_cs[resp_flow[1]]<-resp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract intercompartamental exchanges
T_cs<-matrix(0, nrow=nb_comps, ncol=nb_comps, dimnames=list(names_comps, names_comps))
repeat{
  flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  flow<-flow[!is.na(flow)]
  if(flow[1]<0){  
    break
  }else{
    T_cs[flow[1], flow[2]]<-flow[3]
    materialflows<-materialflows[-1]} }

converted.SCOR<-list(comps, biomass, Z_cs, E_cs, R_cs, T_cs)
names(converted.SCOR)<-c("Number of compartments", "Compartment biomass stock", "Import", "Export",
                         "Energy dissipation", "Intercompartmental exchanges")

return(converted.SCOR)

}#END OF FUNCTION



#Test function:
ConeSpring<-SCOR.convert("SCOR_exampleNetwork.txt")
ConeSpring<-SCOR.convert("SCOR_exampleNetwork2.txt", headline=FALSE)
