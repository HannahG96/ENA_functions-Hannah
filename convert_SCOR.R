#setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/SCOR format")

###CONVERSION OF SCOR-FORMAT NETWORK INTO:
#####-->information of living/total number of compartments
#####-->vector of compartment biomass
#####-->vector of import flows (Z_cs)
#####-->vector of export flows (E_cs)
#####-->vector of energy dissipation flows=respiration (R_cs)
#####-->matrix of intercompartamental exchanges (T_cs) 

#NOTE: some foodwebs have 3 elements as Info of non-living/living compartments
# (!!!)FUNCTION DOES NOT WORK ON eg. file=SCOR format/cypwet.dat"!!!!!

SCOR.convert<-function(file, header=1){

errorline<-header #object to store the line that is currently scanned by the function in the SCOR file

#1.Info of non-living/living compartments:
errorline<-errorline+1

comps<-scan(file, skip=header , nlines=1)
nb_comps<-comps[1]#Total number of compartments
names(comps)<-c("total", "living")

#2.vector of compartment names
lines_skip<-header+1
names_comps<-scan(file, skip=lines_skip, what="character", nlines=nb_comps, sep="\n")

errorline<-errorline+length(names_comps)

#3.Import material flows
lines_skip<-lines_skip+nb_comps
materialflows<-scan(file, skip=lines_skip, what="character", sep="\n")

#Extract compartment biomass
biomass<-rep(0, nb_comps)
names(biomass)<-names_comps
repeat{
  errorline<-errorline+1
  
  biom<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  biom<-biom[!is.na(biom)]
  
  if(biom[1]<0){
    materialflows<-materialflows[-1]  
    break}

  ##SCAN FOR ERRORS  
  if(length(biom) != 2 || biom[1] > nb_comps){
    print(paste("ERROR in line", errorline, ": undefined compartments or material flow"))
    materialflows<-materialflows[-1]
    next}
  if(biom[2] < 0){
    print(paste("ERROR in line", errorline, ": negative material flow"))
    materialflows<-materialflows[-1]
    next
  #################
    
  }else {
    biomass[biom[1]]<-biom[2]
    materialflows<-materialflows[-1]} }

#Extract import flows
Z_cs<-rep(0, nb_comps)
names(Z_cs)<-names_comps
repeat{
  errorline<-errorline+1
  
  imp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  imp_flow<-imp_flow[!is.na(imp_flow)]
  if(imp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break}
  
  ##SCAN FOR ERRORS  
  if(length(imp_flow) != 2 || imp_flow[1] > nb_comps){
    print(paste("ERROR in line", errorline, ": undefined compartments or material flow"))
    materialflows<-materialflows[-1]
    next}
  if(imp_flow[2] < 0){
    print(paste("ERROR in line", errorline, ": negative material flow"))
    materialflows<-materialflows[-1]
    next
  #################
    
  }else{
    Z_cs[imp_flow[1]]<-imp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract export flows
E_cs<-rep(0, nb_comps)
names(E_cs)<-names_comps
repeat{
  errorline<-errorline+1
  
  exp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  exp_flow<-exp_flow[!is.na(exp_flow)]
  
  if(exp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break}
  
  ##SCAN FOR ERRORS  
  if(length(exp_flow) != 2 || exp_flow[1] > nb_comps){
    print(paste("ERROR in line", errorline, ": undefined compartments or material flow"))
    materialflows<-materialflows[-1]
    next}
  if(exp_flow[2] < 0){
    print(paste("ERROR in line", errorline, ": negative material flow"))
    materialflows<-materialflows[-1]
    next
  #################
    
  }else{
    E_cs[exp_flow[1]]<-exp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract respiration flows
R_cs<-rep(0, nb_comps)
names(R_cs)<-names_comps
repeat{
  errorline<-errorline+1
  
  resp_flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  resp_flow<-resp_flow[!is.na(resp_flow)]
  
  if(resp_flow[1]<0){
    materialflows<-materialflows[-1]  
    break}
  
  ##SCAN FOR ERRORS  
  if(length(resp_flow) != 2 || resp_flow[1] > nb_comps){
    print(paste("ERROR in line", errorline, ": undefined compartments or material flow"))
    materialflows<-materialflows[-1]
    next}
  if(resp_flow[2] < 0){
    print(paste("ERROR in line", errorline, ": negative material flow"))
    materialflows<-materialflows[-1]
    next
  #################
    
  }else{
    R_cs[resp_flow[1]]<-resp_flow[2]
    materialflows<-materialflows[-1]} }

#Extract intercompartamental exchanges
T_cs<-matrix(0, nrow=nb_comps, ncol=nb_comps, dimnames=list(names_comps, names_comps))
repeat{
  errorline<-errorline+1
  
  flow<-as.numeric(strsplit(materialflows[1], split=" ")[[1]])
  flow<-flow[!is.na(flow)]
  
  if(flow[1]<0){  
    break}
  
  ##SCAN FOR ERRORS  
  if(length(flow) != 3 || flow[1] > nb_comps || flow[2] > nb_comps){
    print(paste("ERROR in line", errorline, ": undefined compartments or material flow"))
    materialflows<-materialflows[-1]
    next}
  if(flow[3] < 0){
    print(paste("ERROR in line", errorline, ": negative material flow"))
    materialflows<-materialflows[-1]
    next
  #################
    
  }else{
    T_cs[flow[1], flow[2]]<-flow[3]
    materialflows<-materialflows[-1]} }

converted.SCOR<-list(errorline, comps, biomass, Z_cs, E_cs, R_cs, T_cs)
names(converted.SCOR)<-c("Imported lines of file","Number of compartments", "Compartment biomass stock", 
                         "Import", "Export", "Energy dissipation", "Intercompartmental exchanges")

#print(errorline)#CONTROL

return(converted.SCOR)

}#END OF FUNCTION





