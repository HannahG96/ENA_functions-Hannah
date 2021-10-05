setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")
######
###Test SCOR conversion#######
######

#Test function on simple network:
ConeSpring<-SCOR.convert("Network data/SCOR_exampleNetwork.txt")

#Test function on a file with ERRORS:
cypwet_ERRORS<-SCOR.convert("Network data/cypwet_withERRORS.dat") #7 introduced errors

#Test function on 48 foodwebs stored as SCORs.
#All foodwebs are stored in the same file "DATALL.dat", so the function is wrapped in a FOR-loop that
#jumps from 1 foodweb to the next within the .dat-file

##################
SCOR.convert.multiple<- function(file, nwebs){
  
newheader<-1
foodwebs<-vector("list", length=nwebs)#list to store function output for each food web
NAMES<-rep("", nwebs)#vector to store name of each foodweb

for(i in 1:nwebs){

#convert SCOR format of foodweb:
conv.fw<-SCOR.convert(file, header=newheader)

#store converted foodweb in list "foodwebs":
foodwebs[[i]]<-conv.fw

#store foodweb name:
foodwebname<-scan(file, skip=newheader-1, what="character", nmax=1, sep="\n")
NAMES[i]<-foodwebname

#define line where starts the next foodweb:
#linesSCOR<-sum(2 , #header+info of living/non-living comps
               #foodwebs[[i]][[1]][1] * 2 , #compartment names+compartment biomass
              # 1 ,
              # length(which(foodwebs[[i]][[3]] !=0)) , #import flows
              # 1 ,
              #length(which(foodwebs[[i]][[4]] !=0)) , #export flows
              # 1 ,
               #length(which(foodwebs[[i]][[5]] !=0)) , #respiration flows
               #1 ,
               #length(which(foodwebs[[i]][[6]] !=0)) , #intercompartmental exchanges
               #1 )
newheader<-conv.fw[[1]]+1

} #END OF FOR LOOP

names(foodwebs)<-NAMES #attribute a name to each foodweb

return(foodwebs)

} #END OF FUNCTION
################

foodwebs<-SCOR.convert.multiple("Network data/DATALL.dat", nwebs=48)
##Encountered errors in SCOR files:
#-->no space between nod IDs of intercompartmental exchange eg. (10124 instead of 10 124)
#-->length of biomass vector unequals number of compartments (biomass stock of 1 compartement forgotten?)
#-->some flows of value zero are noted down in file
#-->flow value exponent is not well defined eg. .3480000+04 instead of .3480000E+04

#Re-test the function on file WITHOUT errors:
foodwebs.corrected<-SCOR.convert.multiple("Network data/DATALL_corrected.dat", nwebs=48)


