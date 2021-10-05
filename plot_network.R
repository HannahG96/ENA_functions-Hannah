###########################FUNCTION TO VISUALIZE NETWORKS######################################
#-->The visualization requires download of the package `DiagrammeR´.

###
#FUNCTION INPUT
#-->Input data is the matrix of intercompartmental exchanges (T_cs), Import-(Z_cs),
#   Export-(E_cs) and Respiration-(R_cs) vectors
#-->number of non-living compartments needs to be specified (nl)
##Function Options
#-->include.detritus: if TRUE, non-living compartments are plotted on lowest y-position of graph
                     #if FALSE, non-living compartments are not shown in graph (less messy)
#-->selfloop: TRUE/FALSE to indicate wether intra-species trophic interactions should be shown
#-->TP.method: Numeric value (1;2) to define the method of trophic position calculation for nod
              #ordering along y-axis. 
              #If 1= TPs are calculated based on trophic analysis (reduces function speed in case 
              #of large networks; only applicable if nl>0)
              #If 2= TPs are calculated based on identity and partial-feeding matrix -> solve(I-G), 
              #then sum columns (fast method; qualitatively weaker)
#-->tSim: numeric value between 0 and 1. Indicates the lower  trophic similarity-threshold that defines
         #nods groups. These are jointly plotted on x-axis. Default value is 0.6 (60%).
         #ONLY TO BE INDICATED IF XY-coordinates are defined automatically. 
#-->select.xy: optional definition of nod coordinates by user. Input is a list of length 2. First
              #list element is the vector of x-coordinates, second element the vector of y-coordinates
              #(vectors of same length than numbers of compartments to be shown in graph, 
              #X-coordinate order should follow the compartment order as in Transfer matrix)
#-->select.colors: optional definition of nod and edge colors. Input is a vector of same length than
                  #numbers of compartments to be shown in graph. Color order should follow the 
                  #compartment order as in Transfer matrix
#-->select.Attr: optional definition of nod/edge/graph attributes. Input is a dataframe of 3
                #columns and n rows as the number of attributes to be defined. 
                #First column indicates the attribute(s) (in strings), 2nd column indicates
                #the attribute value(s), 3rd column indicates the attribute type (edge/nod/graph,
                #in strings).IMPORTANT: Set column names of dataframe a c("attr","value","attr_type")
                #For further infos, see:
#-->labels: optional definition of nod labels by user. Input is a vector of same length than
           #numbers of compartments to be shown in graph. Color order should follow the 
           #compartment order as in Transfer matrix. 
           #"Default" option - nod labels are compartment names (labels exceeding nod borders, and
                              #label overlaps can be expected)
           #"numeric" option - nods are numbered from 1 to number of compartments to be shown following
                              #the order as in transfer matrix
#-->nodsize: optional definition of nod dimensions by user. Numeric value between 0 and 1. Nodsize defines
            #width and height of nod. If biomass proportionality is specialized nodsize is multiplicated
            #by biomass and a conversion coefficient.
            #=option helps to adjust overlaping nods.
#-->biomass: if nod biomasses are indicated (vector of same length than number of compartments to
            #to be shown in graph, biomass order should follow compartment order as in Transfer matrix)
            #nod dimensions are plotted proportional to their biomass stock

###
#NOTE: Thickness of edges is always plotted proportional to the edge flow value.

###
#AUTOMATED ATTRIBUTION OF XY COORDINATES
#-->y-coordinate [0,9]: nods are ordered along y axis according to their trophic position
                       #Primary Producer are always on y position 1
                       #Non living compartments are always on y position 0    
#-->x-coordinate[0,12]: nods are ordered along x axis according to trophic similarity
   #(i)definition of high trophic similar nod groups=we select the nod of highest median
   #trophic similarity, we select all compartments that share a trophic similarity of tSim (=0.6)
   #with given nod. For each of these compartments we then search across all unselected compartments
   #and select those that share a trophic similarity of tSim with respective compartment. We 
   #define all selected compartments as a group. We remove them from the set of unselected 
   #compartments. We start from beginning to form the next group.
   #(ii)Nod groups are ordered along x-axis according to average between-group trophic similarity.
   #=we select the group with largest number of edges, we search for the 2 nod groups that share the
   #highest average trophic similarity with given nod group. We position them left and right on
   #x axis around given group. For each of these we search for the nod group that is highest trophic 
   #similar and position it next to it....
   #(iii)we define x and y coordinates for all members of each nod group and attribute members of a
   #same group a same color

###
#FUNCTION OUTPUT
#-->output is a diagrammeR graph-object that can be visualized via render_graph(graph object)
################################################################################################

plot.network<-function(Z_cs,T_cs,E_cs,R_cs,nl,include.detritus=FALSE, selfloop=FALSE, TP.method=1,
                        tSim=0.6, select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                        nodsize=NA, biomass=NA){
  
library("DiagrammeR")

#number of nods
comps<-ncol(T_cs)
living<-comps-nl#living compartments

#intercompartmental transfers of living compartments
T_cs.living<-T_cs[c(1:living),c(1:living)]

#nodnames
nod.names<-rownames(T_cs)
#names.living<-nod.names[1:living]#living compartments
if(nl>0)names.nl<-nod.names[(living+1):comps]#names of nonliving compartments

####################################################################################
###A.NOD POSITION AND NOD COLORGROUPS####
####################################################################################

#####################
#I-COLORS##
#####################
#-->we predefine 20 colors for single nods and 20 for nod groups
#-->if the graph consists of more groups: 
#remaining single nods get "azure3" color
#remaining nod groups get "azure4" color
#Info message to tell user that color limit is reached is printed

#SINGLE NODS=pastel, light colors
colS<-c("lightcoral","khaki","sandybrown","lightsteelblue","palegreen","thistle",
        "yellow3","tan1","plum","grey","salmon","lightskyblue3","olivedrab1","orange",
        "mediumaquamarine","seashell3","rosybrown","burlywood3","darkolivegreen4","cornflowerblue")
#NOD GROUPS=bright colors
colL<-c("gold","turquoise3","forestgreen","royalblue","mediumorchid","red","yellow","pink",
        "lawngreen","chocolate3","magenta","limegreen","aquamarine4","darkgoldenrod4",
        "firebrick3","lightseagreen","mediumvioletred","chartreuse","darkred","darkslateblue")


###########################
#II-X/Y COORDINATES##
###########################
#-->x-axis: [0,12]
#-->y-axis: [0,9]

###OPTION1: automated attribution of XY coordinates 
if(length(select.XY)==1){
  
#Create empty dataframes to store nod information needed for diagrammeR plotting:
nod_info<-data.frame(Nod_ID=nod.names, x.position=NA, y.position=NA, color=NA)

#Define nod information for nonliving compartments:
#-->Position Y of 0, right-skewed
#-->color: "coral4"
if(nl>0){
for(i in 1:length(names.nl)){
  nodID<-names.nl[i]
  if(i==1){posX<-6}
  else{posX<-posX+6/nl}
  nod_info[which(nod_info[,"Nod_ID"]==nodID),c(2:ncol(nod_info))]<-c(posX,0,"coral4")}}

#Define nod information for Primary Producer:
#-->Position Y of 1
#-->color: "seagreen"
PP<-nod.names[which(apply(T_cs,2,sum)==0)]
if(include.detritus==FALSE){PPX<-12}else{PPX<-6}
for(i in 1:length(PP)){
  nodID<-PP[i]
  if(i==1){posX<-0}
  else{posX<-posX+PPX/length(PP)}
  nod_info[which(nod_info[,"Nod_ID"]==nodID),c(2:ncol(nod_info))]<-c(posX,1,"seagreen")}


####Define nod information of all other living compartments:
#-->Positions of Y defined between 2 to 9 (with max of 9 y-axis levels=nTP.itv)

#1.Calculate trophic properties to facilitate grouping of nods
INTERMs<-which(is.na(nod_info[,"x.position"])==TRUE)
interms<-length(INTERMs)

##Trophic position
if(nl==0)TP.method<-2
if(TP.method==1){
#via TROPHIC ANALYSIS
#-->slow if network is large; only applicable if nl>0
TP<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)[[2]]
}else{
#DIRTY METHOD: TP=solve(I-G) and then sum columns
#-->faster, BUT(!!!) qualitatively weaker
#-->TO REVIEW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#(i)Produce identity matrix (I)
I.matrix<-matrix(rep(0,comps^2),ncol=comps,dimnames=list(nod.names,nod.names))
pos<-c(1:(comps^2))
pos.diag<-which(pos%%(comps+1)==1)
I.matrix[pos.diag]<-1
#(ii)Produce partial feeding matrix (G)
G.matrix<-matrix(rep(0,comps^2),ncol=comps,dimnames=list(nod.names,nod.names))
IN<-apply(T_cs,2,sum)-Z_cs
for (i in 1:comps){
  for (j in 1:comps)if(IN[j]!=0)G.matrix[i,j] <- T_cs[i,j]/IN[j]}
#(iii)Invert I-G
inv.matrix<-solve(I.matrix-G.matrix)
#(iv)Calculate trophic positions
TP<-apply(inv.matrix,2,sum)}
##########NOTE: Results of these 2 methods are VERY DIFFERENT!!!

TP.INTERMs<-TP[INTERMs]
names(TP.INTERMs)<-nod.names[INTERMs]

##Trophic similarity
T_cs.INTERMs<-T_cs[INTERMs,INTERMs]
#TS.INTERMs<-TS_function(T_cs.INTERMs)

  sm <- dim(T_cs.INTERMs)[1]
  TS.INTERMs <- matrix(rep(0,sm^2), nrow = sm)
  for(i in 1:nrow(T_cs.INTERMs)){   #remove intraspecies interactions
    T_cs.INTERMs[i,i]<-0 }
  ##
  for(i in 1:sm){
    for(j in 1:sm){
      prey_A <- which(T_cs.INTERMs[,i]!=0)
      predator_A <- which(T_cs.INTERMs[i,]!=0)
      trophic_A <- unique(c(prey_A,predator_A))
      ##
      prey_B <- which(T_cs.INTERMs[,j]!=0)
      predator_B <- which(T_cs.INTERMs[j,]!=0)
      trophic_B <- unique(c(prey_B,predator_B))
      ##
      aa <- length(which(
        is.na(match(trophic_A,trophic_B))==TRUE))
      bb <- length(which(
        is.na(match(trophic_B,trophic_A))==TRUE))
      cc <- length(which(
        is.na(match(trophic_A,trophic_B))==FALSE))
      TS.INTERMs[i,j] <- cc/(aa + bb + cc)}}
  ##
  colnames(TS.INTERMs) <- colnames(T_cs.INTERMs)
  rownames(TS.INTERMs) <- colnames(T_cs.INTERMs)

#Convert TS matrix in edgelist+remove self-loops:
#CLUST<-TSM_to_elW(TS.INTERMs)
    ll<-length(which(TS.INTERMs!=0))
    CLUST<-data.frame(cbind(rep(NA,ll),rep(NA,ll)))
    k<-1
    for(i in 1:nrow(TS.INTERMs)){
      for(j in 1:ncol(TS.INTERMs)){
        if(rownames(TS.INTERMs)[i]!=colnames(TS.INTERMs)[j]){    #helps to identify the interactions
          CLUST[k,1]<-rownames(TS.INTERMs)[i]
          CLUST[k,2]<-colnames(TS.INTERMs)[j]
          CLUST[k,3]<-as.numeric(TS.INTERMs[i,j])
          k<-k + 1 } } }
colnames(CLUST)<-c("Nod1","Nod2","Trophic.Sim")

#remove NA rows (if there are):
if(length(which(is.na(CLUST[,"Trophic.Sim"])==TRUE))>0){
  nas<-which(is.na(CLUST[,"Trophic.Sim"])==TRUE)
  CLUST<-CLUST[-nas,]}

##Subset CLUST according to trophic similarity of nods:

ALLGROUPS<-ALLGROUPS.TS<-list()#lists to store nod IDs and TS of each group
nMEMBER<-vector() #vector to store sizes of all groups
nods.rm<-vector()

#Create nod groups of high trophic similarity:
#-->in each repeat loop we pick the nod with the highest median TS and
#   identify all nods that share a TS of >=0.6
#-->we store these nods as a group and remove them from CLUST

CLUST.rm<-CLUST #all identified nod groups are removed from the list
z<-0

repeat{
  
  if(nrow(CLUST.rm)==0){
    #We check wether the number of groupmembers and single nods corresponds the number
    #of nods in the trophic level:
    #-->if not this means that the nod of lowest median within-TP interval TS was not
    #   added in the list of groups (=occurs only if this nod is not integrated in a group)
    #-->we add this nod as single in the group list
    if(sum(nMEMBER)<interms){
      z<-z+1
      #We identify the NOD_ID of the lacking nod
      absentnod<-unique(CLUST[,"Nod1"])
      for(j in 1:length(nods.rm)){
        absentnod<-absentnod[-which(absentnod==nods.rm[j])]}
      #3.We add the nod to ALLGROUPS and ALLGROUPS.TS
      ALLGROUPS[[z]]<-absentnod
      ALLGROUPS.TS[[z]]<-"NULL"
      nMEMBER[z]<-1}
    
    break
    
  }else{
    z<-z+1
    #Pick nod with the highest median trophic similarity (OR RANDOMLY??):
    medianTS<-aggregate(x= CLUST.rm$Trophic.Sim,by = list(CLUST.rm$Nod1),
                        FUN = median)
    mynod<-medianTS[which(medianTS[,2]==max(medianTS[,2])),1]
    if(length(mynod)>1)mynod<-mynod[1]#if there are several canditates, pick the first
    
    #Add it to the group:
    GROUPMEMBERS<-mynod
    #subset all trophic similarites of this nod with other nods:
    mynod.allTS<-CLUST.rm[which(CLUST.rm[,"Nod1"]==mynod),c("Nod2","Trophic.Sim")]
    #identify all trophic similar other nods (TS >= 60%):
    similars<-mynod.allTS[which(mynod.allTS[,"Trophic.Sim"]>=tSim),]
    
    if(nrow(similars)!=0){ #if there are similar other nods
      Sim.similars<-vector()  
      for(u in 1:nrow(similars)){
        Sim<-similars[u,"Nod2"]
        #We add these similar nods to the group:
        GROUPMEMBERS<-c(GROUPMEMBERS,Sim)
        #We collect nods that are trophically similar to these nods:
        Sim.allTS<-CLUST.rm[which(CLUST.rm[,"Nod1"]==Sim),c("Nod2","Trophic.Sim")]
        Sim.sims<-as.vector(Sim.allTS[which(Sim.allTS[,"Trophic.Sim"]>=tSim),"Nod2"])
        #And store them in a vector:
        Sim.similars<-c(Sim.similars,Sim.sims)}#END OF FOR LOOP
      #We add them to the group:
      GROUPMEMBERS<-c(GROUPMEMBERS,Sim.similars)
      GROUPMEMBERS<-unique(GROUPMEMBERS)
      #We remove all groupmembers from CLUST.rm 
      #and store Nod IDs, TS and TP of group in a list:
      keep<-vector()
      remove<-vector()
      for(u in 1:length(GROUPMEMBERS)){
        member<-GROUPMEMBERS[u]
        othermembers<-GROUPMEMBERS[which(GROUPMEMBERS!=member)]
        out<-which(CLUST.rm[,"Nod1"]==member)
        out.out<-which(CLUST.rm[,"Nod2"]==member)
        remove<-c(remove,out,out.out)
        for(j in 1:length(othermembers)){
          ins<-which(CLUST.rm[out,"Nod2"]==othermembers[j])
          keep<-c(keep,out[ins])}}
      ALLGROUPS[[z]]<-GROUPMEMBERS
      nMEMBER<-c(nMEMBER,length(GROUPMEMBERS))
      nods.rm<-c(nods.rm,GROUPMEMBERS)
      remove<-unique(remove)
      keep<-unique(keep)
      ALLGROUPS.TS[[z]]<-CLUST.rm[keep,c("Nod1","Nod2","Trophic.Sim")]
      CLUST.rm<-CLUST.rm[-remove,]
      ####
      
    }else{ #if there are no similar nods
      
      #We remove the nod from CLUST and store nod ID, TP and TS in a list
      out<-which(CLUST.rm[,"Nod1"]==mynod)
      out.out<-which(CLUST.rm[,"Nod2"]==mynod)
      remove<-c(out,out.out)
      remove<-unique(remove)
      ALLGROUPS[[z]]<-mynod
      nMEMBER<-c(nMEMBER,1)
      nods.rm<-c(nods.rm,mynod)
      ALLGROUPS.TS[[z]]<-"NULL"
      CLUST.rm<-CLUST.rm[-remove,] }
  }#END OF OUTER ELSE CLAUSE
}#END OF REPEAT LOOP

###Identify trophic position of the members of each group und store this info
###in a large dataframe:
GROUP_INFO<-as.data.frame(matrix(NA,nrow=0,ncol=4,dimnames=list(NULL,
             c("Member_ID","GROUP","nMEMBER","TP"))))
z<-1

for(i in 1:length(nMEMBER)){
  members<-ALLGROUPS[[i]]
  group<-rep(z,nMEMBER[i])
  size<-rep(nMEMBER[i],nMEMBER[i])
  TPs<-vector()
  for(j in 1:nMEMBER[i]){
    member<-ALLGROUPS[[i]][j]
    member.TP<-TP.INTERMs[which(names(TP.INTERMs)==member)]
    TPs<-c(TPs,member.TP)}
  mygroup<-data.frame(Member_ID=members,GROUP=group,nMEMBER=size,TP=TPs,row.names=NULL)
  GROUP_INFO<-rbind(GROUP_INFO,mygroup)
  z<-z+1}

#################################################################################################
###Attribute x and y axis positions

#Y AXIS
#-->for intermediate compartments defined as [ 2 , 8.975 ]

#1.determine trophic position range of intermediate compartments
TP.min<-round(min(TP.INTERMs),digits=2)
range.TP<-round(max(TP.INTERMs),digits=2)-round(min(TP.INTERMs),digits=2)
#2.convert trophic positions in Y coordinates by multiplying values by a 
#conversion coefficient=Y range/TP range
GROUP_INFO$posY<-2+(GROUP_INFO$TP-TP.min)*(6.975/range.TP)
GROUP_INFO$posY<-round(GROUP_INFO$posY,digits=2)

#X AXIS
#-->X-axis: [ 0 , 12 ]

#1.divide y-axis in 9 levels: 
#   2.000; 2.775; 3.550; 4.325; 5.100; 5.875; 6.650; 7.425; 8.200; 8.975
#=nods located in a same Ylevel should have different x-axis positions
GRID.row<-round(seq(min(GROUP_INFO$posY),max(GROUP_INFO$posY),
                    length.out = 10),digits=2)

#2.estimate the number of nods to be positioned in each y level
#-->interval of ( y1 , y2 ]
GROUP_INFO$GRID.row<-cut(GROUP_INFO[,"posY"], breaks = GRID.row,
                         labels = c(1:9),include.lowest=TRUE)
GROUP_INFO$GRID.row<-as.numeric(GROUP_INFO$GRID.row)
GROUP_INFO[which(GROUP_INFO[,"posY"]==max(GROUP_INFO[,"posY"])),"GRID.row"]<-9
nNODSperRow<-vector(length=9)
for(i in 1:length(nNODSperRow)){
  nNODSperRow[i]<-length(which(GROUP_INFO[,"GRID.row"]==i))}


#3.Produce a grid that indicates possible x positions of intermediate nods 
#in each Ylevel(=9)
#-->x intervals of grid are based on the x coordinates for the row
#   with maximum number of nods to fit
#-->number of x intervals is defined as twice or 3x the max nb of nods/row
#     =allows to interspace nod groups
#-->empty x interval columns are subsequently removed
#-->grid matrix serves to mark x positions that are occupied by nodgroups (with groupID)
X.itv<-max(nNODSperRow)*2
if(X.itv<=length(ALLGROUPS))X.itv<-max(nNODSperRow)*3
GRIDX.free<-matrix(NA,nrow=9,ncol=X.itv,
                   byrow=TRUE,dimnames=list(c(1:9),NULL))
#vector to indicate number of free spaces if less nods than max nods
#need to be fitted in a row
PosX.free<-rep(ncol(GRIDX.free),length(ncol(GRIDX.free)))-nNODSperRow
#vector to indicate nb of unique groups per row
nGROUP.row<-rep(0,length(nNODSperRow))
groups.row<-unname(tapply(GROUP_INFO[,"GROUP"],INDEX=GROUP_INFO[,"GRID.row"],
                          FUN=function(x)length(unique(x))))
nGROUP.row[which(nNODSperRow!=0)]<-groups.row
#matrix that indicates the number of nods to fit per row and group
GROUPmembers.row<-matrix(0,nrow=9,ncol=length(ALLGROUPS),
                         dimnames=list(c(1:9),c(1:length(ALLGROUPS))))
for(i in 1:length(ALLGROUPS)){
  subset.info<-GROUP_INFO[which(GROUP_INFO[,"GROUP"]==i),]
  nMember.row<-tapply(subset.info[,"Member_ID"],
                      INDEX=as.factor(subset.info[,"GRID.row"]),
                      FUN=function(x)length(unique(x)))
  GROUPmembers.row[names(nMember.row),i]<-unname(nMember.row)}

#vector to indicate max interval between nod groups/row
#if spacing should be equal among groups
spacingX.pot<-rep(NA,length(PosX.free))
spacingX.pot[which(nGROUP.row>1)]<-floor(PosX.free[which(nGROUP.row>1)]/
                                           (nGROUP.row[which(nGROUP.row>1)]-1))
spacingX.pot[which(nGROUP.row<=1)]<-PosX.free[which(nGROUP.row<=1)]

##############

###################
#Attribute x-axis positions according to trophic similarity of nods:

#(I)Order groups along X-axis based on between-group trophic similarity

#1.Identify nb of edges (IN+OUTGOING) per nod group
nEDGES.group<-rep(NA,length(ALLGROUPS))
for(i in 1:length(ALLGROUPS)){
  members<-ALLGROUPS[[i]]
  nEDGES.group[i]<-length(which(T_cs[members,]!=0))+length(which(T_cs[,members]!=0))}

#2.Select group with highest nb of edges and order other groups in function
#  to between-group TS relative to this group
myhub<-which(nEDGES.group==max(nEDGES.group))[1]
meanTS.myhub<-rep(NA,length(ALLGROUPS))
pos<-c(1:(ncol(TS.INTERMs)^2))
pos.diag<-which(pos%%(ncol(TS.INTERMs)+1)==1)
allTS.myhub<-TS.INTERMs
allTS.myhub[pos.diag]<-0#remove TS of same compartments
a<-ALLGROUPS[[myhub]]
if(length(a)>1){allTS.myhub<-as.matrix(allTS.myhub[a,])
}else{allTS.myhub<-matrix(allTS.myhub[a,],ncol=ncol(allTS.myhub),
                          dimnames=list(a,colnames(allTS.myhub)))}
for(i in 1:length(ALLGROUPS)){
  if(i==myhub)next
  mygroup<-ALLGROUPS[[i]]
  meanTS.myhub[i]<-mean(allTS.myhub[,mygroup])}

#3.Pick the 2 groups of highest between-group TS and order them right & 
#left from myhub.Focus on neighbour groups(=groups left-side and right-side located 
#from myhub) and find for each the group that is highest trophically similar to them. 
#Repeat (3.)...
XorderGROUPS<-myhub
names(ALLGROUPS)<-c(1:length(ALLGROUPS))
nods.positioned<-ALLGROUPS[[myhub]]
which.positioned<-c()
for(i in 1:length(nods.positioned)){
  which.nod<-which(colnames(T_cs.INTERMs)==nods.positioned[i]) 
  which.positioned<-c(which.positioned,which.nod)}
myhub.l<-which(meanTS.myhub==max(meanTS.myhub,na.rm=TRUE))[1]
meanTS.myhub[myhub.l]<-NA
myhub.r<-which(meanTS.myhub==max(meanTS.myhub,na.rm=TRUE))[1]

repeat{
  
  #exit command:
  #if(length(XorderGROUPS)==length(ALLGROUPS))break
  
  #(i)Calculate avg TS of right and left side located hubs with other groups:
  meanTS.myhub.l<-meanTS.myhub.r<-rep(NA,(length(ALLGROUPS)-length(XorderGROUPS)))
  allTS.myhubs<-as.matrix(TS.INTERMs[-which.positioned,-which.positioned])
  #remove TS of same compartments:
  pos<-c(1:(ncol(allTS.myhubs)^2))
  pos.diag<-which(pos%%(ncol(allTS.myhubs)+1)==1)
  allTS.myhubs[pos.diag]<-0
  #remove groups that are already ordered in XorderGROUP from set of candidates:
  to.position<-as.list(ALLGROUPS[-XorderGROUPS])
  names<-names(ALLGROUPS)[-XorderGROUPS]
  names(to.position)<-as.character(names)
  a<-unname(unlist(to.position[which(names(to.position)==as.character(myhub.l))]))
  b<-unname(unlist(to.position[which(names(to.position)==as.character(myhub.r))]))
  if(length(a)>1){allTS.myhub.l<-as.matrix(allTS.myhubs[a,])
  }else{allTS.myhub.l<-matrix(allTS.myhubs[a,],ncol=ncol(allTS.myhubs),
                              dimnames=list(a,colnames(allTS.myhubs)))}
  if(length(b)>1){allTS.myhub.r<-as.matrix(allTS.myhubs[b,])
  }else{allTS.myhub.r<-matrix(allTS.myhubs[b,],ncol=ncol(allTS.myhubs),
                              dimnames=list(b,colnames(allTS.myhubs)))}
  out.l<-which(names(to.position)==myhub.l)
  out.r<-which(names(to.position)==myhub.r)
  for(i in 1:length(to.position)){
    if(i==out.l||i==out.r)next
    mygroup<-to.position[[i]]
    meanTS.myhub.l[i]<-mean(allTS.myhub.l[,mygroup])
    meanTS.myhub.r[i]<-mean(allTS.myhub.r[,mygroup])}
  #(ii)Pick for left and right hub the group with highest TS relative to them
  #-->if both share same group:look for whom TS is higher and pick for other the
  #   2nd highest TS group:
  neighbour.l<-which(meanTS.myhub.l==max(meanTS.myhub.l,na.rm=TRUE))[1]
  neighbour.r<-which(meanTS.myhub.r==max(meanTS.myhub.r,na.rm=TRUE))[1]
  if(neighbour.l==neighbour.r){
    if(meanTS.myhub.l[neighbour.l]>meanTS.myhub.r[neighbour.r]){
      meanTS.myhub.r[neighbour.r]<-NA 
      if(length(which(is.na(meanTS.myhub.r)==FALSE))>0){
        neighbour.r<-which(meanTS.myhub.r==max(meanTS.myhub.r,na.rm=TRUE))[1]
      }else{neighbour.r<-NULL}
    }else{
      meanTS.myhub.l[neighbour.l]<-NA 
      if(length(which(is.na(meanTS.myhub.l)==FALSE))>0){
        neighbour.l<-which(meanTS.myhub.l==max(meanTS.myhub.l,na.rm=TRUE))[1]
      }else{neighbour.l<-NULL}  }}
  #(iii)remove positioned nods from set of candidates
  nods.positioned<-c(ALLGROUPS[[myhub.l]],ALLGROUPS[[myhub.r]])
  which.nods<-c()
  for(i in 1:length(nods.positioned)){
    which.nod<-which(colnames(T_cs.INTERMs)==nods.positioned[i]) 
    which.nods<-c(which.nods,which.nod)}
  which.positioned<-c(which.positioned,which.nods)
  #(iv)Add myhub.l and myhub.r to XorderGROUPS and redefine myhub.l and myhub.r:
  XorderGROUPS<-c(myhub.l,XorderGROUPS,myhub.r)
  myhub.l<-as.numeric(names(to.position)[neighbour.l])
  myhub.r<-as.numeric(names(to.position)[neighbour.r])
  #exit command:
  if(length(myhub.l)==0){
    XorderGROUPS<-c(XorderGROUPS,myhub.r)
    break
  }else if(length(myhub.r)==0){
    XorderGROUPS<-c(myhub.l,XorderGROUPS)
    break
  }else if(length(to.position)==4){
    XorderGROUPS<-c(myhub.l,XorderGROUPS,myhub.r)
    break}
  
}#END OF REPEAT LOOP

###
#(II)
#(A)Occupy x positions for all nod groups

#identify single nods and nod groups within the Xorder
whichS_L<-nMEMBER[XorderGROUPS]
whichS<-which(whichS_L==1)
whichL<-which(whichS_L>1)
whichS_L[whichS]<-"S"
whichS_L[whichL]<-"L"

#matrix to indicate X-axis space occupied by each group:
Xfield<-matrix(NA,nrow=2,ncol=length(ALLGROUPS),
               dimnames=list(c("Xfrom","Xto"),c(1:length(ALLGROUPS))))

#vectors to indicate which groups have been already fit in GRIDX.free
fittedGROUPS<-vector()
z<-0
repeat{
  
  #-->exit loop:
  if(length(fittedGROUPS)==length(XorderGROUPS))break 
  
  z<-z+1
  #if(z==10)break
  #-->pick groups following the order of XorderGROUPS and select row with highest nb 
  #   of members
  #=if several rows with same amount of members:pick first one
  mygroup<-XorderGROUPS[z]
  nMAX.row<-unname(which(GROUPmembers.row[,mygroup]==
                           max(GROUPmembers.row[,mygroup])))[1]
  members.tofit<-GROUPmembers.row[nMAX.row,mygroup]
  
  #-->occupy x positions for all members of this row  in GRIDX.free
  
  #(i)check number of groups that have been already fit
  # + occupy x positions (if needed) in GRIDX.free for othergroup 
  #   members in row
  
  if(length(fittedGROUPS)>0){
    for(i in 1:length(fittedGROUPS)){
      #if(i==8)break
      if(GROUPmembers.row[nMAX.row,fittedGROUPS[i]]==0){
        next
        
      }else{
        #check Xfield that has been occupied by the group
        if(i>1){
          start.pot<-Xfield["Xfrom",fittedGROUPS[i]]
          startsX.TRUE<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))
          if(startsX.TRUE==start.pot){sX<-0
          }else if(startsX.TRUE<start.pot){
            sX<-length(GRIDX.free[nMAX.row,c(startsX.TRUE:(start.pot-1))])
          }else{sX<-spacingX.pot[nMAX.row]+1}
          
          if(spacingX.pot[nMAX.row]>=sX){
            start<-start.pot
          }else if(startsX.TRUE>Xfield["Xto",fittedGROUPS[i]]){
            start<-startsX.TRUE
            sX<-0
            #PosX.free[nMAX.row]<-PosX.free[nMAX.row]+GROUPmembers.row[nMAX.row,fittedGROUPS[i]]
          }else{start<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))+
            spacingX.pot[nMAX.row]} 
        }else{start<-1}
        end<-start+GROUPmembers.row[nMAX.row,fittedGROUPS[i]]-1
        #occupy x positions for the row.members of the other group in 
        #GRIDX.free+indicate in GROUPmembers.row and nGROUPS.row 
        #that groupmembers have been fit
        GRIDX.free[nMAX.row,c(start:end)]<-fittedGROUPS[i]
        GROUPmembers.row[nMAX.row,fittedGROUPS[i]]<-0
        nGROUP.row[nMAX.row]<-nGROUP.row[nMAX.row]-1
        
        #-->only if other groups have been already fit and spacing is possible (or wanted):
        if(i>1){
          #mark added spaces in GRIDX.free(with zeros)
          if(sX>0){
            sX.start<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))
            sX.end<-start-1
            if(sX.start>sX.end)print(paste("PROBLEM with sX start+end definition! (1) in z=",z))
            GRIDX.free[nMAX.row,c(sX.start:sX.end)]<-0
            #recalculate PosX free+potential spacing
            sX.fit<-length(GRIDX.free[nMAX.row,c(sX.start:sX.end)])
            PosX.free[nMAX.row]<-PosX.free[nMAX.row]-sX.fit
            spacingX.pot[which(nGROUP.row>1)]<-floor(PosX.free[which(nGROUP.row>1)]/
                                                       (nGROUP.row[which(nGROUP.row>1)]-1))
            spacingX.pot[which(nGROUP.row<=1)]<-PosX.free[which(nGROUP.row<=1)] }
          
        } } }
  }#END OF IF CLAUSE=there are other groups that have been already fit    
  
  #(ii)Mark X positions for members of mygroup in GRIDX.free
  
  if(length(fittedGROUPS)>0){
    if(whichS_L[z]=="L"){ #if mygroup has >1 member
      #-->check groupID of the group that most extends along X axis
      #   + determine max X position occupied by this group
      group.left<-which(Xfield["Xto",]==max(Xfield["Xto",],na.rm=TRUE))[1]
      Xoccupied.left<-Xfield["Xto",group.left]
      start.pot<-Xoccupied.left+1 
      startsX.TRUE<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))
      if(startsX.TRUE==start.pot || start.pot>ncol(GRIDX.free)){sX<-0
      }else if(startsX.TRUE<start.pot){
        sX<-length(GRIDX.free[nMAX.row,c(startsX.TRUE:(start.pot-1))])
      }else{sX<-spacingX.pot[nMAX.row]+1}
      
      if(spacingX.pot[nMAX.row]>=sX && start.pot<=ncol(GRIDX.free)){
        start<-start.pot
      }else{start<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))+
        spacingX.pot[nMAX.row]
        sX<-spacingX.pot[nMAX.row]}
      
      
    }else if(whichS_L[z]=="S"){ #if mygroup is a single nod
      #determine the most left available space that does not overlap with other nod groups 
      #of >1member
      most.leftX<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE)) 
      if(length(which(whichS_L[c(1:(z-1))]=="L"))>0){
        min.Xspace.free<-max(Xfield["Xto",which(nMEMBER>1)],na.rm=TRUE)+1
      }else{min.Xspace.free<-most.leftX}
      if(most.leftX>=min.Xspace.free  || min.Xspace.free>ncol(GRIDX.free)){
        start<-most.leftX
        sX<-0
      }else{
        start.pot<-min.Xspace.free
        #check if introduced spacing is possible
        sX<-length(seq(most.leftX,(start.pot-1),by=1))
        if(spacingX.pot[nMAX.row]>=sX && start.pot<ncol(GRIDX.free)){start<-start.pot
        }else{
        start<-most.leftX+spacingX.pot[nMAX.row]
        sX<-spacingX.pot[nMAX.row]} }
    }#END OF NOD GROUP/SINGLE NOD DEFINITION      
    
  }else{start<-1}
  
  end<-start+GROUPmembers.row[nMAX.row,mygroup]-1
  
  #occupy x positions for the row.members of mygroup in 
  #GRIDX.free+indicate in GROUPmembers.row and nGROUPS.row 
  #that groupmembers have been fit
  GRIDX.free[nMAX.row,c(start:end)]<-mygroup
  GROUPmembers.row[nMAX.row,mygroup]<-0
  nGROUP.row[nMAX.row]<-nGROUP.row[nMAX.row]-1
  #indicate X space occupied by mygroup in Xfield
  Xfield["Xfrom",mygroup]<-start
  Xfield["Xto",mygroup]<-end
  
  #-->only if other groups have been already fit and spacing is possible:
  if(length(fittedGROUPS)>0 && sX>0){
    
    #for nod groups:      
    if(whichS_L[z]=="L"){ 
      #mark added spaces in GRIDX.free(with zeros)
      sX.start<-min(which(is.na(GRIDX.free[nMAX.row,])==TRUE))
      sX.end<-start-1
      if(sX.start>sX.end)print(paste("PROBLEM with sX start+end definition! (2) in z=",z))
      GRIDX.free[nMAX.row,c(sX.start:sX.end)]<-0
      #recalculate PosX free+potential spacing
      sX.fit<-length(GRIDX.free[nMAX.row,c(sX.start:sX.end)])
      PosX.free[nMAX.row]<-PosX.free[nMAX.row]-sX.fit
      spacingX.pot[which(nGROUP.row>1)]<-floor(PosX.free[which(nGROUP.row>1)]/
                                                 (nGROUP.row[which(nGROUP.row>1)]-1))
      spacingX.pot[which(nGROUP.row<=1)]<-PosX.free[which(nGROUP.row<=1)] 
      
      #for single nods:
    }else if(whichS_L[z]=="S"){
      #mark added spaces in GRIDX.free(with zeros)
      sX.start<-most.leftX
      sX.end<-start-1
      if(sX.start>sX.end)print(paste("PROBLEM with sX start+end definition! (2) in z=",z))
      GRIDX.free[nMAX.row,c(sX.start:sX.end)]<-0
      #recalculate PosX free+potential spacing
      sX.fit<-length(GRIDX.free[nMAX.row,c(sX.start:sX.end)])
      PosX.free[nMAX.row]<-PosX.free[nMAX.row]-sX.fit
      spacingX.pot[which(nGROUP.row>1)]<-floor(PosX.free[which(nGROUP.row>1)]/
                                                 (nGROUP.row[which(nGROUP.row>1)]-1))
      spacingX.pot[which(nGROUP.row<=1)]<-PosX.free[which(nGROUP.row<=1)] }   
  }#END OF IF CLAUSES=there are other groups that were already fit and spacing is possible
  
  #Add mygroup to fitted groups
  fittedGROUPS<-c(fittedGROUPS,mygroup)
  
}#END OF REPEAT LOOP


#(B)Occupy for all groups of >1 member X positions for all remaining members
nNODS.remaining<-apply(GROUPmembers.row,1,sum)
rows.tofit<-which(nNODS.remaining>0)

#1.Occupy row by row X positions for all remaining groupmembers
#-->group positions are occupied following the order how groups have been allocated
#   from left to right in GRIDX.free
for(i in 1:length(rows.tofit)){
  myrow<-rows.tofit[i]
  GROUPS.tofit<-which(GROUPmembers.row[myrow,]>0)
  nods.tofit<-nNODS.remaining[myrow]
  myrow.free<-rep(0,ncol(GRIDX.free))
  myrow.free[which(is.na(GRIDX.free[myrow,])==TRUE)]<-1
  
  
  for(j in 1:length(fittedGROUPS)){
    if(length(which(GROUPS.tofit==fittedGROUPS[j]))==0){next 
      
    }else{
      mygroup<-fittedGROUPS[j]
      members.tofit<-GROUPmembers.row[myrow,mygroup]
      Xfield.mygroup<-seq(Xfield["Xfrom",mygroup],Xfield["Xto",mygroup],by=1)
      fit.OK<-c()
      
      for(k in 1:length(Xfield.mygroup)){
        if(myrow.free[Xfield.mygroup[k]]==1){
          fit.OK<-c(fit.OK,Xfield.mygroup[k])}}
      #check if X positions are available in GRIDX.free
      if(length(fit.OK)>=members.tofit){fit.OK<-fit.OK[c(1:members.tofit)]
      }else if(length(fit.OK)>0){
        allfree.left<-myrow.free[c(1:fit.OK[length(fit.OK)])]
        if(sum(allfree.left)>=members.tofit){
          fit.mygroup<-seq(to=length(which(allfree.left==1)),by=1,length.out = members.tofit)
          fit.OK<-allfree.left[which(allfree.left==1)][fit.mygroup] 
        }else{fit.OK<-myrow.free[which(myrow.free==1)][c(1:members.tofit)]}
      }else if(length(fit.OK)==0){
        #Look  for closest available Xpositions in GRIDX.free
        free1<-max(which(is.na(GRIDX.free[myrow,])==TRUE))
        free2<-min(which(is.na(GRIDX.free[myrow,])==TRUE))
        free1<-abs(free1-Xfield["Xfrom",mygroup])
        free2<-abs(free2-Xfield["Xfrom",mygroup])
        if(free1<free2){
          fit.OK<-seq(to=max(which(is.na(GRIDX.free[myrow,])==TRUE)),by=1,
                      length.out = members.tofit)}  
        if(free2<=free1){
          fit.OK<-seq(from=min(which(is.na(GRIDX.free[myrow,])==TRUE)),by=1,
                      length.out = members.tofit)} }
      #check if occupying these X positions guarantees enough space for all nods 
      #that remain to be fit
      X.free<-myrow.free[c((fit.OK[length(fit.OK)]+1):length(myrow.free))]
      free.OK<-sum(X.free)
      if(free.OK<(nods.tofit-members.tofit)){
        #redefine fit.OK
        X.free<-which(myrow.free==1)
        free.OK<-seq(to=length(X.free),by=1,length.out = (nods.tofit-members.tofit))
        X.occupied<-seq(to=(length(X.free)-1),by=1,length.out = members.tofit)
        fit.OK<-X.free[X.occupied]}
      #occupy X positions for mygroup
      GRIDX.free[myrow,fit.OK]<-mygroup
      #if there was allocated additional space on left side, occupy space by zeros
      if(min(which(is.na(GRIDX.free[myrow,])))<fit.OK[1]){
        zeros<-seq(from=min(which(is.na(GRIDX.free[myrow,]))),to=(fit.OK[1]-1),by=1)
        GRIDX.free[myrow,zeros]<-0
      }else{zeros<-NULL}
      
      #mark X positions as occupied in myrow.free
      myrow.free[c(zeros,fit.OK)]<-0
      #remove fitted groupmembers from GROUPmembers.row
      GROUPmembers.row[myrow,mygroup]<-0
      #remove fitted groupmembers from nods.tofit
      nods.tofit<-nods.tofit-members.tofit
    }#END OF ELSE CLAUSE=my group has members to fit in myrow  
  }#END OF INNER FOR LOOP=groups to fit per row
}#END OF OUTER FOR LOOP=row by row

#Fill NA values of GRIDX.free with zeros(=spaces)
GRIDX.free[which(is.na(GRIDX.free))]<-0

#(C)Improve Group allocation in GRIDX.free 

#1.Estimate realized X field of each group+add Xfield of single nods
for(i in 1:length(ALLGROUPS)){
  Xfrom<-rep(NA,nrow(GRIDX.free))
  Xto<-rep(NA,nrow(GRIDX.free))
  if(nMEMBER[i]==1){
    Xfield["Xfrom",i]<-Xfield["Xto",i]<-which(GRIDX.free==i,arr.ind = TRUE)[2]
  }else{
    for(k in 1:nrow(GRIDX.free)){
      if(length(which(GRIDX.free[k,]==i))>0){
        Xfrom[k]<-min(which(GRIDX.free[k,]==i)) 
        Xto[k]<-max(which(GRIDX.free[k,]==i))}}
    Xfield["Xfrom",i]<-Xfrom[which(Xfrom==min(Xfrom,na.rm=TRUE))[1]]
    Xfield["Xto",i]<-Xto[which(Xto==max(Xto,na.rm=TRUE))[1]]}}

#2.Rearrange nods:
#-->starting from with group of min Xfrom till group of max Xfrom
#identify overlaps in Xfield with other groups
#-->distinguish between overlaps with other nod groups and single nods
#-->if yes:
#orientate nods along Xaxis so that they least possible overlap with other group
#-->if no:
#orientate nods around a common center

orderGROUP<-order(Xfield["Xfrom",],decreasing=FALSE)
whichsingles<-colnames(Xfield)
whichsingles[which(nMEMBER==1)]<-NA
shift<-rep(NA,length(orderGROUP))

for(i in 1:length(orderGROUP)){
  if(is.na(whichsingles[orderGROUP[i]])==TRUE)next
  Xfield.mygroup<-seq(Xfield["Xfrom",orderGROUP[i]],Xfield["Xto",orderGROUP[i]],by=1)
  overlaps<-rep(0,length(Xfield.mygroup))
  for(k in 1:length(Xfield.mygroup)){
    nGROUPS.col<-unique(GRIDX.free[,Xfield.mygroup[k]])
    nGROUPS.col<-nGROUPS.col[which(nGROUPS.col!=0)]
    if(length(nGROUPS.col)>1){overlaps[k]<-1 }}
  if(sum(overlaps)>0)shift[i]<-"center" #TO REVIEW!!!!!!!
  if(i==1&&sum(overlaps)>0)shift[i]<-"left"
  if(i>1&&overlaps[1]==1)shift[i]<-"right"
  if(i>1&&overlaps[length(overlaps)]==1)shift[i]<-"left"
  if(i>1&&overlaps[1]==1&&overlaps[length(overlaps)]==1)shift[i]<-"center"
  if(sum(overlaps)==0)shift[i]<-"center" 
}#END OF FOR LOOP


for(i in 1:length(orderGROUP)){
  if(i==9)break
  mygroup<-orderGROUP[i]
  if(is.na(whichsingles[mygroup])==TRUE)next
  Xfield.mygroup<-seq(Xfield["Xfrom",mygroup],Xfield["Xto",mygroup],by=1)
  XGRID.mygroup<-as.matrix(GRIDX.free[,Xfield.mygroup])
  rows.mygroup<-unique(which(XGRID.mygroup==mygroup,arr.ind = TRUE)[,"row"])
  
  for(k in 1:length(rows.mygroup)){
    myrow<-XGRID.mygroup[rows.mygroup[k],]
    groupmembers.myrow<-which(myrow==mygroup)
    space.myrow<-which(myrow==0)
    myrow[c(groupmembers.myrow,space.myrow)]<-NA
    if(shift[i]=="left"){myrow[which(is.na(myrow==TRUE))]<-
      c(rep(mygroup,length(groupmembers.myrow)),rep(0,length(space.myrow)))}
    if(shift[i]=="right"){myrow[which(is.na(myrow==TRUE))]<-
      c(rep(0,length(space.myrow)),rep(mygroup,length(groupmembers.myrow)))}
    if(shift[i]=="center"){
      if(length(which(is.na(myrow==TRUE)))==length(myrow)&&length(space.myrow)>0){
        impair<-k%%2!=0
        if(impair==TRUE){
          mid1<-floor(length(space.myrow)/2)
          mid2<-ceiling(length(space.myrow)/2)}
        if(impair==FALSE){
          mid1<-ceiling(length(space.myrow)/2)
          mid2<-floor(length(space.myrow)/2)}
        myrow<-c(rep(0,mid1),rep(mygroup,length(groupmembers.myrow)),
                 rep(0,mid2))
      }else{
        myrow<-XGRID.mygroup[rows.mygroup[k],]}}#else keep nod positions as they are
    XGRID.mygroup[rows.mygroup[k],]<-myrow}#END OF INNER FOR LOOP
  GRIDX.free[,Xfield.mygroup]<-XGRID.mygroup
}#END OF OUTER FOR LOOP

#Remove empty columns (only zeros)from GRIDX.free
colSum<-apply(GRIDX.free,2,sum)
zeros<-which(colSum==0)
if(length(zeros)>0){
GRIDX.free<-GRIDX.free[,-zeros]}

#(III)Attribute colors to groups according to group size
allCOL<-rep(NA,length(XorderGROUPS))
#SINGLE NODS:
if(length(whichS)<=length(colS)){
  COLS<-colS[1:length(whichS)]
}else{
  nNOCOL<-length(whichS)-length(colS)
  COLS<-c(colS,rep("azure3",nNOCOL))
  print("Reached limit of unique colors.")}
allCOL[whichS]<-COLS
#NOD GROUPS
if(length(whichL)<=length(colL)){
  COLL<-colL[1:length(whichL)]
}else{
  nNOCOL<-length(whichL)-length(colL)
  COLL<-c(colL,rep("azure4",nNOCOL))
  print("Reached limit of unique colors.")}
allCOL[whichL]<-COLL
allCOL<-allCOL[order(XorderGROUPS,decreasing=FALSE)]

#(IV)Attribute all nods an x-axis and y-axis position
#-->x-coordinates are based on the number of columns in GRIDX.free
#for row 1,3,5,9: we first attribute the nod with lowest Y position an X position
#for row 2,4,6,8: we first attribute the nod with highest Y position an X position
allPosX<-matrix(c(1:(ncol(GRIDX.free)*nrow(GRIDX.free))),ncol=ncol(GRIDX.free),
                byrow=TRUE)
x.coord<-round(seq(0,12,length.out = ncol(GRIDX.free)),digits=2)
GRIDX<-matrix(rep(x.coord,9),nrow=9,ncol=ncol(GRIDX.free),
              byrow=TRUE,dimnames=list(c(1:9),NULL))

for(i in 1:length(orderGROUP)){
  mygroup<-orderGROUP[i]
  infos<-as.data.frame(GROUP_INFO[which(GROUP_INFO[,"GROUP"]==mygroup),])
  colgroup<-allCOL[mygroup]
  nrows<-unique(infos[,"GRID.row"])[order(unique(infos[,"GRID.row"]),decreasing=FALSE)]
  #define all X positions of my group:requires byrow=TRUE!!!!
  PosX.byrow<-allPosX[which(GRIDX.free==mygroup)]
  PosX.bycol<-which(GRIDX.free==mygroup)
  PosX.mygroup<-PosX.bycol[order(PosX.byrow,decreasing=FALSE)]
  impair<-nrows%%2!=0
  nMember.row<-rep(NA,length(nrows))
  for(p in 1:length(nrows)){
    count<-length(which(infos[,"GRID.row"]==nrows[p]))
    nMember.row[p]<-count}
  for(p in 1:length(nrows)){
    nods.tofit<-as.data.frame(infos[which(infos[,"GRID.row"]==nrows[p]),])
    if(nrow(nods.tofit)>1){
      if(impair[p]==TRUE){
        nods.tofit<-nods.tofit[order(nods.tofit[,"posY"],decreasing=FALSE),]}
      if(impair[p]==FALSE){
        nods.tofit<-nods.tofit[order(nods.tofit[,"posY"],decreasing=TRUE),]} }
    for(j in 1:nrow(nods.tofit)){
      nodID<-nods.tofit[j,"Member_ID"]
      nod_info[which(nod_info[,"Nod_ID"]==nodID),"y.position"]<-nods.tofit[j,"posY"]
      nod_info[which(nod_info[,"Nod_ID"]==nodID),"color"]<-colgroup
      if(p>1){PosX.nod<-sum(nMember.row[c(1:(p-1))])+j}else{PosX.nod<-j}
      CoordX.nod<-GRIDX[PosX.mygroup[PosX.nod]]
      nod_info[which(nod_info[,"Nod_ID"]==nodID),"x.position"]<-CoordX.nod} }}
##################END OF AUTOMATED XY-COORDINATE ATTRIBUTION

#OPTION2: user definition of XY coordinates
}else{
  if(length(select.cols)==1)select.cols<-"azure4"
  coordX<-select.XY[[1]]
  coordY<-select.XY[[2]]
  nodnames<-nod.names[c(1:length(coordX))]
  nod_info<-data.frame(Nod_ID=nodnames, x.position=coordX, y.position=coordY, color=select.cols)
}#END OF ELSE CLAUSE (Option2)

###################################################################################
###B.NOD LABELS####
###################################################################################
if(length(labels)>1){LABS<-labels 
}else if(labels=="Default"){
labs<-strsplit(nod_info[,"Nod_ID"], split=" ")
LABS<-rep(NA,nrow(nod_info))
for(i in 1:length(LABS)){
  mylab<-labs[[i]]
  mylab<-mylab[which(mylab!="")]
  for(j in 1:length(mylab)){
    if(j==1){myLAB<-mylab[j]
    }else{myLAB<-paste(myLAB,"\n",mylab[j],sep="")}}
  LABS[i]<-myLAB}
}else if(labels=="numeric"){LABS<-c(1:comps)}


####################################################################################
###C.ATTRIBUTES####
####################################################################################

#1.Set nod attributes
nodAttr<-data.frame(attr=c("shape","color","width","height","fixedsize"),
                    val=NA, attr_type="node")
nodAttr[1,"val"]<-"circle"
nodAttr[2,"val"]<-"black"
if(is.na(nodsize)==TRUE){
if(length(select.XY)==1){nodsize<-0.5*unname(GRIDX[1,2])
}else{nodsize<-0.5}}
if(length(biomass)==1){
nodAttr[c(3,4),"val"]<-nodsize
nodAttr[5,"val"]<-TRUE}

#2.Set nod labels and label attributes
labelAttr<-data.frame(attr=c("fontsize","fontcolor","fontname","forcelabels"),
                      val=NA, attr_type="node")
labelAttr[1,"val"]<-12
labelAttr[2,"val"]<-"black"
labelAttr[3,"val"]<-"Futura"
labelAttr[4,"val"]<-FALSE

#3.Set edge attributes
edgeAttr<-data.frame(attr=c("arrowsize","arrowhead","arrowtail"),
                     val=NA, attr_type="edge")
edgeAttr[1,"val"]<-0.75
edgeAttr[2,"val"]<-"vee"
edgeAttr[3,"val"]<-"vee"


#4.Set graph attributes
grAttr<-data.frame(attr=c("layout","outputorder","splines"),
                   val=NA,attr_type="graph")
grAttr[1,"val"]<-"neato"
grAttr[2,"val"]<-"edgesfirst"
grAttr[3,"val"]<-"spline"

ALLATTR<-rbind(grAttr,labelAttr,nodAttr,edgeAttr)  

if(length(select.Attr)>1){
  ALLATTR<-merge(select.Attr,ALLATTR,by=c("attr","attr_type"),all=TRUE)
  for(i in 1:nrow(ALLATTR)){
    if(is.na(ALLATTR[i,"val"])==TRUE)ALLATTR[i,"val"]<-ALLATTR[i,"value"]}
  ALLATTR<-ALLATTR[,c("attr","val","attr_type")] }

####################################################################################
###D.EDGES AND EDGECOLORS####
####################################################################################
#-->edge color is defined by the nod from where it goes out

#1.Produce a weighted edgelist for all nods
if(include.detritus==TRUE){
  nNODS<-living+nl
ITM.nodID<-T_cs
}else{
  ITM.nodID<-T_cs.living
  nNODS<-living}
colnames(ITM.nodID)<-rownames(ITM.nodID)<-c(1:ncol(ITM.nodID))

if(selfloop==FALSE){
    pos<-c(1:(ncol(ITM.nodID)^2))
    pos.diag<-which(pos%%(ncol(ITM.nodID)+1)==1)
    ITM.nodID[pos.diag]<-0}
  ll<-length(which(ITM.nodID!=0))
  edge_info<-data.frame(cbind(rep(NA,ll),rep(NA,ll)))
  k<-1
  for(i in 1:nrow(ITM.nodID)){
    for(j in 1:ncol(ITM.nodID)){
      if(ITM.nodID[i,j]!=0){    #helps to identify the interactions
        edge_info[k,1]<-rownames(ITM.nodID)[i]
        edge_info[k,2]<-colnames(ITM.nodID)[j]
        edge_info[k,3]<-as.numeric(ITM.nodID[i,j])
        k<-k + 1} } }
colnames(edge_info)<-c("FROM","TO","WEIGHT")

#2.Add edge colors
edgecols<-rep(NA,nrow(edge_info))
allnods.out<-unique(edge_info[,"FROM"])
for(i in 1:length(allnods.out)){
  mynod<-as.numeric(allnods.out[i])
  mycol<-nod_info[which(nod_info[,"Nod_ID"]==colnames(T_cs)[mynod]),"color"]
  edgecols[which(edge_info[,"FROM"]==mynod)]<-mycol}
if(include.detritus==TRUE){
  for(i in 1:nl){
    det<-as.character(living+i)
    whichdet<-which(edge_info[,"TO"]==det)
  edgecols[whichdet]<-"coral4"}}
edge_info$color<-edgecols

##############################################################################################
###E.CREATE GRAPH####
##############################################################################################
nod_df <- create_node_df(n = nNODS, style = "filled",
                          fillcolor = nod_info[c(1:nNODS),"color"],
                          label = LABS[c(1:nNODS)], 
                          x = nod_info[c(1:nNODS),"x.position"], 
                          y = nod_info[c(1:nNODS),"y.position"],
                          fontcolor="black",fontname="Futura")
if(length(biomass)>1){
  nodsize.biom<-nodsize*abs(log10(unname(biomass[c(1:nNODS)])* 10^6)/4)
  nod_df$width<-nod_df$height<-nodsize.biom}

edge_df <- create_edge_df(from = edge_info[,"FROM"], to = edge_info[,"TO"],
                           penwidth = abs(log10(edge_info[,"WEIGHT"] * 10^6)/4), 
                           color = edge_info[,"color"])

#Create graph:
gr<-create_graph(nodes_df = nod_df, edges_df = edge_df) 
default.attr<-gr[["global_attrs"]]
ALLATTR<-merge(default.attr,ALLATTR,by=c("attr","attr_type"),all=TRUE)
for(i in 1:nrow(ALLATTR)){
  if(is.na(ALLATTR[i,"val"])==TRUE)ALLATTR[i,"val"]<-ALLATTR[i,"value"]}
ALLATTR<-ALLATTR[,c("attr","val","attr_type")]
colnames(ALLATTR)<-c("attr","value","attr_type")
gr[["global_attrs"]]<-ALLATTR

###############################################################################################
###F.FUNCTION OUTPUT####

render_graph(gr)

return(gr)

} #END OF FUNCTION
###############################################################################################
