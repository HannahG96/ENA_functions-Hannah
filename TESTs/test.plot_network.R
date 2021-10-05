setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")

#####IMPORT NETWORKS FOR TESTING############################################

#1.GRAMDRY=66 nods
Gramdry<-SCOR.convert("Network data/gramdry.dat")
Z_cs1<-Gramdry[[4]]#Import
E_cs1<-Gramdry[[5]]#Export
R_cs1<-Gramdry[[6]]#Respiration
T_cs1<-Gramdry[[7]]#intercompartmental exchange
biomass1<-Gramdry[[3]]#stock sizes 
nl1<-Gramdry[[2]][1]-Gramdry[[2]][2]#number of non-living compartments 

#2.GRAMWET=66 nods
Gramwet<-SCOR.convert("Network data/gramwet.dat")
Z_cs2<-Gramwet[[4]]#Import
E_cs2<-Gramwet[[5]]#Export
R_cs2<-Gramwet[[6]]#Respiration
T_cs2<-Gramwet[[7]]#intercompartmental exchange
biomass2<-Gramwet[[3]]#stock sizes 
nl2<-Gramwet[[2]][1]-Gramwet[[2]][2]#number of non-living compartments 

#3.CRYSTAL RIVER CREEK=21 nods
CrystalCreek<-SCOR.convert("Network data/CrystalRiverControl.txt")
Z_cs3<-CrystalCreek[[4]]#Import
E_cs3<-CrystalCreek[[5]]#Export
R_cs3<-CrystalCreek[[6]]#Respiration
T_cs3<-CrystalCreek[[7]]#intercompartmental exchange
biomass3<-CrystalCreek[[3]]#stock sizes 
nl3<-CrystalCreek[[2]][1]-CrystalCreek[[2]][2]#number of non-living compartments 

#4.CYPRESS WET SEASON=68 nods
Cypwet<-SCOR.convert("Network data/cypwet.dat")
Z_cs4<-Cypwet[[4]]#Import
E_cs4<-Cypwet[[5]]#Export
R_cs4<-Cypwet[[6]]#Respiration
T_cs4<-Cypwet[[7]]#intercompartmental exchange
biomass4<-Cypwet[[3]]#stock sizes 
nl4<-Cypwet[[2]][1]-Cypwet[[2]][2]#number of non-living compartments 

#5.MONDEGO ESTUARY=43 nods
Mondest<-SCOR.convert("Network data/MondegoEstuary.txt")
Z_cs5<-Mondest[[4]]#Import
E_cs5<-Mondest[[5]]#Export
R_cs5<-Mondest[[6]]#Respiration
T_cs5<-Mondest[[7]]#intercompartmental exchange
biomass5<-Mondest[[3]]#stock sizes 
nl5<-Mondest[[2]][1]-Mondest[[2]][2]#number of non-living compartments 

#6.ST MARKS RIVER=51 nods
Stmarks<-SCOR.convert("Network data/StMarks.txt")
Z_cs6<-Stmarks[[4]]#Import
E_cs6<-Stmarks[[5]]#Export
R_cs6<-Stmarks[[6]]#Respiration
T_cs6<-Stmarks[[7]]#intercompartmental exchange
biomass6<-Stmarks[[3]]#stock sizes 
nl6<-Stmarks[[2]][1]-Stmarks[[2]][2]#number of non-living compartments 

#############################TEST OPTIONS################################

######################
##LARGE networks:
######################
#-->66 nods
#DEFAULT SETTINGS
Gramdry.plot<-plot.network(Z_cs=Z_cs1,T_cs=T_cs1,E_cs=E_cs1,R_cs=R_cs1,nl=nl1,
                           include.detritus=FALSE, selfloop=FALSE,tSim=0.6,
               select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
               nodsize=NA, biomass=NA)
render_graph(Gramdry.plot)
#INCLUDE BIOMASS PROPORTIONALITY
Gramdry.plotI<-plot.network(Z_cs=Z_cs1,T_cs=T_cs1,E_cs=E_cs1,R_cs=R_cs1,nl=nl1,
                           include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.6,
                           select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                           nodsize=NA, biomass=biomass1)
render_graph(Gramdry.plotI)
####################
#-->68 nods
#DEFAULT SETTINGS
Cypwet.plot<-plot.network(Z_cs=Z_cs4,T_cs=T_cs4,E_cs=E_cs4,R_cs=R_cs4,nl=nl4,
                          include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.6,
                  select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                  nodsize=NA, biomass=NA)
render_graph(Cypwet.plot)
#INCLUDE DETRITUS
Cypwet.plotI<-plot.network(Z_cs=Z_cs4,T_cs=T_cs4,E_cs=E_cs4,R_cs=R_cs4,nl=nl4,
                            include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                            select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                            nodsize=NA, biomass=NA)
render_graph(Cypwet.plotI)
#+INCLUDE BIOMASS PROPORTIONALITY
Cypwet.plotII<-plot.network(Z_cs=Z_cs4,T_cs=T_cs4,E_cs=E_cs4,R_cs=R_cs4,nl=nl4,
                            include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                            select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                            nodsize=NA, biomass=biomass4)
render_graph(Cypwet.plotII)
#INCLUDE SELFLOOPS
Cypwet.plotIII<-plot.network(Z_cs=Z_cs4,T_cs=T_cs4,E_cs=E_cs4,R_cs=R_cs4,nl=nl4,
                          include.detritus=FALSE, selfloop=TRUE, TP.method=1, tSim=0.6,
                          select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                          nodsize=NA, biomass=NA)
render_graph(Cypwet.plotIII)#no selfloops(?)
######################
##MEDIUM networks:
######################
#-->43 nods
#DEFAULT SETTINGS
Mondest.plot<-plot.network(Z_cs=Z_cs5,T_cs=T_cs5,E_cs=E_cs5,R_cs=R_cs5,nl=nl5,
               include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.6,
               select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
               nodsize=NA, biomass=NA)
render_graph(Mondest.plot)
#INCLUDE DETRITUS
Mondest.plotI<-plot.network(Z_cs=Z_cs5,T_cs=T_cs5,E_cs=E_cs5,R_cs=R_cs5,nl=nl5,
                           include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                           select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                           nodsize=NA, biomass=NA)
render_graph(Mondest.plotI)
#+INCLUDE BIOMASS PROPORTIONALITY
Mondest.plotII<-plot.network(Z_cs=Z_cs5,T_cs=T_cs5,E_cs=E_cs5,R_cs=R_cs5,nl=nl5,
                            include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                            select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                            nodsize=NA, biomass=biomass5)
render_graph(Mondest.plotII)

####################
#-->51 nods
Stmarks.plot<-plot.network(Z_cs=Z_cs6,T_cs=T_cs6,E_cs=E_cs6,R_cs=R_cs6,nl=nl6,
                           include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.4,
                           select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                           nodsize=NA, biomass=biomass6)
render_graph(Stmarks.plot)

####################
##SMALL networks:
####################
#-->21 nods
#DEFAULT SETTINGS
CrystalCreek.plot<-plot.network(Z_cs=Z_cs3,T_cs=T_cs3,E_cs=E_cs3,R_cs=R_cs3,nl=nl3,
                     include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.6,
                     select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                     nodsize=NA, biomass=NA)
render_graph(CrystalCreek.plot)
#INCLUDE DETRITUS
CrystalCreek.plotI<-plot.network(Z_cs=Z_cs3,T_cs=T_cs3,E_cs=E_cs3,R_cs=R_cs3,nl=nl3,
                                include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                                select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                                nodsize=NA, biomass=NA)
render_graph(CrystalCreek.plotI)
#+INCLUDE BIOMASS PROPORTIONALITY
CrystalCreek.plotII<-plot.network(Z_cs=Z_cs3,T_cs=T_cs3,E_cs=E_cs3,R_cs=R_cs3,nl=nl3,
                            include.detritus=TRUE, selfloop=FALSE, TP.method=1, tSim=0.6,
                            select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                            nodsize=NA, biomass=biomass3)
render_graph(CrystalCreek.plotII)
#+INCLUDE SELFLOOPS
CrystalCreek.plotIII<-plot.network(Z_cs=Z_cs3,T_cs=T_cs3,E_cs=E_cs3,R_cs=R_cs3,nl=nl3,
                                  include.detritus=TRUE, selfloop=TRUE, TP.method=1, tSim=0.6,
                                  select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                                  nodsize=NA, biomass=biomass3)
render_graph(CrystalCreek.plotIII) #(there are no selfloops in this network)

#########################TEST OWN SELECTION OF COORDINATES/COLORS/LABELS############
###LARGE NETWORK
#-->66 nods
#produce xy coordinates:
mycoordX<-as.numeric(Gramdry.plot[["nodes_df"]][,"x"])
mycoordX<-sample(mycoordX,replace=FALSE)
mycoordY<-as.numeric(Gramdry.plot[["nodes_df"]][,"y"])
mycoordY<-sample(mycoordY,replace=FALSE)
mycoordXY<-list(mycoordX,mycoordY)
#produce color vector:
mycols<-c(rep("red",6),rep("blue",6),rep("green",6),rep("pink",6),
          rep("gold",6),rep("purple",6),rep("thistle",6),rep("azure4",6),
          rep("orange",6),rep("olivedrab",6),rep("cyan",3))
mycols<-sample(mycols,replace=FALSE)
#produce labels:
mylabels<-mycols
Gramdry.plot.XY<-plot.network(Z_cs=Z_cs1,T_cs=T_cs1,E_cs=E_cs1,R_cs=R_cs1,nl=nl1,
                           include.detritus=FALSE, selfloop=FALSE, TP.method=1, tSim=0.6,
                           select.XY=mycoordXY, select.cols=mycols, select.Attr=NA,
                           labels=mylabels,
                           nodsize=NA, biomass=NA)
render_graph(Gramdry.plot.XY)

########################TEST CASE: nl=0#############################################
Mondest.plot.nl0<-plot.network(Z_cs=Z_cs5,T_cs=T_cs5,E_cs=E_cs5,R_cs=R_cs5,nl=0,
                           include.detritus=FALSE, selfloop=FALSE, TP.method=2, tSim=0.6,
                           select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                           nodsize=NA, biomass=NA)
render_graph(Mondest.plot.nl0)
