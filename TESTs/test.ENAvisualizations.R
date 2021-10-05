setwd("C:/Hannah/Biological Oceanography/HiWi Scotti/ENA PACKAGE/TESTs")
############################TEST ENA VISUALIZATIONS##########################

#1.GRAMDRY=66 nods
#Import network:
Gramdry<-SCOR.convert("Network data/gramdry.dat")
Z_cs<-Gramdry[[4]]#Import
E_cs<-Gramdry[[5]]#Export
R_cs<-Gramdry[[6]]#Respiration
T_cs<-Gramdry[[7]]#intercompartmental exchange
biomass<-Gramdry[[3]]#stock sizes 
nl<-Gramdry[[2]][1]-Gramdry[[2]][2]#number of non-living compartments 
#Produce cycling output:
###default:
cyc.Gramdry1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            select="Default", nmax=10000)
###weighted:
#cyc.Gramdry2<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            #select="weighted", nmax=10000)

#I-Test cycling visualization:
PLOTcyc.Gramdry<-plot.cycling(cyc.list=cyc.Gramdry1,Z_cs,T_cs,E_cs,R_cs,nl)
plot(PLOTcyc.Gramdry[[1]])
plot(PLOTcyc.Gramdry[[2]])
plot(PLOTcyc.Gramdry[[3]])
render_graph(PLOTcyc.Gramdry[[4]])
plot(PLOTcyc.Gramdry[[5]])
plot(PLOTcyc.Gramdry[[6]])
render_graph(PLOTcyc.Gramdry[[7]])

#Perform trophic analysis:
troph.Gramdry<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)

#II-Test trophic visualization:
PLOTtroph.Gramdry<-plot.trophic(TROPHIC=troph.Gramdry,biomass=biomass)
render_graph(PLOTtroph.Gramdry[[1]])
plot(PLOTtroph.Gramdry[[2]])
render_graph(PLOTtroph.Gramdry[[3]])

#III-Test indices visualization
PLOTinds.Gramdry<-plot.allindices(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced")
plot(PLOTinds.Gramdry)

################################################################################################
#2.MONDEGO ESTUARY=43 nods
#Import network:
Mondest<-SCOR.convert("Network data/MondegoEstuary.txt")
Z_cs<-Mondest[[4]]#Import
E_cs<-Mondest[[5]]#Export
R_cs<-Mondest[[6]]#Respiration
T_cs<-Mondest[[7]]#intercompartmental exchange
biomass<-Mondest[[3]]#stock sizes 
nl<-Mondest[[2]][1]-Mondest[[2]][2]#number of non-living compartments

#Produce cycling output:
###default:
cyc.Mondest1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            select="Default", nmax=10000)
###weighted:
#cyc.Mondest2<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            #select="weighted", nmax=10000)

#I-Test cycling visualization:
PLOTcyc.Mondest<-plot.cycling(cyc.list=cyc.Mondest1,Z_cs,T_cs,E_cs,R_cs,nl)
plot(PLOTcyc.Mondest[[1]])
plot(PLOTcyc.Mondest[[2]])
plot(PLOTcyc.Mondest[[3]])
render_graph(PLOTcyc.Mondest[[4]])
plot(PLOTcyc.Mondest[[5]])
plot(PLOTcyc.Mondest[[6]])
render_graph(PLOTcyc.Mondest[[7]])

#Perform trophic analysis:
troph.Mondest<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)

#II-Test trophic visualization:
PLOTtroph.Mondest<-plot.trophic(TROPHIC=troph.Mondest,biomass=biomass)
render_graph(PLOTtroph.Mondest[[1]])
plot(PLOTtroph.Mondest[[2]])
render_graph(PLOTtroph.Mondest[[3]])

#III-Test indices visualization
PLOTinds.Mondest<-plot.allindices(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced")
plot(PLOTinds.Mondest)

###################################################################################################
#3.CRYSTAL RIVER CREEK=21 nods
CrystalCreek<-SCOR.convert("Network data/CrystalRiverControl.txt")
Z_cs<-CrystalCreek[[4]]#Import
E_cs<-CrystalCreek[[5]]#Export
R_cs<-CrystalCreek[[6]]#Respiration
T_cs<-CrystalCreek[[7]]#intercompartmental exchange
biomass<-CrystalCreek[[3]]#stock sizes 
nl<-CrystalCreek[[2]][1]-CrystalCreek[[2]][2]#number of non-living compartments 

#Produce cycling output:
###default:
cyc.CrystalCreek1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            select="Default", nmax=10000)
###weighted:
#cyc.CrystalCreek2<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                           #select="weighted", nmax=10000)

#I-Test cycling visualization:
PLOTcyc.CrystalCreek<-plot.cycling(cyc.list=cyc.CrystalCreek1,Z_cs,T_cs,E_cs,R_cs,nl)
plot(PLOTcyc.CrystalCreek[[1]])
plot(PLOTcyc.CrystalCreek[[2]])
plot(PLOTcyc.CrystalCreek[[3]])
render_graph(PLOTcyc.CrystalCreek[[4]])
plot(PLOTcyc.CrystalCreek[[5]])
plot(PLOTcyc.CrystalCreek[[6]])
render_graph(PLOTcyc.CrystalCreek[[7]])

#Perform trophic analysis:
troph.CrystalCreek<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)

#II-Test trophic visualization:
PLOTtroph.CrystalCreek<-plot.trophic(TROPHIC=troph.CrystalCreek,biomass=biomass)
render_graph(PLOTtroph.CrystalCreek[[1]])
plot(PLOTtroph.CrystalCreek[[2]])
render_graph(PLOTtroph.CrystalCreek[[3]])

#III-Test indices visualization
PLOTinds.CrystalCreek<-plot.allindices(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced")
plot(PLOTinds.CrystalCreek)

######################################################################################################
#4.CYPRESS WET SEASON=68 nods
Cypwet<-SCOR.convert("Network data/cypwet.dat")
Z_cs<-Cypwet[[4]]#Import
E_cs<-Cypwet[[5]]#Export
R_cs<-Cypwet[[6]]#Respiration
T_cs<-Cypwet[[7]]#intercompartmental exchange
biomass<-Cypwet[[3]]#stock sizes 
nl<-Cypwet[[2]][1]-Cypwet[[2]][2]#number of non-living compartments 

#Produce cycling output:
###default:
cyc.Cypwet1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                            select="Default", nmax=10000)
###weighted:
#cyc.Cypwet2<-test.function(Z_cs,T_cs,E_cs,R_cs, 
#select="weighted", nmax=10000)

#I-Test cycling visualization:
PLOTcyc.Cypwet<-cycling.analysis(cyc.list=cyc.Cypwet1,Z_cs,T_cs,E_cs,R_cs,nl)
plot(PLOTcyc.Cypwet[[1]])
plot(PLOTcyc.Cypwet[[2]])
plot(PLOTcyc.Cypwet[[3]])
render_graph(PLOTcyc.Cypwet[[4]])
plot(PLOTcyc.Cypwet[[5]])
plot(PLOTcyc.Cypwet[[6]])
render_graph(PLOTcyc.Cypwet[[7]])

#Perform trophic analysis:
troph.Cypwet<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)

#II-Test trophic visualization:
PLOTtroph.Cypwet<-plot.trophic(TROPHIC=troph.Cypwet,biomass=biomass)
render_graph(PLOTtroph.Cypwet[[1]])
plot(PLOTtroph.Cypwet[[2]])
render_graph(PLOTtroph.Cypwet[[3]])

#III-Test indices visualization
PLOTinds.Cypwet<-plot.allindices(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced")
plot(PLOTinds.Cypwet)

#######################################################################################################
#5.ST MARKS RIVER=51 nods
Stmarks<-SCOR.convert("Network data/StMarks.txt")
Z_cs<-Stmarks[[4]]#Import
E_cs<-Stmarks[[5]]#Export
R_cs<-Stmarks[[6]]#Respiration
T_cs<-Stmarks[[7]]#intercompartmental exchange
biomass<-Stmarks[[3]]#stock sizes 
nl<-Stmarks[[2]][1]-Stmarks[[2]][2]#number of non-living compartments 

#Produce cycling output:
###default:
cyc.Stmarks1<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
                           select="Default", nmax=10000)
###weighted:
#cyc.Stmarks2<-cycling.analysis(Z_cs,T_cs,E_cs,R_cs, 
#select="weighted", nmax=10000)

#I-Test cycling visualization:
PLOTcyc.Stmarks<-plot.cycling(cyc.list=cyc.Stmarks1,Z_cs,T_cs,E_cs,R_cs,nl)
plot(PLOTcyc.Stmarks[[1]])
plot(PLOTcyc.Stmarks[[2]])
plot(PLOTcyc.Stmarks[[3]])
render_graph(PLOTcyc.Stmarks[[4]])
plot(PLOTcyc.Stmarks[[5]])
plot(PLOTcyc.Stmarks[[6]])
render_graph(PLOTcyc.Stmarks[[7]])

#Perform trophic analysis:
troph.Stmarks<-TP.CTA(inp=Z_cs,inter=T_cs,outt=E_cs,diss=R_cs,nl=nl)

#II-Test trophic visualization:
PLOTtroph.Stmarks<-plot.trophic(TROPHIC=troph.Stmarks,biomass=biomass)
render_graph(PLOTtroph.Stmarks[[1]])
plot(PLOTtroph.Stmarks[[2]])
render_graph(PLOTtroph.Stmarks[[3]])

#III-Test indices visualization
PLOTinds.Stmarks<-plot.allindices(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced")
plot(PLOTinds.Stmarks)