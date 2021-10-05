###################VISUALIZE RESULTS OF CYCLING ANALYSIS#########################

###FUNCTION################################################################
plot.cycling<-function(cyc.list,Z_cs,T_cs,E_cs,R_cs,nl){

library("DiagrammeR")
library("ggplot2")
###########################################################################
##A-CYCLE LENGTH DISTRIBUTION
#(i)Get cycle length distribution
CYCDIST<-as.data.frame(cyc.list[[2]])
#(ii)Produce barplot
cycdist.levels<-CYCDIST[order(CYCDIST[,"length"],decreasing=FALSE),"length"]
CYCDIST$length<-factor(CYCDIST$length,levels=cycdist.levels)
plotA<-ggplot(CYCDIST,aes(x=length,y=Nb.cycles))+
  theme_minimal()+
  geom_bar(width = 0.7, fill="darkgoldenrod3", stat = "identity") +
  ylab("Number of cycles")+
  xlab("Cycle length")+
  theme(axis.title.y=element_text(size=10,face="bold",color="gray21"),
        axis.text.y = element_text(size=7,color="azure4"),
        axis.title.x=element_text(size=10,face="bold",color="gray21"),
        axis.text.x = element_text(size=10,face="bold",color="azure4"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###########################################################################
##B-WEAK ARC BARPLOTs
#-->REVIEW THIcKNESS OF BARS
#Get weak arc information
WEAK<-cyc.list[[4]]
#(I)dodge barplot of flows: 2 bars=removed flow & true flow
#-->LOG10-SCALE
#(i)Re-organize data for plotting
labs.weak<-paste(WEAK[,"Outgoing"],"-->",WEAK[,"Ingoing"],sep="")
labs.weak<-c(labs.weak,labs.weak)
value.weak<-c(log10(WEAK[,"Flow"]),log10(WEAK[,"Flow.removed"]))
value.weak[which(is.infinite(value.weak)==TRUE)]<-0
flowtype<-c(rep("true",nrow(WEAK)),rep("removed",nrow(WEAK)))
plot.weakflows<-data.frame(WeakID=labs.weak,Flow=value.weak,Type=flowtype)
#(ii)Define barwidth based on number of weak arks
if(nrow(WEAK)>=70){barsize.weak<-0.7;labelsize<-6
}else if(nrow(WEAK)>=50){barsize.weak<-0.75;labelsize<-6.5
}else if(nrow(WEAK>=30)){barsize.weak<-0.8;labelsize<-6.5
}else if(nrow(WEAK>=20)){barsize.weak<-0.9;labelsize<-7
}else if(nrow(WEAK>=10)){barsize.weak<-1;labelsize<-7
}else{barsize.weak<-1.25;labelsize<-8}
#(iii)Produce barplot
plotB1<-ggplot(plot.weakflows,aes(x=WeakID,y=Flow,fill=Type))+
  theme_minimal()+
  geom_bar(position="dodge", stat = "identity",width = barsize.weak) +
  scale_fill_manual("",values=c("gold","darkmagenta"),
                    breaks=c("removed","true"), 
                    labels=c("removed flow","true flow"))+
  ylab("Log10-transformed flows")+
  theme(legend.text =element_text(size=7, face="bold", color="azure4"),
        legend.position = "top",
        axis.title.y=element_text(size=8, face="bold", color="gray21"),
        axis.text.y = element_text(size=9,color="azure4"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=labelsize,angle=90,color="gray21",
                                   face="bold",hjust=0,vjust=0),
        panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line( size=.01, color="azure4"),
        panel.grid.minor.x = element_blank())


#(II)nexus size barplot: number of nexus member cycles~weak arc ID
#(i)Re-organize data for plotting
labs.weak<-labs.weak[c(1:nrow(WEAK))]
Nnex<-WEAK[,"Nb.cycles"]
plot.nex<-data.frame(WeakID=labs.weak,Ncycles=Nnex)
#(ii)Produce barplot
plotB2<-ggplot(plot.nex,aes(x=WeakID, y=Ncycles))+
  theme_minimal()+
  geom_bar(width = barsize.weak, fill="darkgoldenrod3", stat = "identity") +
  ylab("Number of nexus member cycles")+
  theme(axis.title.y=element_text(size=8, face="bold", color="gray21"),
        axis.text.y = element_text(size=9,color="azure4"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=labelsize,angle=90,color="gray21",
                                   face="bold",hjust=0,vjust=0),
        panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line( size=.01, color="azure4"),
        panel.grid.minor.x = element_blank())

###########################################################################
##C-WEAK ARC NETWORK GRAPH
#-->ADJUST THIKNESS OF WEAK ARC EDGES!!!!!!!!!!!!!!!!!!
#-->reqires Input of T_cs,Z_cs,E_cs,R_cs,nl!!!!!!!!!!!!!!!
#-->weak arcs are highlighted in colors, edge thickness is proportional
#   to nexus size
#-->other network interactions as grey thin lines
#(i)extract node and edge information and graph attributes from plot.network()
fullgraph<-plot.network(Z_cs=Z_cs,T_cs=T_cs,E_cs=E_cs,R_cs=R_cs,nl=nl,
                        include.detritus=TRUE, selfloop=TRUE,tSim=0.6,
                        select.XY=NA, select.cols=NA, select.Attr=NA,labels="Default",
                        nodsize=NA, biomass=NA)
edge_info<-fullgraph[["edges_df"]]#edge information
#(ii)identify weak arc rows in edge_info
weakedges<-c()
names_comps<-colnames(T_cs)
for(i in 1:nrow(WEAK)){
  nodIN<-WEAK[i,"Ingoing"]
  nodOUT<-WEAK[i,"Outgoing"]
  numIN<-which(names_comps==nodIN)
  numOUT<-which(names_comps==nodOUT)
  a<-which(edge_info[,"to"]==numIN)
  aa<-which(edge_info[a,"from"]==numOUT)
  myrow<-a[aa]
  weakedges<-c(weakedges,myrow)}
#(iii)define penwidth of weak arcs=proportional to nexus size
#    and other egdges=thin line
#-->we therefore predefine size classes:
#   >=10 000 ; >=7 500 ; >=5000 ; >=2500 ; >=1000 ; >=500 ; >=250 ; >=100 ;
#   >=50 ; >=20 ; >=10 ; <10
#SIZE CLASSES TO BE REVIEWED!

if(nrow(WEAK)>=70){coef<-1
}else if(nrow(WEAK)>=50){coef<-1.25
}else if(nrow(WEAK>=30)){coef<-1.5
}else if(nrow(WEAK>=20)){coef<-1.75
}else if(nrow(WEAK>=10)){coef<-2
}else{coef<-3}
nex<-WEAK[,"Nb.cycles"]
prop.weak<-rep(NA,length(nex))
prop.weak[which(nex<10)]<-0+coef
prop.weak[which(nex>=10)]<-0.1+coef
prop.weak[which(nex>=20)]<-0.2+coef
prop.weak[which(nex>=50)]<-0.3+coef
prop.weak[which(nex>=100)]<-0.5+coef
prop.weak[which(nex>=250)]<-0.75+coef
prop.weak[which(nex>=500)]<-1+coef
prop.weak[which(nex>=750)]<-1.25+coef
prop.weak[which(nex>=1000)]<-1.5+coef
prop.weak[which(nex>=2500)]<-1.75+coef
prop.weak[which(nex>=5000)]<-2+coef
prop.weak[which(nex>=7500)]<-2.25+coef
prop.weak[which(nex>=10000)]<-2.5+coef
edge_info[weakedges,"penwidth"]<-prop.weak
edge_info[-weakedges,"penwidth"]<-0.25
#(iv)Set the color of all non-weak arc interactions as grey
edge_info[-weakedges,"color"]<-"azure4"
#(v)Create graph
weakgraph<-fullgraph
weakgraph[["edges_df"]]<-edge_info
plotC<-weakgraph

###########################################################################
##D-RECYCLED MATTER ~ CYCLE LENGTH
#-->LOG10-SCALE
#(I)absolute values
#(i)Get recycled matter information
RECYCL<-cyc.list[[5]]
#(ii)Organize data for plotting
plot.recycl<-data.frame(Value=log10(RECYCL),Length=names(RECYCL))
plot.recycl[which(is.infinite(plot.recycl[,"Value"])==TRUE),"Value"]<-0
#(iii)Produce barplot
plot.recycl$Length<-factor(plot.recycl$Length,levels=cycdist.levels)
plotD1<-ggplot(plot.recycl,aes(x=Length,y=Value))+
  theme_minimal()+
  geom_bar(width = 0.7, fill="darkgoldenrod3", stat = "identity") +
  ylab("Log10-transformed matter recycled")+
  xlab("Cycle length")+
  theme(axis.title.y=element_text(size=10, face="bold", color="gray21"),
        axis.text.y = element_text(size=9,color="azure4"),
        axis.title.x=element_text(size=10,face="bold",color="gray21"),
        axis.text.x = element_text(size=10,face="bold",color="azure4"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#(II)values normalized by TST
#(i)Get normalized recycled matter information
RECYCLnorm<-cyc.list[[6]]
#(ii)Organize data for plotting
plot.recyclnorm<-data.frame(Value=log10(RECYCLnorm),Length=names(RECYCLnorm))
plot.recyclnorm[which(is.infinite(plot.recyclnorm[,"Value"])==TRUE),"Value"]<-0
#(iii)Produce barplot
plot.recyclnorm$Length<-factor(plot.recyclnorm$Length,levels=cycdist.levels)
plotD2<-ggplot(plot.recyclnorm,aes(x=Length,y=Value))+
  theme_minimal()+
  geom_bar(width = 0.7, fill="darkgoldenrod3", stat = "identity") +
  ylab("Normalized\nlog10-transformed matter recycled")+
  xlab("Cycle length")+
  theme(axis.title.y=element_text(size=9, face="bold", color="gray21"),
        axis.text.y = element_text(size=9,color="azure4"),
        axis.title.x=element_text(size=10,face="bold",color="gray21"),
        axis.text.x = element_text(size=10,face="bold",color="azure4"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###########################################################################
##E-ACYCLIC NETWORK GRAPH
#-->reqires Input of T_cs,Z_cs,E_cs,R_cs,nl!!!!!!!!!!!!!!
#-->NOT INFORMATIVE!!!
#(i)extract node and edge information and graph attributes from plot.network()
grAttr<-fullgraph[["global_attrs"]]#attributes
nod_info<-fullgraph[["nodes_df"]]#node information
#(ii)produce edgelist of acyclic network
acyc.matrix<-cyc.list[[7]]
acyc.matrix<-round(acyc.matrix, digits=7)#-->remove near-zero values
nNODS<-ncol(acyc.matrix)
nod.names<-rownames(acyc.matrix)
living<-nNODS-nl
colnames(acyc.matrix)<-rownames(acyc.matrix)<-c(1:nNODS)
ll<-length(which(acyc.matrix!=0))
edge_info<-data.frame(cbind(rep(NA,ll),rep(NA,ll)))
k<-1
for(i in 1:nrow(acyc.matrix)){
  for(j in 1:ncol(acyc.matrix)){
    if(acyc.matrix[i,j]!=0){    #helps to identify the interactions
      edge_info[k,1]<-rownames(acyc.matrix)[i]
      edge_info[k,2]<-colnames(acyc.matrix)[j]
      edge_info[k,3]<-as.numeric(acyc.matrix[i,j])
      k<-k + 1} } }
colnames(edge_info)<-c("FROM","TO","WEIGHT")
#(iii)attribute colors to edges
edgecols<-rep(NA,nrow(edge_info))
allnods.out<-unique(edge_info[,"FROM"])
for(i in 1:length(allnods.out)){
  mynod<-as.numeric(allnods.out[i])
  mycol<-nod_info[which(nod.names==colnames(T_cs)[mynod]),"fillcolor"]
  edgecols[which(edge_info[,"FROM"]==mynod)]<-mycol}
for(i in 1:nl){
    det<-as.character(living+i)
    whichdet<-which(edge_info[,"TO"]==det)
    edgecols[whichdet]<-"coral4"}
edge_info$color<-edgecols
#(vi)adjust proportionality coefficient of edge thickness to flow value
prop.flow<-abs(log10(edge_info[,"WEIGHT"])) 
prop.flow<-round(prop.flow,digits=2)
prop.flow<-(prop.flow/max(prop.flow))*2
#(v)create edge information df
edge_df <- create_edge_df(from = edge_info[,"FROM"], to = edge_info[,"TO"],
                          penwidth = prop.flow,
                          fontcolor= edge_info[,"color"],fontsize=12,                                       
                          color = edge_info[,"color"])

#Create graph:
gr<-create_graph(nodes_df = nod_info, edges_df = edge_df) 
grAttr<-rbind(grAttr,c("fontcolor","black","edge"))
gr[["global_attrs"]]<-grAttr
plotE<-gr

###########################################################################
##FUNCTION OUTPUT
CYCPLOTS<-list(plotA,plotB1,plotB2,plotC,plotD1,plotD2,plotE)
names(CYCPLOTS)<-c("Cycle Length Distribution","Nexus Analysis I","Nexus Analysis II",
                   "Weak arc network","Cycle Distribution of recycled matter",
                   "Normalized Cycle Distribution of recycled matter","Acyclic network")
return(CYCPLOTS)
}#END OF PLOT:CYCLING FUNCTION
