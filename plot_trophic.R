########################VISUALIZE OUTPUTS OF TROPHIC ANALYSIS####################################

#####FUNCTION
plot.trophic<-function(TROPHIC,biomass=NA){

library("DiagrammeR")
library("ggplot2")
#########################LINDEMAN SPINE#########################
#-->select different fontcolors for edgelabels!!!! ie. respiration in azure4
#-->remove arrowhead of edge-->Edge1
#-->increase edge label size(???)
####
#1.Create all necessary objects to produce Lindeman Spine
####
#Trophic Transformation matrix
A.matrix<-TROPHIC[[1]]
#Canonical Exports
CEXP<-TROPHIC[[3]]
#Canonical Respiration
CRESP<-TROPHIC[[4]]
#Canonical Imports-->only if there occurs migration
if(length(TROPHIC)==13){migr<-TRUE}else{migr<-FALSE}
if(migr==TRUE){
CINP<-TROPHIC[[13]]}
#Returns to Detrital pool
backDET<-as.vector(TROPHIC[[6]])
#Input to detrital pool
inpDET<-TROPHIC[[8]]
#Lindeman Spine
LindS<-TROPHIC[[10]]
#Stock sizes per discrete trophic level
if(length(biomass)>1){
  biom<-A.matrix
  for(i in 1:ncol(A.matrix))biom[,i]<-biom[,i]*biomass[i]
  BIOM<-apply(biom,1,sum)}

####
#2.Define number of discrete trophic levels to be shown in diagramme
####
nTL<-length(which(round(LindS,digits=7)>0))
if(nTL>12)nTL<-12
A.matrix<-A.matrix[c(1:nTL),]
CEXP<-CEXP[c(1:nTL)]
CRESP<-CRESP[c(1:nTL)]
if(migr==TRUE)CINP<-CINP[c(1:nTL)]
backDET<-backDET[c(1:nTL)]
LindS<-LindS[c(1:nTL)]
if(length(biomass)>1)BIOM<-BIOM[c(1:nTL)]
####
#3.Produce edgelists for:
#Lindman Spine=Tlevel-->Tlevel+1
#Export= Tlevel-->OUT
#Import= IN-->Tlevel
#Respiration= Tlevel-->LOSS
#Returns to detrital pool= Tlevel-->Dummy-->Edge1
####
cexp<-cresp<-cinp<-linds<-backdet<-as.data.frame(matrix(NA,nrow=0,ncol=3))
z<-0
for(i in 1:nTL){
  CEXP.outcoming<-i
  CEXP.incoming<-paste("OUT",i)
  CEXP.value<-CEXP[i]
  cexp<-rbind(cexp,c(CEXP.outcoming,CEXP.incoming,CEXP.value))
  CRESP.outcoming<-i
  CRESP.incoming<-paste("LOSS",i)
  CRESP.value<-CRESP[i]
  cresp<-rbind(cresp,c(CRESP.outcoming,CRESP.incoming,CRESP.value))
  if(migr==TRUE){
    CINP.outcoming<-paste("IN",i)
    CINP.incoming<-i
    CINP.value<-CINP[i]
    cinp<-rbind(cinp,c(CINP.outcoming,CINP.incoming,CINP.value))}
  LindS.outcoming<-z
  z<-z+1
  LindS.incoming<-i
  if(i==1){LindS.value<-LindS[i]+inpDET #add input to detrital pool
  }else{LindS.value<-LindS[i]}
  linds<-rbind(linds,c(LindS.outcoming,LindS.incoming,LindS.value))
  backDET.outcoming<-i
  backDET.passby<-paste("Dummy",i,sep="")
  backDET.value<-backDET[i]
  backdet<-rbind(backdet,c(backDET.outcoming,backDET.passby,backDET.value))}
backbackdet<-data.frame(FROM=backdet[,2], TO=rep("eDummy",nTL), WEIGHT=NA)
colnames(backdet)<-c("FROM","TO","WEIGHT")
backdet<-rbind(backdet,backbackdet,c("eDummy","Edge1",NA))
colnames(cexp)<-colnames(cresp)<-colnames(cinp)<-colnames(linds)<-
             c("FROM","TO","WEIGHT")
edge_info<-rbind(linds,cexp,cresp,backdet)
whichFROM<-c(c(0:(nTL-1)),c(1:nTL),c(1:nTL),c(1:nTL),rep("Dummy",nTL),"eDummy")
whichTO<-c(c(1:nTL),rep("OUT",nTL),rep("LOSS",nTL),rep("Dummy",nTL),rep("eDummy",nTL),"Edge1")
if(migr==TRUE){
  edge_info<-rbind(edge_info,cinp)
  whichFROM<-c(whichFROM,rep("IN",nTL))
  whichTO<-c(whichTO,c(1:nTL))}

edge_info$WEIGHT<-as.numeric(edge_info$WEIGHT)

####
#4.Produce nod information and specify attributes
      #Note: attribute vectors work ONLY if node IDs follow the order of
            # 0:nTL ; Dummy ; IN1:INnTL ; OUT1:OUTnTL ; LOSS1:LOSSnTL ; Edge1
            #-->TO BE CHECKED
####
#Node IDs
NOD_ID<-c(edge_info[,"FROM"], edge_info[,"TO"])
NOD_ID<-unique(NOD_ID)
#Labels
labelnames<-c("","D+I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII")

if(length(biomass)>1){
  BIOMlab<-BIOM
logBIOM<-floor(log10(BIOM))
logBIOM[which(logBIOM>=0)]<-0
logBIOM[which(is.infinite(logBIOM)==TRUE)]<-0
BIOMlab[which(logBIOM>=0)]<-round(BIOMlab[which(logBIOM>=0)],digits=2)
for(i in 1:max(abs(logBIOM),na.rm=TRUE)){
  selectLABS<-which(abs(logBIOM)==i)
  BIOMlab[selectLABS]<-round(BIOMlab[selectLABS],digits=(i+1))}
labelnames[c(2:(nTL+1))]<-paste(labelnames[c(2:(nTL+1))],BIOMlab,sep="\n")}

LABS<-rep("",length(NOD_ID))
LABS[c(1:(nTL+1))]<-labelnames[c(1:(nTL+1))]
#Fillcolor
FILL<-rep("white",length(NOD_ID))
FILL[c(2:(nTL+1))]<-"darkgoldenrod3"
#Color
COL<-rep("white",length(NOD_ID))
COL[c(2:(nTL+1))]<-"black"
#Coordinates
#Y-COORD:
if(migr==TRUE){
  COORDY<-c(rep(6,(nTL+1)),rep(5,(nTL+1)),rep(7.5,(nTL*3)),6)
}else{COORDY<-c(rep(6,(nTL+1)),rep(5,(nTL+1)),rep(7.5,(nTL*2)),6)}
#X-COORD:
X.LindS<-seq(from=0,to=12,length.out=(nTL+1))
if(migr==TRUE){
  inps<-X.LindS[-1]-(0.25*X.LindS[2])
  exps<-X.LindS[-1]+(0.25*X.LindS[2])
  X.ENV<-c(inps,exps)
}else{X.ENV<-X.LindS[-1]}
X.RESP<-X.LindS[-1]+(0.5*X.LindS[2])
X.Edge1<-X.LindS[1]+(0.6*X.LindS[2])
X.Dummy<-c(X.LindS[c(2:(nTL+1))],X.Edge1)
COORDX<-c(X.LindS,X.Dummy,X.ENV,X.RESP,X.Edge1)
#Node size
if(nTL<12){
WIDTH<-rep(0.75,length(NOD_ID))
}else{WIDTH<-rep(0.6,length(NOD_ID))}
HEIGHT<-rep(0.75,length(NOD_ID))
for(i in 1:nTL){
  mydummy<-paste("Dummy",i,sep="")
WIDTH[which(NOD_ID==mydummy)]<-HEIGHT[which(NOD_ID==mydummy)]<-0.1}
WIDTH[which(NOD_ID=="eDummy")]<-HEIGHT[which(NOD_ID=="eDummy")]<-0.1
WIDTH[which(NOD_ID=="Edge1")]<-HEIGHT[which(NOD_ID=="Edge1")]<-0.1
#Node style
STYLE<-rep("invis",length(NOD_ID))
STYLE[c(2:(nTL+1))]<-"filled"

####
#5.Add edge attributes to edge list
####
#Label
eLABS<-edge_info[,"WEIGHT"]
logLABS<-floor(log10(eLABS))
logLABS[which(logLABS>=0)]<-0
logLABS[which(is.infinite(logLABS)==TRUE)]<-0
eLABS[which(logLABS>=0)]<-round(eLABS[which(logLABS>=0)],digits=0)
for(i in 1:max(abs(logLABS),na.rm=TRUE)){
  selectLABS<-which(abs(logLABS)==i)
  eLABS[selectLABS]<-round(eLABS[selectLABS],digits=(i+1))}
eLABS[which(is.na(eLABS)==TRUE)]<-""
#Labelfontcolor
LABSCOL<-rep("black",nrow(edge_info))
LABSCOL[which(whichTO=="LOSS")]<-"azure4"
#Color
eCOL<-rep("black",nrow(edge_info))
eCOL[which(whichTO=="LOSS")]<-"azure4"
#Style
eSTYLE<-rep("solid",nrow(edge_info))#or "bold"
eSTYLE[which(whichTO=="LOSS")]<-"dashed"
#Headport/Tailport
#-->"n","ne","e","se","s","sw","w","nw","c","_"
HEADPORT<-rep(NA,nrow(edge_info))
HEADPORT[which(whichTO=="Dummy")]<-"e" #Return to detritus
HEADPORT[which(whichTO=="eDummy")]<-"w" #Return to detritus
HEADPORT[which(whichTO=="Edge1")]<-"w" #Return to detritus
HEADPORT[which(whichTO=="LOSS")]<-"s"  #Respiration
HEADPORT[which(whichFROM=="IN")]<-"nw" #Import
HEADPORT[which(whichTO=="OUT")]<-"s"  #Export
HEADPORT[which(is.na(HEADPORT)==TRUE)]<-"w" #Lindeman Spine

TAILPORT<-rep(NA,nrow(edge_info))
TAILPORT[which(whichTO=="Dummy")]<-"se" #Return to detritus
TAILPORT[which(whichTO=="eDummy")]<-"e" #Return to detritus
TAILPORT[which(whichTO=="Edge1")]<-"w" #Return to detritus
TAILPORT[which(whichTO=="LOSS")]<-"ne"  #Respiration
TAILPORT[which(whichFROM=="IN")]<-"s" #Import
TAILPORT[which(whichTO=="OUT")]<-"n"  #Export
TAILPORT[which(is.na(TAILPORT)==TRUE)]<-"e" #Lindeman Spine

#Arrowhead
ARROWHEAD<-rep("normal",nrow(edge_info))
ARROWHEAD[which(whichTO=="Dummy")]<-"half-open" #dot;odot
ARROWHEAD[which(whichTO=="eDummy")]<-"normal" 
ARROWHEAD[which(whichTO=="Edge1")]<-"none"

#Numeric Nod IDs
edge_info$ID<-c(1:nrow(edge_info))
numFROM<-data.frame(FROM=NOD_ID,num.FROM=c(1:length(NOD_ID)))
numTO<-data.frame(TO=NOD_ID,num.TO=c(1:length(NOD_ID)))
edge_info<-merge(edge_info,numFROM,by="FROM")
edge_info<-merge(edge_info,numTO,by="TO")
edge_info<-edge_info[order(edge_info[,"ID"],decreasing=FALSE),]
####
#6.Create graph
####
nod_df <- create_node_df(n = length(NOD_ID),
                         fillcolor = FILL, style=STYLE,
                         color = COL,width=WIDTH,height=HEIGHT,
                         label = LABS,x = COORDX,y = COORDY,
                         fontcolor="black",fontname="times-bold",shape="rectangle")
edge_df <- create_edge_df(from = edge_info[,"num.FROM"], to = edge_info[,"num.TO"],
                          penwidth = 1.25, color=eCOL,style=eSTYLE,
                          label=eLABS,fontcolor=LABSCOL,fontsize=9,
                          arrowsize=0.85, arrowhead=ARROWHEAD,tailport=TAILPORT,
                          headport=HEADPORT,labelfloat=TRUE,fontname="times-bold")


plotI<-create_graph(nodes_df = nod_df, edges_df = edge_df)%>%
  add_global_graph_attrs(attr="splines",
                         value="spline",
                         attr_type="graph")
#render_graph(plotI)

###########TROPHIC TRANSFORMATION MATRIX#################################
#-->stacked barplot of feeding partitioning among discrete trophic levels
#   for each compartment
#-->requires input of living compartment names


#1.Re-organize data for plotting
Trophic.level<-Compartment<-Feeding.coef<-c()
names.living<-colnames(A.matrix)
for(i in 1:ncol(A.matrix)){
  Trophic.level<-c(Trophic.level,c(1:nTL))
  Compartment<-c(Compartment,rep(names.living[i],nTL))
  Feeding.coef<-c(Feeding.coef,A.matrix[,i])}
plotA<-data.frame(Trophic.level=Trophic.level,Compartment=Compartment,
                  Feeding.coef=Feeding.coef)
#2.Create color and label vectors to produce plot
#12 predefined colors
A.col<-c("seagreen","gold","magenta","royalblue","red","turquoise3","darkgoldenrod4",
         "chocolate3","grey","yellow","pink","aquamarine4")
col.plotA<-A.col[c(1:nTL)]
#12 predefined labels
A.labs<-c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII")
labs.plotA<-A.labs[c(1:nTL)]
#Width of bars
if(ncol(A.matrix)>=70){A.barsize<-0.5
}else if(ncol(A.matrix)>=50){A.barsize<-0.7
}else if(ncol(A.matrix>=30)){A.barsize<-1
}else{A.barsize<-1.25}
#3.Produce barplot
plotA$Trophic.level<-factor(plotA$Trophic.level,levels=c(nTL:1))
plotII<-ggplot(plotA,aes(x=Compartment,y=Feeding.coef,fill=Trophic.level))+
  theme_minimal()+
  geom_bar(position = "stack", stat="identity", width = A.barsize, color="gray21")+
  scale_fill_manual("Trophic level",values=col.plotA[c(nTL:1)],
                    breaks=c(nTL:1),labels=labs.plotA[c(nTL:1)])+
  ylab("Feeding activity")+
  theme(legend.title=element_text(size=9, face="bold", color="gray21"),
        legend.text =element_text(size=10, face="bold", color="azure4"),
        axis.title.y=element_text(size=10, face="bold", color="gray21"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=7,angle=90,color="azure4"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###########TROPHIC EFFICIENCY PYRAMID####################################
#1.Adjust trophic efficiency data
tEFF<-TROPHIC[[11]][c(1:(nTL-1))]
tEFF<-round(tEFF,digits=3)
#2.Produce all relevant objects to create pyramid
label.tEFF<-A.labs[c(1:nTL)]
elabel.tEFF<-paste(tEFF*100,"%",sep="")
col.tEFF<-A.col[c(1:nTL)]
shape.tEFF<-c(rep("trapezium",(nTL-1)),"triangle")
y.tEFF<-seq(from=0,to=9,length.out = nTL)
x.tEFF<-rep(6,nTL)
height.tEFF<-y.tEFF[2]*0.55
width.tEFF<-c(20,20*tEFF)
#3.Produce trophic pyramid
nod_df <- create_node_df(n = nTL,
                         fillcolor = col.tEFF, style="filled",
                         color = "black",width=width.tEFF,height=height.tEFF,
                         label = label.tEFF,x = x.tEFF,y = y.tEFF, fontsize=12,
                         fontcolor="black",fontname="times-bold",shape=shape.tEFF)

edge_df <-create_edge_df(from = c(1:(nTL-1)), to = c(2:nTL),
                         penwidth = 30*tEFF, label=elabel.tEFF,
                         fontcolor= "black", fontsize=12, fontname="times-bold",                                      
                         color = col.tEFF)

plotIII<-create_graph(nodes_df = nod_df, edges_df = edge_df)
#render_graph(plotIII)
###########GRAZING CHAIN#################################################
#(?)
###########FUNCTION OUTPUT#################################################
TROPHPLOTS<-list(plotI,plotII,plotIII)
names(TROPHPLOTS)<-c("Lindeman Spine","Feeding activity of living compartments",
                     "Trophic Efficiency")
return(TROPHPLOTS)
}#END OF PLOT.TROPHIC FUNCTION
