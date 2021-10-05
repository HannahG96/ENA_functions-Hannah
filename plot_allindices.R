###############################VISUALIZE NETWORK INDICES########################################

####################FUNCTION#####################################
#-->requires input of Z_cs,T_cs,E_cs,R_cs,nl
#-->balance=c("unbalanced","inp","out","avg","io","oi","avg2")
plot.allindices<-function(Z_cs,T_cs,E_cs,R_cs,nl,balance="unbalanced"){

  library("ggplot2")

#(i)Create all indice-relevant objects
comps<-ncol(T_cs)#number of compartments
living<-comps-nl#number of living compartments
#(ii)Balance network (OR NOT="unbalanced")
if(balance!="unbalanced"){
  BAL<-network.balance(Z_cs=Z_cs,E_cs=E_cs,R_cs=R_cs,T_cs=T_cs,
                       method=balance)
#Recalculate transfer matrix and vectors:
  Z_f1<-BAL[[3]]
  T_f1<-BAL[[6]]
  E_f1<-BAL[[4]]
  R_f1<-BAL[[5]]
  if(length(Z_f1)==0){
    print("ERROR in Network balancing procedure: det(X)=0. Abort function. 
          Use `balance=unbalanced´ to avoid this error.")
    invokeRestart("abort")}
}else{
  Z_f1<-Z_cs
  T_f1<-T_cs
  E_f1<-E_cs
  R_f1<-R_cs }
#(iii)Calculate all information indices
#=c("TST","Development Capacity (DC)","Ascendency (A)", "A/DC",
#   "Average Mutual Information","Overhead on Imports","Overhead on Exports",
#   "Dissipative Overhead","Redundancy","Total Overhead (OV)","OV/TST",
#   "Internal Capacity","Internal Ascendency","Overall Connectance",
#   "Intercompartmental Connectance","Food web Connectance")
indices<-c("TST", "DC", "H", "A", "A(%)", "AMI",									
           "OI", "OE", "OD", "R", "O", "OI(%)", "OE(%)", "OD(%)", "R(%)", 
           "O(%)", "Hc", "int_Hc","IC", "IA", "IR","OC", "ICC", "FWC") #24
## --- (1) total system throughput
i1 <- Tst(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
##
## --- (2) development capacity
i2 <- dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
##
## --- (3) uncertainty (Shannon's index of diversity, H = DC/TST)
i3 <- i2/i1
##
#Ascendency
output <- Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
## --- (4) ascendency (absolute value)
i4<-unname(output[1])
##
## --- (5) ascendency (ratio)
i5 <- unname(output[2])
##
## --- (6) average mutual information
i6 <- Ami(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
##
##
#Overhead
output<-Overhead(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE, show.components=TRUE)
##
## --- (7) overhead on imports (absolute value)
i7 <- output[3,1]
##
## --- (8) overhead on exports (absolute value)
i8 <- output[4,1]
##
## --- (9) overhead on dissipations (absolute value)
i9 <- output[5,1]
##
## --- (10) redundancy (absolute value)
i10 <- output[2,1]
##
## --- (11) total overhead (absolute value)
i11 <- output[1,1]
##
## --- (12) overhead on imports (relative value)
i12 <- output[3,2]
##
## --- (13) overhead on exports (relative value)
i13 <- output[4,2]
##
## --- (14) overhead on dissipations (relative value)
i14 <- output[5,2]
##
## --- (15) redundancy (relative value)
i15 <- output[2,2]
##
## --- (16) total overhead (relative value)
i16 <- output[1,2]
##
## --- (17) residual diversity (Hc = O/TST)
i17 <- i11/i1
##
## --- (18) internal diversity (int_Hc = R/TST)
i18 <- i10/i1
##
##
##
##
## --- (19) internal capacity
i19 <- Internal.dC(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
##
## --- (20) internal ascendency (relative value)
i20 <- Internal.Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)[1]
##
## --- (21) internal redundancy (relative value)
i21 <- i10/i20
##
## --- (22) overall connectance
i22 <- connectance(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, type="whole")
##
## --- (23) intercompartmental connectance
i23 <- connectance(T_cs=T_f1, type="intercompartmental")
##
## --- (24) food web connectance
i24 <- connectance(T_cs=T_f1,nl=nl, type="foodweb")

inds<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,
        i21,i22,i23,i24)
plot.indices<-data.frame(Value=log10(inds),Indice=indices)
#(iv)Define value labels of indices
labs.inds<-inds
logLABS<-floor(log10(labs.inds))
logLABS[which(logLABS>=0)]<-0
logLABS[which(is.infinite(logLABS)==TRUE)]<-0
labs.inds[which(logLABS>=0)]<-round(labs.inds[which(logLABS>=0)],digits=2)
for(i in 1:max(abs(logLABS),na.rm=TRUE)){
  selectLABS<-which(abs(logLABS)==i)
  labs.inds[selectLABS]<-round(labs.inds[selectLABS],digits=(i+1))}
plot.indices$Label<-labs.inds
#Produce indice barplot
plot.indices$Indice<-factor(plot.indices$Indice,levels=indices)
INDSPLOT<-ggplot(plot.indices, aes(x = Value, y = Indice,label=Label)) + 
  theme_minimal()+
  geom_bar(stat = "identity",width=0.7,fill="darkgoldenrod3",color="azure4") +
  xlab("Log10(Indice)")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.01, color="azure4"),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=9,color="azure4"),
        axis.title.x = element_text(face="bold", size=8, color="grey21"),
        axis.text.y = element_text(size=10,face="bold",color="grey21"))+
  geom_text(aes(fontface=2),size=3,color="darkmagenta",vjust=0.25,hjust=0.5)

###
#FUNCTION OUTPUT
return(INDSPLOT)
}#END OF PLOT.ALLINDICES FUNCTION