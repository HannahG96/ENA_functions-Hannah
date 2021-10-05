####1####
#-->RELEVANT???
######################################################################################################
###												   												   ###
### Functions to convert dataframe input data (in a similar SCOR format) into vectors and matrices ###
###												   												   ###
######################################################################################################


## Conversion of SCOR vectors (available as input dataframes) into whole vectors. It is useful for (B), (Z), (E) and (R).

#vec.ena <- function(df,x){
  #V <- rep(0,x)
  #for(i in 1:nrow(df))V[df[i,1]] <- df[i,2]
  #return(V)
#}


## Conversion of SCOR matrices (available as input dataframes) into whole matrix. 
## Mainly useful for the matrix of inter-compartmental[T] exchanges.

#matr.ena <- function(df,x){
  #tot_el <- x^2 
  #M <- matrix(rep(0,tot_el), nrow=x)
  #for(i in 1:nrow(df))M[df[i,1],df[i,2]] <- df[i,3]
  #return(M)
#}


## Conversion of the dataframe input (.csv) file into vectors and matrices for ENA.
## Filname should be "quoted". Empty cells are converted into zeros.

df.ena <- function(filename){
  ch <- nchar(filename)
  n.ch <- ch - 2
  nch.frag <- substring(filename,n.ch)
  #
  {
    if(nch.frag == "txt")AA <- read.table(filename,header=TRUE)
    #
    else AA <- read.csv2(filename,header=TRUE)
  }
  BB <- as.matrix(AA)
  #
  if(length(which(is.na(AA)==TRUE))!=0){
    aa <- which(is.na(AA)==TRUE)
    BB[aa] <- 0
  }
  #			
  # We create a list including information on: 1) total number of nodes; 2) number of living nodes;
  # 3) number of non-living nodes; 4) name of the nodes; 5) biomasses; 6) imports; 7) exports;
  # 8) respirations; 9) intercompartmental exchanges; 10) exchanges between living compartments
  #
  BB.L <- as.list(rep(NA,10))
  #
  BB.L[[1]] <- nrow(BB)
  tt <- BB.L[[1]]
  BB.L[[4]] <- t(t(colnames(BB)[6:(5+tt)]))
  {
    if(length(which(BB[,1]>0))!=0){
      BB.L[[2]] <- length(which(BB[,1]>0))
    }
    #
    else BB.L[[2]] <- 0
  }
  BB.L[[5]] <- BB[,2]
  BB.L[[6]] <- BB[,3]
  BB.L[[7]] <- BB[,4]
  BB.L[[8]] <- BB[,5]
  BB.L[[9]] <- BB[c(1:tt),c(6:(5+tt))]
  #
  {
    if(length(which(BB[,1]<0))!=0){
      BB.L[[3]] <- length(which(BB[,1]<0))
      neg <- which(BB[,1]<0)
      #
      BB.L[[10]] <- BB.L[[9]][-c(neg),-c(neg)]
    }
    #
    else {
      BB.L[[3]] <- 0
      BB.L[[10]] <- BB.L[[9]]
    }
  }
  #
  return(BB.L)
}

####2####
########################################################################################################################################
###																     ###
### Functions to depict all the flows (import, export, dissipations and inter-compartmentmental exchanges) as one big matrix [T_tot] ###
### In this case the input considered are vector and matrices (not dataframes)							     ###
### In the first function ("bigmatr.ena") the output will also show the name of each compartment, while this addional information    ###
### is not available in the case of the raw result of "bigmatrNOnames.ena"							     ###
###																     ###
########################################################################################################################################

bigmatr.ena <- function(cnames,inp,inter,outt,diss){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  zeros_up <- rep(0, tot_row)
  core <- cbind(zeros_up, inter, outt, diss)
  inputs <- c(0, inp, 0, 0)
  zeros_down <- rep(0, tot_row+3)
  final <- rbind(inputs, core, zeros_down, zeros_down)
  #
  names_tot <- c("IMPORT", cnames, "EXPORT", "RESP")
  rownames(final) <- names_tot
  colnames(final) <- names_tot
  return(final)
}


#bigmatrNOnames.ena <- function(inp,inter,outt,diss){
  #tot_row <- length(inp)
  #tot_el <- tot_row^2
  #zeros_up <- rep(0, tot_row)
  #core <- cbind(zeros_up, inter, outt, diss)
  #inputs <- c(0, inp, 0, 0)
  #zeros_down <- rep(0, tot_row+3)
  #final <- rbind(inputs, core, zeros_down, zeros_down)
  #colnames(final) <- NULL
  #rownames(final) <- NULL
  #
  #return(final)
#}

####3####
#####################################################################################
###										  										  ###
### Functions to calculate Partial Feeding Matrix [G] and Partial Host Matrix [F] ###
###																				  ###
#####################################################################################

## The function "partial.feeding" convert data from (Z) and [T] into partial feeding matrix [G].

partial.feeding <- function(inp,inter){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  G <- matrix(rep(0,tot_el), nrow = tot_row)
  IN <- rep(0,tot_row)	
  #
  for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + sum(inp[j]))
  #
  for (i in 1:tot_row)
    for (j in 1:tot_row)if(sum(IN[j])!=0)G[i,j] <- inter[i,j]/sum(IN[j])
  return(G)
}


## The function "partial.host" convert data from [T], (E) and (R) into partial host matrix [F].

partial.host <- function(inter,outt,diss){
  tot_row <- length(outt)
  tot_el <- tot_row^2
  FF <- matrix(rep(NA,tot_el), nrow = tot_row)
  OUT <- rep(0,tot_row)
  #
  for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + sum(outt[i]) + sum(diss[i]))
  #
  for (i in 1:tot_row)
    for (j in 1:tot_row)if(sum(OUT[i])!=0)FF[i,j] <- inter[i,j]/sum(OUT[i])
  return(FF)
}


####4####
#######################################################################
###																    ###
### Functions to normalize vectors (that could be (Z), (E) and (R)) ###
###																    ###
#######################################################################

## The function to normalize input vector is called "normalize.invector" and it is used for (Z).

normalize.invector <- function(inp,inter){
  tot_row <- length(inp)
  NZ <- rep(0,tot_row)
  IN <- rep(0,tot_row)
  #
  for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + sum(inp[j]))
  #
  for (j in 1:tot_row)if(sum(IN[j])!=0)NZ[j] <- inp[j]/sum(IN[j])
  return(NZ)
}


## The function to normalize both export and dissipation vectors is called "normalize.outvector" and it is used for (E) and (R).

normalize.outvector <- function(inter,outt,diss){
  tot_row <- length(outt)
  NO <- matrix(rep(0,tot_row*2), nrow = tot_row)
  OUT <- rep(0,tot_row)
  #
  for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + sum(outt[i]) + sum(diss[i]))
  #
  for (i in 1:tot_row)if(sum(OUT[i])!=0){NO[i,1] <- outt[i]/sum(OUT[i])
  NO[i,2] <- diss[i]/sum(OUT[i])
  }
  colnames(NO) <- c("NE", "NR")
  return(NO)
}

####5####
########################################################################################################
###																								     ###
### Functions to calculate Inverse Matrix of Leontief [LEO] and Inverse Matrix of Augustinovic [AUG] ###
###												    												 ###
########################################################################################################

## The function "leontief" estimates [LEO] when (Z) and [T] are available.

leontief <- function(inp,inter){
  tot_row <- length(inp)
  G <- partial.feeding(inp,inter)
  #
  ddd <- det(diag(1,tot_row) - G)
  if(ddd == 0)return(NA)
  #
  LEO <- solve(diag(1,tot_row) - G)
  return(LEO)
}


## The function "augustinovic" estimates [AUG] when [T], (E) and (Z) are available.

augustinovic <- function(inter,outt,diss){
  tot_row <- length(outt)
  FF <- partial.host(inter,outt,diss)
  #
  ddd <- det(diag(1,tot_row) - t(FF))
  if(ddd == 0)return(NA)
  #
  AUG <- solve(diag(1,tot_row) - t(FF))
  return(AUG)
}

####6####
#####################################################################################
###																				  ###
### Functions to build Total Dependency [IG] and Total Contribution [IF] matrices ###
###										  										  ###
#####################################################################################

## The function "total.dependency" convert data from (Z), [T], (E) and (R) into [IG].

total.dependency <- function(inp,inter,outt,diss){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  #
  IN <- rep(0,tot_row)	
  for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + inp[j])
  OUT <- rep(0,tot_row)
  for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + outt[i] + diss[i])
  #
  inv.leontief <- leontief(inp,inter)
  if(length(inv.leontief) == 1)return(NA)
  #
  leo <- (inv.leontief - diag(1,tot_row))
  interm1 <- rep(NA,tot_row)
  for(i in 1:tot_row)interm1[i] <- IN[i]/(inv.leontief[i,i] * OUT[i])
  #
  DEP <- matrix(rep(NA,tot_el), nrow = tot_row)
  for (i in 1:tot_row)
    for(j in 1:tot_row)DEP[i,j] <- leo[i,j] * interm1[i]
  return(DEP)
}


## The function "total.contribution" convert data from (Z), [T], (E) and (R) into [IF].

total.contribution <- function(inp,inter,outt,diss){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  #
  IN <- rep(0,tot_row)	
  for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + inp[j])
  OUT <- rep(0,tot_row)
  for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + outt[i] + diss[i])
  #
  inv.augustinovic <- augustinovic(inter,outt,diss)
  if(length(inv.augustinovic) ==1)return(NA)
  #
  aug <- (inv.augustinovic - diag(1,tot_row))
  interm1 <- rep(NA,tot_row)
  for(i in 1:tot_row)interm1[i] <- IN[i]/(inv.augustinovic[i,i] * OUT[i])
  #
  CONTR <- matrix(rep(NA,tot_el), nrow = tot_row)
  for (i in 1:tot_row)
    for(j in 1:tot_row)CONTR[i,j] <- aug[i,j] * interm1[i]
  return(t(CONTR))
}

####7####
########################################################
###												     ###
### Functions to carry out INPUT and OUTPUT analyses ###
###						    						 ###
########################################################

## The function "input.analysis".

input.analysis <- function(inp,inter,outt,diss){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  F <- partial.host(inter,outt,diss)
  NE <- normalize.outvector(inter,outt,diss)[,1]
  NR <- normalize.outvector(inter,outt,diss)[,2]
  #
  ddd <- det(diag(1,tot_row) - t(F))
  if(ddd == 0)return(NA)
  IFT <- solve(diag(1,tot_row) - t(F))						
  #
  # non-zero elements in the input vector (Z)
  #
  count<-0
  for (i in 1:tot_row)if(inp[i]!=0)count <- count+1
  #
  # list for each input different from zero (Z)_i != 0
  #
  final_input <- vector("list", count)
  for(i in 1:count){final_input[[i]] <- vector("list", 4)
  final_input[[i]][[1]] <- rep(0,tot_row)
  final_input[[i]][[2]] <- matrix(rep(0,tot_el), nrow = tot_row)
  final_input[[i]][[3]] <- rep(0,tot_row)
  final_input[[i]][[4]] <- rep(0,tot_row)
  }
  #
  step <- 1
  for(i in 1:tot_row)if(inp[i]!=0){final_input[[step]][[1]][i] <- 1
  step <- step + 1
  }
  #
  intermediate_step <- vector("list", count)
  for (i in 1:count)intermediate_step[[i]] <- IFT %*% final_input[[i]][[1]]
  #
  # defining normalized inter-compartmental exchanges, export and dissipations
  #
  for(k in 1:count){
    for(i in 1:tot_row){
      final_input[[k]][[3]][i] <- intermediate_step[[k]][i] * NE[i]
      final_input[[k]][[4]][i] <- intermediate_step[[k]][i] * NR[i]
      for(j in 1:tot_row)final_input[[k]][[2]][i,j] <- intermediate_step[[k]][i] * F[i,j]
    }
    k <- k+1
  }
  #
  return(final_input)			
}


## The function "output.analysis".

output.analysis <- function(inp,inter,outt){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  G <- partial.feeding(inp,inter)
  NZ <- normalize.invector(inp,inter)
  #
  ddd <- det(diag(1,tot_row) - G)
  if(ddd == 0)return(NA)
  IGT <- solve(diag(1,tot_row) - G)
  #
  # non-zero elements in the export vector (E)
  #
  count<-0
  for (i in 1:tot_row)if(outt[i]!=0)count <- count+1
  #
  # list for each output different from zero (E)_i != 0
  #
  final_output <- vector("list", count)
  for(i in 1:count){final_output[[i]] <- vector("list", 4)
  final_output[[i]][[1]] <- rep(0,tot_row)
  final_output[[i]][[2]] <- matrix(rep(0,tot_el), nrow = tot_row)
  final_output[[i]][[3]] <- rep(0,tot_row)
  final_output[[i]][[4]] <- rep(0,tot_row)
  }
  #
  step <- 1
  for(i in 1:tot_row)if(outt[i]!=0){final_output[[step]][[3]][i] <- 1
  step <- step + 1
  }
  #
  intermediate_step <- vector("list", count)
  for (i in 1:count)intermediate_step[[i]] <- IGT %*% final_output[[i]][[3]]
  #
  # defining normalized import and inter-compartmental exchanges
  #
  for(k in 1:count){
    for(i in 1:tot_row){
      final_output[[k]][[1]][i] <- intermediate_step[[k]][i] * NZ[i]
      for(j in 1:tot_row)final_output[[k]][[2]][j,i] <- intermediate_step[[k]][i] * G[j,i]
    }
    k <- k+1
  }
  #
  return(final_output)
}

####8####
## The Finn Cycling Index (FCI) is computed as the ratio between cycled throughflows and total system throughput.
## Using the "FCI" function, three different versions of the index are calculated. In the first case ("finn80"),
## diagonal elements of the Leontief matrix are used to assess the amount of material cycling within the ecosystem,
## - (Sii - 1)/Sii - and TST stands for the sum of internal transfers plus imports - TST* = sum(Z) + sum([T]).
## In the Ulanowicz's (1986) version ("ulan86"), TST* is replaced by the more widespread formulation (with total system
## throughput corresponding to the total amount of flows occurring in the system) - TST = sum(Z) + sum([T]) + sum(E) + sum(R).
## Finally, "szyr87" is measured using the diagonal elements of the total dependency matrix, with the TST computed as the sum
## of all the ecosystem flows (Szyrmer & Ulanowicz, 1987).

FCI <- function(inp,inter,outt,diss){
  tot_flows <- sum(inp)+sum(inter)+sum(outt)+sum(diss)
  tot_finn <- sum(inter) + sum(inp)
  xx <- ncol(inter)
  #
  IN <- apply(inter,2,sum) + inp
  #
  LEONTIEF <- leontief(inp,inter)
  if(length(LEONTIEF) == 1)return(NA)
  LEONT <- diag(LEONTIEF)
  #
  DEPENDENCY <- total.dependency(inp,inter,outt,diss)
  if(length(DEPENDENCY) == 1)return(NA)
  #
  DEPEN <- diag(DEPENDENCY)
  #
  finn80 <- rep(0,xx)
  szyr87 <- rep(0,xx)
  #
  for(i in 1:xx){
    if(LEONT[i]!=0)finn80[i] <- IN[i] * (LEONT[i]-1)/LEONT[i]
    #
    szyr87[i] <- IN[i] * DEPEN[i]
  }
  #
  ulan86 <- finn80	
  #
  finn.OK <- rep(0,3) 
  finn.OK[1] <- sum(finn80)/tot_finn		## Finn, 1980
  finn.OK[2] <- sum(ulan86)/tot_flows		## Ulanowicz, 1986
  finn.OK[3] <- sum(szyr87)/tot_flows		## Szyrmer & Ulanowicz, 1987
  #
  return(finn.OK)
}


## Comprehensive Cycling Index (CCI), the updated and improved version of the FCI. It contains the contribution to cysling given by:
## (a) simple cycles; (b)compound paths and (c) and compound cycles. As CCI is cumbersome to compute, but correlated with FCI, we used
## a constant coefficient (1.142) for its rough estimation from FCI.

CCI <- function(inp,inter,outt,diss){
  k <- 1.142
  finn <- FCI(inp,inter,outt,diss)[2]
  if(is.na(finn) == TRUE)return(NA)
  allesina <- finn * k
  return(allesina)
}

####9####
###########################################################################################################################################
###																	###
### The following functions are inspired to the pubblication of Puccia and Ulanowicz (1990) that analyzed indirect effects in ecosystem	###
### This is not a typical tool of network analysis and is also available at the Ulanowicz's website as a program called IMPACTS		###
### The function "IMA" was introduced by Vasas and Jordan (2006).									###
###																	###
###########################################################################################################################################

## The output of the function "Q.TI" is the Matrix of Direct Trophic Interactions [Q].

Q.TI <- function(inp,inter,outt,nl){
  pf <- partial.feeding(inp,inter)
  #
  tot.l <- length(inp)
  living <- tot.l - nl
  #
  out.l <- apply(inter,1,sum)
  OOO.l <- out.l + outt
  f.ast <- matrix(rep(0,tot.l^2),nrow=tot.l)
  #
  for(i in 1:tot.l)
    for(j in 1:living){
      if(OOO.l[i]!=0)f.ast[i,j] <- inter[i,j]/OOO.l[i]
    }
  #
  pht <- t(f.ast)
  Q <- (pf - pht)
  #
  return(Q)
}	


## The output of the function "M.TI" is the Matrix of Total Trophic Impacts [M].

M.TI <- function(inp,inter,outt,nl){
  xx <- ncol(inter)
  living <- xx - nl
  #
  I <- diag(1,xx)
  Q <- Q.TI(inp,inter,outt,nl)
  #
  ddd <- det(I-Q)
  if(ddd == 0)return(NA)
  #
  IQinv <- solve(I-Q)
  M <- (IQinv - I)
  #
  return(M)
}


## The "IMA" function is used for summing the absolute values of effects as measured by partial feeding and partial host matrices.
## In this way, we investigate the absolute effects of species and tropho-species (both positive and negative), estimating their
## global interacting power (key species). Using absolute values does not lead to the loss of information on sign structure.

IMA <- function(inp,inter,outt,nl){
  if(length(M.TI(inp,inter,outt,nl)) == 1)return(NA)
  TTI.m <- M.TI(inp,inter,outt,nl)
  P.TI <- abs(TTI.m)
  IMA.v <- apply(P.TI,1,sum)
  #
  return(IMA.v)
}

####10####
####################################################################################################################################
###																 																 ###
### Trophic analysis both in the classical framework of Canonical Trophic Aggregation (CTA) and in the "extended" version of CTA ###
###																																 ###
####################################################################################################################################
#What is returned?
#A list of 13 elements with:
#(1)Lindeman Transformation Matrix
#(2)Vector of effective trophic levels of each species
#(3)Canonical Exports
#(4)Canonical Respiration
#(5)Grazing Chain
#(6)Returns to Detrital Pool
#(7)Detrivory
#(8)Input to Detrital Pool
#(9)Circulation within Detrital Pool
#(10)Lindeman Spine
#(11)Trophic Efficiencies
#(12)Detection of Migratory Imports 
#-->REQUIRES EXTENDED TROPHIC ANALYSIS:
#(13)Canonical Distribution of Inputs (only if migration)

#names(EL)<-c("Lindeman Transformation Matrix", "Effective Trophic Levels of each species", 
#"Canonical Exports", "Canonical Respirations", "Grazing Chain", 
#"Returns to Detrital Pool", "Detrivory", "Input to Detrital Pool",
#"Circulation within Detrital Pool", "Lindeman Spine", "Trophic Efficiency",
#"Detection of migration", "Canonical Imports")[c(1:length(EL))]

#colnames(EL[[1]])<-colnames(T_cs)[c(1:living)]#names of living compartments


## "Extended" version of the CTA. In this case the effect of imports with TPs far from 0 (i.e. prey migration) can be considered.
## It is also allowed to modulate the TL of "non-living imports" (TL = 0 or TL = 1).

TP.extendedCTA <- function(inp,inter,outt,diss,nl,TPimports,TPnonliving){
  xx <- ncol(inter)
  #
  # defining "living" imports (Z_living) with internal flows (T_living) 
  # and computing input from non_living nodes to living (K)
  #
  living <- xx - nl
  Z_living <- inp[c(1:living)]
  TTT <- inter[c(1:living),c(1:living)]
  E_living <- outt[c(1:living)]
  R_living <- diss[c(1:living)]
  #
  # exctracting the matrix of acyclic flows between living compartments
  #
  TTT.s1 <- cycling.analysis(Z_cs=Z_living,T_cs=TTT,E_cs=E_living,R_cs=R_living, 
                             select="Default", nmax=10000)
  #
  {
    if(TTT.s1[[1]] == 0) T_living <- TTT.s1[[2]]
    #
    else T_living <- TTT.s1[[7]]
  }
  #
  {
    if(nl!=0){
      LVSNL <- matrix(inter[c((living+1):xx),c(1:living)],nrow=nl,byrow=FALSE)
    }
    #
    else LVSNL <- matrix(rep(0,living),nrow=1)
  }
  #
  K <- apply(LVSNL,2,sum)
  #
  # vector for the non_living compartments TLs (TK). In this version all the TLs are set to 1
  #
  TK <- rep(0,living)
  for(i in 1:living){
    if(K[i]!=0)TK[i]<-1
  }
  #
  # checking non-living input with TLs equal to 0 (flows from non-living to primary producers);
  # the vector (TPnonliving) is composed of "living" elements
  #
  for(i in 1:living){
    if(K[i]!=0 & TPnonliving[i]==0){
      Z_living[i] <- Z_living[i] + K[i]
      K[i] <- 0
    }
  }
  #
  # building partial feeding matrix G_living [G] and normalized import vectors N_living (N) and NK_living (NK)
  #
  IN <- rep(NA,living)
  G <- matrix(rep(NA,living^2), nrow=living)
  N <- rep(NA,living)
  NK <- rep(NA,living)
  for (j in 1:living) IN[j] <- (sum(T_living[,j]) + Z_living[j] + K[j])
  #
  for (i in 1:living){
    for (j in 1:living){
      if(IN[j]==0) G[i,j] <- 0
      #
      else G[i,j] <- T_living[i,j]/IN[j]
    }
  }
  #
  for (j in 1:living){
    if(IN[j]==0) N[j] <- 0
    #
    else N[j] <- Z_living[j]/IN[j]
  }
  #
  for (j in 1:living){
    if(IN[j]==0) NK[j] <- 0
    #     	
    else NK[j] <- K[j]/IN[j]
  }
  #
  # calculation of [G] matrix powers (computation interrupted after "living" steps)
  #
  G.lista <- as.list(rep(NA,living))
  #
  G.lista[[1]] <- diag(rep(1,living))
  G.lista[[2]] <- G
  for(q in 3:living)G.lista[[q]] <- G%*%G.lista[[q-1]]
  #
  # checking external migratory imports (with TPs far from 0);
  # the vector (TPimports) is composed of "living" elements
  #
  levels <- 0
  for(i in 1:living){
    if((Z_living[i]!=0) & (TPimports[i]!=0)) levels <- levels+1
  }
  #
  # vectors containing TPs of the migratory imports only
  #
  EFFECTIVE <- rep(NA,levels)
  pat <- 1
  for(i in 1:living){
    if(TPimports[i]!=0){
      EFFECTIVE[pat] <- TPimports[i]
      pat <- pat+1
    }
  }
  #
  k <- 1
  initialize <- rep(0,living)
  NP <- rep(0,living)
  TO <- rep(0,living)
  NPC <- rep(0,living)
  TONE <- rep(0,living)
  ooo <- rep(NA,levels)
  eee <- rep(NA,levels)
  normLIST <- as.list(ooo)
  trophicLIST <- as.list(eee)
  #
  if(levels!=0){
    for(i in 1:levels){
      normLIST[[i]] <- initialize
      trophicLIST[[i]] <- initialize
    }
    for(i in 1:living){
      if(Z_living[i]!=0){
        if(TPimports[i]==0)NP[i] <- N[i]
        #
        else{
          normLIST[[k]][i] <- N[i]
          trophicLIST[[k]][i] <- TPimports[i]
          k <- k + 1
        }
      }
    }
  }
  #
  # list containing Transformation Matrices with a Normalized Import Vectors for each migratory import
  #
  TEN <- matrix(rep(0,living^2),nrow = living)
  tot_lev <- levels+2
  el <- rep(NA,tot_lev)
  WHOLE <- as.list(el)
  INPUTS <- as.list(el)		
  #		
  for (i in 1:tot_lev) WHOLE[[i]] <- TEN
  #
  for (i in 1:tot_lev){
    if(i==1)INPUTS[[i]] <- t(NP)
    #
    else{
      if(i==2)INPUTS[[i]] <- t(NK)
      #
      else INPUTS[[i]] <- t(normLIST[[i-2]])
    }
  }		
  #	
  # preparation of a Transformation Matrix for each migratory import
  # and its storage into the list called WHOLE
  #		
  for(i in 1:tot_lev){
    for(q in 1:living)WHOLE[[i]][q,] <- INPUTS[[i]]%*%G.lista[[q]]	
  }
  #
  # preparation of vectors used to calculate TPs of living nodes
  #
  TROPHICVEC <- as.list(el)
  for(i in 1:tot_lev){
    if(i==1)TROPHICVEC[[i]] <- c(1:living)
    #
    else{
      if(i==2)TROPHICVEC[[i]] <- c(2:(living+1))
      #
      else TROPHICVEC[[i]]<- EFFECTIVE[i-2] * rep(1,living) + c(1:living)
    }
  }		
  #
  # computation of partial TPs and their storage into a list called PARTIAL.list
  #
  PARTIAL.list <- as.list(el)	
  for(i in 1:tot_lev)PARTIAL.list[[i]] <- TROPHICVEC[[i]]%*%WHOLE[[i]]
  #
  # TPs of living nodes
  #
  tr.vec <- rep(0,living)
  for(i in 1:tot_lev)tr.vec <- tr.vec + PARTIAL.list[[i]]
  #
  final.effectiveTPs <- t(tr.vec)
  #
  spl <- pat + 1
  AAA <- as.list(rep(NA,pat))
  DDD <- as.list(rep(NA,pat))
  #
  for(i in 2:spl){
    infe <- trunc(TROPHICVEC[[i]][1]-1)
    supe <- trunc(TROPHICVEC[[i]][1])
    supe.el <- TROPHICVEC[[i]][1] - supe
    #
    {
      if(supe.el==0){
        ZEROs <- matrix(rep(0,(living*infe)),nrow = infe)
        cuts <- c((living + 1):(living+infe))
        AAA.s1 <- rbind(ZEROs,WHOLE[[i]])
        AAA[[i-1]] <- AAA.s1[-cuts,]
        DDD[[i-1]] <- matrix(rep(0,living^2),nrow = living)
      }
      else	{
        infe.el <- 1 - supe.el
        #
        M.infe <- WHOLE[[i]] * infe.el
        M.supe <- WHOLE[[i]] * supe.el
        #
        ZEROs.infe <- matrix(rep(0,(living * infe)),nrow = infe)
        cuts.infe <- c((living + 1):(living + infe))
        AAA.s1 <- rbind(ZEROs.infe,M.infe)
        AAA[[i-1]] <- AAA.s1[-cuts.infe,]
        #
        ZEROs.supe <- matrix(rep(0,(living * supe)),nrow = supe)
        cuts.supe <- c((living + 1):(living + supe))
        DDD.s1 <- rbind(ZEROs.supe,M.supe)
        DDD[[i-1]] <- DDD.s1[-cuts.supe,]
      }
    }
  }
  #
  ONE.LIND <- WHOLE[[1]]
  for(i in 1:pat)ONE.LIND <- ONE.LIND + AAA[[i]] + DDD[[i]]
  #
  EL <- as.list(rep(NA,12))
  EL[[1]] <- ONE.LIND														## Lindeman Transformation Matrix
  #
  if(nl!=0)final.effectiveTPs.w <- c(final.effectiveTPs,rep(1,nl))
  #
  EL[[2]] <- final.effectiveTPs.w											## Effective Trophic Levels of each species
  #																		## in presence of migratory imports
  {
    if(nl!=0){
      EL[[3]] <- c(EL[[1]] %*% t(t(E_living)),outt[c(living+1):xx])		## Canonical Exports
      R.det <- sum(diss[c(living+1):xx])
      R.det.v <- rep(0,nl)
      R.det.v[nl] <- R.det
      EL[[4]] <- c(EL[[1]] %*% t(t(R_living)),R.det.v)					## Canonical Respirations
      #
      DetPool <- matrix(inter[c(1:living),c((living+1):xx)],nrow=living,byrow=FALSE)
      EL6.s1 <- apply(DetPool,1,sum)
      EL[[6]] <- EL[[1]] %*% t(t(EL6.s1))									## Returns to Detrital Pool
      #
      EL[[8]] <- sum(inp[c((living+1):xx)])								## Input to Detrital Pool
      EL[[9]] <- sum(inter[c((living+1):xx),c((living+1):xx)])			## Circulation within Detrital Pool
    }
    else	{
      DetPool <- 0
      EL[[3]] <- EL[[1]] %*% t(t(E_living))
      EL[[4]] <- EL[[1]] %*% t(t(R_living))
      EL[[6]] <- 0
      EL[[8]] <- 0
      EL[[9]] <- 0
    }
  }
  #
  detritivory.s1 <- apply(LVSNL,2,sum)
  EL[[7]] <- sum(detritivory.s1[which(TPnonliving==1)])					## Detritivory
  #
  EL5.s1 <- apply(T_living,1,sum)
  EL5.s2 <- EL[[1]] %*% t(t(EL5.s1))
  lll <- living + 1
  IN.liv <- apply(T_living,2,sum)
  prim <- which(Z_living != 0 & IN.liv == 0)
  Z.prim <- sum(Z_living[prim])
  EL5.s3 <- c(Z.prim,EL5.s2)
  detr.to.pp <- sum(detritivory.s1[which(TPnonliving==0)])
  EL5.s3 <- EL5.s3[-lll]
  EL5.s3[1] <- EL5.s3[1] + detr.to.pp
  EL[[5]] <- EL5.s3									## The Grazing Chain
  #
  Det.not.p <- DetPool
  if(nl!=0)Det.not.p[prim,] <- rep(0,nl)
  lind.s1 <- EL[[5]]
  pr.pro <- which(detritivory.s1 != 0 & TPnonliving == 0)
  {
    if(length(pr.pro)==0)lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + EL[[8]]
    #
    else lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + sum(detritivory.s1[pr.pro])
  }
  lind.s1[2] <- lind.s1[2] + EL[[7]]
  EL[[10]] <- lind.s1									## Lindeman Spine
  #
  lind.eff <- rep(0,(living-1))
  for(q in 2:living){
    if(EL[[10]][q-1] > 0.0001){
      lind.eff[q-1] <- EL[[10]][q]/EL[[10]][q-1]
    }
  }
  EL[[11]] <- lind.eff								## Trophic Efficiencies
  #
  mig <- which(Z_living != 0 & IN.liv != 0)
  if(length(mig)==0)mig <- 0
  EL[[12]] <- mig										## Detecting Migratory Imports
  
colnames(EL[[1]])<-colnames(T_cs)[c(1:living)]#names of living compartments  
names(EL)<-c("Lindeman Transformation Matrix", "Effective Trophic Levels of each species", 
             "Canonical Exports", "Canonical Respirations", "Grazing Chain", 
             "Returns to Detrital Pool", "Detrivory", "Input to Detrital Pool",
             "Circulation within Detrital Pool", "Lindeman Spine", "Trophic Efficiency",
             "Detection of migration", "Canonical Imports")[c(1:length(EL))]
  
  return(EL)
}


## Canonical Trophic Aggregation (CTA) with all the imports from outside with Hypothetical TP = 0 and non living nodes set to TP = 1 
## (as primary producers). No migratory imports are considered and only a warning message is supplied in this case.

TP.CTA <- function(inp,inter,outt,diss,nl){
  xx <- ncol(inter)
  #
  # defining "living" imports (Z_living) with internal flows (T_living) 
  # and computing input from non_living nodes to living (K)
  #
  living <- xx - nl
  Z_living <- inp[c(1:living)]
  TTT <- inter[c(1:living),c(1:living)]
  E_living <- outt[c(1:living)]
  R_living <- diss[c(1:living)]
  #
  # exctracting the matrix of acyclic flows between living compartments
  #
  TTT.s1 <- cycling.analysis(Z_cs=Z_living,T_cs=TTT,E_cs=E_living,R_cs=R_living, 
                             select="Default", nmax=10000)
  #
  {
    if(TTT.s1[[1]] == 0) T_living <- TTT.s1[[2]]
    #
    else T_living <- TTT.s1[[7]]
  }
  #
  {
    if(nl!=0){
      LVSNL <- matrix(inter[c((living+1):xx),c(1:living)],nrow=nl,byrow=FALSE)
    }
    #
    else LVSNL <- matrix(rep(0,living),nrow=1)
  }
  #
  K <- apply(LVSNL,2,sum)
  #
  # vector for the non_living compartments TLs (TK). In this version all the TLs are set to 1
  #
  TK <- rep(0,living)
  for(i in 1:living){
    if(K[i]!=0)TK[i]<-1
  }
  #
  # building partial feeding matrix G_living [G] and normalized import vectors N_living (N) and NK_living (NK)
  #
  IN <- rep(NA,living)
  G <- matrix(rep(NA,living^2), nrow=living)
  N <- rep(NA,living)
  NK <- rep(NA,living)
  for (j in 1:living) IN[j] <- (sum(T_living[,j]) + Z_living[j] + K[j])
  #
  for (i in 1:living){
    for (j in 1:living){
      if(IN[j]==0) G[i,j] <- 0
      #
      else G[i,j] <- T_living[i,j]/IN[j]
    }
  }
  #
  for (j in 1:living){
    if(IN[j]==0) N[j] <- 0
    #
    else N[j] <- Z_living[j]/IN[j]
  }
  #
  for (j in 1:living){
    if(IN[j]==0) NK[j] <- 0
    #     	
    else NK[j] <- K[j]/IN[j]
  }
  #
  G.lista <- as.list(rep(NA,living))
  #
  G.lista[[1]] <- diag(rep(1,living))
  G.lista[[2]] <- G
  for(q in 3:living)G.lista[[q]] <- G%*%G.lista[[q-1]]
  #
  # Trophic Transformation Matrices: [A] for input set to TL = 0 
  # and [B] for non-living inputs with TL = 1 (called AAA e BBB)
  #
  EL <- as.list(rep(NA,12))
  TROPHIC <- c(1:living)
  ONE <- rep(1,living)
  #
  # calculation of [A]
  #
  AAA <- matrix(rep(0,living^2),nrow=living)
  for(q in 1:living)AAA[q,] <- N%*%G.lista[[q]]
  #
  # calculation of [B]
  #
  BBB <- matrix(rep(0,living^2),nrow=living)
  for(q in 1:living)BBB[q,] <- NK%*%G.lista[[q]]
  #
  # multiplication of [A] and [B] matrices with the corresponding (TROPHIC) and (ONE) vectors of TLs
  #
  AAAPOS <- TROPHIC%*%AAA
  BBBPOS <- (ONE+TROPHIC)%*%BBB
  #
  # saving [A] and [B]
  #
  BBB <- rbind(rep(0,living),BBB)
  BBB <- BBB[-(living+1),]
  EL[[1]] <- AAA + BBB									## Lindeman Transformation Matrix
  #
  # final trophic positions vector for the living compartments (finalTPs)
  #
  intermTPs <- AAAPOS+BBBPOS
  finalTPs <- t(intermTPs)
  EL[[2]] <- c(finalTPs,rep(1,nl))							## Effective Trophic Levels of each species
  #
  {
    if(nl!=0){
      EL[[3]] <- c(EL[[1]] %*% t(t(E_living)),outt[c(living+1):xx])			## Canonical Exports
      R.det <- sum(diss[c(living+1):xx])
      R.det.v <- rep(0,nl)
      R.det.v[nl] <- R.det
      EL[[4]] <- c(EL[[1]] %*% t(t(R_living)),R.det.v)				## Canonical Respirations
      #
      DetPool <- matrix(inter[c(1:living),c((living+1):xx)],nrow=living,byrow=FALSE)
      EL6.s1 <- apply(DetPool,1,sum)
      EL[[6]] <- EL[[1]] %*% t(t(EL6.s1))					## Returns to Detrital Pool
      #
      EL[[8]] <- sum(inp[c((living+1):xx)])						## Input to Detrital Pool
      EL[[9]] <- sum(inter[c((living+1):xx),c((living+1):xx)])			## Circulation within Detrital Pool
    }
    else	{
      DetPool <- 0
      EL[[3]] <- EL[[1]] %*% t(t(E_living))
      EL[[4]] <- EL[[1]] %*% t(t(R_living))
      EL[[6]] <- 0
      EL[[8]] <- 0
      EL[[9]] <- 0
    }
  }
  #
  EL[[7]] <- sum(LVSNL)									## Detritivory
  #
  EL5.s1 <- apply(T_living,1,sum)
  EL5.s2 <- EL[[1]] %*% t(t(EL5.s1))
  lll <- living + 1
  IN.liv <- apply(T_living,2,sum)
  prim <- which(Z_living != 0 & IN.liv == 0)
  Z.prim <- sum(Z_living[prim])
  EL5.s3 <- c(Z.prim,EL5.s2)
  EL5.s3 <- EL5.s3[-lll]
  EL[[5]] <- EL5.s3												## The Grazing Chain
  #
  Det.not.p <- DetPool
  if(nl!=0)Det.not.p[prim,] <- rep(0,nl)
  lind.s1 <- EL[[5]]
  lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + EL[[8]]
  lind.s1[2] <- lind.s1[2] + EL[[7]]
  EL[[10]] <- lind.s1												## Lindeman Spine
  #
  lind.eff <- rep(0,(living-1))
  for(q in 2:living){
    if(EL[[10]][q-1] > 0.0001){
      lind.eff[q-1] <- EL[[10]][q]/EL[[10]][q-1]
    }
  }
  EL[[11]] <- lind.eff											## Trophic Efficiencies
  #
  mig <- which(Z_living != 0 & IN.liv != 0)
  {
    if(length(mig)==0)mig <- 0
    #
    else	{
      TP.migr <- rep(0,living)
      TP.detr <- rep(0,living)
      #
      TP.migr[mig] <- trunc(EL[[2]])[mig]
      ONEs <- which(apply(LVSNL,2,sum)!=0)
      TP.detr[ONEs] <- 1
      EL <- TP.extendedCTA(inp,inter,outt,diss,nl,TP.migr,TP.detr)
      #
      EL[[13]] <- EL[[1]]%*%t(t(inp[1:living]))					## Canonical Distribution of inputs
      EL[[13]] <- c(EL[[13]],inp[(living+1):xx])					## (only in case of migrations)
      names(EL[[13]]) <- NULL
    }
  }
  EL[[12]] <- mig													## Detecting Migratory Imports
  #
colnames(EL[[1]])<-colnames(T_cs)[c(1:living)]#names of living compartments  
names(EL)<-c("Lindeman Transformation Matrix", "Effective Trophic Levels of each species", 
               "Canonical Exports", "Canonical Respirations", "Grazing Chain", 
               "Returns to Detrital Pool", "Detrivory", "Input to Detrital Pool",
               "Circulation within Detrital Pool", "Lindeman Spine", "Trophic Efficiency",
               "Detection of migration", "Canonical Imports")[c(1:length(EL))]
  
  return(EL)
}

## System Omnivory Index (SOI)

SOI <- function(inp,inter,outt,diss,nl){
  GG <- partial.feeding(inp,inter)
  npp <- max(which(apply(GG,2,sum)==0))
  TP <- TP.CTA(inp,inter,outt,diss,nl)[[2]]
  nr <- nc <- nrow(GG)
  live <- nr - nl
  OI_m <- matrix(rep(0,nr^2), nrow = nr)
  for(i in 1:nr){
    for(j in 1:nc){
      if(GG[i,j]!=0)OI_m[i,j] <- ((TP[i] - (TP[j]-1))^2) * GG[i,j]
    }
  }
  OI <- apply(OI_m,2,sum)
  TI <- inp + apply(inter,2,sum)
  s1 <- s2 <- rep(0,(live-npp))
  for(i in (npp+1):live){
    s1[i] <- OI[i] * log(TI[i])
    s2[i] <- log(TI[i])
  }
  sN <- sum(s1)
  sQ <- sum(s2)
  SOI <- sN/sQ
  return(SOI)
}

