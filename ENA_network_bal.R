#####################################################
###						  						  ###
### Balancing procedure with different algorithms ###
###												  ###
#####################################################

## "ETM" function represents an intermediate step for performing all the balancing procedure listed below.

ETM <- function(inp,inter,outt,diss){
  tot_row <- length(inp)
  tot_el <- tot_row^2
  zeros_side <- rep(0, tot_row)
  core <- cbind(inter, outt, diss, zeros_side)
  #
  inputs <- c(inp, 0, 0, 0)
  zeros_down <- rep(0, tot_row+3)
  #
  final <- rbind(core, zeros_down, zeros_down, inputs)
  #
  return(final)
}


## Input Based Approach.

IN.balance <- function(inp,inter,outt,diss){
  ##
  ## data input as: import vector (inp); intercompartmental exchange matrix [inter]; export (outt) 
  ## and respiration vectors (diss); total number of nodes x
  ## STEP 1
  ##
  x <- nrow(inter)
  matrice.step1 <- ETM(inp,inter,outt,diss)
  INP <- apply(matrice.step1,2,sum)
  OUT <- apply(matrice.step1,1,sum)
  #
  DIFF <- rep(0,x)
  N_INP <- rep(0,x)
  N_OUT <- rep(0,x)
  #
  for(i in 1:x){
    N_INP[i] <- INP[i]
    N_OUT[i] <- OUT[i]
    DIFF[i] <- INP[i] - OUT[i]
  }
  #
  ## DIFF_IN <- (DIFF*100)/N_INP
  ## diff.df <- data.frame(cbind(N_INP,N_OUT,DIFF,DIFF_IN))
  ##
  ## STEP 2 - Partial Host Matrix [F]
  ##
  ele <- x + 3
  F <- matrix(rep(0,ele^2),nrow=ele)
  #
  for(i in 1:ele)
    for(j in 1:ele){
      if(OUT[i]!=0) F[i,j] <- matrice.step1[i,j]/OUT[i]
    }
  ##
  ## STEP 3 - Matrix [R] = t[cF] - [I]
  ## Transpose the N x N part of matrix F and substract the identity matrix to get matrix [R]
  ##
  cF <- matrix(rep(0,x^2),nrow=x)
  #
  for(i in 1:x)
    for(j in 1:x){
      cF[i,j] <- F[i,j]
    }
  #
  I <- diag(1,x)
  R <- t(cF) - I
  ##
  ## STEP 4 - Inversion of matrix [R]
  ##
  ddd <- det(R)
  if(ddd == 0)return(NA)
  RI <- solve(R)
  ##
  ## STEP 5 - Matrix [R] manipulation
  ## Multiply every r_ij coefficient by the corresponding jth input in the [ETM] matrix, and change its sign
  ## r_ij <- -r_ij * (t_N+1,j + t_N+2,j + t_N+3,j)
  ##
  R.man <- matrix(rep(0,x^2),nrow=x)
  #
  for(i in 1:x)
    for(j in 1:x){
      R.man[i,j] <- (-RI[i,j]*inp[j])
    }
  #
  ##
  ## STEP 6 - Vector (U)
  ## Sum i's row of the [R.man] matrix to build vector U(u_i)
  ##
  U <- apply(R.man,1,sum)
  ##
  ## STEP 7 - [T.balIN]
  ## Multiply each f_ij coefficient by the corresponding u_i, to obtain a balanced form [T.balIN] of the [inter] matrix
  ##
  T.balIN <- matrix(rep(0,ele^2),nrow=ele)
  #
  for(i in 1:x)
    for(j in 1:ele){
      T.balIN[i,j] <- F[i,j] * U[i]
    }
  
  for(i in (x+1):ele)
    for(j in 1:ele){
      T.balIN[i,j] <- matrice.step1[i,j]
    }
  #
  return(T.balIN)
}


## Output Based Approach.

OUT.balance <- function(inp,inter,outt,diss){
  ##
  ## data input as: import vector (inp); intercompartmental exchange matrix [inter]; export (outt) 
  ## and respiration vectors (diss); total number of nodes x
  ## STEP 1
  ##
  x <- nrow(inter)
  matrice.step1 <- ETM(inp,inter,outt,diss)
  matrice.step1 <- t(matrice.step1)
  #
  INP <- apply(matrice.step1,2,sum)
  OUT <- apply(matrice.step1,1,sum)
  ##
  ## STEP 2 - Partial Host Matrix [F]
  ##
  ele <- x + 3
  F <- matrix(rep(0,ele^2),nrow=ele)
  #
  for(i in 1:ele)
    for(j in 1:ele){
      if(sum(matrice.step1[i,])!=0) F[i,j] <- matrice.step1[i,j]/sum(matrice.step1[i,])
    }
  ##
  ## STEP 3 - Matrix [R] = t[cF] - [I]
  ## Transpose the N x N part of matrix F and substract the identity matrix to get matrix [R]
  ##
  cF <- matrix(rep(0,x^2),nrow=x)
  #
  for(i in 1:x)
    for(j in 1:x){
      cF[i,j] <- F[i,j]
    }
  #
  I <- diag(1,x)
  R <- t(cF) - I
  ##
  ## STEP 4 - Inversion of matrix [R]
  ##
  ddd <- det(R)
  if(ddd == 0)return(NA)
  RI <- solve(R)
  ##
  ## STEP 5 - Matrix [R] manipulation
  ## Multiply every r_ij coefficient by the corresponding jth input in the [ETM] matrix, and change its sign
  ## r_ij <- -r_ij * (t_N+1,j + t_N+2,j + t_N+3,j)
  ##
  R.man <- matrix(rep(0,x^2),nrow=x)
  ER <- outt + diss
  #
  for(i in 1:x)
    for(j in 1:x){
      R.man[i,j] <- (-RI[i,j]*ER[j])
    }
  ##
  ## STEP 6 - Vector (U)
  ## Sum i's row of the [R.man] matrix to build vector U(u_i)
  ##
  U <- apply(R.man,1,sum)
  ##
  ## STEP 7 - [T.balOUT]
  ## Multiply each f_ij coefficient by the corresponding u_i, to obtain a balanced form [T.balIN]
  ## of the [inter] matrix and finally transpose it to calculate [T.balOUT]
  ##
  T.balIN <- matrix(rep(0,ele^2),nrow=ele)
  #
  for(i in 1:x)
    for(j in 1:ele){
      T.balIN[i,j] <- F[i,j] * U[i]
    }
  #				
  for(i in (x+1):ele)
    for(j in 1:ele){
      T.balIN[i,j] <- matrice.step1[i,j]
    }
  #
  T.balOUT <- t(T.balIN)
  #
  return(T.balOUT)
}


## Average Approach (AVG derived algorithm).

AVG.balance <- function(inp,inter,outt,diss){
  ##
  ## it calculates average coefficients using the corresponding values obtained 
  ## with input-based and output-based approach
  ##
  if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
  if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
  #
  matrix.one <- IN.balance(inp,inter,outt,diss)
  matrix.two <- OUT.balance(inp,inter,outt,diss)
  #
  matrix.sum <- matrix.one + matrix.two
  #
  matrix.end <- 0.5 * matrix.sum
  #
  return(matrix.end)
}


## Input - Output Based Approach.

IO.balance <- function(inp,inter,outt,diss){
  ## 
  ## starting from unbalanced matrix [inter], the input based approach is applied 
  ## to (1/2)[inter], the result is summed with (1/2)[inter]. The resulting matrix is
  ## completely balanced using the output-based algorithm
  ##
  x <- nrow(inter)
  if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
  #
  matrix.one <- 0.5 * IN.balance(inp,inter,outt,diss)
  matrix.two <- 0.5 * (ETM(inp,inter,outt,diss))
  #
  matrix.interm <- matrix.one + matrix.two
  #
  scambi<-matrix(rep(0,x^2),nrow=x)
  inpor <- rep(0,x)
  expor <- rep(0,x)
  res <- rep(0,x)
  #
  for(i in 1:x)
    for(j in 1:x)scambi[i,j] <- matrix.interm[i,j]
  #
  for(i in 1:x){ 
    inpor[i] <- matrix.interm[(x+3),i]
    expor[i] <- matrix.interm[i,(x+1)]
    res[i] <- matrix.interm[i,(x+2)]
  }
  #
  if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
  matrix.final <- OUT.balance(inpor,scambi,expor,res)
  #
  return(matrix.final)
}


## Output - Input Based Approach.

OI.balance <- function(inp,inter,outt,diss){
  ## 
  ## is the same of IO but first output-based approach is implemented, and then input-based
  ##
  x <- nrow(inter)
  if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
  #
  matrix.one <- 0.5 * OUT.balance(inp,inter,outt,diss)
  matrix.two <- 0.5 * (ETM(inp,inter,outt,diss))
  #
  matrix.interm <- matrix.one + matrix.two
  #
  scambi<-matrix(rep(0,x^2),nrow=x)
  inpor <- rep(0,x)
  expor <- rep(0,x)
  res <- rep(0,x)
  #
  for(i in 1:x)
    for(j in 1:x)scambi[i,j] <- matrix.interm[i,j]
  #
  for(i in 1:x){ 
    inpor[i] <- matrix.interm[(x+3),i]
    expor[i] <- matrix.interm[i,(x+1)]
    res[i] <- matrix.interm[i,(x+2)]
  }
  #
  if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
  matrix.final <- IN.balance(inpor,scambi,expor,res)
  #
  return(matrix.final)
}


## Average-2 Approach (AVG2 derived algorithm).
## This is the same of AVG, but using "Input - Output" and "Output - Input" based approaches.

AVG2.balance <- function(inp,inter,outt,diss){
  ##
  ## the same of AVG but with [T]_IO*(bal) and [T]_OI*(bal) matrices
  ##
  if(length(IO.balance(inp,inter,outt,diss)) == 1)return(NA)
  if(length(OI.balance(inp,inter,outt,diss)) == 1)return(NA)
  #
  matrix.one <- IO.balance(inp,inter,outt,diss)
  matrix.two <- OI.balance(inp,inter,outt,diss)
  #
  matrix.sum <- matrix.one + matrix.two
  #
  matrix.end <- 0.5 * matrix.sum
  #
  return(matrix.end)
}


