######################################################################################################
###												   												   ###
### Functions to convert dataframe input data (in a similar SCOR format) into vectors and matrices ###
###												   												   ###
######################################################################################################


## Conversion of SCOR vectors (available as input dataframes) into whole vectors. It is useful for (B), (Z), (E) and (R).

vec.ena <- function(df,x){
		V <- rep(0,x)
		for(i in 1:nrow(df))V[df[i,1]] <- df[i,2]
		return(V)
}


## Conversion of SCOR matrices (available as input dataframes) into whole matrix. 
## Mainly useful for the matrix of inter-compartmental[T] exchanges.

matr.ena <- function(df,x){
			tot_el <- x^2 
			M <- matrix(rep(0,tot_el), nrow=x)
			for(i in 1:nrow(df))M[df[i,1],df[i,2]] <- df[i,3]
			return(M)
}

########################################################################################################################################
###																     ###
### Functions to depict all the flows (import, export, dissipations and inter-compartmentmental exchanges) as one big matrix [T_tot] ###
### In this case the input considered are vector and matrices (not dataframes)							     ###
### In the first function ("bigmatr.ena") the output will also show the name of each compartment, while this addional information    ###
### is not available in the case of the raw result of "bigmatrNOnames.ena"							     ###
###																     ###
########################################################################################################################################

bigmatrNOnames.ena <- function(inp,inter,outt,diss){
					tot_row <- length(inp)
					tot_el <- tot_row^2
					zeros_up <- rep(0, tot_row)
					core <- cbind(zeros_up, inter, outt, diss)
					inputs <- c(0, inp, 0, 0)
					zeros_down <- rep(0, tot_row+3)
					final <- rbind(inputs, core, zeros_down, zeros_down)
					colnames(final) <- NULL
					rownames(final) <- NULL
					#
					return(final)
}





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

############################################################################################
###																						 ###
### System level index for the assessment of the whole ecosystem functioning		 	 ###
### We also included the analysis of some sinle-species features that are extracted from ###
### global indices (e.g. AMI on input and output)					 					 ###
###											 											 ###
############################################################################################

## The function "TST" is used to calculate the Total System Throughput. The inputs are vectors of import (Z), 
## export (E) and dissipations (R) with the matrix of intercompartmental exchanges [T].

TST <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				B <- sum(A)
				return(B)
}


## The function "AMI" calculates the Average Mutual Information of the ecosystem. Also in this case, as in the following,
## the only needed inputs are given by a matrix of internal exchanges [T] with vectors of import (Z), export (E) and respirations (R).

AMI <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(A)
				ami <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)ami <- (ami+0)
							else if (sum(A[,j])==0)ami <- (ami+0)
								else if (sum(A[i,])==0)ami <- (ami+0)
									else ami <- ami+(A[i,j]/tot)*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				return(ami)
}


## The function "DC" calculates the Development Capacity (the upper limit to the development of an ecosystem).

DC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(A)
				dc <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)dc <- (dc+0)
							else dc <- dc+A[i,j]*log2(A[i,j]/tot)
							}
						}
				dev.cap<-(-1*dc)
				return(dev.cap)
}


## The function "AMI.max" estimates the AMI corresponding to the Development Capacity.

AMI.max <- function(inp,inter,outt,diss){
					dev <- DC(inp,inter,outt,diss)
					tot <- TST(inp,inter,outt,diss)
					a.m <- dev/tot
					return(a.m)
}


## The function "ASC" to compute the Ascendency of the system (TST * AMI).

ASC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				xx <- ncol(A)
				asc <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)asc <- (asc+0)
							else if (sum(A[,j])==0)asc <- (asc+0)
								else if (sum(A[i,])==0)asc <- (asc+0)
									else asc <- asc+A[i,j]*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				perc <- asc/develop
				end <- matrix(c(asc,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Ascendency")
				colnames(end) <- c("Value", "Ratio")				
				return(end)
}


## The function "OVERHEAD.imports" determines the Overhead on Imports.

OVERHEAD.imports <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				io <- 0
				for (j in 1:kk){
					if(A[1,j]==0)io <- (io+0)
						else io <- io+A[1,j]*log2((A[1,j]^2)/(sum(A[1,])*sum(A[,j])))
						}
				#
				{
				if(sum(io!=0))import.o <- (-1*io)
				else import.o <- 0
				}
				#
				perc <- import.o/develop
				end <- matrix(c(import.o,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Overhead on Imports")
				colnames(end) <- c("Value", "Ratio")
				return(end)
}


## The function "OVERHEAD.exports" determines the Overhead on Exports.

OVERHEAD.exports <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				uu <- ncol(A)-1
				eo <- 0
				for (i in 1:kk){
					if(A[i,uu]==0)eo <- (eo+0)
						else eo <- eo+A[i,uu]*log2((A[i,uu]^2)/(sum(A[i,])*sum(A[,uu])))
						}
				#
				{
				if(sum(eo!=0))export.o <- (-1*eo)
				else export.o <- 0
				}
				#
				perc <- export.o/develop
				end <- matrix(c(export.o,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Overhead on Exports")
				colnames(end) <- c("Value", "Ratio")
				return(end)
}


## The function "OVERHEAD.resp" determines the Dissipative Overhead.

OVERHEAD.resp <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				xx <- ncol(A)
				kk <- ncol(A)-2
				do <- 0
				for (i in 1:kk){
					if(A[i,xx]==0)do <- (do+0)
						else do <- do+A[i,xx]*log2((A[i,xx]^2)/(sum(A[i,])*sum(A[,xx])))
						}
				#
				{
				if(sum(do!=0))dissipative.o <- (-1*do)
				else dissipative.o <- 0
				}
				#
				perc <- dissipative.o/develop
				end <- matrix(c(dissipative.o,perc),nrow=1,byrow=TRUE)
				rownames (end) <- c("Dissipative Overhead")
				colnames(end) <- c("Value", "Ratio")
				return(end)
}


## The function "REDUNDANCY" determines the Redundancy (Overhead on internal flows).

REDUNDANCY <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				r <- 0
				for (i in 2:kk){
					for (j in 2:kk){
						if(A[i,j]==0)r <- (r+0)
							else r <- r+A[i,j]*log2((A[i,j]^2)/(sum(A[i,])*sum(A[,j])))
							}
						}
				#
				{
				if(sum(r)!=0)redun <- (-1*r)
				else redun <- 0
				}
				#
				perc <- redun/develop
				end <- matrix(c(redun,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Redundancy")
				colnames(end) <- c("Value", "Ratio")
				return(end)
}


## The function "OVERHEAD" determines the total amount of System Overhead.

OVERHEAD <- function(inp,inter,outt,diss){
					o1 <- OVERHEAD.imports(inp,inter,outt,diss)[1,1]
					o2 <- OVERHEAD.exports(inp,inter,outt,diss)[1,1]
					o3 <- OVERHEAD.resp(inp,inter,outt,diss)[1,1]
					o4 <- REDUNDANCY(inp,inter,outt,diss)[1,1]
					develop <- DC(inp,inter,outt,diss)
					o <- o1 + o2 + o3 + o4
					perc <- o/develop
					end <- matrix(c(o,perc),nrow=1,byrow=TRUE)
					rownames(end) <- c("Total Overhead")
					colnames(end) <- c("Value", "Ratio")
					return(end)
}


## The function "internal.DC" is used to calculate the Internal Development Capacity.

internal.DC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				kk <- ncol(A)-2
				ic <- 0
				for (i in 2:kk){
					for (j in 2:kk){
						if(A[i,j]==0)ic <- (ic+0)
							else ic <- ic+A[i,j]*log2(A[i,j]/tot)
							}
						}
				#
				{
				if(sum(ic)!=0)intc <- (-1*ic)
				else intc <- 0
				}
				#
				return(intc)
}


## The function "internal.ASC" allows the computation of the Internal Ascendency.

internal.ASC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				red <- REDUNDANCY(inp,inter,outt,diss)[1,1]
				idc <- internal.DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				iasc <- 0
				for (i in 2:kk){
					for(j in 2:kk){
						if(A[i,j]!=0)iasc <- iasc+A[i,j]*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				perc <- iasc/idc
				perc1 <- red/idc
				end <- matrix(c(iasc,perc,red,perc1),nrow=2,byrow=TRUE)
				rownames(end) <- c("Internal Ascendency", "Redundancy")
				colnames(end) <- c("Value", "Ratio")
				return(end)
}


## The function for computing the overall connectance also includes the effects of exogenous transfers.

O.C <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				el <- length(inp)
				#
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				#
				C.M <- matrix(rep(0,el^2),nrow=el)
				V.I <- rep(0,el)
				V.E <- rep(0,el)
				V.R <- rep(0,el)
				#
				for(i in 1:el)
					for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]/tot) * log(inter[i,j]^2/(in.v[j]*out.v[i]))
				#
				for(i in 1:el){
					if(inp[i]!=0)V.I[i] <- (inp[i]/tot) * log(inp[i]/in.v[i])
					#
					if(outt[i]!=0)V.E[i] <- (outt[i]/tot) * log(outt[i]/out.v[i])
					#
					if(diss[i]!=0)V.R[i] <- (diss[i]/tot) * log(diss[i]/out.v[i])
				}
				#
				TT <- sum(-1 * sum(C.M))
				II <- sum(-2 * sum(V.I))
				EE <- sum(-2 * sum(V.E))
				RR <- sum(-2 * sum(V.R))
				#
				fi <- sum(TT + II + EE + RR)
				#
				m <- exp(fi/2)
				#
				return(m)
}


## The intercompartmental connectance characterizes only the endogenous exchanges.

I.C <- function(inter){
		tot <- sum(inter)
		el <- nrow(inter)
		#
		in.v <- apply(inter,2,sum)
		out.v <- apply(inter,1,sum)
		#
		C.M <- matrix(rep(0,el^2),nrow=el)
		#
		for(i in 1:el)
			for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]/tot) * log(inter[i,j]^2/(in.v[j]*out.v[i]))
		#
		fi <- sum(-1 * sum(C.M))
		#
		m <- exp(fi/2)
		#
		return(m)
}


## The foodweb connectance pertains only to transfers among the living compartments.

FW.C <- function(inter,nl){
		el <- nrow(inter) - nl
		inter.l <- inter[c(1:el),c(1:el)]
		tot <- sum(inter.l)
		#
		in.v <- apply(inter.l,2,sum)
		out.v <- apply(inter.l,1,sum)
		#
		C.M <- matrix(rep(0,el^2),nrow=el)
		#
		for(i in 1:el)
			for(j in 1:el)if(inter.l[i,j]!=0)C.M[i,j] <- (inter.l[i,j]/tot) * log(inter.l[i,j]^2/(in.v[j]*out.v[i]))
		#
		fi <- sum(-1 * sum(C.M))
		#
		m <- exp(fi/2)
		#
		return(m)
}


## Weighted connectance computed with the formula described by Zorach and Ulanowicz (2003).
## It is equivalent to the effective connectance (m) introduced by Ulanowicz and Wolff (1991).

effconn.m <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				el <- nrow(inter)
				#
				C.M <- matrix(rep(1,el^2),nrow=el)
				#
				for(i in 1:el)
					for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]^2/(in.v[j]*out.v[i]))^(-0.5*inter[i,j]/tot)
				#
				m <- prod(C.M)
				#
				return(m)
}


## Link density (LD) estimated with the Bersier et al. (2002) formulation.
## It is slightly different from the Ulanowicz's effective connectance (m).
## With LD, the weighting applies directly to the taxa's equivalent numbers of prey and consumers.

LD.B <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				el <- nrow(inter)
				#
				H.N.M <- matrix(rep(0,el^2),nrow=el)
				H.P.M <- matrix(rep(0,el^2),nrow=el)
				N.v <- rep(0,el)
				P.v <- rep(0,el)
				#
				for(i in 1:el)
					for(j in 1:el){
						if(inter[i,j]!=0){
							H.N.M[i,j] <- (-inter[i,j]/in.v[j]) * log2(inter[i,j]/in.v[j])
							#
							H.P.M[i,j] <- (-inter[i,j]/out.v[i]) * log2(inter[i,j]/out.v[i])
							}
						}
				#
				H.N <- apply(H.N.M,2,sum)
				H.P <- apply(H.P.M,1,sum)
				#
				for(i in 1:el){
					P.v[i] <- (out.v[i]/tot) * 2^H.P[i]
					#
					N.v[i] <- (in.v[i]/tot) * 2^H.N[i]
					}
				#
				LD <- 0.5 * (sum(P.v) + sum(N.v))
				#
				return(LD)
}


## Biomass inclusive Ascendency (Ulanowicz and Abarca-Arenas, 1997).
## It is with intercompartmental exchanges only (equivalent to the Internal Ascendency).

ASC.B <- function(inp,inter,outt,diss,bi){
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(inter)
				B <- sum(bi)
				#
				asc <- 0
				for (i in 1:xx){
					for (j in 1:xx)if(inter[i,j]!=0)asc <- asc + (inter[i,j] * log2((inter[i,j] * B^2)/(tot * bi[i] * bi[j])))
						}
				#
				return(asc)
}


## AMI on input for each compartment. This function estimates two versions of the index:
## 1) the first takes into account all the flows (intercompartmental exchanges both for living and non-living compartments
## as well as vectors of import, export and respiration); 2) the second includes living exchanges only. 
## These indices are described in Scotti et al. (2009).

AMI.in <- function(inp,inter,outt,diss,nl){
				tot <- TST(inp,inter,outt,diss)
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				A.in <- apply(A,2,sum)
				A.ou <- apply(A,1,sum)
				ww <- nrow(A)
				ll <- nrow(inter)
				#
				A.w <- matrix(rep(0,ww^2),nrow=ww)
				#
				for(i in 1:ww)
					for(j in 1:ww)if(A[i,j]!=0)A.w[i,j] <- A[i,j] * log2((A[i,j] * tot)/(A.in[j] * A.ou[i]))
				#
				A.w.c <- apply(A.w,2,sum)
				AMI.w.c <- A.w.c[2:(ll+1)]/A.in[2:(ll+1)]
				#
				#
				live <- ll - nl
				L <- inter[c(1:live),c(1:live)]
				tot.l <- sum(L)
				A.l.in <- apply(L,2,sum)
				A.l.ou <- apply(L,1,sum)
				#
				A.l <- matrix(rep(0,live^2),nrow=live)
				#
				for(i in 1:live)
					for(j in 1:live)if(L[i,j]!=0)A.l[i,j] <- L[i,j] * log2((L[i,j] * tot.l)/(A.l.in[j] * A.l.ou[i]))
				#
				A.l.c <- apply(A.l,2,sum)
				#
				AMI.l.c <- rep(0,live)
				#
				for(i in 1:live){
					if(A.l.c[i]!=0)AMI.l.c[i] <- A.l.c[i]/A.l.in[i]
					}
				#
				# this function returns two outcomes for the AMI on inflows:
				# 1) the first is estimated using the "classical" ENA approach (including
				# non-living nodes and exchanges with surroundings (i.e., import,
				# export and respiration vectors);
				# 2) the second is computed by taking into account the flows between
				# living compartments only
				#
				AMI.IN.FINAL <- as.list(rep(NA,2))
				names(AMI.IN.FINAL) <- c("ENA","living")
				AMI.IN.FINAL[[1]] <- AMI.w.c
				AMI.IN.FINAL[[2]] <- AMI.l.c
				#
				return(AMI.IN.FINAL)
}


## AMI on export for each compartment. This index is measured:
## 1) for the whole flows, with the "classical" network analysis approach; 2) for living exchanges only.
## Further details are provided with the description of the function "AMI.in".

AMI.out <- function(inp,inter,outt,diss,nl){
				tot <- TST(inp,inter,outt,diss)
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				A.in <- apply(A,2,sum)
				A.ou <- apply(A,1,sum)
				ww <- nrow(A)
				ll <- nrow(inter)
				#
				A.w <- matrix(rep(0,ww^2),nrow=ww)
				#
				for(i in 1:ww)
					for(j in 1:ww)if(A[i,j]!=0)A.w[i,j] <- A[i,j] * log2((A[i,j] * tot)/(A.in[j] * A.ou[i]))
				#
				A.w.r <- apply(A.w,1,sum)
				AMI.w.r <- A.w.r[2:(ll+1)]/A.ou[2:(ll+1)]
				#
				#
				live <- ll - nl
				L <- inter[c(1:live),c(1:live)]
				tot.l <- sum(L)
				A.l.in <- apply(L,2,sum)
				A.l.ou <- apply(L,1,sum)
				#
				A.l <- matrix(rep(0,live^2),nrow=live)
				#
				for(i in 1:live)
					for(j in 1:live)if(L[i,j]!=0)A.l[i,j] <- L[i,j] * log2((L[i,j] * tot.l)/(A.l.in[j] * A.l.ou[i]))
				#
				A.l.r <- apply(A.l,1,sum)
				#
				AMI.l.r <- rep(0,live)
				#
				for(i in 1:live){
					if(A.l.r[i]!=0) AMI.l.r[i] <- A.l.r[i]/A.l.ou[i]
					}
				#
				# this function returns two outcomes for the AMI on outflows:
				# 1) AMI on outflows using the "classical" ENA approach (including
				# non-living nodes and exchanges with surroundings (i.e., import,
				# export and respiration vectors);
				# 2) AMI on outflows computed by taking into account the flows between
				# living compartments only
				#
				AMI.OU.FINAL <- as.list(rep(NA,2))
				names(AMI.OU.FINAL) <- c("ENA","living")
				AMI.OU.FINAL[[1]] <- AMI.w.r
				AMI.OU.FINAL[[2]] <- AMI.l.r
				#
				return(AMI.OU.FINAL)
}


##########################################################
###		     				       ###
### Cycling Analysis - The structure of network cycles ###
###		     				       ###
##########################################################

## These functions ("coppia.bis") and ("LINK.to.NODES") are used into the algorithm for decyclization.
## They do not return any ENA outcome, although implemented within the cycling analysis framework.

coppia.bis <- function(V,M){
		rr <- nrow(M)
		A <- ((V[2]-1)*rr) + V[1]
		for(i in 3:length(V)) A <- c(A,(((V[i]-1)*rr)+V[i-1]))
		#
		return(A)
		}


LINK.to.NODES <- function(VET){
	nnn <- length(VET)
	NOD <- c(VET,NA) %% nnn
	NOD[which(NOD==0)] <- nnn
	mmm <- min(which(is.na(NOD)==TRUE))
	NOD[mmm] <- NOD[1]
	#
	return(NOD)
	}


## Ulanowicz (1983, 1986) devised a procedure to extract cycles from the network and compare their activities with the remaining unidirectional flows.
## He adopted a backtracking algorithm to find all the simple cycles (Mateti and Deo, 1976), increasing its efficiency with a suitable pruning method
## (Knuth, 1973). The whole process can be summarized as follows: (a) a preliminary depth-first search to count the number of cycle arcs incident
## to each node; (b) ordering nodes by decreasing number of incoming cycle arcs (compartments with no cycle arcs are eliminated from further analysis);
## (c) scrutinizing each cycle in search for the smallest cycle arc; (d) identifying if this cycle link is common to other cycles and grouping cycles
## sharing it (nexus); (e) starting backtracking algorithm to delete the nexus and distribute its magnitude among member cycles (in proportion to circuit
## probabilities and magnitudes of links). Finally, the acyclic residual oncethrough flow web (tree, in graph theory) and the aggregated network of
## cycled medium result as two separated parts of the whole network.

cyc.analysis <- function(inp,inter,outt,diss){
	#
	nt.system <- nrow(inter)
	CRCC.FULL <- inter
	CRGG <- partial.feeding(inp,inter)
	CRFF <- partial.host(inter,outt,diss)
	CRCC.ONES <- matrix(rep(0,nt.system^2),nrow=nt.system)
	#
	for(i in 1:nt.system)
		for(j in 1:nt.system)if(inter[i,j]!=0)CRCC.ONES[i,j] <- 1
	#
	ROW.CRCC <- apply(CRCC.ONES,1,sum)
	COL.CRCC <- apply(CRCC.ONES,2,sum)
	#
	SUM.CRCC <- ROW.CRCC + COL.CRCC
	names(SUM.CRCC) <- c(1:nt.system)
	#
	FROM.ROW <- as.list(rep(NA,nt.system))
	for(i in 1:nt.system)FROM.ROW[[i]] <- which(CRCC.ONES[i,]!=0)
	#
	#
	## A) Defining nodes showing both inflows and outflows. The others are excluded
	## (nodes with only input or output are excluded) --> "BLACK" nodes
	## B) Ordering nodes for their total number of links (IN + OUT)
	#
	color <- rep("WHITE",nt.system)
	#
	for(i in 1:nt.system)if(ROW.CRCC[i]==0 | COL.CRCC[i]==0){
							color[i] <- "BLACK"
							SUM.CRCC[i] <- 0
							}
	#
	SRT <- sort(SUM.CRCC,decreasing=TRUE)
	EFF.SRT <- as.integer(names(which(SRT>0)))
	tr <- length(EFF.SRT)
	p <- c(1:nt.system)
	MASSI <- rep(0,nt.system)
	#
	for(i in 1:nt.system){
			if(length(FROM.ROW[[i]])!=0){
				gt <- length(FROM.ROW[[i]])
				#
				W <- which(EFF.SRT==FROM.ROW[[i]][1])
				#
				if(gt > 1){
					for(j in 2:gt)W <- c(W,which(EFF.SRT==FROM.ROW[[i]][j]))
					}
				if(length(W) > 0)MASSI[i] <- max(W)
				}
			}
	#
	#
	## Looking for self loops and deleting them, when found
	#
	self.cyc <- 0
	inds <- 1
	#
	for(i in 1:nt.system)if(CRCC.ONES[i,i]==1){
				self.cyc <- self.cyc + 1
				CRCC.ONES[i,i] <- 0
				inter[i,i] <- 0
				{
				if(inds==1){
					SC <- c(i,i,1,inter[i,i])
					inds <- 2
					}
				else SC <- rbind(SC,c(i,i,1,inter[i,i]))
				}
				#
				}
	#
	#
	## Selecting the more connected node as the starting point for the
	## DFS algorithm.
	#
	cycles.count <- rep(0,nt.system)
	cycles <- 0
	index <- 1
	smile <- 1
	STOPPA <- 0
	one <- 1
	WEAK <- 0
	#
	#
	## Backtracking algorithm including procedure to save and collect weakest arcs
	#
	{
	if(inds!=1)LOO <- 1
	else LOO <- 0
	}
	#
	while(length(color=="WHITE")>1 & sum(CRCC.ONES)!=LOO){
	#
	starting <- EFF.SRT[one]
	color[starting] <- "GRAY"									## Discover vertex starting <- u
	#
	u <- starting
	CIC <- u
	#
	entra <- 0
	#
	#
	counter <- 1
	STOPPA <- 0
	#
	#
	while(STOPPA != 1){
		for(v in counter:tr){
			#
			if(STOPPA == 1)break
			#
			if(entra==2){
				entra <- 0
				aaa <- length(CIC)
				bbb <- length(CIC) - 1
				nott <- CIC[aaa]
				meas <- CIC[bbb]
				#
				if(bbb==0){
					STOPPA <- 1
					break
					}
				#
				{
				if(which(EFF.SRT==nott)<MASSI[meas]){
						 counter <- which(EFF.SRT==nott) + 1
						 CIC <- CIC[-aaa]
						 }
				#
				else {
					CIC <- CIC[-aaa]
					aaa <- length(CIC)
					bbb <- length(CIC) - 1
					#
					if(bbb==0){
						STOPPA <- 1
						break
						}
					#
					nott <- CIC[aaa]
					meas <- CIC[bbb]
					#
					while(which(EFF.SRT==nott)==MASSI[meas]){
									CIC <- CIC[-aaa]
									aaa <- length(CIC)
									bbb <- length(CIC) - 1
									#
									if(bbb==0){
										STOPPA <- 1
										break
										}
									#
									nott <- CIC[aaa]
									meas <- CIC[bbb]
									}
					#
					counter <- which(EFF.SRT==nott) + 1
					CIC <- CIC[-aaa]
					}
				}
				#
				u <- CIC[length(CIC)]
				break
			}
			#
			if(CRCC.ONES[u,EFF.SRT[v]]!=0){							## Examine edge (u,EFF.SRT[v])
				{
				if(color[EFF.SRT[v]]=="WHITE" & length(CIC[CIC==EFF.SRT[v]])==0){	## (u,EFF.SRT[v]) is a tree edge
					#
					CIC <- c(CIC,EFF.SRT[v])
					u <- EFF.SRT[v]
					counter <- 1
					entra <- 0
					break
					}
				else {
					if(color[EFF.SRT[v]]=="GRAY"){					## (u,EFF.SRT[v]) is a back edge
						#
						entra <- 0
						CIC <- c(CIC,EFF.SRT[v])
						names(CIC) <- NULL
						cycles <- cycles + 1
						#
						{
						if(smile==1){				
							LI <- coppia.bis(CIC,inter)
							LL <- length(LI)
							#
							{
							if(index == 1)G.prod <- prod(CRGG[LI])
							#
							else G.prod[index] <- prod(CRGG[LI])
							}
							#
							cycles.count[LL+1] <- cycles.count[LL+1] + 1
							{
							if(LL < nt.system)WEAK <- c(LI,rep(0,nt.system - LL))
							#
							else WEAK <- LI
							}
							index <- index + 1
							}
						else{
							LI <- coppia.bis(CIC,inter)
							LL <- length(LI)
							#
							{
							if(index == 1)G.prod <- prod(CRGG[LI])
							#
							else G.prod[index] <- prod(CRGG[LI])
							}
							#
							cycles.count[LL+1] <- cycles.count[LL+1] + 1
							#
							{
							if(LL < nt.system){
								WEAK <- cbind(WEAK,c(LI,rep(0,nt.system - LL)))
								index <- index + 1
								}
							else	{
								WEAK <- cbind(WEAK,LI)
								index <- index + 1
								}
						     	}
						     }
						}
						smile <- smile + 1
						#
						aaa <- length(CIC)
						bbb <- length(CIC) - 1
						nott <- CIC[aaa]
						meas <- CIC[bbb]
						#
						{
						if(which(EFF.SRT==nott)<MASSI[meas]){
								counter <- which(EFF.SRT==nott) + 1
								CIC <- CIC[-aaa]
								}
						#
						else {	CIC <- CIC[-aaa]
							aaa <- length(CIC)
							bbb <- length(CIC) - 1
							nott <- CIC[aaa]
							meas <- CIC[bbb]
							#					
							while(which(EFF.SRT==nott)==MASSI[meas]){
											CIC <- CIC[-aaa]
											aaa <- length(CIC)
											bbb <- length(CIC) - 1
											#
											if(bbb==0){
												STOPPA <- 1
												break
												}
											#
											nott <- CIC[aaa]
											meas <- CIC[bbb]
											}
							#
							counter <- which(EFF.SRT==nott) + 1
							CIC <- CIC[-aaa]
							#						
							}
						}
						#
						u <- CIC[length(CIC)]
						break
						#
						}
					}
				}
			}
		}	
	#
	entra <- entra + 1
	#
	}
	#
	color[starting] <- "BLACK"
	one <- one + 1
	#
	for(h in 1:nt.system){
		CRCC.ONES[h,starting] <- 0
		CRCC.ONES[starting,h] <- 0
		}
	#
	if(one==length(EFF.SRT))break
	}
	#
	##
	##
	##
	#
	br <- 0
	if(length(which(diag(inter)!=0)) > 0){
				l.sc <- length(which(diag(inter)!=0))
				#
				{
				if(sum(WEAK)==0)
					#
					{
					if(l.sc == 1)
						{
						WEAK <- c(which(diag(inter)!=0)[1],rep(0,nt.system - 1))
						break
						}
					#
					else	
						{
						WEAK <- c(which(diag(inter)!=0)[1],rep(0,nt.system - 1))
						for(i in 2:l.sc)WEAK <- cbind(WEAK,c(which(diag(inter)!=0)[i],rep(0, nt.system - 1)))
						break
						}
					}
				#
				else {
					for(i in 1:l.sc)WEAK <- cbind(WEAK,c(which(diag(inter)!=0)[i],rep(0, nt.system - 1)))
					}
				}
	}
	#
	if(is.matrix(WEAK)==FALSE & sum(WEAK)!=0){
		WEAK <- matrix(WEAK,ncol=1,byrow=TRUE)
		}
	#
	rownames(WEAK) <- NULL
	FULL.CYCLES <- ncol(WEAK)
	#
	##
	##
	##
	#
	CRCC.NA <- matrix(rep(0,nrow(inter)^2),nrow=nrow(inter))
	for(i in 1:nrow(inter))
		for(j in 1:nrow(inter)){
			if(inter[i,j]>0)CRCC.NA[i,j] <- inter[i,j]
			}
	#
	CLASS.cyc <- as.list(NA)
	WEAK.copy <- WEAK
	indaco <- 1
	#
	#
	while(sum(WEAK)!=0){
			#
			CRCC.vec <- as.vector(CRCC.NA)
			CRCC.min <- CRCC.vec
			WEAK.vec <- as.vector(WEAK)
			#
			least <- which.min(CRCC.vec)
			nexus.elements <- which(WEAK.vec==least)
			#
			while(length(nexus.elements)==0){
				CRCC.min[least] <- NA
				least <- which.min(CRCC.min)
				nexus.elements <- which(WEAK.vec==least)
				}
			#
			nexus.vec <- trunc((nexus.elements-1)/nt.system) + 1
			lun <- length(nexus.vec)
			#
			CLASS.cyc[[indaco]] <- matrix(WEAK[,c(nexus.vec)], nrow=lun, byrow=TRUE)
			#
			#
			PRD <- G.prod[c(nexus.vec)]/sum(G.prod[c(nexus.vec)])
			WEAK[,nexus.vec] <- rep(0,nt.system)
			#
			for(i in 1:lun){
				abc <- CLASS.cyc[[indaco]][i,]
				abc.p <- abc[abc > 0]
				CRCC.NA[abc.p] <- CRCC.NA[abc.p] - (CRCC.vec[least] * PRD[i])
				}
			#
				{
				if(indaco==1) nxs.val <- CRCC.NA[least]
				#
				else nxs.val <- c(nxs.val,CRCC.NA[least])
				}
			#
			CRCC.NA[least] <- 0
			indaco <- indaco + 1
		#
	}
	#
	#
	#
	inter <- CRCC.FULL
	#
	if(self.cyc!=0){
		cycles.count[2] <- self.cyc
		FULL.CYCLES <- sum(cycles.count)
		ZER <- matrix(rep(0,nt.system^2),nrow=nt.system)
		diag(ZER) <- diag(inter)
		dia <- which(ZER!=0)
		dia.l <- length(dia)
		WEAK.sl <- matrix(rep(0,(dia.l * nt.system)), nrow = nt.system)
		#
		for(t in 1:dia.l){
			WEAK.sl[,t] <- c(dia[t],rep(0,(nt.system) - 1))
		}
		#
		{
		if(sum(WEAK.copy)!=0)WEAK.copy <- cbind(WEAK.copy,WEAK.sl)
		#
		else WEAK.copy <- WEAK.sl
		}
	}	
	#
	#
	#
	{
	if(sum(FULL.CYCLES)!=0){
		deci <- 1
		#
		WEAK.copy.NAs <- WEAK.copy
		#
		CRCC.LINE <- inter
		WEAK.copy.NAs[which(WEAK.copy.NAs==0)] <- NA
		jump <- 1
		CYCLE.distr <- rep(0,nt.system)
		#
		while(sum(na.omit(WEAK.copy.NAs))!=0){
			#
			if(jump == 1){
				CRCC.NAs <- CRCC.LINE	
				CRCC.NAs[which(CRCC.NAs <= 0)] <- NA
				}
			#
			un.elem <- unique(c(WEAK.copy.NAs))
			un.elem.n <- un.elem[which(is.na(un.elem)==FALSE)]
			#
			nex.min <- which.min(CRCC.NAs)
			nex.an <- which(WEAK.copy.NAs==nex.min)
			#
			{
			if(length(nex.an)!=0){
					rig <- trunc(nex.an/nt.system,1) + 1
					{
					if(deci == 1){
						WHOLE.cyc <- as.list(c(NA))
						WHOLE.FFF <- as.list(c(NA))
						WHOLE.FFF.v <- as.list(c(NA))
						WHOLE.nodes <- as.list(c(NA))
						#
						WHOLE.cyc[[deci]] <- t((WEAK.copy.NAs[,rig]))
						WHOLE.FFF[[deci]] <- matrix(CRFF[WHOLE.cyc[[deci]]], byrow = FALSE, nrow = length(nex.an))
						WHOLE.FFF[[deci]][which(is.na(WHOLE.FFF[[deci]])==TRUE)] <- 1
						WHOLE.FFF.v <- apply(WHOLE.FFF[[deci]],1,prod)
						WHOLE.FFF.s <- sum(WHOLE.FFF.v)
						WHOLE.FFF.r <- WHOLE.FFF.v/WHOLE.FFF.s
						#
						WHOLE.nodes[[deci]] <- matrix(rep(NA,((nt.system + 1)*length(nex.an))),nrow = length(nex.an))
						#
						for(q in 1:length(nex.an)){
							ve <- WHOLE.cyc[[deci]][q,which(is.na(WHOLE.cyc[[deci]][q,])==FALSE)]
							nexi <- CRCC.NAs[nex.min]
							rela <- WHOLE.FFF.r[q]
							dele <- nexi * rela
							CRCC.LINE[ve] <- CRCC.LINE[ve] - (dele)
							#
							WHOLE.nodes[[deci]][q,] <- LINK.to.NODES(WHOLE.cyc[[deci]][q,])
							#
							ll.cc <- length(which(is.na(WHOLE.cyc[[deci]][q,])==FALSE))
							CYCLE.distr[ll.cc] <- CYCLE.distr[ll.cc] + (dele * ll.cc)
							#
						}
						#
						WEAK.ARCS <- inter[nex.min]
						#
						aj <- nex.min%%nt.system
						if(aj == 0) aj <- nt.system
						WHOLE.nexus <- c(aj,trunc(nex.min/nt.system)+1)
						#
						deci <- deci + 1
					}
					#
					else	{
						WHOLE.cyc[[deci]] <- t((WEAK.copy.NAs[,rig]))
						WHOLE.FFF[[deci]] <- matrix(CRFF[WHOLE.cyc[[deci]]], byrow = FALSE, nrow = length(nex.an))
						WHOLE.FFF[[deci]][which(is.na(WHOLE.FFF[[deci]])==TRUE)] <- 1
						WHOLE.FFF.v <- apply(WHOLE.FFF[[deci]],1,prod)
						WHOLE.FFF.s <- sum(WHOLE.FFF.v)
						WHOLE.FFF.r <- WHOLE.FFF.v/WHOLE.FFF.s
						#
						WHOLE.nodes[[deci]] <- matrix(rep(NA,((nt.system + 1)*length(nex.an))),nrow = length(nex.an))
						#
						for(q in 1:length(nex.an)){
							ve <- WHOLE.cyc[[deci]][q,which(is.na(WHOLE.cyc[[deci]][q,])==FALSE)]
							nexi <- CRCC.NAs[nex.min]
							rela <- WHOLE.FFF.r[q]
							dele <- nexi * rela
							CRCC.LINE[ve] <- CRCC.LINE[ve] - (dele)
							#
							WHOLE.nodes[[deci]][q,] <- LINK.to.NODES(WHOLE.cyc[[deci]][q,])
							#
							ll.cc <- length(which(is.na(WHOLE.cyc[[deci]][q,])==FALSE))
							CYCLE.distr[ll.cc] <- CYCLE.distr[ll.cc] + (dele * ll.cc)
							#
						}
						#
						WEAK.ARCS <- c(WEAK.ARCS,inter[nex.min])
						#
						aj <- nex.min%%nt.system
						if(aj == 0) aj <- nt.system
						WHOLE.nexus <- rbind(WHOLE.nexus,c(aj,trunc(nex.min/nt.system)+1))
						#
						deci <- deci + 1
						#
						}
					}
					jump <- 1
					WEAK.copy.NAs <- as.matrix(WEAK.copy.NAs[,-rig])
					CRCC.NAs[nex.min] <- NA
				}
			#
			else 	{
				CRCC.NAs[nex.min] <- NA
				jump <- 2
				}
			}
		}
		#
		#
		#
		M.LINEAR <- CRCC.LINE
		diag(M.LINEAR) <- 0
		M.LINEAR[which(M.LINEAR < 1e-12)] <- 0
		M.CYCLES <- inter - CRCC.LINE
		names(WEAK.ARCS) <- NULL
		rownames(WHOLE.nexus) <- NULL
		N.CYCLE.distr <- CYCLE.distr/TST(inp,inter,outt,diss)
		cycles.count <- cycles.count[-1]
		#
		CY <- as.list(rep(NA,9))
		CY[[1]] <- FULL.CYCLES					## Total number of cycles removed
		CY[[2]] <- cycles.count					## Cycle analysis (number of cycles classified by cycle length)
		CY[[3]] <- WHOLE.nodes					## Full cycle analysis (classified by nexus)
		CY[[4]] <- WEAK.ARCS					## Link strength of the removed nexus
		CY[[5]] <- WHOLE.nexus					## Weak arcs representing nexus (ranked by removal order)
		CY[[6]] <- CYCLE.distr					## Cycle distribution
		CY[[7]] <- N.CYCLE.distr				## Normalized distribution
		CY[[8]] <- M.LINEAR						## Residual flows
		CY[[9]] <- M.CYCLES						## Aggregated cycles
		#
		return(CY)
		}
	#
	else{
		CY <- as.list(rep(NA,2))
		CY[[1]] <- 0
		CY[[2]] <- inter
		return(CY)
		}
	}
}

