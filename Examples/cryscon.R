
cryscon <- as.matrix(read.table(file = "Examples/cryscon.txt", header = TRUE, sep = "\t"))
rownames(cryscon) <- colnames(cryscon)

## total number of nodes
tn <- ncol(cryscon) - 3

## non-living nodes
nl <- 1

## import vector (Z)
Z <- cryscon["IMPORT",2:(tn+1)]

## export vector (E)
E <- cryscon[2:(tn+1),"EXPORT"]

## respiration vector (R)
R <- cryscon[2:(tn+1),"RESP."]

## matrix of intercompartmental exchanges [T]
T <- cryscon[2:(tn+1),2:(tn+1)]
