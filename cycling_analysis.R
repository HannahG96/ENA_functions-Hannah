##############################CYCLING ANALYSIS#########################
#NOTE: steady-state of system is not checked
#######################################################################

cycling.analysis<-function(Z_cs,T_cs,E_cs,R_cs, select="Default", nmax=10000){
  
  ###Multiple weak arc selection: select weak arc going out from most connected nod OR weak arc that is heavier weighted
  if(select=="weighted"){select<-"weighted"}
  else{select<-"Default"}
  
  ####################################################################
  
  #STEP1
  #Identify most connected compartment as start point and order compartments downwards according to total 
  #links:
  restore.order<-colnames(T_cs)#store the initial order in a vector to rearrange matrix at the end of analysis
  
  connections<-matrix(rep(NA,nrow(T_cs)), nrow=nrow(T_cs), ncol=1, dimnames=list(c(1:nrow(T_cs)), NULL))
  for(i in 1:nrow(T_cs)){
    OUT<-length(which(T_cs[i,]!=0))
    IN<-length(which(T_cs[,i]!=0))
    connections[i,1]<-IN+OUT}
  connections<-connections[order(connections[,1], decreasing=TRUE),]
  nod_DFS<-colnames(T_cs)[as.numeric(names(connections))]#That is the order how nods will be analysed by DFS algorithm
  T_cs<-T_cs[as.numeric(names(connections)),as.numeric(names(connections))]#reorganize IA matrix according to order
  Z_cs<-Z_cs[as.numeric(names(connections))] #reorder Input
  E_cs<-E_cs[as.numeric(names(connections))] #reorder Output
  R_cs<-R_cs[as.numeric(names(connections))] #reorder Respiration
  
  
  
  #STEP2
  #Detect self-loops, store them as cycles+flow values and remove them
  #=if there are self-loops the DFS search algorithm may over/underestimate the number of cycle
  CYCLE1<-as.data.frame(matrix(rep(NA), nrow=0, ncol=length(nod_DFS)+1)) #store all cyclepaths of the foodweb

  ###Multiple weak arc selection based on cycle weights###
  if(select=="weighted"){
    multiWeak<-list() #list to store nod IDs and flow values of multiple weak arcs in cycles (stored as dataframe)
    multiW.element<-0 #indicates the number of elements in list
    CYCLE2<-as.data.frame(matrix(rep(NA), nrow=0, ncol=1, dimnames=list(NULL, "Flow.sum"))) } #store all cycleflows of the foodweb 
  ########################################################
  
  Weak<-vector() #vector to store weakest arc of each cycle
  WeakOUT<-vector()#vector to store outgoing nods of weak arcs
  WeakIN<-vector()#vector to store ingoing nods of weak arcs
  restore.selfloop<-vector()#vector to store selfloop position and restore flow value (selfloops are removed for DFS search)
  
  for(i in 1:nrow(T_cs)){
    if(rownames(T_cs)[i]==colnames(T_cs)[i] && T_cs[i,i]!=0){
      selfloop<-c(rownames(T_cs)[i], rownames(T_cs)[i], rep(NA, ncol(CYCLE1)-2))
      weakarc<-T_cs[i,i]
      CYCLE1<-rbind(CYCLE1, selfloop)
      
      ########################################################    
      ###Multiple weak arc selection based on cycle weights###    
      if(select=="weighted"){
        selfloopflow<-T_cs[i,i]
        CYCLE2<-rbind(CYCLE2, selfloopflow)}
      ########################################################
      
      WOut<-rownames(T_cs)[i]#store nod IDs of selfloop weak arc
      WIn<-rownames(T_cs)[i]
      WeakIN<-c(WeakIN,WIn)
      WeakOUT<-c(WeakOUT,WOut)
      Weak<-c(Weak, weakarc)#store selfloop weak arc flow 
      restore.selfloop<-c(restore.selfloop, i)#store position of selfloop
      T_cs[i,i]<-0} #remove selfloop from IA matrix
    else{next}}
  
  #STEP3
  #DFS algorithm using colors to enumerate+extract cycles from a grah
  #1.we visit first nod in the list (=grey)
  #2.we detect all cycles that include this nod+extract nod sequence & flow values
  #3.we mark this nod as visited (=black) and remove it from our search
  #method based on Ulanowicz (1983, 1986)
  
  color<-rep("WHITE", length(nod_DFS))#vector to store unvisited(=white)/being visited (=grey)/visited(=black) nods
  
  for(u in 1:length(nod_DFS)){
    
    if(length(Weak) >= nmax){  #LIMIT COMPUTING TIME: set a maximal number of cycles to be detected
      print("Reached limit of nmax removed cycles.")
      break}
    #1)visit nod following the list order
    mynod<-nod_DFS[min(which(color=="WHITE"))]
    path<-rep(NA, length(which(color=="WHITE")))#start a new DFS search pathway to detect all (remaining) cycles that include mynod
    flow<-path<-rep(NA, length(which(color=="WHITE")))#start a new flow vector to store flow values and extract weakest arcs
    DFS_path.rm <- rep(0, length(nod_DFS))#vector to store nods that are part of current pathway=nods can not appear twice in current pathway
    path[1]<-mynod
    DFS_path.rm[which(nod_DFS==mynod)]<- 1 #mynod now is part of current pathway
    color[which(nod_DFS==mynod)]<-"GREY" # mark this nod as being visited
    
    #2)create an edgelist to use as template=visualize pathways of DFS algorithm for mynod
    #-->edges are ordered from weakest to strongest
    #-->remove all edges towards mynod (this will be separatedly checked in each step of the search)
    DFS_path<-matrix(rep(NA), nrow=length(which(color!="BLACK")), ncol=length(which(color=="WHITE")))
    rownames(DFS_path)<-nod_DFS[which(color!="BLACK")]
    for(i in 1:nrow(DFS_path)){ 
      edges<-as.matrix(T_cs[which(color!="BLACK"), which(color=="WHITE")])
      if(nrow(edges)==0){
        DFS_path<-NA
        break }
      a<-unname(edges[i,])
      b<-order(a, decreasing = FALSE)[-c(1:length(which(a==0)))]
      DFS_path[i,]<-c(colnames(edges)[b], rep(NA,ncol(DFS_path)-length(b)))}
    
    if(length(which(is.na(DFS_path)==TRUE))==length(DFS_path)){    #if there are no existing connections to mynod, mark it as visited and move to next nod
      color[which(nod_DFS == mynod)] <- "BLACK"
      next }  
    if(length(which(is.na(DFS_path[which(rownames(DFS_path)==mynod),])==TRUE))==ncol(DFS_path)  ){   #if mynod has no (remaining) outgoing links, mark it as visited and move to next nod
      color[which(nod_DFS == mynod)] <- "BLACK"
      next }
    
    #3)Follow this template for the DFS search and store all (remaining) cycles that include mynod
    template <- DFS_path #copy of our path for restoring edges that need to be re-examined
    visited<-mynod #we leave starting "mynod" 
    visit<-unname(DFS_path[which(rownames(DFS_path)==mynod), 
                           min(which(is.na(DFS_path[which(rownames(DFS_path)==mynod),])==FALSE))]) #we visit the first nod linked to mynod
    DFS_path[which(rownames(DFS_path)==mynod), 
             min(which(is.na(DFS_path[which(rownames(DFS_path)==mynod),])==FALSE))] <-NA #remove this edge from DFS search to not repeat it
    backedges <-rep(NA, length(nod_DFS[which(color=="WHITE")])) #empty vector to store all detected backedges to mynod-->helps to detect only once each cycle
    
    repeat{
      
      if(length(Weak) >= nmax) break #LIMIT COMPUTING TIME: set a maximal number of cycles to be detected
      
      path[min(which(is.na(path)==TRUE))] <- visit  #add the currently visited nod to the current pathway
      DFS_path.rm[which(nod_DFS==visit)] <- 1 #mark this nod as part of the current pathway to not add it >1 in current pathway
      
      if(visited!=visit){
        flow[min(which(is.na(flow)==TRUE))] <- T_cs[which(rownames(T_cs)==visited), which(colnames(T_cs)==visit)] }#add corresponding flow of the edge
      
      ##SEARCH FOR A BACKEDGE AND IF FOUND STORE CYCLE PATTERN AND WEAKEST ARC    
      if(T_cs[which(rownames(T_cs)==visit), which(colnames(T_cs)==mynod)]!=0 &&
         length(which(backedges==visit))==0){ #check if there exists a backedge to mynod+this backedge has not been recorded yet
        
        mycycle<-c(path[which(is.na(path)==FALSE)], mynod, rep(NA, ncol(CYCLE1)-length(which(is.na(path)==FALSE))-1))
        CYCLE1<-rbind(CYCLE1, mycycle) #if YES, store cycle
        colnames(CYCLE1)<-NULL #rbind attributes some weird columnnames to the dataframe...
        
        cycleflow<-c(flow[which(is.na(flow)==FALSE)], 
                     T_cs[which(rownames(T_cs)==visit), which(colnames(T_cs)==mynod)]) #extract cycle flows
        
        weakarc<-min(cycleflow)
        WOut<-mycycle[which(cycleflow==weakarc)]#indicates nod of outgoing weak arc
        WIn<-mycycle[which(cycleflow==weakarc)+1]#indicates nod of ingoing weak arc 
        
        #######################################################################################################   
        ##WEAK ARC SELECTION (IF MULTIPLE WEAK ARCS IN ONE CYCLE)
        #######################################################################################################
        
        ###Multiple weak arc selection based on cycle weights###
        if(select=="weighted"){
          cycleflow.sum<-sum(cycleflow) #calculate sum of flows transported through cycle
          CYCLE2<-rbind(CYCLE2, cycleflow.sum)
          colnames(CYCLE2)<-"Flow.sum" #rbind attributes some weird columnnames to the dataframe...
          
          if(length(WOut)>1){
            multiW.element<-multiW.element+1
            multiW<-data.frame(Out=WOut, In=WIn, Flow=weakarc)
            multiWeak[[multiW.element]]<-multiW #define v!!!
            WOut<-NA
            WIn<-NA  }   }
        
        ###Default Option:pick the weak arc going out of the most connected nod
        else{
          if(length(WOut)>1 ){
            rank<-rep(NA, length(WOut))
            for(l in 1:length(WOut))rank[l]<-which(nod_DFS==WOut[l])
            weakSelect<-which(rank==min(rank))
            if(length(weakSelect)>1)print("PROBLEM")
            WOut<-WOut[weakSelect]
            WIn<-WIn[weakSelect]  } } #store weak arc flow in matrix
        ###################################################################################
        
        Weak<-c(Weak,weakarc)#store weakest arc of the cycle
        WeakOUT<-c(WeakOUT, WOut)#store corresponding outgoing nod
        WeakIN<-c(WeakIN, WIn)#store corresponding ingoing nod
        backedges[min(which(is.na(backedges)==TRUE))] <- visit  #note down nodname from which goes out the backedge to not repeat cycle detection while backtracking
        path[max(which(is.na(path)==FALSE))] <- NA #we remove this nod from current pathway (as it will be added again at begin of loop)
        flow[max(which(is.na(flow)==FALSE))] <-NA #same approach for flow vector
        next }
      
      ###IF NO CYCLE IS FOUND, VISIT NEXT NOD
      visited <- visit #we leave the currently visited nod
      
      #####
      #PROBLEM=a nod CAN NOT APPEAR >1 IN PATH 
      #SOLUTION=only visit nods which are not part of current pathway
      for(p in 1:ncol(DFS_path)){
        visit<-unname(DFS_path[which(rownames(DFS_path)==visited), p])  #attention shifts to next connected nod
        if(is.na(visit)==TRUE){success<-0 ; next}
        if(DFS_path.rm[which(nod_DFS==visit)]==1){success<-0 ; next} #if this nod is part of current pathway we look for other connected nods
        else{success<-1 ; break} }
      #######
      
      if(length(which(is.na(DFS_path[which(rownames(DFS_path)==visited),])==FALSE)) != 0 &&
         success == 1){ #check if there remain unvisited edges going out from current visited nod that are NOT part of current pathway
        
        DFS_path[which(rownames(DFS_path)==visited), 
                 which(DFS_path[which(rownames(DFS_path)==visited),]==visit)] <-NA } #if YES we remove the examined edge from DFS search to not repeat it
      
      ###IF THERE DO NOT REMAIN UNVISITED EDGES THAT ARE NOT PART OF CURRENT PATHWAY: BACKTRACK    
      else{ 
        visit<-visited #we go back to the lastly visited nod
        path[which(path==visit)]<-NA #remove the lastly visited nod from current pathway
        backedges[which(backedges==visit)] <-NA #remove nod also (if present) from the backedges list
        DFS_path.rm[which(nod_DFS==visit)] <- 0 #we mark this nod as not part of current pathway
        DFS_path[which(rownames(DFS_path)==visit), ] <- template[which(rownames(DFS_path)==visit), ]#restore all edges of DFS path going out from this nod--> need to be re-examined
        
        if(length(which(is.na(path)==TRUE))==length(path)) {
          color[which(nod_DFS == mynod)] <- "BLACK" #if there remain no nods in the current pathway we are done with mynod and can remove it from search
          break} # EXIT INNER LOOP
        
        else{ 
          flow[max(which(is.na(flow)==FALSE))] <-NA #remove corresponding flow value
          visit<-path[max(which(is.na(path)==FALSE))]  #attention shifts again to previously visited nod
          if(visit==mynod){ #current pathway consists only of 1 nod 
            visited<-mynod } #we define the before visited nod as the same as "visit"-->there is no flow value=means we do not need to store a flow at the start of next inner loop
          else{
            visited<-path[max(which(is.na(path)==FALSE))-1] #redefine the before visited nod (as it is needed to define flow value at start of inner loop)
            flow[max(which(is.na(flow)==FALSE))] <- NA } #remove 2nd flow values-->this one will be restored at begin of loop
          
          path[max(which(is.na(path)==FALSE))] <- NA #remove this nod from current pathway (as it will be added again at begin of next loop)
          next } }
      
    } #END OF INNER LOOP
  } #END OF OUTER LOOP
  
if(length(Weak)==0 ||
   length(which(round(Weak, digits=9)==0))==length(Weak)){ #DFS SEARCH DID NOT DETECT CYCLES IN NETWORK-->WE ARE FINISHED
    cyc_cs<-list()
    cyc_cs[[1]]<-0 #Nb of cycles=0
    names(cyc_cs)<-"Number of removed cycles"
}else{
  
  #STEP4
  #Identify nods involved in weak arcs + store nod IDs of weakest arcs ordered from weakest to strongest
  
  #######################################################################################################
  ###Multiple weak arc selection based on cycle weights###
  #REVIEW!!!!-->evtl simplify???
  
  if(select=="weighted" && length(multiWeak)!=0){
    
    #1.identify multiple weak arcs and calculate corresponding weights
    #=sum of cycleflows that share this weak arc
    
    multis<-as.data.frame(matrix(rep(NA), nrow=0, ncol=3, dimnames=list(NULL, c("Out", "In","Flow"))))
    for(p in 1:length(multiWeak))multis<-rbind(multis, multiWeak[[p]])
    multis<-unique(multis)
    multis$Weight<-NA
    multis.membercycles<-vector("list", length=nrow(multis))  
    
    for(i in 1:nrow(multis)){
      members<-vector()
      
      for(p in 1:length(Weak)){ #PART 1: check which cycles with a unique weak arc share a multi weak arc
        if(is.na(WOut[p]==TRUE)){next}
        if(Weak[p]==multis[i,"Flow"] && WeakOUT[p]==multis[i,"Out"] && WeakIN[p]==multis[i,"In"]){
          members<-c(members, p) } 
        else{next} }
      
      undefinedW<-which(is.na(WeakOUT)==TRUE) #PART 2: check which cycles with multiple weak arcs share a multi weak arc
      
      for(p in 1:length(undefinedW)){
        success<-0
        
        for(j in 1:nrow(multiWeak[[p]])){
          if(multiWeak[[p]][j,"Flow"]==multis[i,"Flow"] && multiWeak[[p]][j,"Out"]==multis[i,"Out"] && 
             multiWeak[[p]][j,"In"]==multis[i,"In"]){success<-1}else{next} }
        
        if(success==1){members<-c(members, undefinedW[p])} } 
      
      multis.membercycles[[i]] <- members 
      multis[i, "Weight"]<-sum(CYCLE2[members,"Flow.sum"])}
    
    multis.membercycles<-multis.membercycles[order(multis[,"Weight"], decreasing=FALSE)]
    multis<-multis[order(multis[,"Weight"], decreasing=FALSE),]
    
    #2.Define weak arc IDs for cycles with multiple weakarcs based on the heaviest weighted weak arc
    for(i in 1:nrow(multis)){
      WeakOUT[multis.membercycles[[i]]]<-multis[i, "Out"]
      WeakIN[multis.membercycles[[i]]]<-multis[i, "In"] }
    
  }#END OF OUTER LOOP
  #####################################################################################
  ###(Default Option: multiple weak arc selection based on most connected outgoing nod
  ###-->defined during DFS search)
  
  
  #STEP5 
  #restore selfloops in T_cs
  for(i in 1:length(restore.selfloop)) T_cs[restore.selfloop[i], restore.selfloop[i]]<-Weak[i]
  
  #STEP6
  #####NEXUS REMOVAL
  #Note: ordering vectors and matrix according to weak arc flow values leads to ERRORS!
  #Note: counting of matrix elements is done by columns
  
  CYCLE1_rm<-as.matrix(CYCLE1)#matrix indicating all cycles to be removed=0 rows at end
  if(select=="weighted")CYCLE2_rm<-as.matrix(CYCLE2)#+corresponding total cycleflows (only needed for multiple weak arc selection based on weights)
  WeakOUT_rm<-WeakOUT#+corresponding outgoing nods
  WeakIN_rm<-WeakIN#+corresponding ingoing nods
  Weak_rm <- Weak #+corresponding weak arcs of cycles from
  
  rank_nod<-nod_DFS #vector to order nods based on number of based on number of out/in- links
  
  if(select=="weighted"){ 
    multiWeak<-list()#list to store multiple weak arcs of cycles if selection is based on cycle weights
    whichcycle<-1} #object to define list element
  
  nex<-list()#list to store membercycles of removed nexus=stored as matrix
  nexMEMBER<-1 #object defines list element of nex
  
  WEAK<-vector()
  WEAK_out<-vector()
  WEAK_in<-vector()
  Membercycles<-vector()
  if(select=="weighted")WEIGHT<-vector()
  
  fraction_WEAK<-vector() #vector to store fraction of weak arc attributed to each cycle of the nexus=needed to calculate cycledistribution
  cyclelength<-vector() #+vector to store corresponding length of the cycle
  
  T_cs_acyclic <- T_cs
  
  repeat{ 
    
    #1.Define weak arc flow of the nexus to be removed
    #=smallest flow + most connected outgoing (+ingoing) nod 
    #                      OR involved in heaviest cycles
    
    myNex<-min(Weak_rm)
    
    #Check for multiple weak arc candidates:
    a<-which(Weak_rm==myNex)
    if(length(unique(WeakOUT_rm[a]))>1 || length(unique(WeakIN_rm[a]))>1){
      
      ###Multi-weak arc selection based on cycle weights:      
      if(select=="weighted"){
        myNexOUT<-unique(WeakOUT_rm[a])
        myNexIN<-unique(WeakIN_rm[a])
        nexcandidate<-vector("list", length=length(myNexOUT)*length(myNexIN))#list to store weak arc IDs of candidates
        Weights<-rep(NA, length(myNexOUT)*length(myNexIN))#empty vector to store weight of each weak arc candidate
        candidate.ele<-1 #defines corresponding element of list and weightvector
        for(p in 1:length(myNexOUT)){
          for(l in 1:length(myNexIN)){
            nexcandidate[[candidate.ele]]<-c(myNexOUT[p], myNexIN[l])#define the weak arc candidate
            aa<-a[which(WeakOUT_rm[a]==nexcandidate[[candidate.ele]][1])]
            aaa<-aa[which(WeakIN_rm[aa]==nexcandidate[[candidate.ele]][2])]#define membercycles of this weak arc
            Weights[candidate.ele]<-sum(CYCLE2_rm[aaa, "Flow.sum"])#store sum of weights of these membercycles
            candidate.ele<-candidate.ele+1 }  }
        if(length(which(Weights==max(Weights)))>1)print("PROBLEM1") #CONTROL
        myNexOUT<-nexcandidate[[which(Weights==max(Weights))]][1]
        myNexIN<-nexcandidate[[which(Weights==max(Weights))]][2]
      }#END OF WEIGHTED SELECTION CLAUSE
      
      ###Default: Multi-weak arc selection based on most connected outgoing (+ingoing) nod:      
      else{
      myNexOUT<-unique(WeakOUT_rm[a])
      rank<-rep(NA, length(myNexOUT))
      for(l in 1:length(myNexOUT)){
        rank[l]<-which(rank_nod==myNexOUT[l])}
      myNexOUT<-myNexOUT[which(rank==min(rank))]
      if(length(myNexOUT)>1)print("PROBLEM2") #CONTROL
      
      aa<-a[which(WeakOUT_rm[a]==myNexOUT)]
      myNexIN<-unique(WeakIN_rm[aa])
      rank<-rep(NA, length(myNexIN))
      for(l in 1:length(myNexIN)){
        rank[l]<-which(rank_nod==myNexIN[l])}
      myNexIN<-myNexIN[which(rank==min(rank))]
      if(length(myNexIN)>1)print("PROBLEM3") #CONTROL
       }#END OF DEFAULT SELECTION
    } #END OF MULTIPLE WEAK ARC CANDITATE SELECTION
    
    else{ #only 1 weak arc candidate
      myNexOUT<-unique(WeakOUT_rm[a])
      myNexIN<-unique(WeakIN_rm[a])  
      if(select=="weighted")Weights<-NA}
    
    
    
    
    #2.Detect membercycles of the weak arc nexus to be removed
    
    membercycles<-vector()
    
    for(p in 1:nrow(CYCLE1_rm)){
      if(Weak_rm[p]==myNex && WeakOUT_rm[p]==myNexOUT && WeakIN_rm[p]==myNexIN){
        membercycles<-c(membercycles, p)}
      else{next} }
    
    #3.Remove nexus from foodweb
    
    if(length(membercycles)>1){ #several cycles are part of nexus
      
      #PART 1: calculate circuit probabilities of nexus membercycles    
      circuitProb<-rep(NA, length(membercycles))
      
      for(j in 1:length(membercycles)){ 
        TotOut<-rep(NA, length(which(is.na(CYCLE1_rm[membercycles[j],])==FALSE))-1)
        cycleflows<-rep(NA, length(which(is.na(CYCLE1_rm[membercycles[j],])==FALSE))-1)
        for(p in 1:length(TotOut)){
          TotOut[p]<-sum(T_cs_acyclic[which(rownames(T_cs_acyclic)==CYCLE1_rm[membercycles[j],p]),]) +
            E_cs[which(names(E_cs)==CYCLE1_rm[membercycles[j],p])] + #export flows
            R_cs[which(names(R_cs)==CYCLE1_rm[membercycles[j],p])]  #respiration flows
          cycleflows[p]<-T_cs_acyclic[which(rownames(T_cs_acyclic)==CYCLE1_rm[membercycles[j],p]),
                                      which(colnames(T_cs_acyclic)==CYCLE1_rm[membercycles[j],p+1])]}
        fracOut<-cycleflows/TotOut
        circuitProb[j]<-prod(fracOut)}
      
      #PART 2: remove cycle in proportion to circuit probability  
      for(j in 1:length(membercycles)){ 
        cycle_rm<-matrix(0,nrow=length(nod_DFS), ncol=length(nod_DFS), 
                         dimnames=list(rownames(T_cs), colnames(T_cs)))#matrix indicating flows values to be removed from each cycle link 
        
          for(p in 1:length(which(is.na(CYCLE1_rm[membercycles[j],])==FALSE))-1){
            cycle_rm[which(rownames(cycle_rm)==CYCLE1_rm[membercycles[j],p]),
                     which(colnames(cycle_rm)==CYCLE1_rm[membercycles[j],p+1])] <- myNex *
                                                                        circuitProb[j]/sum(circuitProb) }
          
          fraction_WEAK<-c(fraction_WEAK, myNex * circuitProb[j]/sum(circuitProb)) #store removed weak arc fraction=needed to calculate cycle distribution
          cyclelength<-c(cyclelength, length(which(is.na(CYCLE1_rm[membercycles[j],])==FALSE))-1)#store corresponding length of the cycle
        
        T_cs_acyclic <- T_cs_acyclic - cycle_rm } #redefine IA_matrix after cycle removal
         #END OF REMOVAL CONDITION 1
      
    } else { 
      #weak arc is unique to the cycle=remove (remaining) weakest arc flow value from each cycle link
      #=if(length(membercycles)==1)
      cycle_rm<-matrix(0,nrow=length(nod_DFS), ncol=length(nod_DFS), 
                       dimnames=list(rownames(T_cs), colnames(T_cs)))#matrix indicating flows values to be removed from each cycle link 
      
      for(p in 1:length(which(is.na(CYCLE1_rm[membercycles,])==FALSE))-1){
        cycle_rm[which(rownames(cycle_rm)==CYCLE1_rm[membercycles,p]),
                 which(colnames(cycle_rm)==CYCLE1_rm[membercycles,p+1])] <- myNex  }
      
      fraction_WEAK<-c(fraction_WEAK, myNex) #store removed flow fraction=needed to calculate cycle distribution
      cyclelength<-c(cyclelength, length(which(is.na(CYCLE1_rm[membercycles,])==FALSE))-1) #store corresponding length of the cycle
      
      T_cs_acyclic <- T_cs_acyclic - cycle_rm #redefine IA_matrix after cycle removal
    } #END OF REMOVAL CONDITION 2
    
    
    #4. Store membercycles of the removed nexus in list + corresponding flow fraction of each cycle
    #-->needed for function output of full cycle analysis
    nexfullcycles<-list(myNex, myNexOUT, myNexIN, CYCLE1_rm[membercycles,])
    if(length(membercycles)==1){
      if(nrow(CYCLE1_rm)==1){
      nexfullcycles[[4]]<-as.data.frame(CYCLE1_rm[membercycles,])}
      else{nexfullcycles[[4]]<-as.data.frame(t(CYCLE1_rm[membercycles,]))} }
    nex[[nexMEMBER]]<-nexfullcycles
    nexMEMBER<-nexMEMBER+1
    
    #5. Reorder nods based on total number of out/in-going links
    #-->cycleflows are rounded till 9th digit=this value is picked randomly!
    totlinks<-rep(NA,length(rank_nod))
    for(i in 1:ncol(T_cs_acyclic)){
      outlinks<-unname(T_cs_acyclic[i,])
      outlinks<-round(outlinks, digits=9)
      inlinks<-unname(T_cs_acyclic[,i])
      inlinks<-round(inlinks, digits=9)
      OUT<-length(which(outlinks!=0))
      IN<-length(which(inlinks!=0))
      totlinks[i]<-IN+OUT}
    rank_nod<-colnames(T_cs_acyclic)[order(totlinks, decreasing=TRUE)]
    
    #6.Remove the removed cycles from CYCLE1_rm, (CYCLE2_rm), Weak_rm, WeakOUT_rm and Weak_IN_rm
    CYCLE1_rm<-CYCLE1_rm[-membercycles, ]
    if(select=="weighted")CYCLE2_rm<-matrix(CYCLE2_rm[-membercycles, ], dimnames=list(NULL,"Flow.sum"))
    Weak_rm<-Weak_rm[-membercycles]
    WeakOUT_rm<-WeakOUT_rm[-membercycles]
    WeakIN_rm<-WeakIN_rm[-membercycles]
    
    if(length(Weak_rm)==1){
      CYCLE1_rm<-as.data.frame(t(CYCLE1_rm))}
    
    if(length(Weak_rm)==0){
      WEAK<-c(WEAK, rep(myNex, length(membercycles)))
      WEAK_out<- c(WEAK_out, rep(myNexOUT, length(membercycles)))
      WEAK_in<- c(WEAK_in, rep(myNexIN, length(membercycles)))
      Membercycles<-c(Membercycles, rep(length(membercycles), length(membercycles)))
      if(select=="weighted")WEIGHT<-c(WEIGHT, rep(max(Weights), length(membercycles)))
      break }#EXIT REPEAT LOOP=ALL NEXUS HAVE BEEN REMOVED
    
    #7. Redefine weak arcs after nexus removal
    
    for(p in 1:length(Weak_rm)){
      flow<-rep(NA, length(which(is.na(CYCLE1_rm[p,])==FALSE))-1)
      
      for(j in 1:length(flow)){
        flow[j]<-T_cs_acyclic[which(rownames(T_cs_acyclic)==CYCLE1_rm[p, j]),
                              which(colnames(T_cs_acyclic)==CYCLE1_rm[p, j+1])]}
      
      if(select=="weighted")CYCLE2_rm[p, "Flow.sum"]<-sum(flow) #update cycle weight after weak arc removal
      
      weak_evaluate<-which(flow==min(flow))
      
      ###MULTIPLE WEAK ARC SELECTION#################################################
      if(length(weak_evaluate)>1 ){ #>1 weak arc canditates for update
        
        ###############################################################################
        ###Multiple weak arc selection based on cycle weights###
        #-->weak arcs are defined after all cycle weights have been recalculated
        if(select=="weighted"){
          multiWeakOUT<-unname(CYCLE1_rm[p, weak_evaluate])
          multiWeakIN<-unname(CYCLE1_rm[p, weak_evaluate+1])
          multiWeak[[whichcycle]]<-list(multiWeakOUT, multiWeakIN)
          whichcycle<-whichcycle+1
          weak_evaluate<-"undefined" } 
        
        ###############################################################################
        ###Default: multiple weak arc selection based on most connected outgoing nod
        else{ 
          rank<-rep(NA, length(weak_evaluate))
          for(l in 1:length(weak_evaluate)){
            rank[l]<-which(rank_nod==CYCLE1_rm[p, weak_evaluate[l]])}
          weak_evaluate<-weak_evaluate[which(rank==min(rank))]
          if(length(weak_evaluate)>1)print("PROBLEM4") #Control
        }  #END OF DEFAULT SELECTION CLAUSE
        
      } #END OF MULTIPLE WEAK ARC SELECTION WITHIN LOOP
      ##################################################################################### 
      
      #Update weak arcs in Weak_rm, WeakOUT_rm, WeakIN_rm:
      
      ###Part 1: Multiple weak arc selection based on cycle weights##########################
      #-->mark weak arcs as undefined=weak arcs are defined after all cycle weights have been recalculated        
      if(select=="weighted" && weak_evaluate=="undefined"){
        Weak_rm[p]<-NA 
        WeakOUT_rm[p]<-NA
        WeakIN_rm[p]<-NA }
      ###############################################################################
      
      else if(Weak_rm[p]!=flow[weak_evaluate] || 
         WeakOUT_rm[p]!=CYCLE1_rm[p, weak_evaluate] ||
         WeakIN_rm[p]!=CYCLE1_rm[p, weak_evaluate+1]){ #NEW WEAK ARC!  
        
        Weak_rm[p]<-flow[weak_evaluate] #update weak arc flow
        WeakOUT_rm[p]<-CYCLE1_rm[p, weak_evaluate] #update corresponding nod IDs
        WeakIN_rm[p]<-CYCLE1_rm[p, weak_evaluate+1]}
      
      else{next}
    }#END OF WEAK ARC UPDATE CLAUSE (OUTER FOR LOOP)
    
    ###########################################################################################
    ###Part 2: Multiple weak arc selection based on cycle weights##############################
    #-->define missing weak arcs based on updated cycle weights 

    if(select=="weighted" && length(which(is.na(Weak_rm)==TRUE))!=0) {
      undefinedW<-which(is.na(Weak_rm)==TRUE)
      multiWeak_rm<-multiWeak #cycles are removed step by step from this list as weak arc of cycle gets defined
      
      for(p in 1:length(undefinedW)) {
        candidateWOUT<-multiWeak[[p]][[1]] #identify weak arc candidates of each cycle with undefined weak arc 
        candidateWIN<-multiWeak[[p]][[2]]
        members<-vector("list", length=length(candidateWOUT))#list to store which cycles share these weak arc candidates
        
        for(l in 1:length(candidateWOUT)){
          #1: 
          #Search for all cycles with unique weak arc that share this candidate
          a<-which(WeakOUT_rm==candidateWOUT[l])
          members[[l]]<-a[which(WeakIN_rm[a]==candidateWIN[l])]
          #2: 
          #Search for all cycles with undefined weak arc that share this candidate
          undefinedW_rm<-which(is.na(Weak_rm)==TRUE)
          for(j in 1:length(multiWeak_rm)){
            aa<-which(multiWeak_rm[[j]][[1]]==candidateWOUT[l])
            if(length(aa)>1){print("PROBLEM5")}
            if(length(aa)==0){next}
            if(multiWeak_rm[[j]][[2]][aa]==candidateWIN[l]){members[[l]]<-c(members[[l]], undefinedW_rm[j])}
            else{next}  } }
        
        #3: 
        #Compare sums of weight of membercycles of weak arc candidates and define weak arc for cycle
        rank<-rep(NA, length(candidateWOUT))
        
        for(l in 1:length(rank)){
          rank[l]<-sum(CYCLE2_rm[members[[l]], "Flow.sum"])}
        WeakOUT_rm[undefinedW[p]]<-multiWeak[[p]][[1]][which(rank==max(rank))]
        WeakIN_rm[undefinedW[p]]<-multiWeak[[p]][[2]][which(rank==max(rank))]
        Weak_rm[undefinedW[p]]<-T_cs_acyclic[which(rownames(T_cs_acyclic)==WeakOUT_rm[undefinedW[p]]),
                                             which(colnames(T_cs_acyclic)==WeakIN_rm[undefinedW[p]])]
        #4:
        #Remove the cycle from multiWeak_rm as weak arc has been defined
        multiWeak_rm[[1]] <- NULL 
      } #END OF FOR LOOP-->Switch to next cycle with undefined weak arc
      
      #Restore storage space for new multi weak arc selection:
      multiWeak<-list()
      whichcycle<-1  
      
   }  #END OF MULTIPLE WEAK ARC SELECTION BASED ON CYCLE WEIGHTS (Part 2)
    ######################################################################################    
    
    
    #8.store the NEXUS ID+number of membercycles (+Nexus weight)
    WEAK<-c(WEAK, rep(myNex, length(membercycles)))
    WEAK_out<- c(WEAK_out, rep(myNexOUT, length(membercycles)))
    WEAK_in<- c(WEAK_in, rep(myNexIN, length(membercycles)))
    Membercycles<-c(Membercycles, rep(length(membercycles), length(membercycles)))
    if(select=="weighted")WEIGHT<-c(WEIGHT, rep(max(Weights), length(membercycles)))
    
  }#END OF REPEAT LOOP
  
  #STEP7
  #1.Create summary table of removed nexus:
  if(select=="weighted"){
  Weak_ID<-data.frame(Outgoing=WEAK_out,Ingoing=WEAK_in, Flow.removed=WEAK, Flow=NA, 
                      Weight=WEIGHT, Nb.cycles=Membercycles)
  }else{
  Weak_ID<-data.frame(Outgoing=WEAK_out,Ingoing=WEAK_in, Flow.removed=WEAK, Flow=NA, 
                      Nb.cycles=Membercycles) }
  Weak_ID<-unique(Weak_ID) 
  
  #2.Add the flow values of the weak arcs, as they occur in the T_cs matrix
  for(i in 1:nrow(Weak_ID)){
    Weak_ID[i, "Flow"]<-T_cs[which(rownames(T_cs)==Weak_ID[i, "Outgoing"]),
                             which(colnames(T_cs)==Weak_ID[i, "Ingoing"])] }
  
  
  
  ###FUNCTION OUTPUT###
  cyc_cs<-list()
  
  ##1. total number of cycles removed
  cyc_cs[[1]]<-nrow(CYCLE1)
  
  ##2. Cycle analysis (number of cycles classified by cycle length):
  cyc_length<-matrix(0, nrow=length(c(1:max(unique(cyclelength)))), ncol=2, dimnames=list(NULL, c("length", "Nb.cycles")))
  cyc_length[,"length"]<-c(1:max(unique(cyclelength)))
  cyc_length[,"Nb.cycles"]<-sapply(cyc_length[,"length"], FUN=function(x)length(which(cyclelength==x)))
  cyc_cs[[2]]<-cyc_length
  #What happens if cyc_length consists only of 1 row??? Should work...-->TO CHECK!!
  
  ##3. full cycle analysis (classified by nexus) ---> list of compartment IDs, see Fig. 2(b) at page 222
  ## the order in which cycles are listed is the order of removal
  #(starts from those with the weakest nexus (nexus = weakest link of each cycle)??)
  fullcycles<-vector("list", length=length(nex))
  namesnex<-rep(NA, length(nex))
  for(i in 1:length(nex)){
   nexcycles<-vector("list", length=nrow(nex[[i]][[4]]))
  for(j in 1:length(nexcycles))nexcycles[[j]]<-unname(nex[[i]][[4]][j, which(is.na(nex[[i]][[4]][j, ])==FALSE)])
   namesnex[i]<-paste(nex[[i]][[1]], ";", nex[[i]][[2]], "->", nex[[i]][[3]], sep=" ") 
   fullcycles[[i]]<-nexcycles}
   names(fullcycles)<-namesnex
  
  cyc_cs[[3]]<-fullcycles #Looks a bit ugly
  
  
  ##4.Nexus analysis
  #     link strength of the removed nexus: removed flow value & flow value of T_cs matrix
  #         + compartments' IDs of removed nexus
  #         + number of cycles that share the nexus
  #ordered by the order of removal
  cyc_cs[[4]]<-Weak_ID
  
  
  ##6. cycle distribution ---> amount of matter recycled by each compartment-->NEED FOR A MORE ELEGANT SOLUTION!!!
  cyc_distribution<-rep(NA, nrow(cyc_length))
  for(i in 1:nrow(cyc_length)){
    lengthmembers<-which(cyclelength==cyc_length[i, "length"])
    cyc_distribution[i]<-sum(fraction_WEAK[lengthmembers])*cyc_length[i, "length"]}
  cyc_dist<-cyc_distribution
  names(cyc_dist)<-cyc_length[,"length"]
  cyc_cs[[5]]<-cyc_dist
  
  ##7. normalized distribution ---> normalized contribution to cycling of each compartment
  TST<-sum(T_cs)+sum(R_cs)+sum(E_cs)+sum(Z_cs)
  cyc_dist.norm<-cyc_distribution / TST
  names(cyc_dist.norm)<-cyc_length[,"length"]
  cyc_cs[[6]]<-cyc_dist.norm
  
  ##8. residual flows ---> residual acyclic network after the removal of all cycles, see Fig. 2(c) at page 222
  resflows<-T_cs_acyclic[restore.order, restore.order]
  cyc_cs[[7]]<-resflows
  
  ## aggregated cycles ---> sum of all material flowing through the five cycles, see Fig. 2(b) at page 222
  aggrcycles<-(T_cs - T_cs_acyclic)[restore.order, restore.order]
  cyc_cs[[8]]<-aggrcycles

names(cyc_cs)<-c("Number of removed cycles","Cycle length distribution",
                 "Full cycle analysis","Nexus analysis","Cycle distribution",
                 "Normalized cycle distribution","Residual flows",
                 "Aggregated cycles")
}#END OF OUTER ELSE CLAUSE=cycles were detected in network
  
  return(cyc_cs)
  
} #END OF FUNCTION




