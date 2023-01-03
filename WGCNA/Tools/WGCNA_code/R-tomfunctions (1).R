#****************************************************************************************
#
#                       Gene co-Expression Network Analysis Package
#         
#
#        Objective: 
#                   Automatic construction of gene co-expression networks &
#                   Decomposition of a network into tightly co-expressed modules
#
#        Input:     DNA microarray data file
#
#        Output:    Topological overlap matrix (TOM) and gene co-expressed modules
#
#        Authors:   Bin Zhang and Steve Horvath
#
#                   Array Data Analysis Group (ADAG)
#                   Departments of Human Genetics and Biostatistics
#                   University of California at Los Angeles                
#
#        Contact:   {binzhang,shorvath}@mednet.ucla.edu
#
#        Copyright  @2004-2008 ADAG,UCLA
#
#        Date:      Aug. 18, 2004
#
#        
#        REV            DATE            BY           DESCRIPTION
#
#        
#****************************************************************************************

#------------------------------------------------------------------------------------------------------------------
#Function: cluster idnetification based on dynamic programming
#
#
#
#Input:
#
#    adjmatrix           ~ adjacency matrix
#    hierclust           ~ hierarchical clustering structure
#    minModuleSize       ~  min size of module
#
#    Choices of clustering coefficient (L= # of links): 
#       1) L/n,   favour bigger modules 
#       2) L/n^2, favour small clusters as small modules are likely to have high CC
#       3) L/(n*1.5+const), no bias to big or small modules
#
#       4) foldchange_cc is used to terminate the search when the module coherence is already small
#
#       5) apower is used for defining module efficiency CC*n^apower, so smaller apower favour small modules
#           and so bigger apower favour big modules
#
#Output:  module colors for every node
# 
#------------------------------------------------------------------------------------------------------------------
identifyModuleByDynamicProgram = function(adjmatrixSorted, hierclust, minModuleSize=5, npower=1.5, nconst=50, maxModules=169, heicutoff=0.99, useNumberAsLabel=F, startlabel=1, foldchange_cc=100.0)
{
 
  nodes = dim(adjmatrixSorted)[1]
  cc    = matrix(0, nodes, nodes)

  # find break points to avoid false modules
  #
  #staticClu  = cutTreeStatic(hiercluster=hierclust, heightcutoff=heicutoff, minsize1=0)

  # make sure that there is at least two clusters after the static cut
  curheicutoff=heicutoff
  heistep     = 0.02
  no.bounds   = 0
  while(no.bounds <=1 & curheicutoff>0){
    staticClu  = cutTreeStatic(hiercluster=hierclust, heightcutoff=curheicutoff, minsize1=0)
    ostaticClu = staticClu[ h1row$order] #ordered
    bounds     = findBoundary(ostaticClu)
    bounds     =c(bounds, nodes+1)
    no.bounds  = length(bounds)
    curheicutoff = curheicutoff - heistep
  }
  if (curheicutoff<=0 | no.bounds <=1){     
     return(rep("grey", nodes))
  }

  breakpoints= rep(0, nodes)
  for(i in c(1:(no.bounds-1) )){
     # index of elements belong to the current cluster i
     # as bounds[i+1] is the start point of next cluster, 
     # bounds[i+1]-1 is the end point of the cluster i
     idx = c(bounds[i]:(bounds[i+1]-1))
     breakpoints[idx] = bounds[i+1]-1
  }
  breakpoints

  #assign colors for modules
  module.assign = rep(0, nodes)
  module.cnt=1

  # 1) number of links of a network with nodes from 1 to j, use an accumulated ways to save time
  i= 1
  for (j in c((i+1):nodes) ){
     cc[i,j] = cc[i, j-1] + sum(adjmatrixSorted[j, c(i:j)] )*2
  }

  # 2) number of links of a network with nodes from  i to j
  for (i in c(2:nodes) ){
    n = cumsum( adjmatrixSorted[i-1, ] )
    n = n-n[i-1]
    idx = c((i+1):nodes)
    if (i+1>nodes){
       break
    }

    cc[i, idx ] = cc[i-1, idx] - n[idx]*2
  }

  # 2.1) enforce breakpoints
  for (i in c(1:(nodes-1)) ){
      if(breakpoints[i] <nodes){
          bidx         = c((breakpoints[i]+1):nodes)
      }else{
          bidx         = c(breakpoints[i]:nodes)
          break
      }
      cc[i, bidx ] = 0 
  }

  # 3) normalization by the total number of nodes in each network
  for (i in c(1:nodes) ){
    n       = rep(1, nodes)
    idx     = c(i:nodes)
    n[idx]  = idx-i+1

    apower = n**npower + nconst
    #apower = n + nconst

    #cc[i, ] = cc[i,]/(apower)
    #cc[i, ] = cc[i,]/n

    #cc[i,]  = ifelse(n<minModuleSize, 0, cc[i,]/n)
    #cc[i, ] = cc[i,]/(apower)    

    #cc[i, ] = cc[i,]/(n-1)/((n)^0.5)
    cc[i, ] = cc[i,]/(n-1)/((n)^(1-apower))

    #cc[i,]  = ifelse(n<minModuleSize, 0, cc[i,])
  }

  if(F) {
  # module identification
  cc1d = -as.numeric(cc) # convert a matrix into a vector
  maxIdx     = order(cc1d)
  Mc = cc1d[ maxIdx[1] ]
  # location of the max in the original matrix
  Ii  = maxIdx[1] %% nodes  #Ii as row of max
  if(Ii==0){#last column
    Ii=nodes
  }
  i = as.integer(maxIdx[1]/nodes) #i as column of max
  fi = maxIdx[1]/nodes
  if (fi>i){
      i=i+1
  }
  }

  # MCI[1, ]: max value and MCI[2, ] is index
  #
  mcI = apply(cc, 2, maxValueIndex)
  mc  = mcI[1, ]
  I   = mcI[2, ]
  Mci = maxValueIndex(mc)
  Mc  = Mci[1]
  i   = Mci[2]
  Ii  = I[i]

  maxMc=max(mc)
  cutMc= maxMc/foldchange_cc

  mresult = NULL
  #while(Mc > 0 & module.cnt <maxModules){ #original definition
  while(Mc > cutMc & module.cnt <maxModules){
      
     nres = c(Ii, i, Mc, i-Ii)
     mresult = rbind(mresult, nres)

     mstr = paste("max_cc=", as.character(Mc), " size=", as.character(i-Ii+1),"\n",sep="")
     cat(mstr)

     #assign modules
     module.size = i-Ii+1
     if(module.size >=minModuleSize){
        #assign module lable
        module.assign[c(Ii:i)] = rep(module.cnt, module.size)
        module.cnt = module.cnt + 1
     }

     # reset the nodes already in the module
     cc[1:i, c(Ii:nodes) ] = 0 #original

     #cc[1:i, c(Ii:nodes) ] = 0
     #cc[Ii:i, c(Ii:nodes) ] = 0
     #cc[1:Ii, c(Ii:i) ] = 0

     if(F) {
     cc1d  = -as.numeric(cc)
     maxIdx= order(cc1d)
     Mc    = cc1d[ maxIdx[1] ]
     # location of the max in the original matrix
     Ii  = maxIdx[1] %% nodes  #Ii as row of max
     if(Ii==0){#last column
       Ii=nodes
     }

     i = as.integer(maxIdx[1]/nodes) #i as column of max
     fi = maxIdx[1]/nodes
     if (fi>i){
        i=i+1
     }
     }
 
     # MCI[1, ]: max value and MCI[2, ] is index
     #
     mcI = apply(cc, 2, maxValueIndex)
     mc  = mcI[1, ]
     I   = mcI[2, ]
     Mci = maxValueIndex(mc)
     Mc  = Mci[1]
     i   = Mci[2]
     Ii  = I[i]
  }

  #mcol     = dim(mresult)[2]
  #selected = mresult[,mcol]>=minModuleSize
  #sel.results = mresult[selected, ]

  colcode.reduced.order = assignModuleColor(module.assign, minsize1=minModuleSize-1,  
               anameallmodules=F, auseblackwhite=F, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)

  #re-order to the normal one with sequential signleton index
  recov.order     = order( hierclust$order)
  colcode.reduced = colcode.reduced.order[recov.order]

  rm(cc)
  collect_garbage()

  return (colcode.reduced)

}

maxValueIndex = function(mvect){
   index  = c(1:length(mvect))
   maxval = max(mvect, na.rm=T)
   sel    = maxval==mvect
   sel    = ifelse(is.na(sel), F, sel)
   mindex = index[sel]
   return ( c(maxval, mindex[1]) )
}
#a=cbind(c(1:3), c(5,6,7), c(9,10,8))
#b=apply(a, 1, maxValueIndex)

# a= 1 1 2 2 3 3 3
# findBoundary(a) => [1] 1 3 5
#
findBoundary=function(vect){
   no.elements = length(vect)
   shifted= c(vect[no.elements], vect[1:(no.elements-1)] )
   sel = shifted != vect

   return ( c(1:no.elements)[sel] )
}

#------------------------------------------------------------------------------------------------------------------
#Function: cut hierarchical clusering tree, base don internal structure of ordered dendrogram 
#          use the information of dendrogram height with the mean height of a module to determine the cut points,
#             if this doesn't split  the module, try the differences with (max_height+mean_height)/2
#          the whole process can be iterated until number of clusters become stable (if deepSplit is set as TRUE)
#Input:
#    hierclust           ~  hierarchical clustering object
#    deepSplit           ~  perform iterative searching of sub clusters if set to be TRUE, otherwise, do
#    minModuleSize       ~  min size of module
#    minAttachModuleSize ~  min size of module to be attached, as a major module with larger size
#                                startpoint of current module 
#Output:  module colors for every node
# 
#------------------------------------------------------------------------------------------------------------------

cutreeDynamic = function(hierclust, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50, minAttachModuleSize=100, nameallmodules=FALSE, useblackwhite=FALSE, useNumberAsLabel=F, startlabel=1)
{
    #dim(hierclust$merge)
    #length(hierclust$height)
    #length(hierclust$order)
    #hierclust$merge[1:10, ]
    #round(hierclust$height[1:20], 2)
    #orderHei   = hierclust$height[hierclust$order]
    #orderMerge = hierclust$merge[, ]
    #plclust(hierclust, labels=F, xlab="",ylab="",main="",sub="")

    if(maxTreeHeight >=1){
      staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=0.99, minsize1=minModuleSize)
    }else{
      staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=maxTreeHeight, minsize1=minModuleSize)
    }      

    #get tree height for every singleton
    #node_index   tree_height
    demdroHeiAll= rbind( cbind(hierclust$merge[,1], hierclust$height), cbind(hierclust$merge[,2], hierclust$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(hierclust$order)

    #> demdroHei.sort[1:10,]
    #       [,1]      [,2]
    #[1,] -3000 0.7943783
    #[2,] -2999 0.7863851
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort    = demdroHeiAll.sort[c(1:no.singletons), ]

    #finally, we got tree height for each of the singleton inorder of 1 to no.singletons
    #     [,1]      [,2]
    #[1,]    1 0.8389184
    #[2,]    2 0.8772433
    #[3,]    3 0.8308482
    demdroHei      = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]  = -demdroHei[,1]

    # combine with prelimilary cluster-cutoff results
    demdroHei  = cbind(demdroHei, as.integer(staticCutCluster))

    # re-order the order based on the dendrogram order hierclust$order
    demdroHei.order = demdroHei[hierclust$order, ]

    # get start and end posiiton of every cluster
    # [1,]  173  315
    # [2,]  676  793
    static.clupos = locateCluster(demdroHei.order[, 3])
    static.no     = dim(static.clupos)[1]

    static.clupos2 =     static.clupos
    static.no2     =     static.no

    #split individual cluster if there are sub clusters embedded
    if(F){
        clusterDemdroHei=mydemdroHei.order;
        cminModuleSize       = minModuleSize;
        cminAttachModuleSize = minAttachModuleSize;
    }

    mcycle=1
    while(1==1){
        clupos = NULL
        for (i in c(1:static.no)){
           mydemdroHei.order = demdroHei.order[ c(static.clupos[i,1]:static.clupos[i,2]), ] #index to [1, clusterSize]
           mydemdroHei.order[, 1] = mydemdroHei.order[, 1] - static.clupos[i, 1] + 1

           #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
           #cat("i=", as.character(i), "\n")

           iclupos = processInvididualCluster(mydemdroHei.order, 
                                              cminModuleSize       = minModuleSize, 
                                              cminAttachModuleSize = minAttachModuleSize)

           iclupos[,1] = iclupos[,1] + static.clupos[i, 1] -1 #recover the original index
           iclupos[,2] = iclupos[,2] + static.clupos[i, 1] -1

           clupos  = rbind(clupos, iclupos) #put in the final output buffer
        }

        if(deepSplit==FALSE){
          break
        }

        if(dim(clupos)[1] != static.no) {
           static.clupos = clupos
           static.no     = dim(static.clupos)[1]
        }else{
          break
        }
       mcycle = mcycle + 1
       #static.clupos
       
    }
        
    final.cnt = dim(clupos)[1]

    #assign colors for modules
    module.assign = rep(0, no.singletons)
    module.cnt=1
    for (i in c(1:final.cnt ))
    {
       sdx = clupos[i, 1] #module start point
       edx = clupos[i, 2] #module end point

       module.size = edx - sdx +1 

       if(module.size <minModuleSize){
         next
       }
       #assign module lable
       module.assign[sdx:edx] = rep(module.cnt, module.size)
       #update module label for the next module
       module.cnt = module.cnt + 1
    }

    colcode.reduced.order = assignModuleColor(module.assign, minsize1=minModuleSize,  anameallmodules=nameallmodules, 
                             auseblackwhite=useblackwhite, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)
    #tapply(demdroHei.order[1:no.singletons, 2], colcode.reduced.order, mean)

    #re-order to the normal one with sequential signleton index
    recov.order     = order( demdroHei.order[,1])
    colcode.reduced = colcode.reduced.order[recov.order]

    if(1==2){
        aveheight = averageSequence(demdroHei.order[,2], 2)
        procHei = demdroHei.order[,2]-mean(demdroHei.order[,2])

        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(hierclust, labels=F, xlab="",ylab="",main="",sub="",axes = F)
        #barplot(demdroHei.order[,2]-min(demdroHei.order[,2]), 
        #        col= "black", space=0,
        #        border=F,main="", axes = F, axisnames = F)
        #barplot(aveheight-mean(aveheight), 
        barplot(procHei, 
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(colcode.reduced)), 
                col= as.character(colcode.reduced[hierclust$order]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

    colcode.reduced
}


#input is the cluster demdrogram of an individual cluster, we want to find its embbeded subclusters
#execution order: mean-height ==> (mean+max)/2 ==> (mean+min)/2
#useMean: =0 ~ use mean-height   as calibation line
#         =1 ~ use (mean+max)/2  as calibation line to detect relatively a small cluster sitting on the head of a bigger one,
#                      so mean-height is too low to detect the two modules.
#         =-1~ use (mean+min)/2  as calibation line to detect relatively a small cluster sitting on the tail of a bigger one,
#                      so mean-height & (mean+max)/2 are too high to detect the two modules

#processInvididualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
processInvididualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    no.cnodes = dim(clusterDemdroHei)[1]
    
    cmaxhei   = max(clusterDemdroHei[, 2])
    cminhei   = min(clusterDemdroHei[, 2])
    
    cmeanhei  = mean(clusterDemdroHei[, 2])
    cmidhei = (cmeanhei + cmaxhei)/2.0
    cdwnhei = (cmeanhei + cminhei)/2.0

    if (useMean==1){
        comphei = cmidhei
    }else if (useMean==-1){
        comphei = cdwnhei
    }else{ #normal case
        comphei = cmeanhei
    }
        
    # compute height diffrence with mean height
    heidiff       = clusterDemdroHei[,2] - comphei
    heidiff.shift = shiftSequence(heidiff, -1)

    # get cut positions
    # detect the end point of a cluster, whose height should be less than meanhei 
    #  and the node behind it is the start point of the next cluster which has a height above meanhei
    cuts.bool = (heidiff<0) & (heidiff.shift > 0)
    cuts.bool[1]         = TRUE
    cuts.bool[no.cnodes] = TRUE
    
    if(sum(cuts.bool)==2){
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }else{
             new.clupos = rbind(c(1, no.cnodes))
          }
          return (new.clupos)
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    cutindex =c(1:no.cnodes)[cuts.bool]
    no.cutps = length(cutindex)
    runlens  = rep(999, no.cutps)
    cuts.bool2 = cuts.bool
    for(i in c(2:(no.cutps-1)) ){
       seq = c( (cutindex[i-1]+1):cutindex[i] )
       runlens[i] = runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)

       if( (runlens[i]<minTailRunlength) & (runlens[i]<cminModuleSize) ){
          #cat("run length=", runlens[i], "\n")
          cuts.bool2[ cutindex[i] ] = FALSE
       }
    }

    #attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
    cuts.bool3=cuts.bool2
    if(sum(cuts.bool2) > 3) {
       curj = 2
       while (1==1){
           cutindex2 =c(1:no.cnodes)[cuts.bool2]
           no.clus = length(cutindex2) -1
           if (curj>no.clus){
              break
           }
           pre.sdx = cutindex2[ curj-1 ]+1 #previous module start point
           pre.edx = cutindex2[ curj ] #previous module end   point
           pre.module.size = pre.edx - pre.sdx +1 
           pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
         
           cur.sdx = cutindex2[ curj ]+1 #previous module start point
           cur.edx = cutindex2[ curj+1 ] #previous module end   point
           cur.module.size = cur.edx - cur.sdx +1 
           cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])

           #merge to the leftside major module, don't change the current index "curj"
           #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
           if( (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize) ){
                cuts.bool2[ cutindex2[curj] ] = FALSE
           }else{ #consider next cluster
                curj = curj + 1
           }
       }#while
   }#if 

   cutindex2 =c(1:no.cnodes)[cuts.bool2]
   no.cutps = length(cutindex2)

    #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
    #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
    if(no.cutps > 2){
      if( (cutindex2[no.cutps] - cutindex2[no.cutps-1]+1) < cminModuleSize ){
        cuts.bool2[ cutindex2[no.cutps-1] ] =FALSE  
      }
    }

    if(1==2){
        myseqnce = c(2300:3000)
        cutdisp = ifelse(cuts.bool2==T, "red","grey" )
        #re-order to the normal one with sequential signleton index
        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(hierclust, labels=F, xlab="",ylab="",main="",sub="",axes = F)
        barplot(heidiff[myseqnce],
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(cutdisp[myseqnce])), 
                col= as.character(cutdisp[myseqnce]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

   cutindex2  = c(1:no.cnodes)[cuts.bool2]
   cutindex2[1]=cutindex2[1]-1 #the first 
   no.cutps2  = length(cutindex2)

   if(no.cutps2 > 2){
     new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
   }else{
     new.clupos = cbind( 1, no.cnodes)
   }

   if ( dim(new.clupos)[1] == 1 ){   
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }   
   }
   new.clupos
}


findClustersSignificant=function(mysequence, modulecolor)
{
  
   modnames= names( table(modulecolor) )
   mysize     = length(modulecolor)
   validseq     = rep(TRUE, mysize)
   for (each in modnames ){
      mybool = (modulecolor==each)
      mymodulesig = mysequence[mybool]
      mydiff      = abs(mymodulesig - mymodulesig[1])
      if(sum(mydiff)==0){
         validseq = ifelse(mybool==TRUE, FALSE, validseq)
      }
   }
  validseq
}


#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
    seqlen = length(mysequence)   
    if(leftOrright<0){
        pseq = rev(mysequence)
    }else{
        pseq = mysequence
    }

    if(mysign<0){ #see where the first POSITIVE number occurs
       nonezero.bool = (pseq > 0)
    }else{ #see where the first NEGATIVE number occur
       nonezero.bool = (pseq < 0)
    }
    if( sum(nonezero.bool) > 0){
      runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
    }else{
      runlength = 0
    }
}


#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}


#no of neighbors behind and before the point used for average
averageSequence=function(mysequence, noneighbors){
   sumseq = mysequence
   for(i in c(1:noneighbors)){
        iseq = shiftSequence(mysequence, i)
        sumseq =    sumseq + iseq
        iseq = shiftSequence(mysequence, -i)
        sumseq =    sumseq + iseq
   }
   sumseq = sumseq/(1+2*noneighbors)
}

#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}

# compute frequence in bins for a given vector
# return a matrix with 1st column as the interval middles and 2nd as the frequences
#
freqenceInBins = function(myvect, minval=0, binsize=20000, cisvect=NULL, elementnames=NULL){
    
    maxpos        = max(myvect)
    minpos        = min(myvect)

    no.elements = length(myvect)
    
    if (!is.null(minval) ){
        minpos    = minval
    }

    if(is.null(cisvect) ){
        mycisvect = rep(F, no.elements)
    }else{
        mycisvect = ifelse(cisvect==1, T, F)
        mycisvect = ifelse(is.na(mycisvect), F, mycisvect)
    }

    if(is.null( elementnames) ){
        myelementnames = as.character( c(1:no.elements) )
    } else{
        myelementnames = elementnames
    }
    
    noBins        = as.integer( (maxpos-minpos)/binsize) + 1
    # intervals
    intervals     = c(0:noBins) * binsize
    
    # middle position of intervals
    intervalsMidd = minpos + c(1:noBins) * binsize - binsize/2

    no.intervals  = length(intervalsMidd)

    freqs         = rep(0,  no.intervals)
    binElements   = rep("", no.intervals)

    binCisFreqs   = rep(0,  no.intervals)
    binCisElements= rep("", no.intervals)

    halfbinsize   = binsize/2

    # count number of elements in each bin
    for (i in c(1:no.intervals)) {
        iselL    = myvect>= intervalsMidd[i] - binsize/2
        iselR    = myvect < intervalsMidd[i] + binsize/2
        isel     = iselL & iselR
        freqs[i] = sum(isel)

        if (freqs[i] == 0){
           next
        }

        # get element names in the bin
        binElements[i] = concatenate(elementnames[isel],'; ')
        
        # look at the cisQTL genes        
        icisSel           = mycisvect & isel
        binCisFreqs[i]    = sum( icisSel )
        if (sum(icisSel) >0) {
            binCisElements[i] = concatenate(elementnames[icisSel],'; ')
        }

    }
    
    retNumerical = cbind(intervalsMidd, freqs)
    retCharacter = cbind(binElements, binCisFreqs, binCisElements)

    colnames(retNumerical) <- c("interval(middle)", "QTL count") 
    colnames(retCharacter) <- c("QTL genes", "cis-QTL count", "cis-QTL genes") 

    return ( list(retNumerical,retCharacter) )
}


#find the middle of each cluster and label the middle position with the corrsponding color
getDisplayColorSequence=function(colordered){
 mylen = length(colordered)
 colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
 colordiff   = (colordered != colordered2)
 colordiff[1] = TRUE
 colordiff[mylen] = TRUE
 #mydispcolor = ifelse(colordiff==TRUE, colordered, "")
 mydispcolor = rep("", mylen)
 mytrueseq = c(1:mylen)[colordiff]
 for (i in c(1:(length(mytrueseq)-1)) ){
    midi =  (mytrueseq[i] + mytrueseq[i+1])/2
    mydispcolor[midi] = colordered[midi]
 }
 fdispcolor = ifelse(mydispcolor=="grey", "", mydispcolor)
 fdispcolor
}

# find the start and end indice of each color
#
getDisplayColorPos=function(colordered){

 mylen       = length(colordered)
 myidx       = c(1:mylen)
 colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
 colordiff   = (colordered != colordered2)
 colordiff[1]= TRUE

 xColor       = colordered[colordiff]
 xcolStartIdx = myidx[colordiff]
 xcolEndIdx   = c(xcolStartIdx[-1]-1, mylen)

 # excluding grey
 selNoGrey = xColor != "grey" 

 zColor = xColor[selNoGrey]
 zStart = xcolStartIdx[selNoGrey]
 zEnd   = xcolEndIdx[selNoGrey]

 zcolpos= data.frame(rbind(zStart, zEnd ))
 colnames(zcolpos) <- zColor

 zcolpos
}




#use height cutoff to remove
cutTreeStatic = function(hiercluster,heightcutoff=0.99, minsize1=50) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree2(hiercluster,h=heightcutoff)
    sort1=-sort(-table(labelpred))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 > minsize1
    no.modules=sum(modulebranch)

    colorhelp = rep(-1, length(labelpred) )
    if ( no.modules==0){
        print("No mudule detected\n")
    }
    else{
        for (i in c(1:no.modules)) {
            colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
        }
    }
    colorhelp
}


cutree2 = function (tree, k = NULL, h = NULL)
{

    if (is.null(n1 <- nrow(tree$merge)) || n1 < 1)
        stop("invalid 'tree' (merge component)")
    n <- n1 + 1
    if (is.null(k) && is.null(h))
        stop("either 'k' or 'h' must be specified")
    if (is.null(k)) {
        ### tree$height must be sorted
        temp= tree$height-c(0, tree$height[-length(tree$height)])

        sorted=T

        #if(sum(sign(temp[-1])<0)>0) sorted=F

        if (sorted==F)
            stop("the 'height' component of 'tree' is not sorted\n(increasingly); consider applying as.hclust() first")

        k <- integer(length(h))

        k <- n + 1 - apply(outer(c(tree$height, Inf), h, ">"), 2, which.max)

        if (getOption("verbose"))
            cat("cutree(): k(h) = ", k, "\n")

    } else {
        k <- as.integer(k)
        if (min(k) < 1 || max(k) > n)
            stop(gettextf("elements of 'k' must be between 1 and %d",n), domain = NA)
    }

    ans <- .Call("R_cutree", tree$merge, k, PACKAGE = "stats")

    if (length(k) == 1) {
        ans <- as.vector(ans)
        names(ans) <- tree$labels
    } else {
        colnames(ans) <- if (!is.null(h)) h else k

        rownames(ans) <- tree$labels
    }

    return(ans)

}


#locate the start/end positions of each cluster in the ordered cluster label sequence 
#where "-1" indicating no cluster
#3-1 -1 1 1 1 1 2 2 2
#3 3 -1-1 1 1 1 1 2 2 2  (shift)
#---------------------------------
#0-4  0 2 0 0 0 1 0 0 0   (difference)
#       *     * @
locateCluster = function(clusterlabels)
{
 no.nodes = length(clusterlabels)
 clusterlabels.shift = c(clusterlabels[1], c(clusterlabels[1:(no.nodes-1)]) )
 
 #a non-zero point is the start point of a cluster and it previous point is the end point of the previous cluster
 label.diff = abs(clusterlabels - clusterlabels.shift)

 #process the first and last positions as start/end points if they belong to a cluster instead of no cluster "-1"
 if(clusterlabels[1]       >0) {label.diff[1]=1} 
 if(clusterlabels[no.nodes]>0) {label.diff[no.nodes]=1} 

 flagpoints.bool = label.diff > 0
 flagpoints = c(1:no.nodes)[flagpoints.bool]
 no.points  = length(flagpoints)

 myclupos=NULL
 for(i in c(1:(no.points-1)) ){
   idx = flagpoints[i]
   if(clusterlabels[idx]>0){
      if(flagpoints[i+1]==no.nodes) {#boundary effect
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]) )
         break
      }else{
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]-1) )
      }
   }
 }
 myclupos
}


#row-wise reorder a matrix based on "orderedList", default key column in the
#matrix is the first column
orderMergedMatrix = function(disorderMatrix, orderedList, keyCol=1){
  no.samples = length(orderedList)
  cnt = 1
  seqc  =c(1:no.samples)
  rightOrder = rep(0, no.samples)
  orderedMatrix = NULL
  for(each in orderedList){
    whichTrue = ( as.character(each)==as.character(disorderMatrix[, keyCol]) )
    idx = sum(whichTrue * seqc)
    #cat(as.character(each), " : ", as.character(disorderMatrix[idx,keyCol]),as.character(idx), "\n" )
    rightOrder[idx] = cnt
    cnt = cnt + 1
    orderedMatrix = rbind(orderedMatrix, disorderMatrix[idx,])
  }
  colnames(orderedMatrix) <- colnames(disorderMatrix)
  orderedMatrix
}

#-------------------------------------------------------------------------
#Function: convert the upper
#Input:   
#         mvector     ~ adjacency vector with first column is 
#                       the ROW index (i) of the vector in the original matrix A=[a(i,j)]
#         mcutoff     ~ color codes (modules) for genes
#         mgraphfname ~ file for storing output tripples of (i, j, a(i,j)), j>i
#         mfudgefactor~ fudge factor to spread (<1) or compress (>1) a network
#Output:  tripples of (i, j, a(i,j)), j>i
#-------------------------------------------------------------------------
vectorToPairs=function(mvector, mcutoff, mgraphfname,mfudgefactor=1)
{
#mvector = kdatEdge[12, ]
#mcutoff =cutoff
#mgraphfname = graphfname
   i = as.integer(mvector[1])

   #actual data vector
   datv = as.numeric(as.vector(mvector[-c(1)]))
   mno.nodes = length(datv)

   #node index 
   index <- c((i+1):mno.nodes)   
   av = datv[index] 
   bv = av > mcutoff

   no.selected = sum(bv)

   ret=c(0,0,0)
   if(no.selected>0){
     c1=rep(i, no.selected)
     c2=index[bv]
     c3=as.numeric(av[bv])^mfudgefactor
     ret = as.matrix(cbind(c1,c2,c3))
     write.table(ret, mgraphfname, append=T, sep=" ", quote=F, col.names=F, row.names=F)
   }
   #cat(i," :",as.character(no.selected),  "\n")  
   ret
}

mergeClusterByPCA=function(mdendro, genecluster, pccluster, pcnames){

    geneclusize    = table(colcode.reduced)
    finalgenecolor = as.character(genecluster)
    pcclunames = names(table(pccluster))

    # get ordered gene cluster assignment
    orderedCluRaw = as.character(genecluster[mdendro$order])
    fdispcolor = getDisplayColorSequence( orderedCluRaw)
    orderedgeneclus= fdispcolor[fdispcolor!=""]
    

    # get the start and end indices for each module excluding the grey module
    #
    geneModulePos = getDisplayColorPos( orderedCluRaw)

    #we need include grey module, but use different position to differentiate it from the ordinary clusters
    # so that the recombination will not consider the greay module
    orderedgeneclus = c(orderedgeneclus, "grey")

    geneclupos     = c(1:(length(orderedgeneclus)-1), 99999)
    names(geneclupos) <- orderedgeneclus

    # get ordered cluster-sizes
    orderedGeneclusize  = NULL
    for (each in orderedgeneclus){
       orderedGeneclusize  = c(orderedGeneclusize, geneclusize[[each]])
    }
    names(orderedGeneclusize)  <- orderedgeneclus

    for (each in pcclunames ){
       if (each =="grey"){
         next
       }
       
       # PCs belongs to the current cluster 
       sel.pcs      = as.character(pccluster) == each
       if (sum(sel.pcs)<=1){
          next
       }
       sel.pcnames  = pcnames[sel.pcs]
       no.selpcs    = length(sel.pcnames) 

       #get the PCs' sequential positions in the original dendrogram
       sel.pcpos    = NULL
       for (eachmerged in sel.pcnames){
           sel.pcpos  = c(sel.pcpos, geneclupos[[eachmerged]])
       }
       sel.order = order(sel.pcpos)

       # get the neighboring segments of the PCs, so we only merge the neighboring PCs
       #> order.selPCpos: 30 31 32 36 37 38 39 40
       #> shift.selPCpos: 29 30 31 32 36 37 38 39
       #>  diff.selPCpos:  1  1  1  4  1  1  1  1

       order.selPCnames= sel.pcnames[sel.order]
       order.selPCpos  = sel.pcpos[sel.order]
       shift.selPCpos  = c(order.selPCpos[1]-2, order.selPCpos[c(1:(no.selpcs-1))])
       diff.selPCpos   = order.selPCpos - shift.selPCpos
       bool.startpos   = (diff.selPCpos !=1)   
       startpos        = c(1:no.selpcs)[bool.startpos]
       endpos          = c(startpos[-1]-1, no.selpcs)
       nosegments      = length(endpos)
    
       #we choose the cluster with maximal size as the module assigment for this SEGMENT of clusters
       for (iseg in c(1:nosegments) ){
           seg = c(startpos[iseg]:endpos[iseg])
           seg.pcnames = order.selPCnames[seg]
           if (length(seg.pcnames)==1){
               next
           }
           mergedcolor= getNameOfMaxElement(orderedGeneclusize, seg.pcnames)
           mergedcolorPos = geneModulePos[mergedcolor]

           cat("merged color=", mergedcolor, "\n")
           for (eachmerged in seg.pcnames){
              epos = geneModulePos[eachmerged]
              # so modules to be merged has to be adjacent to the major module
              if(abs(mergedcolorPos[1,1]-epos[2,1])==1 | abs(mergedcolorPos[2,1]-epos[1,1])==1  ){
                  cat("   ", eachmerged, "\n")
                  finalgenecolor = ifelse(finalgenecolor==eachmerged, mergedcolor, finalgenecolor)                  
              }
           }
       }

    }
   
    # now we need sort the module assignment with new names by using an existing function
    # to do so, we need convert the current colors into numerals which are the input of this 
    # existing function
    retcolors = reassignModuleNames(finalgenecolor)

    retcolors
}

# now we need sort the module assignment with new names by using an existing function
# to do so, we need convert the current colors into numerals which are the input of this 
# existing function
reassignModuleNames=function(mycolorvect, minmodulesize=0){
    fgenecolor = as.character(mycolorvect)
    ztable = table(fgenecolor)
    zfinal = rep(0, length(fgenecolor))
    iclu   = 1 
    for (each in names(ztable)){
       if (each=="grey")
          next
       if (ztable[[each]] < minmodulesize)
          next

       zfinal = ifelse(fgenecolor==each, iclu, zfinal)
       iclu=iclu+1
    }

    retcolors=assignModuleColor(labelpred=zfinal, minsize1=1, anameallmodules=T, auseblackwhite=FALSE)
    retcolors
}


#orderMergedMatrix = function(disorderMatrix, orderedList, keyCol=1){
#  orderM = order(disorderMatrix[keyCol, ])
#  orderL = order(orderedList)
#  orderL2=order(orderL)  
#  orderM2= orderM[orderL2]
#  disorderMatrix[orderM2, ]
#}

#-------------------------------------------------------------------------
#Function: compute whithin-module cluster coefficient for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module cluster coefficients
#-------------------------------------------------------------------------
computeModuleCC = function(adjMatrix, colorcodeC, weighted=T)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   clustercoeff = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if(as.character(each)=="grey" ){
         next
      }
      whichmod = each==colorcodeC
      if( sum(whichmod)>1){
         icc = computeClusterCoefficient(adjMatrix[whichmod,whichmod], weighted)
      }else{
         icc = 1
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the cc's to the
      clustercoeff[idxmod] = icc
   }
   clustercoeff
}

#------------------------------------------------------------------------
#Function: compute whithin-module conformity for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module conformity
#-------------------------------------------------------------------------
computeModuleConformity = function(adjMatrix, colorcodeC)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   conformity = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if(as.character(each)=="grey" ){
         next
      }
      whichmod = each==colorcodeC

      module.size = sum(whichmod)
      if (module.size==1){
        next
      }

      iconform = SDADJ1(adjMatrix[whichmod,whichmod])     

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the conformity to the global list
      conformity[idxmod] = iconform
   }
   conformity
}


#-------------------------------------------------------------------------
#Function: compute total number of connections for each gene excluding 
#          genes in the grey module
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#
#Output:  the number of connections of each gene to the gens 
#          in non-grey modules
#-------------------------------------------------------------------------
computeTotalLinksToNongreygenes = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
   modnames= names( table(colorcodeC) )
   no.genes     = dim(adjMatrix)[1]
   links        = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   nongrey      = colorcodeC != "grey"

   total.nonegrey = sum(nongrey)

   totallinks <- apply(adjMatrix[, nongrey], 1,sum, na.rm=TRUE)
   if(!isAdjacency){
      totallinks <- (total.nonegrey - totallinks)
   }

   #normalize against the module size
   if(normalized==TRUE){
      totallinks = totallinks /total.nonegrey
   }      
    
  #put the links's to the buffer      
  #links[nongrey] = totallinks
  
  totallinks
}


#-------------------------------------------------------------------------
#Function: compute whithin-module number of connections for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module number of connections of each gene
#-------------------------------------------------------------------------
computeModuleLinks = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   links        = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if((usegreymodule==F) & (as.character(each)=="grey" ) ){
         next
      }

      whichmod    = each==colorcodeC
      module.size = sum(whichmod)

      if (module.size==1){
        next
      }

      modk <- apply(adjMatrix[whichmod,whichmod],2,sum, na.rm=TRUE) 
      if(!isAdjacency){
          modk <- (module.size -modk)
      }

      #normalize against the module size
      if(normalized==TRUE){
         modk = modk/module.size
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the links's to the buffer      
      links[idxmod] = modk
   }
   links
}


#---------------------------------------------------------------------------------
#Function: compute whithin-module number of connections for each row-based element
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  a matrix of cross-module numbers of connections of each of row-based elements
#         columns as modules and rows as elements of interest
#
#Notice that here it is different from previous function, here adjMatrix may be
#     asymmetric, columns are color-coded
#
#----------------------------------------------------------------------------------
computeMultiModuleConnectivity = function(adjMatrix, colorcodeC)
{
   modnames   = names( table(colorcodeC) )
   kin.matrix = NULL
   for (each in modnames ){
      whichmod    = each==colorcodeC
      if( sum(whichmod) >1){
         modk <- apply(adjMatrix[, whichmod], 1, sum, na.rm=TRUE)
      }else{
         modk <- as.integer(adjMatrix[, whichmod])
      }
      kin.matrix = cbind(kin.matrix, modk)
   }
   colnames(kin.matrix) <- modnames
   kin.matrix
}

#process markers location into a sequence in order of chromsomes 1 to 21
processChromPosSequence = function(avector, colorcodeC)
{
  #initilaization
  no.elements= length(avector)
  outpos = avector
  idxseq =c(1:no.elements)

   modnames= names( table(colorcodeC) )
   currMax=0
   for (each in modnames ){      
      whichmod = each==colorcodeC
      #indices of elements in the current module
      idxmod = idxseq * whichmod
      outpos[idxmod] = outpos[idxmod] + currMax 
      
      #last marker pos in the current Xmosome
      imax   = max( avector[whichmod] )
      currMax=  currMax + imax + 1
   }
   
   #tried to find the Xmosoome boundaries
   shifted = c(0, colorcodeC[1:(no.elements-2)], 0)
   boundary= (shifted != colorcodeC)
   boundary.colors=ifelse(boundary==TRUE, colorcodeC, "")
   boundary.colors[ length(boundary.colors) ] = "."

   list(outpos, boundary.colors)
}


#perform ttest for two groups (expression values across samples)
#two group are represented by two boolean vectors with the same size as exprvector
#
# if usePseudoSD=T, then, for the group with 1 sample, we use an addiitonal pseudo sample
#  to ensure pvalue is computable, in this way we can know that the standard deviation of 
#   the group with more than one sample
# 
#
#ttestVector (exprvector=c(1:10), tgroupA=c(T, rep(F,9)), tgroupB=c(rep(F,8), T,T), usePseudoSD=F)
ttestVector = function(exprvector, tgroupA, tgroupB, usePseudoSD=F){

    #all the elements of one group are NA
    if( (sum(is.na(exprvector[tgroupA]))==sum(tgroupA)) ||  (sum(is.na(exprvector[tgroupB]))==sum(tgroupB)) )
      return (1)

    #too few elements in one group to be used to perform ttest
    if (!usePseudoSD) {
       if( (sum(!is.na(exprvector[tgroupA]))<2) ||  (sum(!is.na(exprvector[tgroupB]))<2) )
           return (1)
    }

    myexprA=exprvector[tgroupA]
    myexprB=exprvector[tgroupB]

    delta=5
    #too small group-size
    if ( (sum(tgroupA)<2)  || (sum(tgroupB)<2) ){

       if(!usePseudoSD) #no need for further computation of pvalue
           return (1)

       if ( sum(tgroupA)==1)
          myexprA =c(myexprA-delta, myexprA+delta)
       
       if ( sum(tgroupB)==1)
          myexprB =c(myexprB-delta, myexprB+delta)
    }
    
    tt    = t.test(myexprA,myexprB)
    tt$p.value
}



coxSurvModel = function(exprvector, survector){
    exprvector.num = as.numeric(exprvector)
    mpvalue=1
    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num) )

    if( (var(exprvector.num) >1) & (length(levels) >2) ) {
      cox1=coxph(survector ~ exprvector.num, na.action="na.omit")
      mpvalue=1-pchisq(cox1$score, 1)
    }
    mpvalue
}

coxSurvModelComplex = function(exprvector, survector, minvar=1){
    exprvector.num = as.numeric(exprvector)
    mpvalue = 1
    sde     =-1
    lower.95=-1
    upper.95=-1
    hr      =-1

    #print(exprvector)

    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num) )

    if( (var(exprvector.num, na.rm =T) >=minvar) & (length(levels) >2) ) {

       cox1=coxph(survector ~ exprvector.num, na.action="na.omit")
       mpvalue=1-pchisq(cox1$score, 1)

       groups=I(exprvector.num > median(exprvector.num, na.rm =T) )

       cox2 = coxph(survector ~ groups)
       sde= sqrt(cox2$var)

       lower.95=exp(cox2$coef-1.96*sde)
       upper.95=exp(cox2$coef+1.96*sde)
       hr      =exp(cox2$coef)
    }

    outphr = c(mpvalue, hr, lower.95, upper.95, sde)
}


computeTOM= function(adjMatrix) {
   no.singletons <- dim(adjMatrix)[1]
   nolinks.reduced  <- apply(adjMatrix, 2, sum, na.rm=TRUE)
   #Let's calculate topological overlap matrix
   numTOM= adjMatrix %*% adjMatrix + adjMatrix
   dist1 = matrix(NA, no.singletons, no.singletons)
   diag(dist1) <- 0
   for (i in 1:(no.singletons-1) ){
      for (j in (i+1):no.singletons){
          denomTOMij = min(nolinks.reduced[i], nolinks.reduced[j]) + 1 - adjMatrix[i,j]          
          dist1[i,j] = 1- numTOM[i,j]/denomTOMij
          dist1[j,i] = dist1[i,j]
      }
   }
  dist1
}

computeLinksInNeighbors <- function(x, imatrix)
{
  y= x %*% imatrix %*% x
  y
}

computeClusterCoefficient = function(adjMatrix, weighted=F) {

        no.genes <- dim(adjMatrix)[1]
        nolinksNeighbors <- c(rep(-666,no.genes))
        total.edge <- c(rep(-666,no.genes))

        #for (i in 1:no.genes){
        #     nolinksNeighbors[i] <-  adjMatrix[i,] %*% adjMatrix %*% adjMatrix[,i]
        #     #total.edge[i] <-  adjMatrix[i,] %*% Pmax %*% adjMatrix[,i]
        #}
        nolinksNeighbors <- apply(adjMatrix, 1, computeLinksInNeighbors, imatrix=adjMatrix)

        plainsum  <- apply(adjMatrix, 1, sum)
        if(weighted) {
           squaresum <- apply(adjMatrix^2, 1, sum)
           total.edge = plainsum^2 - squaresum
        }else{ # for unweighted network, 1^2 = 1
           total.edge = plainsum^2 - plainsum
        }

        # in case of single node, this will not affect the CC computation
        #
        total.edge = ifelse(total.edge==0, 1, total.edge)

        cluster.coef = nolinksNeighbors/total.edge
        cluster.coef = ifelse(total.edge>0,cluster.coef, 0) 

        cluster.coef
}

pajekColorcode = function(bincolorcode)
{
   colorcodeC=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow","grey")
   colorcodeN=c("28",        "4",   "14",   "1",     "21",   "3",  "13",   "5",   "17",     "8",     "24",         "30", "22",   "0",    "18",           "31",       "33",     "15",         "16",         "38" )

   rcolors=bincolorcode
   clevels = length(colorcodeC)
   for (i in c(1:clevels) ){
      whichcolor = rcolors==colorcodeC[i]
      rcolors = ifelse(whichcolor, colorcodeN[i], rcolors)
   }
   rcolors
}

statistics_edgelinks = function(netmatrix, filename, keyword="", usedirect=F)
{
  no.links = dim(netmatrix)[1]
  mnetmatrix = as.matrix(netmatrix)
  if (usedirect){
    nodenames = c(as.character(mnetmatrix[,1]), as.character(mnetmatrix[,2]))
  }else{
    #nodenames = as.character(mnetmatrix[,1])
    nodenames = c(as.character(mnetmatrix[,1]), as.character(mnetmatrix[,2]))
  }

  # number of links for each unique gene
  ntab = table(nodenames)
  
  #no of unique genes
  no.uniquegenes = length(ntab)

  # average links per gene
  avglinks=  no.links/no.uniquegenes
 
  # output table to a file
  tabfname = paste(filename, "_nolinks.xls",sep="")
  ordertab = order(-ntab )
  write.table(ntab[ordertab], tabfname, sep="\t",quote=FALSE, col.names=F, row.names=T)
  
  # output unique nodeIDs
  idfname = paste(filename, "_nodeIDs.txt",sep="")
  write.table(as.matrix(names(ntab)), idfname, sep="\t",quote=FALSE, col.names=F, row.names=F)

  # draw histogram, compute frequency of no of links
  ntabtab = table(as.numeric(ntab))

  linkfreq = as.numeric(ntabtab)
  linkidx  = as.numeric( names(ntabtab) )

  maxlinks = max(linkidx)
  mylabel = paste(keyword, "(max no of links=", as.character(maxlinks), ")",sep="" )

  outimg = paste(filename, "_nolinksHistogram.png",sep="")
  openImgDev(outimg, iwidth = 600, iheight = 600)
  barplot(linkfreq, names.arg= as.character(linkidx), 
     xlab="number of links of a node", ylab="frequency", main= mylabel)

  #plot(linkidx, linkfreq,xlab="number of links of a node", ylab="frequency", main= keyword, type="h")
  #histogram(as.numeric(ntab), br=max(linkidx)-1, 
  #         xlab="number of links of a node", ylab="frequency", main= keyword)
  dev.off()

  # return
  c(no.links, no.uniquegenes, avglinks)
}




panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * r)
}

cor1=function(x) {
 if(dim(x)[1]==1){
    omatrix = matrix(0, dim(x)[2], dim(x)[2] )
    colnames(omatrix) <- colnames(x)
    rownames(omatrix) <- colnames(x)
    return (omatrix)
 }

 col=cor(x,use="p", method="pearson")
 col=ifelse(is.na(col), 0, col)
 signif(col,2)
}

corSpearman=function(x) {
 #print(dim(x))
 if(dim(x)[1]==1){
    omatrix = matrix(0, dim(x)[2], dim(x)[2])
    colnames(omatrix) <- colnames(x)
    rownames(omatrix) <- colnames(x)
    return (omatrix)
 }

 col=cor(x,use="p", method ="spearman")
 col=ifelse(is.na(col), 0, col)
 signif(col,2)
}


# this function computes the standard error
stderror <- function(x){ sqrt( var(x)/length(x) ) }

# Error bars for barplot
# written by: Uli Flenker
# Institute of Biochemistry
# German Sports University
# Cologne Carl-Diem-Weg 6
# 50933 Cologne / Germany
# Phone 0049/0221/4982-493 
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")
 }
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") 
    }
 }
 MW<-0.25*(max(xval)/length(xval)) 
 ERR1<-daten+error 
 ERR2<-daten-error
 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
} 

#works for binary outcomes
krusktest=function( x ) {
   len1=dim(x)[[2]]-1
   out1=rep(666, len1);

   #if outcomes are the same value, no need to do the test
   totalOnes=sum(x[,1])
   if (totalOnes==0 | totalOnes == dim(x)[[1]]){
      for (i in c(1:len1) ) {out1[i]= 1.0}
   } else{
      for (i in c(1:len1) ) {out1[i]= signif( kruskal.test(x[,i+1], x[,1] )$p.value ,2) }
   }

   data.frame( variable=names(x)[-1] , kruskP=out1)
}

# access the pvalues of correlation between continuous vectors
#
CORtest=function( x ) {
   len1=dim(x)[[2]]-1
   no.observations= dim(x)[[1]]
   out1=rep(666, len1);
   for (i in c(1:len1) ) {
       if(no.observations>2){
         out1[i]= signif( cor.test(x[,i+1], x[,1] ,method="p",use="p")$p.value ,2) 
       }else{
         out1[i]= 1
       }
   }
   data.frame( variable=names(x)[-1] , cor.Pvalue=out1)
}

xcortest<-function(x,y){
   res<-cor.test(x, y)
   return( c(res$p.value, res$estimate) )
}
cor_multitraits <- function(x, ys){
    no.ys = dim(ys)[2]
    mres  = NULL
    for (j in c(1:no.ys)){
        jres = xcortest(x, ys[,j])
        mres = c(mres, jres)
    }
    return (mres)
}


# pvalue = 0.00021 ==>"2.1e-4"
pvalue2scientificrep=function(mypvalue){
   if(mypvalue==0){
      return("<e-22")
   }
   if(mypvalue>=1 | mypvalue<0){
      return(as.character(mypvalue))
   }

   base = 10
   cnt  = 0
   while(T){
     cnt     = cnt + 1
     mybase  = base^cnt
     foldval = mypvalue*mybase
     if(foldval>1.0){
        retvalue = paste(as.character(foldval), "e-", as.character(cnt),sep="")
        return (retvalue)      
     }
  }  
}

simpleCorTest=function(x,y){
 signif( cor.test(x,y,method="p",use="p")$p.value ,2) 
}

# no of rows of amatrix is the same as the length of myvect
corTest4multivects=function(myvect, amatrix){
 pvals = apply(amatrix, 2, simpleCorTest, y=myvect)
 #cat(pvals[], "\n")
 as.numeric(pvals)
}

# get index of elements in a given vector
#
getMatchedIndex=function(cvector, subvect){
  subindex=NULL
  fullindex = c(1:length(cvector) )
  for (each in subvect){
      idx = fullindex[cvector== each]
      subindex=c(subindex, idx)
  }

  return (subindex)
}


# to split "abc|123", use sep="\\|"
splitString =function(mystring, separator="; "){
  splitted = NULL
  for (each in mystring){
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a)
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}


concatenate=function(myvect, mysep="")
{
  noitems = length(myvect)
  if (noitems==0){
    return ("")
  }else if (noitems==1){
    return (as.character(myvect) )
  }

  tmpfn = "tmp.txt"
  write.table(t(as.character(myvect)),tmpfn,sep=mysep,quote=FALSE, col.names=F, row.names=FALSE)
  concatenated <- read.delim(tmpfn, sep="!", header=F)
  return (as.character(as.matrix(concatenated) ))
}


# list files with certain keywords embedded inside
#
dirfiles=function(path, pattern){
   mylist=dir(path, pattern=mypatt)
   matches=NULL
   for (each in mylist){
       splitted = splitString(each, pattern)
       if (length(splitted) > 1){
           matches = c(matches, each)
       }
   }

   return (matches)
}


#get the filename without extension
#
getFileExtension=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    
    if( length(splitted) >1){
      return (splitted[length(splitted)])
    } else{
      return ("")
    }
}

#get the filename without extension
getFileName=function(fullfname){
    ext=getFileExtension(fullfname)
    if(ext ==""){
       return (fullfname)
    }
    extd = paste(".", ext, sep="")
    splitted=splitString(fullfname, extd)

    splitted[1]
}

#get the filename without extension
getFileNameOld=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    mfname=""
    for ( i in c(1:(length(splitted)-1)) ){
       mfname =paste(mfname, splitted[i], sep="")
    }
    mfname
}



#get second part: 31357-31351 ==> 31351 
getSecondPart=function(fullfname, sep="-", whichpart=-1){
    splitted=unlist( strsplit(fullfname, sep) )
    if (whichpart==-1)
       ret=splitted[ length(splitted) ]
    else
       ret=splitted[whichpart]

    ret
}



#get the filename without extension
getOntologyNameFromPath=function(fullfname){
    fn=getFileName(fullfname)
    splitted=unlist( strsplit(fn, "\\_") )
    splitted[length(splitted) ]
}


#get the filename without extension and path information
getFileNameNopath=function(fullfname){
    myfilename = getFileName(fullfname)
    splitted=unlist( strsplit(myfilename, "/") )
    splitted[length(splitted) ]
}

# get a particulr field of splitted string
getAFieldBySplit=function(strvect,separator="_", fieldId=1)
{
  splitted = NULL
  for (each in strvect){
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a[fieldId])
  }
  return(splitted )
}

patchNchars=function(intStr, patched="0", number=3)
{
   zeros     = rep(patched, number)

   digits    = unlist(strsplit(intStr,''))# get single digits for the given number
   ndigits   = length(digits)
   
   # put back digits in the zeros
   for (i in c(1:ndigits) ) {
      zidx = number-i+1
      zeros[zidx]=digits[ndigits-i+1]
   }
   finalstr=""
   for (each in zeros){
      finalstr=paste(finalstr, each, sep="")
   }
   return(finalstr)
}

replaceChars=function(fullfname, oldchar, newchar){
    splitted=strsplit(fullfname,'')
    i= length( splitted[[1]] )
    for (j in c(1:i) ){
	if (splitted[[1]][j] == oldchar){
           splitted[[1]][j] = newchar
        }
    }
    mfname =''
    for (j in c(1:i) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}


#get the filename without extension
reverseString=function(mystring){
    splitted=strsplit(mystring,'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}

#get the filename without extension
DNApair=c("A","C","G","T")
names(DNApair) <- c("C","A","T","G")
complementDNAString=function(mystring, mappingTable=DNApair){
    splitted=strsplit(mystring,'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) ){
       jchar = splitted[[1]][j]
       cchar = as.character(mappingTable[jchar]) #get complement character
       mfname = paste(mfname, cchar,sep='')
    }
    mfname
}


#get the filename without extension
reverseString4apply=function(mystringvect){
    splitted=strsplit(as.character(mystringvect[1]),'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}


appendStringToFile=function(fname, mstring, newline=T){
    fp <- file(fname, "a")
    if(newline){
     cat(mstring, "\n", file=fp)
    }else{
     cat(mstring, file=fp)
    }
    close(fp)    
}

appendListToFile=function(fname, mlist, listTitle=""){
  fp <- file(fname, "a")
  if (length(listTitle)>0){
       cat(as.character(listTitle), "\n", file=fp)
  }  
      no.fields=length(mlist)
      #write column title
      for (z in 1:no.fields ){
         cat(as.character(names(mlist[z])),"\t", file=fp)
      }      
      cat("\n", file=fp)

      for (z in 1:no.fields ){
         itab=mlist[z]
         cat(as.character(itab[[1]]),"\t", file=fp)
      }      
      cat("\n", file=fp)
      close(fp)
 }


appendTableToFile=function(fname, mtable, tableTitle="", myappend=T){
   if ( is.null(mtable) ){
     return
   }

   if(myappend==T){
     fp <- file(fname, "a")    
   }else{
     fp <- file(fname, "w")
   }
   if (length(tableTitle)>0){
        cat(as.character(tableTitle), "\n", file=fp)    
   }
  if ( (!is.na( dim(mtable)[1])) & is.na( dim(mtable)[2]) ) {#only one row in the table
      #write column title
      coltitles=names(mtable)
      for (z in 1:length(coltitles) ){
         cat(as.character(coltitles[z]),"\t", file=fp)
      }
      cat("\n", file=fp)

      for (i in 1:(dim(mtable)[1]) ){
          cat(as.character(mtable[i]), "\t", file=fp)
      }
      cat("\n", file=fp)

   }else{ # normal table
       cat(" \t", file=fp)
       #write column title
       coltitles=colnames(mtable)
       for (z in 1:length(coltitles) ){
          cat(as.character(coltitles[z]),"\t", file=fp)
       }
       cat("\n", file=fp)

       rowsname = rownames(mtable)
       for (i in 1:(dim(mtable)[1]) ){
           cat(as.character(rowsname[i]), "\t", file=fp)
           for(j in 1:(dim(mtable)[2])){
              cat(as.character(mtable[i, j]), "\t", file=fp)
            }
           cat("\n", file=fp)
       }   
   }
   cat("\n", file=fp)
   close(fp)
}

appendMultiTablesToFile=function(fname, multitables){
    #table of tables
    if(is.na( dim(multitables)[2]) ) {
        titles=names(multitables)
        for (i in 1:(dim(multitables)[1]) ){
            if ( is.null(multitables[[i]]) )
                 next
            appendTableToFile(fname, multitables[[i]], as.character(titles[i]))
         }
    }else{#single table
      appendTableToFile(fname, multitables)
    }
}

openImgDev=function(imgname, iwidth = 1024, iheight = 1024, ipointsize = 12)
{
  imgtype = getFileExtension(imgname)
  
  if (imgtype=="ps"){
     postscript(imgname,width=iwidth, height=iheight, pointsize=ipointsize)
  }else if (imgtype=="png"){
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize)
  }else if (imgtype=="jpg" || imgtype=="jpeg"){
     jpeg(imgname, width=iwidth, height=iheight, pointsize=ipointsize,quality =100)
  }else if (imgtype=="pdf" || imgtype=="PDF"){
     pdf(imgname)
     return   
  }else{
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize)
  }
  trellis.device(new = FALSE, col = TRUE) 
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: compute sigmoid functions, single/vector inputs
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#compute sigmoid function
sigmoid <- function(x,alpha,tau)
{
  y=1.0+exp( -alpha*(x-tau) )
  y=1.0/y
  y
}

#compute sigmoid values for a matrix of inputs
sigmoidMatrix <- function(xmatrix,alpha,tau0)
{
 ey = exp( -alpha*(xmatrix - tau0) )
 #ey = 3^( -alpha*(xmatrix - tau0) )
 y= 1.0/(1.0+ey)
 rm(ey)
 y
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation cutoff Computation 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations

#iteratively call gc() untill there is no change in memory hit
#i.e., all freed memory has been released to the system
#Usage: 1. immediately call this function after you call a function or
#       2. rm()
collect_garbage=function(){
    while (gc()[2,4] != gc()[2,4]){}
}

dichotcut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, 
mincutval=0.2, maxcutval=0.95, cutstep=0.01) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        # perform sampling 	
        subsetsize=min(samplingSize, orgSize)
        subset1=sample(c(1:orgSize),subsetsize,replace=F)
        cor1 <- corrlMatrix[subset1, subset1]

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }

        cutvector=seq(mincutval, maxcutval, by=cutstep)

        colname1=c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]

             if(sampled) {
                 dichotCorhelp  <-  I(abs(cor1)>cut1 )+0.0
             }else{
                 dichotCorhelp  <-  I(abs(corrlMatrix)>cut1 )+0.0
             }

             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             print( paste(cut1, ": ", no.genes, "genes, ", sum(nolinkshelp), " links") )

             str_cutoff=paste("tau=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue  = 2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut(nolinkshelp,no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)/length(nolinkshelp)
             binned.k=tapply(nolinkshelp, cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             datout[i,]=signif(c( cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }
        }

        #output data
        #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")
        print(datout);
        write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

        corrlcutoff=NA
        msqrcut = RsquaredCut
        while(is.na(corrlcutoff)){
             ind1   = datout[,3] > msqrcut
             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             msqrcut     = msqrcut - 0.05
        }

        # consider the truncated scalefree index if the scalefree index is not well satisfied
        if (msqrcut+0.05<0.75){
            truncatedSel = (datout[,4] >= 0.9) & (datout[,6]<=60) & (datout[,3] >=0.5)
            if (sum(truncatedSel)>0){
                 colcutoffs = datout[truncatedSel, 1]
                 corrlcutoff= colcutoffs[1]
            }
        }

        corrlcutoff
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation POWER Computation 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations
powercut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, 
mincutval=1,maxcutval=12, bystep=0.5, plotimg=T) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }

        cutvector=seq(mincutval, maxcutval, by=bystep)

        colname1 =c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             if(sampled) {
                 dichotCorhelp      <-  abs(cor1)^cut1
             }else{
                 dichotCorhelp      <-  abs(corrlMatrix)^cut1
             }

             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             print( paste(cut1, ": ", no.genes, "genes, ", sum(nolinkshelp), " links") )

             str_cutoff=paste("beta=", as.character(cut1), sep="")
             #fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             if(plotimg) {
               fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             }else{
               fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="sclfree.png", outputFitness=TRUE)
             }

             pvalue  = -1

             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut((nolinkshelp+1),no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)
             binned.k=tapply(nolinkshelp+1,cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             #pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             pvalue=-1
             datout[i,]=signif(c( cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }
        }

        #in case that the largest R^2 is not bigger than RsquaredCut,
        #we have to reduce it gradually
        corrlcutoff=NA
        msqrcut = RsquaredCut
        while(is.na(corrlcutoff)){
             ind1   = datout[,3] > msqrcut
             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             msqrcut     = msqrcut - 0.05
        }

        # consider the truncated scalefree index if the scalefree index is not well satisfied
        if (msqrcut+0.05<0.75){
            truncatedSel = (datout[,4] >= 0.9) & (datout[,6]<=60) & (datout[,3] >=0.5)
            if (sum(truncatedSel)>0){
                 colcutoffs = datout[truncatedSel, 1]
                 corrlcutoff= colcutoffs[1]
            }
        }

       #output data
       #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")
       print(datout);

       write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

       corrlcutoff
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation Sigmoid Computation, search for alpha
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations
sigmoidcut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, tau=0.5, out=0) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        mincutval=1
        maxcutval=12

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }


        cutvector=seq(mincutval, maxcutval, by=0.2)

        colname1=c("tau", "alpha","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             dichotCorhelp      <-  sigmoidMatrix( abs(cor1), cut1, tau)
             collect_garbage()

             if(sampled) {
                 dichotCorhelp   <-  sigmoidMatrix( abs(cor1), cut1, tau)
             }else{
                 dichotCorhelp   <-  sigmoidMatrix( abs(corrlMatrix), cut1, tau)
             }


             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             str_cutoff=paste("beta=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue  = -1
             datout[i,]=signif(c(tau, cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)



             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut((nolinkshelp+1),no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)
             binned.k=tapply(nolinkshelp+1,cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             #pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             pvalue=-1
             datout[i,]=signif(c(tau, cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }

             rm(dichotCorhelp)
             collect_garbage()
        }

        #in case that the largest R^2 is not bigger than RsquaredCut,
        #we have to reduce it gradually
        corrlcutoff=NA
        msqrcut = RsquaredCut
        while(is.na(corrlcutoff)){
             ind1   = datout[,4] > msqrcut
             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             msqrcut     = msqrcut - 0.05
        }

       #output data
       #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")
       print(datout);

       if (out==0){
          returnval = corrlcutoff
          write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
       }
       else
          returnval = datout
       #return returnval
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##        Function: correlation Sigmoid Computation, search for alpha and tau
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigmoidcutTauAlpha = function(mcorrlMatrix,mnoOFsamples, mmlogFname, msamplingSize=3600, mRsquaredCut=0.8) 
{
   mintau = 0.0
   maxtau = 1.0

   tauvector=seq(mintau, maxtau, by = 0.02)
   
   allData = c()
   for (i in tauvector ){
      idata = sigmoidcut(corrlMatrix=mcorrlMatrix, noOFsamples=mnoOFsamples, mlogFname=mmlogFname, 
                         samplingSize=msamplingSize, RsquaredCut=mRsquaredCut, tau=i, out=1) 
      allData = rbind(allData, idata)
   }
   write.table(allData, mmlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
   
   corrlcutoff=NA
   msqrcut = mRsquaredCut
   while(is.na(corrlcutoff[1])){
     ind1   = allData[,4] > msqrcut
     indcut = NA
     indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)

     corrlcutoff = allData[ indcut[[1]], c(1:2)]
     msqrcut     = msqrcut - 0.05
   }

   #a 2-d vector with [1]: tau and [2]: alpha
   corrlcutoff
}


#------------------------------------------------------------------------------------------------------------------
#Function: get most connected genes based on a scale free network derived from a power adjacency function
#
#Input:
#    minputfname    ~  gene expresion matrix
#    geneinforCols  ~  number of columns as gene information
#    R2Cut          ~  scale free fitness:ouput most connected genes from a the network with scale free fitness>R2Cut OR
#    tructR2Cut     ~  scale free fitness > tructR2Cut
#    topN           ~  top N most connected genes will be output 
#    mincutval      ~  min power value
#    maxcutval      ~  max power value
#
#Output:  1) scale free plot 2) logfile 3) array data with most connected genes
#
#------------------------------------------------------------------------------------------------------------------
getMostConnectedGenesByPowerAdj = function(minputfname, geneinforCols, R2Cut=0.85, tructR2Cut=0.90, topN=3600, mincutval=4, maxcutval=16, maxMeanLinks=30)
{
       fname       =getFileName(minputfname)
       mlogFname   =paste(fname, "_logo.txt",       sep='')
       imgScaleFree=paste(fname, "_imgScaleFree.png",sep='')

       #------- STEP 0: read in gene information, expression data
       allMatrix <- read.delim(minputfname,sep="\t", header=T)
       dim(allMatrix)

       genesInfor <- allMatrix[,c(1:geneinforCols)]
       rowTitles=names(allMatrix)

       #These are the expression values
       datExpr <- t(allMatrix[,-c(1:geneinforCols)])
       no.samples <- dim(datExpr)[1]
       dim(datExpr)

       corhelp <- cor(datExpr, use = "pairwise.complete.obs") 
       no.genes <- dim(corhelp)[[2]]
       dim(corhelp)

       diag(corhelp) <- 0

        cutvector=seq(mincutval, maxcutval, by=1)

        colname1=c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")
        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             
             #get no links for each gene
             nolinkshelp   = rep(0.0, no.genes)
             for (j in c(1:no.genes)){
                sel =  as.numeric(abs(corhelp[j,])) ^ cut1
                nolinkshelp[j]  = sum(sel)
             }

             str_cutoff= paste("beat=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue=-1
             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if( mean(nolinkshelp)<maxMeanLinks & (fitness[1] >= R2Cut || fitness[2]>tructR2Cut) ){
                break
             }
        }
       print(datout);
       write.table(datout[c(1:i),], mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

       str_cutoff=paste("beta=",as.character(cut1), sep="")
       suminfor=ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="")
       ScaleFreePlot(nolinkshelp, no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile=imgScaleFree)

       #*----------------------------- choose the top N genes ----------------------------------------
       orderLink = order(-nolinkshelp)
       nolinks.ordered = nolinkshelp[orderLink]
     
       restFname = paste(fname, "_p", as.character(topN), ".xls",       sep='')
       finalMatrix = allMatrix[orderLink[1:topN], ]
       write.table(finalMatrix, restFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 
       str_cutoff
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: automatically detect & label modules
##
##
##   nameallmodules=FALSE: label modules with all possible colors
##                 =TRUE:  when # of modules exceeds length(colorcode), we use false color 
##                          names to label the reamining modules
##
##   useblackwhite=FALSE: label as normal
##                =TRUE:  label extra modules by black and white alternatively
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#The input is an hclust object.
moduleDetectLabel = function(hiercluster,heightcutoff=0.5,minsize1=20, nameallmodules=FALSE, useblackwhite=FALSE, useNumberAsLabel=F, startlabel=1) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree2(hiercluster,h=heightcutoff)
    sort1    =-sort(-table(labelpred))
    #sort1
    modulename   = as.numeric(names(sort1))
    modulebranch = sort1 > minsize1
    no.modules   = sum(modulebranch)

    if (useNumberAsLabel){
       # now make cluster label
       #
       colorcode=NULL
       for (i in c(startlabel:(startlabel+no.modules-1)) ){
          ipad = patchZeros(i)
          colorcode= c(colorcode, ipad)
       }
    }else{
    # now we assume that there are fewer than 10 modules
    colorcode=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
    }

    # "Grey" means not in any module;
    colorhelp=rep("grey",length(labelpred))
    if ( no.modules==0){
        print("No mudule detected\n")
    }
    else{
        if ( no.modules > length(colorcode)  ){
            print( paste("Too many modules \n", as.character(no.modules)) )
        }

        if ( (nameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
            labeledModules = min(no.modules, length(colorcode) )
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            }
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
            }
        }else{#nameallmodules==TRUE and no.modules >length(colorcode)
            maxcolors=length(colorcode)
            labeledModules = no.modules
            extracolors=NULL
            blackwhite=c("black", "white")
            for(i in c((maxcolors+1):no.modules)){
              if(useblackwhite==FALSE){
                  icolor=paste("module", as.character(i), sep="")
              }else{#use balck white alternatively represent extra colors, for display only
                  icolor=blackwhite[1+(i%%2) ]
              }
              extracolors=c(extracolors, icolor)
            }
            allcolorcode=c(colorcode, extracolors)
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            }
            
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
            }
        }
    }
    colorhelp
}

# make label by patching zeros
patchZeros = function(intnumber, digits=4){
  strnumber = as.character(intnumber)
  splitted  = splitString(strnumber,"")
  mydigits  = length(splitted) 

  patchedDigits = digits - mydigits
  if (patchedDigits<=0){
     return(strnumber)
  }
  myzeros=c("0","00","000","0000", "00000","000000","0000000", "00000000", "000000000")

  patched = paste(myzeros[patchedDigits], strnumber, sep="")

  return (patched)
}




#"0" is for grey module
assignModuleColor = function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE, useNumberAsLabel=F, startlabel=0) {
    # here we define modules by using a height cut-off for the branches
    #labelpred= cutree2(hiercluster,h=heightcutoff)
    #cat(labelpred)

    #"0", grey module doesn't participate color assignment, directly assigned as "grey"
    labelpredNoZero = labelpred[ labelpred >0 ]
    sort1=-sort(-table(labelpredNoZero))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 > minsize1
    no.modules=sum(modulebranch)

    if (useNumberAsLabel){
       # now make cluster label
       #
       colorcode=NULL
       for (i in c(startlabel:(startlabel+no.modules-1)) ){
          ipad = patchZeros(i)
          colorcode= c(colorcode, ipad)
       }
    }else{
    # now we assume that there are fewer than 10 modules
    #colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow")
    # now we assume that there are fewer than 10 modules
    colorcode=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
    }


    #"grey" means not in any module;
    colorhelp=rep("grey",length(labelpred))
    if ( no.modules==0){
        print("No mudule detected\n")
    } else{
        if ( no.modules > length(colorcode)  ){
            print( paste("Too many modules \n", as.character(no.modules)) )
        }

        if ( (anameallmodules==FALSE) | (no.modules <=length(colorcode)) ){
            labeledModules = min(no.modules, length(colorcode) )
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            }
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
            }

        }else{#nameallmodules==TRUE and no.modules >length(colorcode)
            maxcolors=length(colorcode)
            labeledModules = no.modules
            extracolors=NULL
            blackwhite=c("red", "black")
            for(i in c((maxcolors+1):no.modules)){
              if(auseblackwhite==FALSE){
                  icolor=paste("module", as.character(i), sep="")
              }else{#use balck white alternatively represent extra colors, for display only
                  #here we use the ordered label to avoid put the same color for two neighboring clusters
                  icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
              }
              extracolors=c(extracolors, icolor)
            }

            #combine the true-color code and the extra colorcode into a uniform colorcode for 
            #color assignment
            allcolorcode=c(colorcode, extracolors)

            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            }

            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
            }
        }
    }

    colorhelp
}


# label each grey node as a single module until all colors used up
#
turnGreyNodesIntoModules = function(incolorcode) {

    # now we assume that there are fewer than 10 modules
    #colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow")
    # now we assume that there are fewer than 10 modules
    colorcode=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")

    maxcolors = length(colorcode)
    colorcodeIndex = cbind(colorcode, c(1:maxcolors)) 
    
    cincolorcode   = as.character(incolorcode)

    # find the maximum index of the last color used
    #
    merged = merge(incolorcode, colorcodeIndex, by.x=1, by.y=1, all=F)
    maxColIdx = max(as.integer(as.character(merged[,2])))


    # turn each grey into a new color
    #
    no.items = length(cincolorcode)
    selgrey  = cincolorcode=="grey"
    if(sum(selgrey)==0 | maxColIdx>maxcolors){
       return (incolorcode)
    }

    greyIdx  = c(1:no.items)[selgrey]
    no.greys = length(greyIdx)
    
    curcolorIdx = maxColIdx + 1
    for (g in c(1:no.greys)){
       cincolorcode[ greyIdx[g] ] = colorcode[curcolorIdx]
       curcolorIdx = curcolorIdx + 1
       if(curcolorIdx>maxcolors){
           break;
       }
    }

    colorhelp=factor(cincolorcode)

    return (colorhelp)
}





#merge the minor cluster into the major cluster
merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = ifelse(as.character(mycolorcode)==minorclusterColor, mainclusterColor, as.character(mycolorcode) )
  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}

# unnamed modules are represented by like "module12", thus these module names include keyword "module"
# > unlist(strsplit("module21","module"))
# [1] ""   "21"
# > unlist(strsplit("yellow","module"))
# [1] "yellow"
assignPseudoColors = function(mdendro, incolorcode){

  #locate ordered locations of clusters in the dendrogram
  fdispcolor = getDisplayColorSequence(as.character(incolorcode[mdendro$order]))
  labeledColors= fdispcolor[fdispcolor!=""]

  mycolorcode= as.character(incolorcode)

  pseudocolors=c("white","black")
  cnt.unnamedcolor = 0
  for ( each in labeledColors){
    is.unnamedcolor = unlist(strsplit(each,"module"))	
    if (length(is.unnamedcolor) == 1){ #the module is labeled with true color, and keep it
        next
    }
    
    # the module is labeled with true color, so we need replace it with 
    if (cnt.unnamedcolor==1){
        cnt.unnamedcolor =2
    }else{
        cnt.unnamedcolor =1
    }

    ipseduocol = pseudocolors[cnt.unnamedcolor]

    #cat(each, " to ", ipseduocol, "\n")

    mycolorcode = ifelse(mycolorcode==each, ipseduocol, mycolorcode)
  }

  fcolorcode   =factor(mycolorcode)
  fcolorcode
}


exchange2Clusters = function(mycolorcode, colorA, colorB){
  colorC   = "great"
  fcolcode = mycolorcode
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorC, minorclusterColor=colorA) #A->C
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorA, minorclusterColor=colorB) #B->A
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorB, minorclusterColor=colorC) #C->B
  fcolorcode   =factor(fcolcode)
  fcolorcode
}

# return the name of the max element from a set of element names
# 
getNameOfMaxElement=function(mytable, selectednames){
   maxsize=0
   maxname=""
   for (each in selectednames){
       if (mytable[[each]] >maxsize){
           maxsize = mytable[[each]]
           maxname = each
       }
   }
   return (maxname)
}


# cut dendrogram to get the required number of clusters
#
moduleDetectByFixedClusterno = function(ihcluster, nocluster=2) {
   maxhei   = max(ihcluster$height)
   minhei   = min(ihcluster$height)
   curheightcutoff = minhei
   curno = 0
   cnt=0
   while(curno != nocluster){
       cnt = cnt + 1
       if (cnt > 6)
          break
       curheightcutoff = (curheightcutoff + maxhei)/2.0
       colcode  = moduleDetectLabel(hiercluster=ihcluster, curheightcutoff, minsize1=1)
       colcode  = as.character(colcode)
       # we need as.character(colcode) to deal with the case of grey is zero
       #table(colcolors)
       #turquoise      blue      grey 
       #12             6         0 
       colr = as.integer(table(colcode))
       curno    = sum(colr>0)
   }
   colcode   
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                           Get Most Varying Genes  
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get nvars most varying genes, consider only genes with less than 20% missing values
## default use of standard deviation
## other option: coefficient of variation (CV), i.e., the standard deviation divided by mean
getMostVaryingGenes=function(myinputfname, nvars, myHeadCol=8, missingRate=0.2, mysep="\t", vartype=0){
    
    mfname      =getFileName(myinputfname)
    mvgOutFname =paste(mfname, "_", nvars, "mvgenes.xls", sep='')

    allMatrix <- read.delim(myinputfname,sep=mysep, header=T)
    dims=dim(allMatrix)

    genesInfor <- allMatrix[,c(1:myHeadCol)]
    rowTitles=names(allMatrix)
    rowTitles=rowTitles[1:dims[2]]

    dat1 <-t(allMatrix[,-c(1:myHeadCol)])

    ## --------- compute genes with missing values
    naMatrix <- is.na(dat1)
    nasumVector  <- apply(naMatrix,2,sum, na.rm=TRUE)
    naprobVector <- nasumVector/dim(dat1)[[1]]

    dat1 <- dat1[,naprobVector<missingRate]
    dim(dat1)
    if(myHeadCol>1){
      genesInfor_NAclean <- genesInfor[naprobVector<missingRate,]
    }else{#geneinfor is a vector
      genesInfor_NAclean <- genesInfor[naprobVector<missingRate]
    }

    var1 <- apply(dat1,2,sd,   na.rm=TRUE)
    means <- apply(dat1,2,mean, na.rm=TRUE)

    if (vartype !=0 ){
       var1 <- var1/means        
    }

    rk1  <- rank(-var1)
    dat2 <- dat1[,rk1<nvars+1]

    noGenes   <- dim(dat2)[2]
    noSamples <- dim(dat2)[1]

    if(myHeadCol>1){
      filteredGenesInfor <- genesInfor_NAclean[rk1<nvars+1, ]
    }else{
      filteredGenesInfor <- genesInfor_NAclean[rk1<nvars+1]
    }

    newCol = length(rowTitles)
    yfile <- file(mvgOutFname, "w+")
    #first row
    for (z in 1:newCol){
       cat(as.character(rowTitles[z]), file=yfile)
      if (z==newCol){
         cat("\n", file=yfile)
      }
      else{
         cat("\t", file=yfile)
      }
    }

    for (z in 1:noGenes){
      if(myHeadCol>1){
          for (i in 1:dim(filteredGenesInfor)[2] ){
             cat(as.character(filteredGenesInfor[z,i]), file=yfile)
             cat("\t", file=yfile)
          }
      }else{
          cat(as.character(filteredGenesInfor[z]), file=yfile)
          cat("\t", file=yfile)
      }

      for (i in 1:noSamples){
         cat(as.character(dat2[i,z]), file=yfile)
         if (i==noSamples){
            cat("\n", file=yfile)
         }
         else{
            cat("\t", file=yfile)		
         }
      }
    }
    close(yfile)

    mvgOutFname
}

#----------------------------------------------------------------------------------------------------
#write a table with row and colnames into Excel or Latex output
#Input
#     myMatrix: table 
#     fname:    output filename
#     latex:    output Latex(T) or Excel file(F)
#     myappend: append to the file or not(T/F)
#     usesignif: use short repsentation of floating numbers (T/F)
#----------------------------------------------------------------------------------------------------
writeTableIntoLatex=function(myMatrix, fname, myappend=F, latex=T, mycol.names=T, myrow.names=T, usesignif=F){
   
   mycolnames = colnames(myMatrix)
   myrownames = rownames(myMatrix)

   #put colnames at the first row and rownames as the 1st column, leave a blank for the cell [1,1]
   latexMatrix1 = rbind(mycolnames, as.matrix(myMatrix) )
   latexMatrix  = cbind( c(" ", myrownames), as.matrix(latexMatrix1))
   no.rows=dim(latexMatrix)[1]
   no.cols=dim(latexMatrix)[2]

   if(latex==FALSE){
      if(usesignif==TRUE){
         for (i in c(2:no.rows) ){
           for (j in c(2:no.cols) ){
              latexMatrix[i, j] = signif(as.numeric(latexMatrix[i,j]),2)
           } 
         }
      }
      write.table(latexMatrix,fname, append=myappend, sep="\t", quote=FALSE, col.names=F, row.names=F)
      return (1)
   }

   for (i in c(1:no.rows) ){
     #add & to each cell except those in the last column
     for (j in c(1:(no.cols-1)) ){
         if((usesignif==TRUE)&(i>1)&(j>1) ){
           latexMatrix[i, j] = paste(signif(as.numeric(latexMatrix[i,j]),2), "&", sep="")
         }else{
           latexMatrix[i, j] = paste(latexMatrix[i,j], "&", sep="")
         }
     }

     #\\\\\\hline to each row
     if((usesignif==TRUE)&(i>1) ){
         latexMatrix[i, no.cols] = paste(signif(as.numeric(latexMatrix[i, no.cols]),2), "\\\\\\hline", sep="")
     }else{
         latexMatrix[i, no.cols] = paste(latexMatrix[i, no.cols], "\\\\\\hline", sep="")
     }
   }

   write.table(latexMatrix,fname, append=myappend, sep="\t", quote=FALSE, col.names=F, row.names=F)
}



compareTwoModules = function(gifA, gifB, moduleNameA, moduleNameB, uniqueIdCol=1, moduleColA, moduleColB, totalGenes)
{
  restrictA = gifA[, moduleColA]== moduleNameA 
  restrictB = gifB[, moduleColB]== moduleNameB
  noA= sum(restrictA)
  noB= sum(restrictB)

  moduleSetA = gifA[restrictA, ]
  moduleSetB = gifB[restrictB, ]

  mergeMatrix = merge(moduleSetA, moduleSetB, by.x=uniqueIdCol, by.y=uniqueIdCol, sort=F,all=FALSE)
  intersectNo = dim(mergeMatrix)[1]

  #phyper(89,702,4000-702,280, lower.tail=F)
  pval = phyper(intersectNo, noA, totalGenes-noA, noB, lower.tail=F)

  ret=c(intersectNo, pval, noA, noB)
  ret
}

overlapTwoLists = function(listA, listB, totalGenes, retentries=F)
{
  noA = length(listA)
  noB = length(listB)
  mergeMatrix = merge(listA, listB, by.x=1, by.y=1, sort=F,all=FALSE)
  intersectNo = dim(mergeMatrix)[1]
  
  if (intersectNo >0){
    #phyper(89,702,4000-702,280, lower.tail=F)
    pval = phyper(intersectNo, noA, totalGenes-noA, noB, lower.tail=F)
    commongenes = concatenate( as.character(as.matrix(mergeMatrix)),";")
  }else{
    pval = 1
    commongenes = ""
  }

  if (retentries){
    ret=c(intersectNo, pval, commongenes)
  }else{
    ret=c(intersectNo, pval)
  }
  ret
}



#************************************
#* Given d and r, we compute s, s= d/(sin(pi/2 - atan(r))
#  /|
#/  |  /
#   |/
#  /|\ d  /
#/ s|  \/
#   | /
#   /)r
# /-------------
loessfitbounds=function(y, x){
     myseq = as.numeric(seq(min(x),to=max(x), length=20 )) # length=length(x)) )
     length(myseq)
     new <- data.frame( myseq )

     corrl= cor(x,y)
     lm1=lm(y ~ x)
     r= lm1$coefficients[2]
     d=1.0
     if (r >= 5*pi/180){
       s= d/( sin(pi/2 - atan(r)) )
     }else{
       s= d
     }

     loe1=loess(y ~ x, span=1, degree=1)
     predloe1=predict(loe1, new, se=T)
     myfit= predloe1$fit

     #return x_fitted, y_fitted, mean(se)
     #list(myseq, myfit, mean(predloe1$se))
     list(myseq, myfit, predloe1$se)
     #list(myseq, myfit, s)
}

# textPosY as the portion of the range Y to be shifted for text display
#
scatterPlot=function(aa.pvector, bb.pvector, myxlab="", myylab="", imgfilename=NULL, 
                     pointcolor=NULL, textlabels=NULL, textcolor=NULL, shortnote="", 
                     useloess=F, imgwidth=400, imgheight=400, icex=0.6, textPosY=0)
{   
   lm1=lm(bb.pvector ~ aa.pvector)

   if(useloess) #return x_fitted, y_fitted, mean(se)
      fitted=loessfitbounds(y=bb.pvector, x=aa.pvector)

   no.points = length(aa.pvector)

   corrl= cor(aa.pvector, bb.pvector)
   corrl.p = signif( cor.test(aa.pvector, bb.pvector, method="p",use="p")$p.value ,2)
   if (corrl.p>0){
       corrp.string = as.character(corrl.p)
       mytitle=paste(shortnote, "correlation=", as.character(signif(corrl, 2)),", pvalue =", corrp.string, sep="")
   }else{
       corrp.string = " < e-22"
       mytitle=paste(shortnote, "correlation=", as.character(signif(corrl, 2)),", pvalue", corrp.string, sep="")
   }


   #plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab )

   mse=2
   bcolor="blue"

   #ostring=paste(as.character(signif(corrl,2)), "\t", as.character(corrl.p),"\n" )
   #cat(ostring, file=fp)

   if (is.null(pointcolor) ){
        pcolors=rep("black", no.points)
   }else{
        pcolors=pointcolor
   }
   if (is.null(textcolor) ){
        tcolors=rep("black", no.points)
   }else{
        tcolors=textcolor
   }

   if(!is.null(imgfilename)){
     openImgDev(imgfilename, iwidth =imgwidth, iheight =imgheight, ipointsize = 12)
     #openImgDev(imgfilename)
     par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)

     plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab, col=pcolors)

     if(useloess){
       lines(x=fitted[[1]], y=fitted[[2]], col=bcolor)
       lines(x=fitted[[1]], y=fitted[[2]]+mse*fitted[[3]], col=bcolor,lty=2)
       lines(x=fitted[[1]], y=fitted[[2]]-mse*fitted[[3]], col=bcolor,lty=2)
     }else{
       abline(lm1, col=bcolor)
     }

     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(textlabels) ){
        adj=( max(bb.pvector)-min(bb.pvector))*textPosY
        text(aa.pvector, bb.pvector+adj, labels = textlabels, col=tcolors, cex=0.6)
     }
     dev.off()
   }

     #par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
     plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab, col=pcolors)
     

     #par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(textlabels) ){
        #yshift= rep((max(bb.pvector)-min(bb.pvector))*textPosY, length(bb.pvector) )
        adj=( max(bb.pvector)-min(bb.pvector))*textPosY
        text(aa.pvector, bb.pvector+adj, labels = textlabels, col=tcolors, cex=icex)
     }


     if(useloess){
       lines(x=fitted[[1]], y=fitted[[2]], col=bcolor)
       lines(x=fitted[[1]], y=fitted[[2]]+mse*fitted[[3]], col=bcolor,lty=2)
       lines(x=fitted[[1]], y=fitted[[2]]-mse*fitted[[3]], col=bcolor,lty=2)
     }else{
       abline(lm1, col=bcolor)
     }

   mytitle   
}

# here we identify outliers with difference of fitted value and true value is 
#   equal to or bigger than mnoSDsAsOutlier*sd(fitted differences)
identify_linearfit_outliers=function(yy, xx, mnoSDsAsOutlier=2){
   lm1=lm(yy~xx)
   #plot(a, c)
   #abline(lm1)
   fitted= as.numeric(lm1$fitted.values)
   diffs = abs(fitted-c)
   mean.diff =  mean(diffs)
   sd.diff   =  sd(diffs)
   threshold = mean.diff + sd.diff*mnoSDsAsOutlier

   goods  = (diffs < threshold) # points with fitted value less than mean + n*standard deviation
   goods
}

# we consider an option of removing ouliers
singlePlotComplex = function(outputImage=NULL, rawValA, rawValB, plotTitle="", myylab="",myxlab="", dataSet="", userawdata=F, numPoints=20, removeOutliers=F, noSDsAsOutlier=2){

   if(!userawdata){
	#in general, we will have the thing being tested be rawValA, and K be rawValB
	cut_factor = cut(rawValB,numPoints)
	# Now give me the mean proportion of essentiality in each chunk of the essentiality vector 
	ValA = tapply(rawValA, cut_factor, mean)
	# Now give me the mean proportion of k in each chunk of the essentiality vector 
	ValB = tapply(rawValB, cut_factor, mean)
        
        ValA = ValA[!is.na(ValA)]
        ValB = ValB[!is.na(ValB)]

   }else{
        ValA = rawValA
        ValB = rawValB
   }


   if(removeOutliers){
      good.points = identify_linearfit_outliers(ValA, ValB, mnoSDsAsOutlier=noSDsAsOutlier)
      #print(length(good.points) )
      #print(sum(good.points) )

      ValA= ValA[good.points]
      ValB= ValB[good.points]
   }
        #print(ValA)
        #print(ValB)


	#Now we need to just make a line fit
	line =lm(ValA ~ ValB)

	corr_p   = signif( cor(ValB,ValA),2)	#pearson 
        corr_p.p = signif( cor.test(ValB,ValA, method="p",use="p")$p.value ,2)

	corr_s   = signif(cor(ValB,ValA, method = "s"),2)	#spearman 
        corr_s.p = signif( cor.test(ValB,ValA, method="s",use="p")$p.value ,2)

	titleinfor = paste(as.character(dataSet),",",as.character(plotTitle), ",", "p_r=", as.character(corr_p), ",p_p=", as.character(corr_p.p)
	,"(s_r=", as.character(corr_s),",s_p=", as.character(corr_s.p),")",sep="" )

	# do a plot() of these two vectors against each other. 
	if (!is.null(outputImage) )
	    openImgDev(outputImage)

        plot(ValB,ValA,main=titleinfor, xlab=myxlab, ylab=myylab)
     	abline(line, col="red")

	# and then draw the line
	if (!is.null(outputImage) )
             dev.off()
        #titleinfor
    
    return ( c(corr_p, corr_p.p, corr_s, corr_s.p) )
}



#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleLabels=function(mdendro, modulecolors, save2file=NULL, plotLabels=FALSE, secondColorbar=NULL, secondLabel=NULL, 
      gramcolor="black", plotSampleNames=F, pwid=1024,phei=600)
{
   fdispcolor = getDisplayColorSequence(as.character(modulecolors[mdendro$order]))
   labeledColors= fdispcolor[fdispcolor!=""]
   labeledColXcoor= c(1:length(fdispcolor))[fdispcolor!=""]

   pseudocolor = assignPseudoColors(mdendro=mdendro, incolorcode=modulecolors)

   if(!is.null(save2file)){ 
      openImgDev(save2file, iwidth = pwid, iheight = phei)
   }

   if(!is.null(secondColorbar) & plotLabels) {
     par(mfrow=c(3,1), mar=c(5,0,0,0), cex=0.6)
   } else if(plotLabels){ #show color labels
     par(mfrow=c(2,1), mar=c(10, 0, 0, 0), cex=0.6)
   } else if(!is.null(secondColorbar) ) {
     par(mfrow=c(3,1), mar=c(0,0,0,0) )
   } else{
     par(mfrow=c(2,1), mar=c(0,0,0,0) )
   } 

   if (plotSampleNames){
     plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }else{
     plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }

   #plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
   #          ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
   #          lwd=1, axes=F, frame=F)

   plot(c(1:length(modulecolors)), rep(1, length(pseudocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(pseudocolor[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)


   if(plotLabels){ #show color labels
      axis(1, at =labeledColXcoor, labels = labeledColors, las=2)
   }
   #barplot(height=rep(1, length(modulecolors)), names.arg=fdispcolor, 
   #       col= as.character(modulecolors[mdendro$order]), space=0, las=2,
   #      border=F,main="", axes = T, axisnames = T)#, cex.lab=0.5)

   if(!is.null(secondColorbar) ) {
     colorlabel2 = as.character(secondColorbar[mdendro$order])
     plot(c(1:length(secondColorbar)), rep(1, length(secondColorbar)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=colorlabel2, las=2,
             lwd=1, axes=F, frame=F)     

     if(!is.null(secondLabel)) {
        labeledColXcoor2 = c(1:length(colorlabel2))[colorlabel2!=""]
        lablels2         = secondLabel[colorlabel2 !=""]
        axis(1, at =labeledColXcoor2, labels =lablels2 , las=2)
     }

   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

   if(!is.null(save2file)){ 
      dev.off()
   }
}



#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleLabelsHubs=function(mdendro, modulecolors, save2file=NULL, plotLabels=FALSE, secondColorbar=NULL, secondLabel=NULL, 
               gramcolor="black", plotSampleNames=F, pwid=1024,phei=800)
{
   fdispcolor = getDisplayColorSequence(as.character(modulecolors[mdendro$order]))
   labeledColors= fdispcolor[fdispcolor!=""]
   labeledColXcoor= c(1:length(fdispcolor))[fdispcolor!=""]

   pseudocolor = assignPseudoColors(mdendro=mdendro, incolorcode=modulecolors)

   if(!is.null(save2file)){ 
      openImgDev(save2file, iwidth = pwid, iheight = phei)
   }

   if(!is.null(secondColorbar) & plotLabels) {
     par(mfrow=c(1,1), mar=c(10,0,0,0), cex=0.8)
   } else if(plotLabels){ #show color labels
     par(mfrow=c(2,1), mar=c(10, 0, 0, 0), cex=0.6)
   } else if(!is.null(secondColorbar) ) {
     par(mfrow=c(3,1), mar=c(0,0,0,0) )
   } else{
     par(mfrow=c(2,1), mar=c(0,0,0,0) )
   } 


   if (plotSampleNames){
     if(!is.null(secondLabel)) {
       plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor, labels=secondLabel)
       #plclust(mdendro, hang=2, xlab="",ylab="",main="",sub="",axes = T, labels=secondLabel)
     }else{
       plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
     }
   }else{
     plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }


   if(!is.null(secondColorbar) ) {
     colorlabel2 = as.character(secondColorbar[mdendro$order])
     secondLabel2 = secondLabel[mdendro$order]
     #plot(c(1:length(secondColorbar)), rep(1, length(secondColorbar)), ylim=c(0,1),xlab="", main="",
     #        ylab="",  type = "h",col=colorlabel2, las=2,
     #        lwd=1, axes=F, frame=F)     

     if(!is.null(secondLabel)) {
        labeledColXcoor2 = c(1:length(colorlabel2))[colorlabel2!=""]
        lablels2         = secondLabel2[colorlabel2 !=""]
        axis(1, at =labeledColXcoor2, labels =lablels2 , las=2, 
               col = "violet", col.axis="dark violet", lwd = 1)
     }

   }

   if(1==2){
   #plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
   #          ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
   #          lwd=1, axes=F, frame=F)

   plot(c(1:length(modulecolors)), rep(1, length(pseudocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(pseudocolor[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)


   if(plotLabels){ #show color labels
      axis(1, at =labeledColXcoor, labels = labeledColors, las=2)
   }
   #barplot(height=rep(1, length(modulecolors)), names.arg=fdispcolor, 
   #       col= as.character(modulecolors[mdendro$order]), space=0, las=2,
   #      border=F,main="", axes = T, axisnames = T)#, cex.lab=0.5)

   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

   if(!is.null(save2file)){ 
      dev.off()
   }
}


#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleMitosis=function(mdendro, modulecolors, addmitocolor, save2file=NULL)
{
   if(!is.null(save2file)){ 
      openImgDev(save2file)
   }
   par(mfrow=c(3,1), mar=c(0,0,0,0) )
   plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T)
   plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)
   plot(c(1:length(addmitocolor)), rep(1, length(addmitocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(addmitocolor[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)
   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
   if(!is.null(save2file)){ 
      dev.off()
   }
}


#---------------------------------------------------------------------------------------------------------
#ModulePrincComps, which computes the first principal component of each module. 
#Also it provides the variance explained for the first 5 PCs of a module.
#It is based on the singular value decomposition.
#It takes as input datExpr  for which the rows are samples and the columns are genes.
#Here is how you would use it
#PC1=ModulePrinComps2[datExpr,color1]
#Then PC1[[1]] provides a data frame with the first PC of each module
#PC1[[2]] provides a data frame where each column contains the percentages of 
#the variance explained by the first 10 PCs by module.
#---------------------------------------------------------------------------------------------------------
# the output is a list and each element is a matrix whose number of columns equals the number of modules
# the first element of the list contains the 1st principal components of all modules, and 
# the second element of the list contains the 2nd principal components of all modules, etc...
# the last element is the expained variations of top "no.pcs" PCs in each modules
#---------------------------------------------------------------------------------------------------------

ModulePrinComps = function(datexpr,couleur, min_modulesize=10) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels

  for(i in c(1:length(modlevels)) ){
 
    print(paste(i, modlevels[i]) )
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename

    if(modlevels[i]=="grey"){
       next
    }

    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])

    # check whether some samples have missing rate (missing value in column) >80%
    naDat      = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff   = 0.8*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
        print("patch samples with missing rate >=0.8")
        naSampleIdx = c(1:length(nasBySample))[selSamples]
        for(x in naSampleIdx){
            #print(paste("Sample Idx=", x) )
            datModule[,x]=ifelse(is.na( datModule[,x]), 0, datModule[,x])
        }
    }

 
    if(sum(restrict1)<min_modulesize){
       listPCs[[1] ][,i] = datModule[1,]
       next
    }

    datModule=impute.knn(as.matrix(datModule)) #$data for R version < 2.0

    datModule=t(scale(t(datModule)))
    svd1=svd(datModule)
    mtitle=paste("PCs of ", modulename," module", sep="")

    no.samples = dim(datModule)[2]

    actualpcs=min(dim(svd1$v)[2], no.pcs)

    #cat(modulename, as.character(i), as.character(actualpcs), "\n")

    #explained variation (%)
    listPCs[[ no.pcs+1] ][,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)

    # this is the first principal component, the ith column of the j-th element in the list
    for (j in c(1:actualpcs) ){
         #print(j)        
         pcj=svd1$v[,j]

         # detect NAs
         jnas = is.na(pcj)
         if (sum(jnas)>0){
             break;
         }

         signhj=sign(sum(cor(pcj,  t(datModule))))
         if (signhj != 0)  pcj=signhj* pcj
         listPCs[[j] ][,i] = pcj
    }

  }

  listPCs
}

ModuleHubProfiles = function(datexpr, adjmatrix, couleur, min_modulesize=10) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels


    # take the profile of the most connected gene in each module as PC
    #
    colorI = as.character(colcode.reduced)
    kin    = computeModuleLinks(adjmatrix, couleur)

    # find the hub in each module 
    #
    kinIdxed= cbind(kin, c(1:length(couleur)), couleur)
    orderK  = order(-kin)
    kinIdxed= kinIdxed[orderK, ]
    orderK  = order(kinIdxed[,3])
    kinIdxed= kinIdxed[orderK, ]
    
    hubIdx    = rep(0, length(modlevels) )
    for(z in c(1:length(modlevels)) ) {        
        isel      = modlevels[z] == kinIdxed[,3]
        ikinIdxed = kinIdxed[isel, ]

        # extract hubs' profiles
        #
        listPCs[[1] ][,z] = datexpr[,as.integer(ikinIdxed [1,2])]
        hubIdx[z] = ikinIdxed [1,2]
    }

    # extract hubs' profiles
    #
    #listPCs[[1]] = t(datexpr[,as.integer(hubIdx)])
    names(listPCs[[1]]) <- modlevels

    return(listPCs)
}


plotPrinComps=function(pcVariances,fieldname="", save2file=FALSE)
{
  modulename = names(pcVariances)
  imgname = paste(fieldname, "_", modulename, "_PCs.png", sep="")

  if (save2file)
    openImgDev(imgname)
 
  no.dispPSs = min(8, length( as.matrix(pcVariances)) )
  
  dispnames=rep("Comp.", no.dispPSs )
  for (i in c(1:no.dispPSs) ){
     dispnames[i] =paste("Comp.", as.character(i), sep="")
  }

  mytitle = paste("PCs of ", modulename, " module",sep="")

  par(mfrow=c(1,1), mar=c(4, 6, 4, 2) + 0.1, cex=1)
  barplot(as.numeric(as.matrix(pcVariances))[c(1:no.dispPSs)], 
          names.arg=dispnames,ylim=c(0,1), col=heat.colors(no.dispPSs),
          ylab="Variances",main=mytitle )
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

  if (save2file)
     dev.off()

}

#take the PCs_list, Traits Matrix, moduleIndex (which module)  as input and output correlation between each of
# top N PCs and each of 20 traits, into a figure and output the table into the logo file
#columns of mydatTraits are traits and rows are mice
ComputePlotCorrelOfPCsAndTraits_InModule=function(pcs_list, moduleIndex, mydatTraits, TopNpcs = 4, fieldname=NULL, logfile=NULL){

        topNpcs = min(TopNpcs, length(pcs_list)-1 )

        #module name
        modulename=colnames(pcs_list[[1]])[moduleIndex]

        corrlPCtraits = NULL
        pcnames = NULL
        for (i in c(1:topNpcs)){
           corrlvectA=abs( cor(pcs_list[[i]][,moduleIndex], mydatTraits, use = "pairwise.complete.obs") )
           corrlvect =ifelse(is.na(corrlvectA), 0, corrlvectA)
           corrlPCtraits=rbind(corrlPCtraits, corrlvect)
           pcnames =c(pcnames, paste("PC.", as.character(i), sep=""))
        }
        rownames(corrlPCtraits) <- pcnames

        mycolors=heat.colors( dim(corrlPCtraits)[1])
        mytitle =paste(modulename, " module:", " correlations of traits and the top ", 
                        as.character(topNpcs), " PC.s", sep="")

        if( !is.null(fieldname)){
          pcTraitsfname=paste(fieldname, "_", modulename, "_PCs_corrlTraits.png", sep='')         
          openImgDev(pcTraitsfname,iwidth = 1278, iheight = 768)
        }
        par(mfrow=c(1,1), mar=c(12, 6, 4, 2) + 0.1, cex=1)
        mp  <- barplot(corrlPCtraits)
        tot <- colMeans(corrlPCtraits)
        text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
        barplot(corrlPCtraits, beside = TRUE, ylab="correlation", las=2, 
             col =mycolors, 
             main=mytitle,
             legend = rownames(corrlPCtraits) )
        abline(h=0, col="black")
        abline(h=0.3, col="grey")
        #abline(h=-0.3, col="grey")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        if( !is.null(fieldname)){
           dev.off()
        }

        if( !is.null(logfile))
           appendTableToFile(logfile, round(corrlPCtraits,3), tableTitle=mytitle, myappend=T)

        corrlPCtraits
}


ComputePlotPCsQTL_InModule=function(pcs_list, moduleIndex, TopNpcs = 4, fieldname=NULL, logfile=NULL,minLOD=2.0){
#ComputePlotPCsQTL_InModule=function(pcs_list, moduleIndex, mydatGeno, TopNpcs = 4, fieldname=NULL, logfile=NULL){       

        topNpcs = min(TopNpcs, length(pcs_list)-1 )

        #module name
        modulename=colnames(pcs_list[[1]])[moduleIndex]

        qtl_list = as.list( rep(NULL,topNpcs+1) )
        pcnames = NULL
        for (i in c(1:topNpcs)){
           datGeno$pheno[,2] <- as.numeric( pcs_list[[i]][,moduleIndex] )           
           qtl_list[[i]] = scanone(datGeno, pheno.col=2, method="hk")
           pcnames =c(pcnames, paste("eQTL of PC.", as.character(i), sep=""))
        }

        mycolors=heat.colors(topNpcs)
        mytitle =paste(modulename, " module:", "eQTL of the top ", 
                        as.character(topNpcs), " PC.s",sep="")


        if( !is.null(fieldname)){
          pcQTLfname=paste(fieldname, "_", modulename, "_PCs_QTL.png", sep='')
          openImgDev(pcQTLfname,iwidth = 1278, iheight = 1400)
        }

        par(mfrow=c(topNpcs,1), mar=c(4, 5, 4, 2) + 0.1, cex=0.8)
        for (i in c(1:topNpcs)){            
            mytitle = paste(modulename, " module:", pcnames[i], sep="")
            plot(qtl_list[[i]], col = "green", main=mytitle)
            abline(h=minLOD, col="grey")
        }
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

        if( !is.null(fieldname)){
           dev.off()
           plot(qtl_list[[1]], col = "green", main=pcnames[1])
        }
        qtl_list
}

# return the interval of the maximum valley (blank area in the profile)
#--
#  \                 /\             |\
#   ----------------/  \------------| \----
#      max valley
#
findMaxValleyInProfile =function(profile, xpos){   

   nopoints = length(profile)
   index    = c(1:nopoints)

   sel   = profile<Inf
   nonInfprofile = profile[sel]
   pmin = min(nonInfprofile )
   pmax = max(nonInfprofile )
   yrange = pmax-pmin

   threshB= pmin + yrange*0.75
   pminusB= profile  - threshB
   selB   = pminusB >0

   # transition points are those points that are just exactly above the threshold
   #
   transitpoints = c(1, index[selB], nopoints)

   # get the transition points' physical position on x-axis
   #
   transitxpos   = c(xpos[1], xpos[selB])
   notransits    = length(transitxpos)
   shifttransxpos= c(transitxpos[c(2:notransits)], xpos[nopoints])
   distance      = as.vector(shifttransxpos-transitxpos)
   maxDidx       = which.max(abs(distance)) #left side of max valley 

   maxValleyLeftxPos= transitxpos[maxDidx]

   # give a little space between
   #
   if(distance[maxDidx]>0){
      maxValleyLeft = maxValleyLeftxPos + distance[maxDidx]/8
   
   } else{ # descendant order of the index
      maxValleyLeft = maxValleyLeftxPos + distance[maxDidx] - distance[maxDidx]/8
   }

   return (maxValleyLeft)   
}

#######################################################################################################
#                    Plot multiple eQTL profiles in one figure
#
#  eQTL_profiles: row as markers and column as genes
#
#    type: 1-character string giving the type of plot desired.  The
#          following values are possible, for details, see 'plot': '"p"'
#          for points, '"l"' for lines, '"o"' for overplotted points and
#          lines, '"b"', '"c"') for (empty if '"c"') points joined by
#          lines, '"s"' and '"S"' for stair steps and '"h"' for
#          histogram-like vertical lines.  Finally, '"n"' does not
#          produce any points or lines.

plotMultieQTLProfiles=function(eQTL_profiles, smarkerpos, smarkerXlabel, sboundaries, filename=NULL, minLOD=2.5, singleView=F, ylab="LOD score", plottype='l'){
        no.profiles  = dim(eQTL_profiles)[2]
        no.pmarkers   = dim(eQTL_profiles)[1]
        profilenames <- colnames(eQTL_profiles)

        colorcodeC=c("red","black","blue","green","magenta","yellow","brown","pink","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow","grey")
        # make x labels by concatenating chromosome labels
            bool.chroms=smarkerXlabel != ""
            chroms=smarkerXlabel[bool.chroms]
            no.chroms=length(chroms)-1 # last one is "."
            if (no.chroms>=16){
               xlabel="genome "
            }else{
               xlabel="markers on chromosome "
               for (i in c(1:no.chroms) ){
                    if (i==1){
                      if (no.chroms==1)
                         xlabel=paste(xlabel, chroms[i], sep="") 
                      else
                         xlabel=paste(xlabel, chroms[i], ",", sep="")
                    }else if (i==no.chroms){
                      xlabel=paste(xlabel, chroms[i], sep="")
                    }else{
                      xlabel=paste(xlabel, chroms[i], ",", sep="")
                    }
               }
            }

    showaxis=T
    if (no.chroms>=2)
      showaxis=F

    if(singleView){
        myhei = 400
        mywid = 600
    }else{
        myhei = 300*no.profiles
        mywid = 1278
    }

    if( !is.null(filename)){
          openImgDev(filename,iwidth = mywid, iheight = myhei)
    }

    if(singleView){
            #itype=c(1:no.profiles)
            itype=rep(1, no.profiles)
            icolor=colorcodeC[c(1:no.profiles)]

            # simple way to figure out where to put the legend, x position is tricky
            maxprofiles = apply(eQTL_profiles, 1, max)
            if(1==2){
            maxprforders = order(-maxprofiles)
            peakposIdx = maxprforders[1]     # peak index
            peakpos = smarkerpos[peakposIdx] # peak position
            leftmean = mean(maxprofiles[c(1:(peakposIdx-1))]) # mean max-profiles on the left of Peak value
            rightmean= mean(maxprofiles[c(peakposIdx:no.pmarkers)]) # mean max-profiles on the right of Peak value

            minpos=min(smarkerpos)
            maxpos=max(smarkerpos)
            adjpos= maxpos-(maxpos-minpos)/4
            if ( abs(peakpos-mean(smarkerpos)) < (maxpos-minpos)/8 ){
                # in this case, the peak is in the middle, so we need consider 
                # the mean profiles in both sides of the peak
                legendx = ifelse(leftmean<rightmean, minpos, adjpos)
            }else{
               legendx = ifelse(peakpos-minpos>maxpos-peakpos, minpos, adjpos)
            }
            }

            legendx = findMaxValleyInProfile(maxprofiles, xpos=smarkerpos)


            par(mfrow=c(1, 1), mar=c(5, 5, 3, 2) + 0.1, cex=0.9)
            matplot(x=smarkerpos, y=eQTL_profiles, type=plottype, lty=itype, col = icolor,
                xlab=xlabel, 
                ylab=ylab, main="", axes=showaxis)

            if(!is.null(profilenames)) { 
              if ( max(maxprofiles) < Inf){
                 legend(legendx, max(eQTL_profiles), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
              }else{
                 legend(legendx, max(maxprofiles[maxprofiles < Inf]), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
              }
            }

            if(!showaxis){
               axis(1, at =smarkerpos[sboundaries], labels = smarkerXlabel[sboundaries], las=1.2)
               axis(2,las=2)
            }
            abline(h=minLOD, col="grey")

    }else{ # multiple sub plots
        par(mfrow=c(no.profiles, 1), mar=c(4, 5, 4, 2) + 0.1, cex=0.8)
        for (i in c(1:no.profiles)){            
            #plot(eQTL_profiles[,i], col = "green", main=profilenames[i])

            plot(smarkerpos, eQTL_profiles[,i], type=plottype, col = "green",
                xlab=xlabel, 
                ylab=ylab, main=profilenames[i], axes=showaxis)
            if(!showaxis){
               axis(1, at =smarkerpos[sboundaries], labels = smarkerXlabel[sboundaries], las=1)
               axis(2,las=2)
            }
            abline(h=minLOD, col="grey")
        }
     }

    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if( !is.null(filename)){
        dev.off()
    }
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ A simplfied version of multi curves in one plot $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

plotMultiProfiles=function(eQTL_profiles, xindex, xlabel="",ylabel="", mainlabel="", filename=NULL, plottype="l", ltype=1, myhei = 400, mywid = 500){
    no.profiles  = dim(eQTL_profiles)[2]
    no.pmarkers  = dim(eQTL_profiles)[1]
    profilenames = colnames(eQTL_profiles)
 
    colorcodeC=c("black", "red", "blue","green","grey","yellow","pink","magenta", "purple", "greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow", "grey")

    showaxis=T

    if( !is.null(filename)){
          openImgDev(filename,iwidth = mywid, iheight = myhei)
    }

    #itype=c(1:no.profiles)
    itype=rep(1, no.profiles)
    icolor=colorcodeC[c(1:no.profiles)]

    # simple way to figure out where to put the legend, x position is tricky
    maxprofiles = apply(eQTL_profiles, 1, max)
    if(1==2){
     maxprforders = order(-maxprofiles)
     peakposIdx = maxprforders[1]     # peak index
     peakpos = xindex[peakposIdx] # peak position
     leftmean = mean(maxprofiles[c(1:(peakposIdx-1))]) # mean max-profiles on the left of Peak value
     rightmean= mean(maxprofiles[c(peakposIdx:no.pmarkers)]) # mean max-profiles on the right of Peak value

     minpos=min(xindex)
     maxpos=max(xindex)
     adjpos= maxpos-(maxpos-minpos)/4
     if ( abs(peakpos-mean(xindex)) < (maxpos-minpos)/8 ){
         # in this case, the peak is in the middle, so we need consider 
         # the mean profiles in both sides of the peak
         legendx = ifelse(leftmean<rightmean, minpos, adjpos)
     }else{
         legendx = ifelse(peakpos-minpos>maxpos-peakpos, minpos, adjpos)
     }
    }

     legendx = findMaxValleyInProfile(profile=maxprofiles, xpos=xindex) 

     par(mfrow=c(1, 1), mar=c(5, 5, 3, 2) + 0.1, cex=0.9)
     #colnames(eQTL_profiles)<- c("l")
     
     #matplot(x=xindex, y=eQTL_profiles, type=ltype, lty=itype, col = icolor,
     #matpoints(x=xindex, y=eQTL_profiles, type=ltype, lty=itype, col = icolor,

     if(ltype[1]=="b") {
         matplot(x=xindex, y=eQTL_profiles, type=plottype, lty=ltype, col = icolor,
                xlab=xlabel, ylab=ylabel, main=mainlabel, axes=showaxis, pch="*")
     }else{
         matplot(x=xindex, y=eQTL_profiles, type=plottype, lty=ltype, col = icolor,
                xlab=xlabel, ylab=ylabel, main=mainlabel, axes=showaxis, pch="*")
     }


     #  legend(legendx, max(eQTL_profiles), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
     if ( max(maxprofiles) < Inf){
        legend(legendx, max(eQTL_profiles), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
     }else{
        legend(legendx, max(maxprofiles[maxprofiles < Inf]), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
     }

     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
 
     if( !is.null(filename)){
        dev.off()
     }
}


# input is a globale genotype data structure "datGeno" and 
#         a gene's expression vector across all samples
# output is a vector of the QTL LOD scores on the whole genome
# iqtl[1:10,]
#       chr       pos        lod       AA       AB       BB    sigma
#p45178   1  0.205459 0.04779396 4.255317 4.111057 4.560112 5.634021
#p45787   1  1.872961 0.04778489 4.255144 4.111150 4.560071 5.634021
#p44584   1  6.692455 0.04311972 4.166914 4.149415 4.560018 5.634316
#
qtlOnGenome = function(exprvector){
    datGeno$pheno[,2] <- as.numeric(exprvector)
    iqtl = scanone(datGeno, pheno.col=2, method="hk")
    round(as.numeric(iqtl$lod),2)
}



# plot correlation of each trait pf PCs across modules, different from pervious plot which is based on
# module(all correlations are plot in one figure for each module)
# list_corrlModulePCsTraits[[1]], the first module; list_corrlModulePCsTraits[[2]], 2nd module
#        WeightG      AbFat    OtherFat    TotalFat X100xfat.weight      Trigly
#PC.1 -0.3852039  0.2714077 -0.01301723  0.18511831      0.37228074 -0.17832184
#PC.2 -0.2425847 -0.1467176 -0.09857543 -0.15055786     -0.03087848 -0.08506316
#PC.3  0.3856605 -0.4090452  0.01609531 -0.26830352     -0.45573164  0.32946309
#PC.4 -0.2344006  0.1694582 -0.11794072  0.07000631      0.15279737 -0.06884160

plotCorrlPCsTraitsCrossModule=function(myCorrlModulePCsTraits_list,modulenames, whichPC=1, fieldname=NULL){

        no.modules = length(myCorrlModulePCsTraits_list)
        no.pcs   = dim(myCorrlModulePCsTraits_list[[1]])[1]
        no.traits=dim(myCorrlModulePCsTraits_list[[1]])[2]
        subfigures=4 # number of traits in one figure
        no.figures= as.integer(no.traits/4) + (no.traits%%4>0)

        traitnames = colnames(myCorrlModulePCsTraits_list[[2]])
        pcnames    = rownames(myCorrlModulePCsTraits_list[[1]])

        for ( ifig in c(1:no.figures) ){
            start= (ifig-1)*subfigures + 1   #first trait
            end= ifig*subfigures             #last trait
            if(end>no.traits)
                end=no.traits
            
            actualsubfigures = end - start + 1

            #here we re-organize the correlation data into a matrix with rows as traits and 
            # columns as modules
            corrlMatrix = NULL
            mytitles =  NULL
            for (j in c(start:end) ){ #traits as rows
                rvect=NULL
                for(m in c(1:no.modules) ) #module as columns
                   rvect=c(rvect, (myCorrlModulePCsTraits_list[[m]])[whichPC, j] )
                corrlMatrix = rbind(corrlMatrix, rvect)
                mytitle = paste("r(", traitnames[j], " , ", pcnames[whichPC], ")", sep="")
                mytitles = c(mytitles, mytitle)    
            }

            if( !is.null(fieldname)){
              imgfname=paste(fieldname, "_corrl_PC", as.character(whichPC), 
                             "_Traits", as.character(start), "to", as.character(end), ".png", sep='')
              openImgDev(imgfname,iwidth = 1278, iheight = 1400)
            }
            par(mfrow=c(actualsubfigures,1), mar=c(6, 5, 4, 2) + 0.1, cex=0.8)

            #display corelations of PC with trait in each module
            for (j in c(1:actualsubfigures) ){               
               barplot(as.numeric(corrlMatrix[j,]), ylab="correlation",                       
                 col = modulenames, main=mytitles[j], names.arg=modulenames, las=2)
                 abline(h=0.3, col="grey")
                 abline(h=0, col="black")
                 abline(h=-0.3, col="grey")
            }
            par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
            if( !is.null(fieldname)){
              dev.off()
            }
      }
}



# plot Gene Significance across modules,  a single subfigure contains the means of a trait in
# all modules, a figure will have 4 sub figures
#> genesig_means[1:5,1:5]
#         WeightG     AbFat   OtherFat   TotalFat X100xfat.weight
#black 0.33819315 0.2371843 0.04326499 0.16219943       0.3253637
#blue  0.51582538 0.1813136 0.19555568 0.12881716       0.2365684
#brown 0.62463205 0.1791032 0.24054805 0.11426613       0.2576446
#cyan  0.26215453 0.1969778 0.11479211 0.10464036       0.1801414
#green 0.08332405 0.0964498 0.09881134 0.07656996       0.1036234
plotMeanGeneSignifCrossModule=function(means_matrix, sderrs_matrix, fieldname=NULL, signifiance_level=0.3, imgwid=1278){

        no.modules = dim(means_matrix)[1]
        no.traits  = dim(means_matrix)[2]
        subfigures=4 # number of traits in one figure
        no.figures= as.integer(no.traits/4) + (no.traits%%4>0)

        traitnames = colnames(means_matrix)
        modulenames= rownames(means_matrix)

        for ( ifig in c(1:no.figures) ){
            start= (ifig-1)*subfigures + 1   #first trait
            end= ifig*subfigures             #last trait
            if(end>no.traits)
                end=no.traits
            
            actualsubfigures = end - start + 1

            if( !is.null(fieldname)){
              imgfname=paste(fieldname, "_genesignif",
                             "_Traits", as.character(start), "to", as.character(end), ".png", sep='')
              openImgDev(imgfname,iwidth = imgwid, iheight = 1400)
            }
            par(mfrow=c(actualsubfigures,1), mar=c(6, 5, 4, 2) + 0.1, cex=0.8)

            #display trait-based gene significance in each module
            for (j in c(start:end) ){               
               barplot(as.numeric(means_matrix[,j]), ylab="mean(|r|)",                       
                 col = modulenames, main=traitnames[j], names.arg=modulenames, las=2)
                 abline(h=signifiance_level, col="grey")
                 abline(h=0, col="black")
               err.bp(as.vector(means_matrix[,j]), as.vector(sderrs_matrix[,j]), two.side=T)
            }
            par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
            if( !is.null(fieldname))
              dev.off()
      }
}


#********************************************************************************************************
# the corrlmatrix is the correlations between a few network properties (as row) 
#   and a trait (traitname) across modules (column)
# 
#
plotCorrelOfNetiesAndTrait_AcrossModules=function(mycorrlmatrix, mytraitname, fieldname=NULL){

        no_neties   = dim(mycorrlmatrix)
        #corrlmatrix = abs(mycorrlmatrix)
        corrlmatrix = mycorrlmatrix

        traitname = replaceChars(mytraitname, ".", "")

        #mycolors=heat.colors( dim(corrlmatrix)[1])
        mycolors=c("lightgreen","blue", "red", "yellow")
        mytitle =paste("correlations between gene network properties and gene significance based on ", traitname, sep="")

        if( !is.null(fieldname)){
          netiesTraitfname=paste(fieldname, "_genesignif_Neties_", traitname, ".png", sep='')         
          openImgDev(netiesTraitfname,iwidth = 1278, iheight = 768)
        }
        par(mfrow=c(1,1), mar=c(11, 6, 4, 2) + 0.1, cex=1)
        mp  <- barplot(corrlmatrix)
        tot <- colMeans(corrlmatrix)
        text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
        barplot(corrlmatrix, beside = TRUE, ylab="correlation", las=2, 
             col =mycolors, 
             main=mytitle,
             legend = rownames(corrlmatrix) )
        abline(h=0, col="black")
        abline(h=0.3, col="grey")
        abline(h=-0.3, col="grey")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        if( !is.null(fieldname)){
           dev.off()
        }
}



#count LOD scores on genome for all the genes in a module
LODscoreCOUNT=function(TRAITS, LODscoreCut=2.5, title1="eQTL profile"){
   QTLlist=list(NULL)

   maxLOD = 0

   datGeno$pheno=data.frame( datGeno$pheno[,1], TRAITS)
   no.traits= dim(datGeno$pheno)[[2]]

   QTLlist=scanone(datGeno, pheno.col=2, method="hk")

   maxLODmarker=data.frame(matrix(-666, nrow=no.traits-1, ncol= length(QTLlist$pos)  +1) )
   chromo=QTLlist$chr
   names(maxLODmarker)[[1]]="Trait"

   names(maxLODmarker)[-1]=paste(names(maxLODmarker)[-1],"chr",chromo,sep=".") 
   maxLODmarker[,1]=names(TRAITS)
   for (i in c(1:no.traits) ) {
      datGeno$pheno[,2] <- as.numeric( pcs_list[[i]][,moduleIndex] )           
      QTLlist=scanone(datGeno, pheno.col=i, method="hk");
      maxLODmarker[i-1,-1]=QTLlist$lod
   }

   count1=apply( maxLODmarker[,-1]> LODscoreCut, 2,sum)
   myylab=paste("no. of LODs above ", as.character(LODscoreCut), sep="")

   plot(count1,col=as.character(chromo), main=title1,xlab="marker number", ylab=myylab)
   list(maxLODmarker =maxLODmarker, chromo=chromo)

}

#here, we attach the information in the minor Matrix to the major matrix, if some row keys 
#in the major matrix don't match any in the minor matrix, we use missing label for those keys
# keepAllPrimary=T: return all primary row-elements
# keepAllPrimary=F: return all MISSED primary row-elements
# keepAllPrimary=NULL: return the common elements
mergeTwoMatricesByKeepAllPrimary=function(primaryMatrix, minorMatrix, missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
{
  no.promarycols <- dim(primaryMatrix)[2]
  no.mustbegenes <- dim(primaryMatrix)[1]

  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  keyword="mustbeused"
  mustbeGenesMatrix = cbind(primaryMatrix, c(1:no.mustbegenes), rep(keyword, no.mustbegenes) )

  if (is.null(colnames(primaryMatrix)) ){
    colnames(mustbeGenesMatrix) <- c( c(1:no.promarycols), "primorder", keyword)
  }else{
    colnames(mustbeGenesMatrix) <- c( colnames(primaryMatrix), "primorder", keyword)
  }
  dim(mustbeGenesMatrix)

  if(is.null(keepAllPrimary) ){ #normal merge: to have the common elements
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=F,sort=F,all=F)
  }else{
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=T,sort=F,all=T)
  }

  dim(myMatrix)
  nocols.mymatrix <- dim(myMatrix)[2]

  #the mustbeused genes which are not included in minor have NAs in the column $mustbeused
  #so we can use this information to figure out which mustbeused genes missing in minorMatrix
  myMatrix[,nocols.mymatrix] = ifelse( is.na(myMatrix[,nocols.mymatrix]), missinglabel, as.character(myMatrix[,nocols.mymatrix]) )

  orders = order( as.numeric(as.matrix(myMatrix[, no.promarycols+1])))
  if (keepPrimaryOrder)
      myMatrix = myMatrix[orders,]

  if (is.null(keepAllPrimary) ){
     selected = rep(T, dim(myMatrix)[1])
  }else{
     if (keepAllPrimary)
       selected = !(is.na(myMatrix[, no.promarycols+2]))
     else #return the row-elements in minor which are missed in primary
       selected = is.na(myMatrix[, no.promarycols+2])
  }

  sum(selected)

  #keep the primary matrix and remove the mustbeused column
  myMatrix[selected, -c(no.promarycols+1, no.promarycols+2)]
}



# Here we return the boolean vactor which indicates whether each element in the primary list is
# in theminor list, the boolean vector is in the same order as the primary list
#
findListInPrimarysPosition=function(primaryMatrix, minorMatrix)
{
  #no.mustbegenes <- dim(minorMatrix)[1]
  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  #keyword="mustbeused"
  #mustbeGenesMatrix = cbind(minorMatrix, rep(keyword, no.mustbegenes) )
  #colnames(mustbeGenesMatrix) <- c( colnames(minorMatrix), keyword)
  #dim(mustbeGenesMatrix)
  #It is always dangerous to use merge if you need keep the selected.bool.vector in the same order 
  # as primaryMatrix, whatever, the order will be changed after merging
  #
  #myMatrix = merge(mustbeGenesMatrix, primaryMatrix, by.x=1, by.y=1,all.y=T, sort=F,all=F)
  #dim(myMatrix)
  #nocols.mymatrix <- dim(myMatrix)[2]
  #selected = ifelse(is.na(myMatrix$mustbe), FALSE, TRUE)
  #sum(selected)

  minorlist  = as.character( as.matrix(minorMatrix[,1]) )
  primarylist= as.character( as.matrix(primaryMatrix[,1]) )

  no.primary = dim(primaryMatrix)[1]
  selected=rep(FALSE, no.primary)
  primaryseq = c(1:no.primary)

  for (each in minorlist){
      isel = (each==primarylist)
      if (sum(isel) >0){
        selected[ primaryseq[isel] ] = TRUE
	#cat(sum(isel),primaryseq[isel] , "\n")
      }
  }
  
  selected
}

# re-organize the QTL in the gene order 
#
# turn gene name into gene symbols
# 
# 8	chr8 1.3e+08	7	709	25	11466	6.00E-04	NM_018412;NM_024525;NM_024604;NM_015990;Contig2399_RC;NM_022051;NM_021205
#
# called by C:\Zhang\WntSignaling\HypoxiaNew_Oct20\R-decompose_qtl-tables.R
#
expandQTLenrich_by_genes = function(fqtlenrich, modulecolname="LocusBin", 
                                 genecolname="Probe.Set", geneinfo, signat_source, fout)
{
    allMatrix<- read.delim(fqtlenrich, sep="\t", header=T)
    allMatrix<- as.matrix(allMatrix)
    dim(allMatrix)

    no.rows = dim(allMatrix)[1]

    acolnames = colnames(allMatrix)
    aindices  = getMatchedIndex(cvector=acolnames, 
                                subvect=c(modulecolname, genecolname) )
    midx      = aindices[1]
    gidx      = aindices[2]

    # to shorten the gene information matrix
    #
    allgenes  =splitString(allMatrix[, gidx], ";")

    # remove the duplicates
    #
    uniques   = names(table(allgenes))
    geneinfo2 = merge( cbind(uniques),geneinfo, by.x=1, by.y=1, all=F)
    
    for(i in c(1:no.rows)){
        imodule  = allMatrix[i, midx]
        iallgenes= allMatrix[i, gidx]
        igenes  =splitString(iallgenes, ";")

        # remove the duplicates
        #
        iuniques= names(table(igenes))
        imerged = merge( cbind(iuniques),geneinfo2, by.x=1, by.y=1, all=F)
        imerged = as.matrix(imerged)
        irows   = dim(imerged)[1]

        final   = cbind(imerged[,2], imerged[,1], 
                        rep(imodule,irows), rep(signat_source, irows) )

        write.table(final, fout, sep="\t",quote=FALSE, col.names=F, 
                    row.names=FALSE, append=T)
    }

    rm(geneinfo2)        
}



#--------------------------------------- Steve -----------------------------------------------------------
#--------------------------------------- Steve -----------------------------------------------------------

#The function ScaleFreePlot1 creates a plot for checking scale free topology
ScaleFreePlot = function(kk,no.breaks=15, mtitle="",truncated1=TRUE, tofile="", outputFitness=FALSE){
    
    cut1=cut(kk,no.breaks)
    binned.k=tapply(kk,cut1,mean)
    freq1=tapply(kk,cut1,length)/length(kk)

    #remove NAs
    noNAs = !(is.na(freq1) | is.na(binned.k))
    freq.noNA= freq1[noNAs]
    k.noNA  = binned.k[noNAs]

    #remove Zeros
    noZeros  = !(freq.noNA==0 | k.noNA==0) 
    freq.log = as.numeric(log10(freq.noNA[noZeros]))
    k.log    = as.numeric(log10(k.noNA[noZeros]))

    if(tofile!=""){
       openImgDev(tofile, iwidth = 700, iheight = 700, ipointsize = 10)
    }

    par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1, cex=1.4)

    plot(k.log, freq.log, xlab="log10(k)",ylab="log10(p(k))" )
    lm1=lm( freq.log ~ k.log, na.action=na.omit) 
    lines(k.log,predict(lm1),col=1)
    if (truncated1==TRUE) { 
       lm2=lm(freq.log ~ k.log +I(10^k.log) );
       lines(k.log,predict(lm2),col=2);

       if (mtitle !=""){
           ititle=paste(mtitle, ", scale R^2=",
                          as.character(round(summary(lm1)$adj.r.squared,2)) , 
                          ", trunc. R^2=",
                          as.character(round(summary(lm2)$adj.r.squared,2)),
                          ", slope=", as.character(round(lm1$coef[[2]], 2)),
                          sep="")
       }else{
           ititle=paste("scale R^2=",
                          as.character(round(summary(lm1)$adj.r.squared,2)) , 
                          ", trunc. R^2=",
                          as.character(round(summary(lm2)$adj.r.squared,2)),
                          ", slope=", as.character(round(lm1$coef[[2]], 2)),
                          sep="")
       }

    }else{
       if (mtitle !=""){
            ititle = paste(mtitle, ", scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)), 
                     ", slope=", as.character(round(lm1$coef[[2]], 2)), sep="")
       }else{
            ititle = paste("scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)), 
                     ", slope=", as.character(round(lm1$coef[[2]], 2)), sep="")
       }
    }
    title(ititle)

    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if(tofile !=""){
       dev.off()
    }
        
    if(outputFitness==TRUE){
      if(truncated1){
        myoutput = c(summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]])
      }else{
        myoutput = c(summary(lm1)$adj.r.squared, lm1$coef[[2]])
      }
    }else{
      myoutput=ititle
    }

    myoutput
} # end of function

# The function VariableRelationship visualizes the relationships between the columns of a data frame.
# It is an MDS plot on the basis of the Spearman correlations between the variables
# the relationship between variables
VariableRelationshipPlot = function(dath, fname, colorcode,myimageType="png")
{
   globalfname = paste(fname, "_allmodules.",myimageType, sep="")
   openImgDev(globalfname)
   mds=cmdscale(1-abs(cor(dath, use="p",method="s")),k=2)
   plot(mds,type="n" ,main="Across Modules")
   text(mds,label=names(dath))
   dev.off()

   colorlevels = names(table(colorcode))
   for (icolor in  colorlevels){
       if (icolor=="grey"){
          next
       }
       sel = as.character(colorcode)==as.character(icolor)
       modulefname = paste(fname, "_", as.character(icolor), ".",myimageType,    sep="")
       mtitle      = paste("Inside ",  as.character(icolor), " Module", sep="")
       openImgDev(modulefname)
       plot(varclus(as.matrix(dath[sel,])), main=mtitle, xlab=mtitle)
       dev.off()
   }
}

#System        /GeneCategory      /ListHits/ListTotal/PopulationHits/PopulationTotal/EASE score
#BiologicalProc/mitotic cell cycle/54      /166      /305           /10642          /7.23e-042
# a=rbind( c(305-54, 10642-251-112-54), c(54, 112) )
# a
# signif(fisher.test(a)$p.value,2)
#[1] 4.1e-43
# signif(chisq.test(a)$p.value,2)
#
# pop_list_hits: in order of PopulationTotal, PopulationHits, List Total, List Hits
#             Category           Non-Category
#Population     a[1,1]             a[1,2]
#List           a[2,1]             a[2,2]
#
#
#************ example, Fisher exact Test & Hypergeometric test ********************
#
# 196 of 673 genes overlap the 948 genes, total geens 24000
#a=rbind( c(673-196, 24000-(673-196)-(948-196)-196), 
#          c(196,     948-196))
#
#  signif(fisher.test(a)$p.value,2)
# hypergeometric test: 
#    phyper(196, 673, 24000-673, 948, , lower.tail=F)
#
# By chance you expect 196/23756*948 = 8 overlap
#
# if population hits < minPopulationHits, then the test doesn't continue
#
#input vector:=[popilation, population-hit, list, list-hit]
fisherTest = function(poplisthits, minPopulationHits=5)
{
   if ( (poplisthits[2]<minPopulationHits) | (poplisthits[4]==0) )
      return (1)

   #old way, maybe wrong
   #a=rbind( c(poplisthits[2]-poplisthits[4], poplisthits[1]-poplisthits[2]-poplisthits[3]+poplisthits[4]), 
   #         c(poplisthits[4],           poplisthits[3]-poplisthits[4] ) )
   #cat(as.character(poplisthits[1]), as.character(poplisthits[2]),as.character(poplisthits[3]),as.character(poplisthits[4]),"\n")
   #signif(fisher.test(a)$p.value,3)
   q= poplisthits[4] #list hit
   m= poplisthits[2] #population hit
   n= poplisthits[1]-poplisthits[2] #population non-hit
   k= poplisthits[3] #list total

   myp=1-phyper(q-1, m, n, k)   
   signif(myp,3)
}

# 11 out of 14 overlap 38, total 2588
fisherExactTest=function(overlap, poolA, poolB, total){
   a=rbind( c(poolA-overlap, total- (poolA-overlap)-(poolB-overlap)-overlap),
            c(overlap,       poolB-overlap)
     )
   signif(fisher.test(a)$p.value,2)
}

forcenumericSum=function(myvector){
   sum( as.numeric(as.matrix(myvector)), na.rm=T)
}


tapplySumFunc=function(myvector, mycolor){
   restab=tapply(as.numeric(as.matrix(myvector)), mycolor, sum)
   as.matrix(restab)
}

#imput is a boolean vector and another vector whose components are to be selected
# here we use a trick of temporary file
applySelectFunc=function(selvector0or1, myvector, removeduplicates=F){
   selvector= as.integer(as.matrix(selvector0or1)) & TRUE
   if(sum(selvector)<=0)
      return (" ")

   tmpfn="tmp.txt"

   selcomp = as.character(myvector[selvector])
   if (removeduplicates){
     selgenes = sort( as.character(myvector[selvector]) )
     no.s =length(selgenes)
     selgenes.shift = c("0000000", selgenes[c(1:(no.s-1))])
     unique.bool = (selgenes.shift != selgenes)
     selcomp  = selgenes[unique.bool]
   }

   write.table( t(as.matrix(selcomp)), tmpfn, sep=";",quote=FALSE, col.names=F, row.names=FALSE)
   mystring <- read.delim(tmpfn,sep="\t", header=F)
   as.character(as.matrix(mystring))
}



# The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher's exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
fisherPvectorOld=function(couleur, annotation1, minNumberAnnotation=10) {

    levelsannotation1 =levels(annotation1)
    levelscouleur     =levels(couleur)
    no.couleur        =length(levelscouleur)
    restannotation1   =table(annotation1)>minNumberAnnotation
    no.annotation     =sum( restannotation1)

    datout=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )

    names(datout)        =levelscouleur
    restlevelsannotation1= levelsannotation1[restannotation1]
    row.names(datout)    = restlevelsannotation1

    for (i in c(1:no.annotation) ) {
      for (j in c(1:no.couleur) ){
        tab1=table( annotation1==restlevelsannotation1[i], couleur==levelscouleur[j])
        datout[i,j]=signif(fisher.test(tab1)$p.value,2)
    }}
    datout

} # end of function fisherPvector


#outformat="pvalue":     output pvalue table
#outformat="ratio":      output tables of ratios of two proportions
#outformat="both":       output pvalue(proportion) table
fisherPvector =function(couleur, annotation1, outformat="pvalue", minNumberAnnotation=10) {

    levelsannotation1 =levels(annotation1)
    levelscouleur     =levels(couleur)
    no.couleur        =length(levelscouleur)
    restannotation1   =table(annotation1)>minNumberAnnotation
    no.annotation     =sum( restannotation1)

    datout=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )

    names(datout)        =levelscouleur
    restlevelsannotation1= levelsannotation1[restannotation1]
    row.names(datout)    = restlevelsannotation1

    annotation1.chars = as.character(annotation1)

    for (i in c(1:no.annotation) ) {
      for (j in c(1:no.couleur) ){

        bool.annoted= (annotation1.chars==restlevelsannotation1[i])
        bool.colorj = (couleur==levelscouleur[j])

        # some genes are not annoted, so we need remove those from our proportion computation
        bool.annoted.noNAs = bool.annoted[ !is.na(bool.annoted) ]
        bool.colorj.noNas  = bool.colorj[  !is.na(bool.annoted) ]

        tab1=table(bool.annoted.noNAs, bool.colorj.noNas)
        pvalue = signif(fisher.test(tab1)$p.value,2)

        #                    rest module
        #bool.annoted.noNAs FALSE TRUE
        # nofunc       FALSE 3257   292
        #   func       TRUE    29     1
        portion.rest   = tab1[2,1]/(tab1[1,1]+tab1[2,1]) # portion of 
        portion.colorj = tab1[2,2]/(tab1[1,2]+tab1[2,2]) #=length(bool.colorj.noNas)

        #proportion
        ratioOFportions   = portion.colorj/portion.rest
        
        if (outformat=="pvalue"){
            datout[i,j]=pvalue
        }else if( outformat=="ratio" ){
            datout[i,j]=signif(ratioOFportions,2)
        }else{
            myout = paste(as.character(pvalue),"(", 
                          as.character(signif(portion.colorj,2)),",",
                          as.character(signif(portion.rest,2)),  ")", sep="")
            datout[i,j]= myout
        }

    }}
    datout

} # end of function fisherPvector





# this function computes the standard error
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }

choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1
}


####  Function pamPsClNo computes prediction strength and returns the estimated number of clusters  
#### based on PAM clustering
                                        #k=6  at most k clusters
                                        #v=2 fold cross validation
                                        #c=5  no. of cross validations
pamPsClNo <- function(original.data,
                      klimit=5,  # max no. of clusters
                      cvFold=2,  # how many fold cross validation
                      cvCount=5, # how many cross validations
                      m=2,       # number of points
                      double1=FALSE, diss=FALSE,
                      cut1=1){
  if (double1) {original.data <- rbind(original.data,original.data) }
  clustlearn <- list(NULL);
  clusttest  <- list(NULL);
  nData <- nrow(original.data);
  cps1 <- matrix(1,nrow=klimit,ncol=cvCount)
  criterion1 <- matrix(1,nrow=klimit,ncol=cvCount) 
  if (diss) {
    alldist <- as.matrix(original.data)
  } else {
    alldist <- as.matrix(dist(original.data))
  }
  alldist <- signif(alldist,10)
  ## for each cross-validation set
  for (cc in 1:cvCount) {
    ## two matrices used to store prediction strength calculated by
    ## four kinds of metrics
    ps1 <- matrix(1, nrow=klimit, ncol=cvFold)

    ## utility vector indicating split of data set into cvFold sets 
    rest1 <- nData-as.integer(nData/cvFold)*cvFold
    sam1 <-  sample(c(rep(c(1:cvFold), as.integer(nData/cvFold)),
                      sample(1:cvFold, rest1)))

    ## for each possible number of clusters,
    for (kk in 2:klimit){
      ## cvFold fold splitting for cross validation
      for (vv in 1:cvFold){

        ## indices of test and training sets
        test.index <- c(1:nData)[sam1==vv]
        learn.index <- c(1:nData)[sam1!=vv]
        no.test <- length(test.index)
        no.learn <- length(learn.index)

        if (no.test <= kk || no.learn <= kk) {
          ## clustering too few points into too many clusters
          ps1[kk,vv] <- 0
          next
        }
        ## distances between points in test and training sets
        test.dist <- alldist[test.index, test.index]
        learn.dist <- alldist[learn.index, learn.index]

        ## perform clusterings on test and training sets
        ##print(paste("Performing pam with kk=", kk, "no.test", no.test, "no.learn", no.learn))
        clustlearn[[kk]] <- pam(as.dist(learn.dist), kk, diss=T)
        clusttest[[kk]] <-  pam(as.dist(test.dist), kk, diss=T)
        
        ## this assigns a cluster to each test set observation, based on the
        ## clustering of the training set.
        Cl <- rep(0, no.test)
        d <- rep(10000000, kk) #difference matrix for assigning Cl
        
        ## determine which medoid each test set point is closest to/i.e.,
        ## which cluster it belongs to
        index1 <-  clustlearn[[kk]]$medoids # length is kk
        for (i in 1:no.test){
          for (j in 1:kk){
            d[j] <- alldist[index1[[j]], test.index[i]]
            ## note: this assumes that the medoids are in the original dataset
          }
          mincluster <- c(1:kk)[rank(d) == min(rank(d))]
          if (length(mincluster) == 1) {
            Cl[i] <- mincluster
          }
          else if (length(mincluster)>1) {
            Cl[i] <- sample(mincluster, 1)
          }
        }  # end for for over i in 1:no.test
        
        ## now we compute how often m samples are co-clustered
        tab0 <- table(Cl, clusttest[[kk]]$clustering)
        pshelp <- rep(10000000, kk)
        for (l in 1:kk){
          ## marginals
          tab1=tab0
          tab1[tab0<cut1]=0
          nkl <- sum(tab1[,l])
          if (nkl < m)  {
            pshelp[l] <- 0
          } else { 
            pshelp[l] <- sum(choosenew( tab1[,l], m ) )/ choosenew(nkl,m)     
          }
        } # end of  for l in 1:kk
        ps1[kk,vv] <- min(pshelp)              # Min
      } # end of vv in 1:cvFold
    } # end of kk in 2:klimit
    cps1[,cc] <- apply(ps1,1,mean)
    ## gives max over vv for cvFold=2
    criterion1[,cc] <-  cps1[, cc] + apply(ps1, 1, stderr1)
  } # end of for cc in 1:cvCount
  psse <- signif(apply(criterion1, 1, mean), 3)
  ##psse <- signif(apply(cps1, 1, mean), 3)
  psse
  ##print(psse)
  ##thres <- 0.8-(m-2)*0.05 
  ##for(nn in 1:klimit) { 
  ##  if(psse[nn]>=thres) {clusterNo <- nn}
  ##}
  ##clusterNo
} # end of function pamPsClNo


kpredictPS2=function(ps,span1=0.75,title1="prediction strength") {
    kk=c(1:length(ps))
    kout=1
    if (length(kk)>1) { 
        lm1=loess(ps~ kk,span=span1,degree=2)    
        rank1=rank(round(resid( lm1 ),digits=10))
        kouthelp=  kk[   rank1==max(rank1) ]
        kestimate=kouthelp[length(kouthelp)]
        plot(ps,main=title1,xlab="no. of clusters k",ylab="PSE",ylim=c(0,1))
        lines(kk,predict(lm1))
        abline(v=kestimate,col="red")
        abline(h=.8)
        points(kk,resid(lm1) ,pch="*",col="blue")
        kestimate
    }
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##     Function: Make a plot, run a couple tests and output everything to proper files 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scatterCorrelTester=function(rawValA, rawValB, binaryvar=FALSE, PlotsTitle ="", myxlab="",myylab="",outputImage="tmp.png", outputFile=tests_Fname, numPoints=20 ){

        if( (sd(rawValA, na.rm = T)==0) | (sd(rawValB, na.rm = T)==0) ){
            titleinfor = paste(PlotsTitle, "error in scatterCorrelTester, sd()=0", paste="")
            appendListToFile(outputFile, titleinfor) 
            return
        }

	#in general, we will have the thing being tested be rawValA, and K be rawValB
	cut_factor = cut(rawValB,numPoints)

	# Now give me the mean proportion of essentiality in each chunk of the essentiality vector 
	mean_rawValA = tapply(rawValA, cut_factor, mean)

	# Now give me the mean proportion of k in each chunk of the essentiality vector 
	mean_rawValB = tapply(rawValB, cut_factor, mean)

	#now we need a p-value for this:
        corrl   = signif(cor(rawValB,rawValA),2)
        #if (binaryvar==FALSE){
           corrl.p = signif( cor.test(rawValB,rawValA, method="p",use="p")$p.value ,2)
        #}else{
	#   kruskal = kruskal.test(rawValB,rawValA)
        #   corrl.p = signif(kruskal$p.value, 2)
        #}
        titleinfor = paste(PlotsTitle, ": r=", as.character(corrl),",p=", as.character(corrl.p), paste="")

	#Now we need to just make a line fit
	line =lm(mean_rawValA ~ mean_rawValB)
	sum1 = summary(line)
	sum1

	#Now put in a title for this part of the data:
	appendListToFile(outputFile, titleinfor)
	appendListToFile(outputFile, sum1)

	# do a plot() of these two vectors against each other. 
	openImgDev(outputImage)
        par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
	plot(mean_rawValB,mean_rawValA, main=titleinfor, xlab=myxlab, ylab=myylab)
	# and then draw the line
	abline(line, col="red")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
	dev.off()
}

conservationVectorCleaner = function(rawVal, nasVal) {

	# I need a way to always filter out the NAs from the mean_e_val vector...
	no.nas = nasVal
	# variable log_mean_e_val must be defined for use in plotting conservation plots
	mean_e_valNoNA = rawVal[no.nas]
	# 2nd I need to replace the 0 values with extremely small numbers (smaller than you see so e-200)
	mean_e_valNoNANorZero = ifelse(mean_e_valNoNA==0, 10^(-200), mean_e_valNoNA)
	# 3rd I need to transform the e value since its range is just HUGE
	log(mean_e_valNoNANorZero)
	#R just returns the last "thing" listed in a function
}


# the following function computes the Rand index between 2 clusterings
# assesses how similar to clusterings are
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {
        out1[i] <- 0
    }else {
        out1[i] <- choose(n[i], k)
    }
  }
  out1
}

Rand <- function(itab,adjust=T) {
  no.rows <- nrow(itab)
  no.cols <- ncol(itab)
    tab  =  matrix(as.numeric(as.matrix(itab)), no.rows, no.cols)

     a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
     n <- nrow(tab)
     
     for (i in 1:n) {
          for(j in 1:n) {
               a <- a+choosenew(tab[i,j],2)
               nj <- sum(tab[,j])
               c <- c+choosenew(nj,2)
          }
          ni <- sum(tab[i,])
          b <- b+choosenew(ni,2)
          nn <- nn+ni
     }

     if(adjust==T) {
          d <- choosenew(nn,2)
          adrand <- (a-(b*c/n)/d)/(0.5*(b+c/n)-(b*c/n)/d)
          adrand
     } else {
          b <- b-a
          c <- c/n-a
          d <- choosenew(nn,2)-a-b-c
          rand <- (a+d)/(a+b+c+d)
          rand
     }
}

combinBy2=function(a)
{
  a*(a-1)*0.5
}

#based on the paper: http://faculty.washington.edu/kayee/pca/supp.pdf
RandIndex= function(tab)
{
  no.rows <- nrow(tab)
  no.cols <- ncol(tab)
  
  A  =  matrix(as.numeric(as.matrix(tab)), no.rows, no.cols)
  ni <- apply(A, 1, sum) #sum over columns for each row
  nj <- apply(A, 2, sum) #sum over rows for each column
  n  =  sum(ni)

  n.comb2   = combinBy2(n)
  ni.comb2  = combinBy2(ni)
  nj.comb2  = combinBy2(nj)
  nij.comb2 = combinBy2(A)

  sum.ni.comb2 = sum(ni.comb2)
  sum.nj.comb2 = sum(nj.comb2)

  sum.nij = sum( apply(nij.comb2, 1, sum) )
  
  a= sum.nij - sum.ni.comb2*sum.nj.comb2/n.comb2
  b= 0.5*(sum.ni.comb2+sum.nj.comb2) - sum.ni.comb2*sum.nj.comb2/n.comb2
  rand=a/b
  abs(rand)
}


# This function computes the GTOMm DISSIMILARITY
#
# Input:
# - adjmat1, a symmetric adjacency matrix with binary entries
# - m, the order of GTOM
#
# Output:
# - The GTOMm dissimilarity matrix
#
# Andy M. Yip and Steve Horvath
# January 2005

#if(exists("GTOMmdist1")) rm(GTOMmdist1);
GTOMmdist1 = function(adjmat1,m=1){
    if (m!=round(abs(m))){
        stop("m must be a positive integer!!!", call.=TRUE);}
    if (any(adjmat1!=0 & adjmat1!=1)){
        stop("The adjacency matrix must be binary!!!", call.=TRUE);}

    B <- adjmat1;
    if (m>=2) {
        for (i in 2:m) {
            diag(B) <- diag(B) + 1;
            B = B %*% adjmat1;}}   # number of paths with length at most m connecting each pair
    B <- (B>0);                    # m-step reachability matrix
    diag(B) <- 0;                  # exclude each node being its own neighbor
    B <- B %*% B;                  # number of common k-step neighbors that each pair of nodes share

    Nk <- diag(B);                 # number of common k-step neighbors that each node possesses
    B <- B +adjmat1;
    diag(B) <- 1;
    denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
    diag(denomTOM) <- 1;
    1 - B/denomTOM                 # turn the GTOM matrix into a dissimilarity
}

#ModulePrincComps, which computes the first principal component of each module. 
#Also it provides the variance explained for the first 5 PCs of a module.
#It is based on the singular value decomposition.
#It takes as input datExpr  for which the rows are samples and the columns are genes.
#Here is how you would use it
#PC1=ModulePrinComps2[datExpr,color1]
#Then PC1[[1]] provides a data frame with the first PC of each module
#PC1[[2]] provides a data frame where each column contains the percentages of 
#the variance explained by the first 10 PCs by module.

ModulePrinComps.old = function(datexpr,couleur) {
  no.pcs=10

  modlevels=levels(factor(couleur))

  PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  PrinComps2=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  varexplained= data.frame(matrix(666,nrow=no.pcs,ncol= length(modlevels)))

  colnames(PrinComps)    <- modlevels
  colnames(PrinComps2)   <- modlevels
  colnames(varexplained) <- modlevels

  for(i in c(1:length(modlevels)) ){
    #print(i)   
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    datModule=impute.knn(as.matrix(datModule))$data

    datModule=t(scale(t(datModule)))
    svd1=svd(datModule)
    mtitle=paste("PCs of ", modulename," module", sep="")
    varexplained[,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)

    # this is the first principal component
    pc1=svd1$v[,1]
    signh1=sign(sum(cor(pc1,  t(datModule))))
    if (signh1 != 0)  pc1=signh1* pc1

    # this is the second principal component
    pc2=svd1$v[,2]
    signh2=sign(sum(cor(pc2,  t(datModule))))
    if (signh2 != 0)  pc2=signh2* pc2

    PrinComps[,i]= pc1
    PrinComps2[,i]= pc2
  }

  list(PrinComps, PrinComps2, varexplained)
}

#============================================================
# adj is the input adjacency matrix. It should be symmetric and all its
#diagonal elements should be 1. Do NOT assign zeros to the diagonal entries
#as before!
#
# kpred=T will output the predicted within-module connectivity.
#
# Umat=T will output the U matrix in the Spectral Decomposition. 
# CAUTION: U is an N*N matrix which will eat up a lot of memory and CPU when N is big
#  (such as N=3000).
#
# D=T will output the eigenvalues of the matrix adj.
#
# In default, the output includes Module Size, Module Cohesiveness, Weight,
#   Module Conformity and Within-module connectivity.

SDADJ1 = function(inadj, kPred=F, Umat=F, D=F) 
{
    # Check if adj is a valid adjacency matrix:  square matrix, positive
    #entries, symmetric and diagonal elements being 1.
    adj=inadj
    diag(adj) <- 1

    n=dim(adj)
    tol=0.000000001

    if ( n[1] != n[2])               stop("The adjacency matrix is not a square matrix!")
    if ( sum(is.na(adj))>0 )         stop("There are missing values in the adjacency matrix!")
    if ( sum(adj<0)>0 )              stop("There are negative entries in the adjacency matrix!")
    if ( max(abs(adj-t(adj))) > tol) stop("The adjacency matrix is not symmetric!")
    if ( max(abs(diag(adj)-1))> tol) warning("The diagonal elements are not all one!")

    # Calculate kWithin, MCH, MCF, and/or kPred, Umat, D
    n=n[1]
    sd=eigen(adj) # Spectral Decomposition

    kWithin=apply(adj, 2, sum)                # Within Module Connectivities
    MCohesiveness=sd$values[1]/sum(sd$values) # Module Cohesiveness

    MConformity=abs(sd$vectors[,1])*sqrt(n) # Module Conformity

    Weight= sqrt(MCohesiveness*sum(kWithin))

    output=list(MSize=n,MCohesiveness=MCohesiveness,Weight=Weight,MConformity=MConformity,kWithin=kWithin)

    if (kPred) {
      kPred= Weight* MConformity
      output$kPred=kPred
    }

    if (Umat) output$Umat=sd$vectors
    if (D)    output$D=sd$values
    #output

    MConformity
}
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================

# whichway="row" : cluster by "column","row", or "rowcolumn"
# clusterno.row  : <=0 no row-wise cluster detection, >0 identify only the specified number of clusters
# clusterno.col  : <=0 no col-wise cluster detection, >0 identify only the specified number of clusters

#************** Hierarchical Clustering  ********************************************

## whichway = "row"/"column"/"rowcolumn"
hiercluster2way=function(datMatrix, whichway="row",clusterno.row=2, clusterno.col=3, imagename=NULL,showlabels=FALSE, clustmethod="average"){
  mno.cols = dim(datMatrix)[2]
  mno.rows = dim(datMatrix)[1]
  maxlabels= 100
  mcex=0.9
  cormatrix = NULL

  iint = as.integer(mno.cols/50) + 1
  mcex = mcex*(0.9^iint)

  if (!is.null(imagename) ){ #get file names for column and row dendrograms
       generalname    = getFileName(imagename)
       extensionname  = getFileExtension(imagename)
       myname.row = paste(generalname, "_dendrogramRow.", extensionname, sep="")
       myname.col = paste(generalname, "_dendrogramColumn.", extensionname, sep="")
  }
  
  if( (whichway=="column") | (whichway=="rowcolumn") ){
     #--------- compute correlation coefficient 
     cormatrix <- cor(datMatrix, use = "pairwise.complete.obs") 
     #--------- hierarchical clustering
     distmatrix = 1-cormatrix
     hierclu.col <- hclust(as.dist(distmatrix), method= clustmethod) 

     # we expect two or three clusters of samples from  group comparison
     if (clusterno.col>0){
       colcolors = moduleDetectByFixedClusterno(ihcluster=hierclu.col, nocluster=clusterno.col)
     }else{
       colcolors   = rep("white", mno.cols)
     }

     if (!is.null(imagename) )
         openImgDev(myname.col)
     par(mfrow=c(1,1), cex=mcex)
     if ((showlabels) && (mno.cols<= maxlabels)){
        plot(hierclu.col, xlab="",ylab="",main="",sub="")
     }else{
        plot(hierclu.col, labels=F, xlab="",ylab="",main="",sub="")
     }
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(imagename) )
          dev.off()

  }else {
     hierclu.col = NA
     colcolors   = rep("white", mno.cols)
  }

  if(whichway=="row" | (whichway=="rowcolumn") ) {
     cormatrix <- cor(t(datMatrix), use = "pairwise.complete.obs")
     distmatrix = 1-cormatrix

     hierclu.row <- hclust(as.dist(distmatrix), method= clustmethod) 

     # we expect exactly two  clusters of genes from  group comparison
     if (clusterno.row>0){
       rowcolors = moduleDetectByFixedClusterno(ihcluster=hierclu.row, nocluster=clusterno.row)
     }else{
       rowcolors   = rep("white", mno.rows)
     }
   
    if (!is.null(imagename) )
        openImgDev(myname.row)
     par(mfrow=c(1,1), cex=mcex)
     if ( (showlabels) && (mno.rows<=maxlabels)){
        plot(hierclu.row, xlab="",ylab="",main="",sub="")
     }else{
        plot(hierclu.row, labels=F, xlab="",ylab="",main="",sub="")
     }
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    if (!is.null(imagename) ) 
          dev.off()
  }else {
     hierclu.row = NA
     rowcolors   = rep("white", mno.rows)
  }

  #heatmap(as.matrix(datMatrix), Rowv=as.dendrogram(hierclu.row), Colv=as.dendrogram(hierclu.col), 
  #           RowSideColors=rowcolors, ColSideColors=colcolors,
  #           scale="row", revC=F, xlab="", col = rgcolors.func(50))

  #cat(colcolors, "\n")

  plotTwowayDendrogramsOnExprArrays(exprarray=datMatrix, hiercol=hierclu.col, hierrow=hierclu.row, 
                                              columnmodules=colcolors, rowmodules=rowcolors, imagename=imagename)

  if(!is.null(cormatrix)){
     rm(cormatrix)
     rm(distmatrix)
     collect_garbage()
  }

  list(rowcolors, colcolors)
}


plotTwowayDendrogramsOnExprArrays = function(exprarray, hiercol=NA, hierrow=NA, columnmodules=NA, rowmodules=NA, imagename=NA)
{
  max.elements = 2500
  no.cols = dim(exprarray)[2]
  no.rows = dim(exprarray)[1]

  mcex=1
  iint = as.integer(no.cols/50) + 1
  mcex = mcex*(0.9^iint)


  if (!is.na(imagename)){
     openImgDev(imagename)
     par(mfrow=c(1,1), mar=c(7, 4, 4, 2) + 0.1, cex=0.7)
  }

  orderedMatrix = as.matrix(exprarray)
  if( no.rows<=max.elements && !is.na(hierrow) ){ #normal case
     hierrow4disp  = as.dendrogram(hierrow) #have to use as.dendrogram(..)
     rowmodulesZ = rowmodules     
  }else if( no.rows>max.elements && !is.na(hierrow) ){
     #too many genes,so we order the arrays, but don't draw the dendrogram
     orderedMatrix = orderedMatrix[hierrow$order,]
     hierrow4disp  = NA
     rowmodulesZ   = rowmodules[hierrow$order]
  }else{
     hierrow4disp  = NA
     rowmodulesZ = rowmodules
  }

  if( (no.cols <= max.elements) && (!is.na(hiercol)) ){
     hiercol4disp=as.dendrogram(hiercol)
     columnmodulesZ = columnmodules
  }else if( (no.cols >max.elements) && (!is.na(hiercol)) ){
     orderedMatrix = orderedMatrix[, hiercol$order]
     hiercol4disp  = NA
     columnmodulesZ = columnmodules[hiercol$order]
  }else{
     hiercol4disp  = NA
     columnmodulesZ = columnmodules
  }

  #mycexCol = 0.2 + 1/log10(nc)

  heatmap(orderedMatrix, Rowv=hierrow4disp, Colv=hiercol4disp, 
             RowSideColors=as.character(rowmodulesZ), ColSideColors=as.character(columnmodulesZ),
             margins = c(10, 5),cexCol = mcex,
             scale="row", revC=F, xlab="", col = rgcolors.func(50))


  if (!is.na(imagename)){
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()
  }

}

##################################### EASE ANALYSIS ###############################################################
#

EASEOntologyEnrichmentAnalysis=function(genesInfor, ontologyfn, fname, maxSignifLevel=0.05,myMinPopulationHits=3, OntologyType="Ontology", background=0)
{
    #ontologyname = getOntologyNameFromPath(ontologyfn)

    cols.geneinfor <- dim(genesInfor)[2]

    modulecolor = as.character(genesInfor$module)
    ctable = table(modulecolor )
    #ctable

    #*-------------------------------------------------------------------------------------
    #* STEP 1: read in ontology information, columns as ontology category
    #*
    #*FORMAT:
    #*  System	             GeneCategory	                     PopulationHits	PopulationTotal   LocusLinkNumbers
    #* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
    #*
    #*
    #*-------------------------------------------------------------------------------------
    ontologyMatrix.all <- read.delim(ontologyfn, sep="\t", header=T)
    ontologyMatrix.all <- as.matrix(ontologyMatrix.all)

    #get unique annotated genes
    allgogenes = NULL
    for (each in as.character(ontologyMatrix.all[, 5]) ){
         #each_lowercase= tolower(each)
         each_lowercase= toupper(each)
         llids=unlist( strsplit( as.character(each_lowercase), "; ") )
         allgogenes = c(allgogenes, llids)
    }
    length(allgogenes )
    #use old trick to get unique LLIDs based on shift operations
    allgogenes.sorted = sort(allgogenes )
    allgogenes.shift  = c(-1, allgogenes.sorted[1:(length(allgogenes)-1)] )
    uniques = (allgogenes.sorted != allgogenes.shift)

    annotatedgenes = allgogenes.sorted[uniques]
    no.annotatedgenes = length(annotatedgenes)
    length(annotatedgenes)
    no.categories = dim(ontologyMatrix.all)[1]

    # population/background
    #
    if(background>0) {
       ontologyMatrix.all[,4]=as.character(background)
    }else{
       ontologyMatrix.all[,4]=as.character(no.annotatedgenes)
    }

    modulenames      = names(ctable)
    no.modules       = length(modulenames)


    #eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
    #                   "Population Hits", "Population Total", 
    eventualColnames=c("System","GeneCategory", "ModuleOverlap", "GeneSetCount", 
                       "PopulationOverlap", "Population", 
                       "pvalue_corrected", "Fisher_pvalue", "Probe Set", "Gene Symbol") 

    #process module by module
    for (eachmodule in modulenames){

       cat("Module ", eachmodule, "\n")
       modulesel    = modulecolor==eachmodule
       modulesize   = sum(modulesel)
       modgeneinfor = genesInfor[modulesel, ]
       #modgeneinfor[,1] = toupper(as.character(modgeneinfor[,1]))

       dim(modgeneinfor)

       #find module genes annotated 
       mergedmodule = merge(annotatedgenes, modgeneinfor, by.x=1, by.y=1, sort=F,all=FALSE)

       #cat(mergedmodule,"\n")

       moduletotals = dim(mergedmodule )[1]

       #write column names
       outfilename  = paste(fname, "_", OntologyType, "_",  eachmodule, ".xls",  sep='')
       #write.table(t(as.matrix(eventualColnames)), outfilename, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=F)
       
       if(background>0) { # use actual module size and given background
          a=apply(ontologyMatrix.all, 1, EASE4module, 
             modulegeneInfor=modgeneinfor, ModuleListTotals=modulesize, resultfname=outfilename, 
             maxSignifLevel=maxSignifLevel, myMinPopHits=myMinPopulationHits, no.goTerms=no.categories)

       }else{ # use actual matched genes in module and all annotated genes as background

          a=apply(ontologyMatrix.all, 1, EASE4module, 
             modulegeneInfor=modgeneinfor, ModuleListTotals=moduletotals, resultfname=outfilename, 
             maxSignifLevel=maxSignifLevel, myMinPopHits=myMinPopulationHits, no.goTerms=no.categories)
       }

       # sort the table by pvalue in ascendant order
       mymoduleGOMatrix <- read.delim(outfilename, sep="\t", header=T)
       pvalue.order     <- order(as.numeric(as.matrix(mymoduleGOMatrix$Fisher_pvalue)) )
       
       if( dim(mymoduleGOMatrix)[1]==1) {
          mytitle <- colnames(mymoduleGOMatrix)
          mymoduleGOMatrix <- as.matrix(mymoduleGOMatrix )
          mymoduleGOMatrix <- rbind(mymoduleGOMatrix)
          mymoduleGOMatrix[,4] = rep(as.character(moduletotals), dim(mymoduleGOMatrix)[1])
          fm2= cbind( rep(modulesize, dim(fm)[1]), fm)
          colnames(fm2) = c("ModuleSize", colnames(mymoduleGOMatrix))
       }else{
          mymoduleGOMatrix <- as.matrix(mymoduleGOMatrix )
          mymoduleGOMatrix[,4] = rep(as.character(moduletotals), dim(mymoduleGOMatrix)[1])
          fm = mymoduleGOMatrix[pvalue.order,]
          #print(fm)
          fm2= cbind( rep(modulesize, dim(fm)[1]), fm)
          colnames(fm2) = c("ModuleSize", colnames(fm))
       }

       write.table(fm2, outfilename, sep="\t", quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }
}

EASE4module=function(OneGOcategoryInfor, modulegeneInfor, ModuleListTotals, resultfname, 
                     maxSignifLevel=0.05, myMinPopHits=3, no.goTerms=1)
{

    #*-------------------------------------------------------------------------------------
    #*OneGOcategoryInfor FORMAT:
    #*
    #*  System	             GeneCategory	                     PopulationHits	PopulationTotal   LocusLinkNumbers
    #* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
    #*
    #*
    #*-------------------------------------------------------------------------------------
    go.fields = length(OneGOcategoryInfor)

    # ++++++++++ number of annotated genes for this category ++++++++++
    System           <- as.character(as.matrix(OneGOcategoryInfor[1]))
    GeneCategory     <- as.character(as.matrix(OneGOcategoryInfor[2]))
    PopulationHits   <- as.integer(as.matrix(OneGOcategoryInfor[3]))
    population_total <- as.integer(as.matrix(OneGOcategoryInfor[4]))

    #cat(System, GeneCategory, "\n")

    goLLIDsA = as.character(as.matrix(OneGOcategoryInfor[5]))
    #goLLIDs = tolower(goLLIDsA)
    goLLIDs  = toupper(goLLIDsA)

    ontologyGenes = unlist( strsplit(goLLIDs, "; ") )

    #case 1. blank
    if( length(ontologyGenes )<1 ){
        print(OneGOcategoryInfor)
        cat("no genes gtt splitted at all\n")
        return (0)
    }

    #*-------------------------------------------------------------------------------------
    #* STEP 2: EASE ANALYSIS of Ontology 
    #*         by merge module genes with annotated genes
    #*         
    #*        
    #* 
    #*-------------------------------------------------------------------------------------
    # Here, we avoid the merge of the enormous matrix, instead we make a matrix with geneIds and 
    # their corresponding ROW indices in the big matrix, then perform merge operation
    
    mergedAnnotGeneInfor  = merge(modulegeneInfor, ontologyGenes, by.x=1, by.y=1, sort=F,all=FALSE)
    dim(mergedAnnotGeneInfor)

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    ModuleListHits <- dim(mergedAnnotGeneInfor)[1]

    if ( (ModuleListHits ==0) | (PopulationHits<myMinPopHits) )
       return (-1)

    #eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
    #                   "Population Hits", "Population Total", 
    #                   "Pvalue_FisherExactTest","Gene Symbol","Unique ID") 

    # make a big matrix including all number required for computing pvalue
    fishervector = c(population_total, PopulationHits,
                     ModuleListTotals, ModuleListHits)

    #cat(population_total, PopulationHits, ModuleListTotals, ModuleListHits, "\n")

    mypvalue     = fisherTest(fishervector, minPopulationHits=myMinPopHits)

    # The corrected p-values represent the Bonferroni-corrected p-values 
    #    (nominal p-value multiplied by the number of GO categories searched)
    #
    mypvalue.corrected = mypvalue*no.goTerms

    if ( is.na(mypvalue) )
       return (-2)

    if(mypvalue > maxSignifLevel)
       return (-3)

    if(mypvalue.corrected>1) { mypvalue.corrected=1}

    # Concatenate selected probset names, eg, "1367659_s_at;1367659_s_at;1367659_s_at;1367659_s_at"
    # genes are row-wise, singificant function categories are column-wise
    #
    sel.genes       = rep(T, ModuleListHits)
    signif.probsets = applySelectFunc(sel.genes, as.character(mergedAnnotGeneInfor[,1]) )
    signif.symbols  = applySelectFunc(sel.genes, as.character(mergedAnnotGeneInfor[,2]) )
    signif.mypvalues= as.character(signif(mypvalue, 2))
    signif.pvalcorrect= as.character(signif(mypvalue.corrected, 2))

    sel_rslt = c(System,GeneCategory,    ModuleListHits, ModuleListTotals,
                       PopulationHits,   population_total, 
                       signif.pvalcorrect, signif.mypvalues, 
                       signif.symbols,   signif.probsets)

    write.table( t(as.matrix(sel_rslt)), resultfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=T)
 
    return (1)
}






###################################################################################################################

OntologyEnrichmentAnalysis=function(genesInfor, ontologyfn, fname, maxSignifLevel=0.05,myMinPopulationHits=3,OntologyType="Ontology")
{
    ontologyname = getOntologyNameFromPath(ontologyfn)

    cols.geneinfor <- dim(genesInfor)[2]

    color1 = genesInfor$module
    ctable = table(color1)
    #ctable

    #*-------------------------------------------------------------------------------------
    #* STEP 1: read in ontology information, columns as ontology category
    #*
    #*
    #*-------------------------------------------------------------------------------------
    ontologyMatrix.all <- read.delim(ontologyfn, sep="\t", header=T)
    ontologyGenes      <- as.character(ontologyMatrix.all[,1])
    dim(ontologyMatrix.all)

    # total number of categories
    no.categories       <- ( dim(ontologyMatrix.all)[2] - 1 )

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    # population/background
    #
    if(background>0) {
       population_total =  background
    }else{
       population_total =  dim(ontologyMatrix.all)[1]
    }

    PopulationHits   <- apply(ontologyMatrix.all[, c(2:(no.categories+1))], 2, sum ) #vector


    #*-------------------------------------------------------------------------------------
    #* STEP 2: EASE ANALYSIS of Ontology 
    #*           by merge module genes with annotated genes
    #*         
    #*        
    #* 
    #*-------------------------------------------------------------------------------------
    # Here, we avoid the merge of the enormous matrix, instead we make a matrix with geneIds and 
    # their corresponding ROW indices in the big matrix, then perform merge operation
    #
    #mergeMatrix.all = merge(genesInfor, ontologyMatrix.all, by.x=1, by.y=1, sort=F,all=FALSE)
    #dim(mergeMatrix.all)
    no.ontogenes = dim(ontologyMatrix.all)[1]
    tmpMatrix    = cbind(as.matrix(ontologyMatrix.all[,1]), c(1:no.ontogenes) )
    colnames(tmpMatrix) <- c("gene_id", "index")

    tmpmergeMatrix = merge(genesInfor, tmpMatrix, by.x=1, by.y=1, sort=F,all=FALSE)
    dim(tmpmergeMatrix)

    selectedIndices = as.integer(as.matrix(tmpmergeMatrix$index))

    mergeMatrix.all = cbind(as.matrix(tmpmergeMatrix), as.matrix(ontologyMatrix.all[selectedIndices, -c(1)]) )
    mergeMatrix.all = as.data.frame(mergeMatrix.all)
    dim(mergeMatrix.all)

    cat("  finished merge, sizes (rows, cols): ", as.character(dim(mergeMatrix.all)[1]), as.character(dim(mergeMatrix.all)[2]), "\n")

    ordered.genes      = order(as.character(as.matrix(mergeMatrix.all$module) ))
    orderedMergedMatrix <- mergeMatrix.all[ordered.genes,]

    #mergedAnnotGeneInfor <- orderedMergedMatrix[,  c(1, gene_symbolcol,cols.geneinfor, cols.geneinfor+1)]
    mergedAnnotGeneInfor <- orderedMergedMatrix[,  c(1:(cols.geneinfor+1))]
    mergedAnnotMatrix    <- orderedMergedMatrix[, -c(1:(cols.geneinfor+1))]
    dim(mergedAnnotMatrix)

    cols.mergedinfo = dim(mergedAnnotGeneInfor)[2]

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    OverallListTotal = dim(mergedAnnotMatrix)[1]
    OverallListHits <- apply(mergedAnnotMatrix, 2, forcenumericSum) #vector across all categorries

    color1.merged    = as.character(mergedAnnotGeneInfor$module)
    ModuleListTotals = table(color1.merged)
    modulenames      = names(ModuleListTotals)
    no.modules       = length(modulenames)

    # ++++++++++ 
    # we count how many genes in each module have certain function category
    ModuleListHits = apply(mergedAnnotMatrix, 2, tapplySumFunc, color1.merged)

    if (!is.matrix(ModuleListHits)){
       ModuleListHits = t(as.matrix(ModuleListHits))       
    }

    #print(modulenames)
    #print(ModuleListHits)
    #print(dim(ModuleListHits))

    rownames(ModuleListHits) <- modulenames

    eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
                       "Population Hits", "Population Total", 
                       "Pvalue_FisherExactTest","Gene Symbol","Unique ID") 

    for (imod in c(1:no.modules) ){
       
       cat(modulenames[imod],"\n")

       # make a big matrix including all number required for computing pvalue
       fisherMatrix = rbind(rep(population_total, no.categories), 
                                PopulationHits,
                                rep(ModuleListTotals[imod], no.categories),
                                ModuleListHits[imod, ])

       mypvalues     = apply(fisherMatrix, 2, fisherTest, minPopulationHits=myMinPopulationHits)
       mypvalues.sel = mypvalues <= maxSignifLevel #for column selection

       if (sum(mypvalues.sel) ==0)
           next

       # for row selection
       mymodule      = modulenames[imod]
       mod.sel       = as.character(mergedAnnotGeneInfor[,cols.geneinfor])==mymodule

       #annotated genes as a significant group for each function category
       geneSelMatrix = ( (mergedAnnotMatrix==1) * mod.sel)[, mypvalues.sel]
       geneSelMatrix = as.matrix(geneSelMatrix )

       # Concatenate selected probset names, eg, "1367659_s_at;1367659_s_at;1367659_s_at;1367659_s_at"
       # genes are row-wise, singificant function categories are column-wise
       #
       signif.probsets = apply(geneSelMatrix, 2, applySelectFunc, as.character(mergedAnnotGeneInfor[,1]))
       signif.symbols  = apply(geneSelMatrix, 2, applySelectFunc,  as.character(mergedAnnotGeneInfor[,2]))

       signif.ModuleListHits= ModuleListHits[imod, mypvalues.sel]
       signif.PopulationHits=PopulationHits[mypvalues.sel]
       signif.mypvalues     = mypvalues[mypvalues.sel]

       finaltestmatrix=NULL
       categorynames = names(signif.mypvalues)
       for (jc in c(1:length(categorynames)) ){
          sel_rslt = c(ontologyname,categorynames[jc], signif.ModuleListHits[jc], ModuleListTotals[imod],
                       signif.PopulationHits[jc], population_total, 
                       signif.mypvalues[jc], signif.symbols[jc], signif.probsets[jc])
          finaltestmatrix = rbind(finaltestmatrix, sel_rslt)
       }

       colnames(finaltestmatrix ) <- eventualColnames
       orderpvalue = order(signif.mypvalues)

       #zfinaltestmatrix = as.matrix(rbind(eventualColnames, finaltestmatrix))
       if (length(orderpvalue)>1 ){
         zfinaltestmatrix=finaltestmatrix[orderpvalue, ]
       }else{
         zfinaltestmatrix=t(finaltestmatrix[orderpvalue, ])
       }

       testoutputfname = paste(fname, "_", OntologyType, "_",  modulenames[imod], ".xls",  sep='')
       write.table(zfinaltestmatrix, testoutputfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}

##------------------------------------------------------------------------------------------------------------------
## inputfname contains: gene information, expression data, module information (last column) 
## identifier_col     : unique gene identifier column
##                      =1 for  normal case (probeset)
##                      =3 for bxh (locus link number)
##
## gene_symbolcol     : we assume the first col as probset and this column gives the real gene names
## ontologyfnlist     : lists of GeneFisher ontology files 
## maxSignifLevel     : report only the categories with FisherExactTest Pvalue < maxSignifLevel
## useEASEannotation  : use Affymetrix or EASE annotation file
## outdir             : =="", put the results under the same directory as the input, otherwise under the new directory
##
##  OntologyType="Ontology"/"TF"/"QTL"
##
OntologyAnalysisDull=function(inputfname, identifier_col=1, gene_symbolcol=1, ontologyfnlist, signifLevel=0.05, minPopHits=3, 
                              useEASEannotation=T, useAllModules=F,OntologyType="Ontology", outdir="", ctopNrows=2, background=0)
{

    #* STEP 0: read in gene information, expression data, module information

    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    #attach(allMatrix)
    dim(allMatrix)

    if(outdir=="") {
       mfname        = getFileName(inputfname)
    } else{
       mfname        = getFileNameNopath(inputfname)
       mfname        = paste(outdir, "/", mfname, sep="")
    }

    no.cols  = dim(allMatrix)[2]

    if (!useAllModules){ # consider modules
       inforCols = c(identifier_col, gene_symbolcol, no.cols)
       mgenesInfor <- allMatrix[,inforCols]
    }else{ # consider all gene in the list
       inforCols = c(identifier_col, gene_symbolcol)
       mgenesInforA <- allMatrix[,inforCols]
       allinfo = rep("all", dim(allMatrix)[1]) # make a false module with all genes
       allinfo=factor(allinfo,levels=allinfo)
       mgenesInfor <- cbind(mgenesInforA, allinfo)
       colnames(mgenesInfor) <- c( colnames(mgenesInforA), "module")
    }

    #print(mgenesInfor[1:5,])

    #unify upper/lower cases
    mgenesInfor[, 1] = toupper( as.character(mgenesInfor[, 1]))

    rowTitles=names(allMatrix)
    dim(mgenesInfor)

    color1 = mgenesInfor$module
    ctable = table(color1)
    ctable
    
    if(OntologyType=="QTL"){
        top2colnames=c("Chromosome","LocusBin")
    }else if (OntologyType=="TF"){
        top2colnames=c("TranscriptionFactor","TF")
    } else if (OntologyType=="knockout"){
        top2colnames=c("knockout","knockout")
    } else{
        top2colnames=c("System","Gene Category")
    }
    
    #eventualColnames=c(top2colnames, "List Hits", "List Total", 
    #                   "Population Hits", "Population Total",
    eventualColnames=c(top2colnames, "ModuleOverlap", "GeneSetCount", 
                       "PopulationOverlap", "Population",
    #                   "Pvalue_FisherExactTest","Gene Symbol","Unique ID")
                       "pvalue_corrected", "Fisher_pvalue", "Probe Set", "Gene Symbol") 


    #delete exiting reuslts
    firstfile = NULL
    outfiles  =NULL
    for (each in names(ctable) ){
       testoutputfname = paste(mfname, "_",OntologyType, "_",  each, ".xls",  sep='')
       write.table(t(as.matrix(eventualColnames)), testoutputfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
       if ( is.null(firstfile) ){
           firstfile = testoutputfname 
       }
       outfiles = c(outfiles, testoutputfname)
    }

    for (each in ontologyfnlist){
       cat(each, "\n")
       if (useEASEannotation){
           EASEOntologyEnrichmentAnalysis(genesInfor=mgenesInfor, ontologyfn=each, 
                            fname=mfname, maxSignifLevel=signifLevel, myMinPopulationHits=minPopHits,
                            OntologyType=OntologyType, background=background)
       }else{
           OntologyEnrichmentAnalysis(genesInfor=mgenesInfor, ontologyfn=each, 
                                      fname=mfname, maxSignifLevel=signifLevel, myMinPopulationHits=minPopHits,
                                      OntologyType=OntologyType, background=background)
       }
    }

    # combine all
    topMatrix=NULL
    coutputfname = paste(mfname, "_",OntologyType, ".xls",  sep='')
    for(j in c(1:length(outfiles) ) ) {
       efo = outfiles[j]
       print (efo)
       ematrix= retrieveTopNrows(efo, topNrows=ctopNrows,  OntologyType=OntologyType)
       if (!is.null(ematrix) ){
          topMatrix=rbind(topMatrix, ematrix)
       }
    }
    #write.table(topMatrix[,-c(2)], coutputfname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    # use shortnames
    systems= as.character(topMatrix[,3])
    topMatrix[,3]=ifelse(systems=="GO Biological Process", "BP",
                     ifelse(systems=="GO Cellular Component","CC", 
                     ifelse(systems=="GO Molecular Function","MF", systems) ))

    msizes = as.integer(topMatrix[,2])
    szorder= order(-msizes)

    write.table(topMatrix[szorder,], coutputfname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    # return the first Ontology output file
    return (firstfile )
}

# get LLIds of the last column
#* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
unlistLLI = function(easeline){
   llids=unlist( strsplit( as.character(easeline[5]), "; ") )
   llids
}


###########################################################################################################
#
#  for small gene set: in output file, genes as column and markers as row
#
###########################################################################################################

computeQTLmatrix_4smallgeneset = function(inputfname, headCol=8, genoCSVfname, outdir="")
{
    #----------------------------- 1. READ expression -------------------------------------
    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    attach(allMatrix)
    dim(allMatrix)

    qtlFileDir ="./"

    no.cols <- dim(allMatrix)[2]
    rowTitles <- rownames(allMatrix)
    colTitles <- colnames(allMatrix)

    fname      = getFileName(inputfname)
    eqtlFname  = paste(outdir,fname, "_eQTLmatrix.xls", sep='')
    countFname  = paste(outdir,fname, "_eQTLcounter.xls", sep='')

    #These are the expression values, the last three cols are module infomation
    datExpr <- t(allMatrix[, -c(1:headCol)])
    dim(datExpr)

    genesInfor  <- allMatrix[,c(1:headCol)]
    genenames   <- allMatrix[,1]

    #----------------------------- 2. READ geno markers -------------------------------------
    #datGenoTable  <- read.cross("csv", dir=qtlFileDir, file=genoCSVfname)
    #datGeno        = calc.genoprob(datGenoTable, step=0, error=0.01)

    genoMatrix <- read.delim(genoCSVfname,sep=",", header=F)
    dim(genoMatrix)
    markerInfor <- t(genoMatrix[c(1:3),-c(1,2)])
    dim(markerInfor)

    colnames(markerInfor) <- c("SNP","chr","pos")

    #----------------------------- 3. Compute eQTL matrix -----------------------------------
    #
    # eQTL for all genes, rows as SNPs and column as genes
    eQTLmatrix = apply(datExpr, 2, qtlOnGenome) 

    colnames( eQTLmatrix ) <- as.character(genesInfor[,1])
    eQTLmatrix.4save = cbind(markerInfor, eQTLmatrix)

    #here we threshold eQTL matrix to get the total number of genes with eQTL at each marker
    thresholds = seq(2, 6, 0.5)
    countermatrix = NULL
    for (ithresh in thresholds){
	bool.eQTLmatrix = eQTLmatrix >= ithresh
        icount=apply(bool.eQTLmatrix, 1, sum, na.rm =T)
        countermatrix = cbind(countermatrix, icount)
    }
    colnames(countermatrix) <- as.character(thresholds)
    countmatrix.4save = cbind(markerInfor, countermatrix)
    write.table(countmatrix.4save,countFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 

    write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 
    dim(eQTLmatrix)

    rm(eQTLmatrix)
    rm(allMatrix)
}









































###########################################################################################################
#
#  for huge gene set: markers as column and genes as row
#
###########################################################################################################

computeQTLmatrix_4hugegeneset = function(inputfname, headCol=8, genoCSVfname, outdir="")
{
    #----------------------------- 1. READ expression -------------------------------------
    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    attach(allMatrix)
    dim(allMatrix)

    qtlFileDir ="./"

    no.cols <- dim(allMatrix)[2]
    rowTitles <- rownames(allMatrix)
    colTitles <- colnames(allMatrix)

    fname      = getFileName(inputfname)
    eqtlFname  = paste(outdir,fname, "_eQTLmatrix.xls", sep='')
    countFname  = paste(outdir,fname, "_eQTLcounter.xls", sep='')

    #These are the expression values, the last three cols are module infomation
    datExpr <- t(allMatrix[, -c(1:headCol)])
    dim(datExpr)

    no.genes = dim(datExpr)[2]

    genesInfor  <- allMatrix[,c(1:headCol)]
    genenames   <- allMatrix[,1]

    #----------------------------- 2. READ geno markers -------------------------------------
    #datGenoTable  <- read.cross("csv", dir=qtlFileDir, file=genoCSVfname)
    #datGeno        = calc.genoprob(datGenoTable, step=0, error=0.01)

    genoMatrix <- read.delim(genoCSVfname,sep=",", header=F)
    dim(genoMatrix)
    markerInfor <- t(genoMatrix[c(1:3),-c(1,2)])
    dim(markerInfor)

    no.markers = dim(markerInfor)[1]
    colnames(markerInfor) <- c("SNP","chr","pos")

    #----------------------------- 3. Compute eQTL matrix -----------------------------------
    #
    # eQTL for all genes, rows as SNPs and column as genes

    #here we threshold eQTL matrix to get the total number of genes with eQTL at each marker
    thresholds    = seq(2, 6, 0.5)
    no.thresh     = length(thresholds)
    countermatrix = matrix(0, no.markers, no.thresh)
    colnames(countermatrix) <- as.character(thresholds)
    
    # decompose the genes into several small datasets
    intvls    = seq(1, no.genes, 50)
    intvls    = c(intvls, no.genes)
    no.intvls = length(intvls)
    for (i in c(1:(no.intvls-1)) ){
       if (i < no.intvls-1){
          indices        = c(intvls[i]:(intvls[i+1]-1) )
       }else{
          indices        = c(intvls[i]:(intvls[i+1]) )
       }

       # returned matrix with genes as cols and markers as rows
       eQTLmatrix     = apply(datExpr[, indices], 2, qtlOnGenome)
       
       if (i==1){
          #row names for the first subset
          combeQTLmatrix = t(cbind(markerInfor, eQTLmatrix) ) # marker as the top rows
          
          #row names for the first subset
          irownames = c(colnames(markerInfor), as.character(genesInfor[indices,1]) )
          eQTLmatrix.4save = cbind(irownames, combeQTLmatrix)
          write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE) 

       }else {
          irownames = c(as.character(genesInfor[indices,1]) )
          eQTLmatrix.4save = cbind(irownames, t(eQTLmatrix))
          write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T) 
       }


       for ( k in c(1:no.thresh) ){
          kthresh = thresholds[k]
    	  bool.eQTLmatrix = eQTLmatrix >= kthresh
          icount=apply(bool.eQTLmatrix, 1, sum, na.rm =T)
          countermatrix[,k] = countermatrix[,k] + icount
       }

    }

    countmatrix.4save = cbind(markerInfor, countermatrix)
    write.table(countmatrix.4save,countFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
}

############################### NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN #######################################

# input: "literature-network_unique_tommodules_Ontology_blue.xls"
# output: "blue"
getModuleNameFromEASEresfile=function(easeRESfile, OntologyPatt="_Ontology_"){
  a=unlist( strsplit(easeRESfile, OntologyPatt) )
  if (length(a) <=1)
       return ("")
  b=unlist( strsplit(a[2], ".xls") )
    if (length(b) <=0)
       return ("")
  return (b)
}

retrieveTopNrows=function(inputfname, topNrows=3, OntologyType="Ontology"){

  ontopatt = paste("_", OntologyType, "_",sep="")
  
  modulename = getModuleNameFromEASEresfile(inputfname, OntologyPatt=ontopatt)
  if (modulename=="")
     return (NULL)

  allMatrix <- read.delim(inputfname,sep="\t", header=T)
  dim(allMatrix)
  nrows = dim(allMatrix)[1]
  if (nrows<=0)
     return (NULL)
  actualrows = min(topNrows, nrows)

  extracolumn = rep(modulename, actualrows)
  
  outmatrix   = cbind(extracolumn, allMatrix[c(1:actualrows),])
  newcolnames = c("module", colnames(allMatrix) )

  colnames(outmatrix)  <- newcolnames

  return (outmatrix)
}


writeHugeTable = function(mymatrix, outfname, colnames=F, myappend=F)
{

  if (colnames){
     write.table(t(as.matrix( colnames(mymatrix) )), outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=myappend)
  }

  # write into file in many times for saving memory
  norows=dim(mymatrix)[1]
  step=3000
  intervals = seq(from=1, to=norows, by=step )

  # test if the last one equal the no.rows, otherwise, append no.rows to the end
  no.intrervals= length(intervals)
  if (intervals[no.intrervals] != norows){
    intervals = c(intervals, norows)
  }
  no.intrervals= length(intervals)
  intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

  total_matched=0
  for (i in c(2:no.intrervals) ){
   sidx= intervals[i-1]
   eidx= intervals[i]-1
   #cat(sidx,"\n")
   #cat(eidx,"\n")
   irows= c(sidx:eidx)
   imatrix=mymatrix[irows, ]

   write.table(imatrix, outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
  }

}

#################### Histogram -based Functions ###########################################


# Here we compute the histogram of a list of correlation coefficients [0,1]
# but we customize the output as the frequency in fixed bins specified by
#  the input parameter "nobreaks", eg, 20 intervals in the range of [0,1]
#
histogram_customize = function(corvect,nobreaks=20,xmin=NULL,xmax=NULL,vectorname="", xlabel="", imgfilename=NULL)
{
  # draw histogram, compute frequency of no of links
  #numvect = as.numeric(corvect)
  maxsize = 1000000

  nfold = nobreaks/10
  if(length(corvect) <= maxsize ){
      numvect = as.numeric(corvect)+0.05/nfold
  }else{
      tmpvect= sample(corvect, maxsize, replace=F)
      numvect = as.numeric(tmpvect)+0.05/nfold
  }

  ntab = table(round(numvect*nfold, 1) )
  linkfreq = as.numeric(ntab)
  linkidx  = as.numeric( names(ntab) )/nfold
  linkidx

  no.items = length(numvect)

  histMatrix = cbind(as.character(linkidx), as.character(linkfreq))

  if(is.null(xmin)){
     mxmin=min(numvect,na.rm=T)
  }else{
     mxmin=xmin
  }
  if( is.null(xmax) ){
     mxmax=max(numvect,na.rm=T)
  }else{
     mxmax=xmax
  }

  completeIntervals  = seq(from=mxmin, to=mxmax, by=1/nobreaks)


  completeHistMatrix = merge(as.matrix(as.character(completeIntervals)), as.matrix(histMatrix), by.x=1, by.y=1, all.x=T)
  completeHistMatrix <- as.matrix(completeHistMatrix)
  #completeHistMatrix
  completeHistMatrix = ifelse(is.na(completeHistMatrix),0, completeHistMatrix)
  #completeHistMatrix
  completeHistMatrix2 <- matrix(as.numeric(completeHistMatrix), ncol=2)
  #completeHistMatrix2

  completeHistMatrix2[,2] = completeHistMatrix2[,2] /no.items

  #mylabel = vectorname
  #outimg = paste(filename, "-HIST", vectorname, ".png", sep="")

     index=completeHistMatrix2[,1]
     freq =completeHistMatrix2[,2]

     maxval = round(max(corvect,na.rm=T),2)
     meanval= round(mean(corvect,na.rm=T), 2)

     mylabel = paste(vectorname, ": max=", as.character(maxval), 
                           ", mean=", as.character(meanval), sep="" )

  if (!is.null(imgfilename) ){
     openImgDev(imgfilename, iwidth = 600, iheight = 600)
  }

     barplot(freq, names.arg= as.character(index),
       xlab=xlabel, ylab="frequency", main= mylabel, col="blue")

    #plot(c(0,1), c(0,1), xlab=xlabel,ylab="frequency",col='white' )
    #lines(binned.x, linkfreq, col="black");

    #plot(linkidx, linkfreq,xlab="number of links of a node", ylab="frequency", main= keyword, type="h")
    #histogram(as.numeric(corvect), br=seq(from=0,to=1, by=0.05), col="blue",
    #         xlab=xlabel, ylab="frequency", main= mylabel)

  if (!is.null(imgfilename) ){
    dev.off()
  }

  return (completeHistMatrix2)
}

# bar plot of count and percent of integer vector
#
# if the no of unique values >100, the program automatically compute the binned distribution
#  otherwise, no bins are used for distribution computation
#
histogram_4integers =function(ivector,fkeyname=NULL,keyword="",hxlabel="k",no.cuts= 15){

    maxsize = max(ivector)
    minsize = min(ivector)
    if (length(ivector) >1) {
       szrange = c(minsize : maxsize)
    }else{
       szrange = c((ivector-1):(ivector+1) )
    }

    uniques = names(table(ivector))
    if (length(uniques) < 100) {
        no.cuts = length(szrange)
        binned.k= szrange
    } else{
        intervalSize = (maxsize -minsize)/no.cuts
        binned.k     = rep(minsize-intervalSize/2, no.cuts)
        for (kk in c(1:no.cuts) ) {
           binned.k[kk] = as.integer( (binned.k[kk] + kk*intervalSize) )
        }
    }
    
    cut1 = cut(ivector, no.cuts )
    frequence= tapply(ivector,cut1,length)
    frequence= as.numeric(frequence)
    frequence= ifelse( is.na(frequence),0, frequence )
    percent  = frequence/length(ivector)

    imgCount  =paste( fkeyname, "-Count.png",sep='')  #png only
    imgPercent=paste( fkeyname, "-Percent.png",sep='')  #png only

    # save statistics into log file
    #
    counts = c("count",   length(ivector), frequence)
    freqs  = c("percent(%)", 100,         signif(percent,3)*100)
    title  = c(keyword, "total", as.character(binned.k) )
    logMatrix = rbind(title, counts, freqs)

    if(is.null(fkeyname)){
        return (logMatrix)
    }

    maintitle =""
    if (keyword !=""){
        #maintitle = paste("distribution of number of ", keyword, sep="") 
        maintitle = paste("distribution of ", keyword, sep="") 
    }

    openImgDev(imgCount,  iwidth = 600, iheight = 600)
    barplot(frequence, names.arg= as.character(binned.k),  col="green",
     ylab="count", xlab="k", main= maintitle)
    dev.off()

    openImgDev(imgPercent,  iwidth = 600, iheight = 600)
    #histogram(kcliquesSize, type="percent",xlab=hxlabel)
    barplot(percent*100, names.arg= as.character(binned.k), col="green",
     ylab="percent", xlab=hxlabel, main= maintitle )
    dev.off()

    return(logMatrix)
}



#
# wil.cox test is very slow for huge dataset, so here we use sampling method to solve this dilemma
# if a vector is too big (>maxsize), then we need sample the dataset into a set with length
#   < maxsize
#
wilcoxtest4Hugedata = function(numericA, numericB, maxsize = 50000)
{
   #maxsize = 50000
   if ( length(numericA)>maxsize ){
       snumericA = sample(numericA, maxsize, replace=F)
   }else{
       snumericA = numericA
   }
   if ( length(numericB)>maxsize ){
        snumericB = sample(numericB, maxsize, replace=F)
   }else{
        snumericB = numericB
   }

   wtest = wilcox.test(snumericA, snumericB)
   wtest$p.value
}


# perform wil.cox ranking test on the column based vectors
# return a table of pvalue table
wilcoxtest4sets = function(datalist, usewilcox=T)
{
  no.vects = length(datalist)
  datout    = matrix('-',nrow=no.vects,ncol=no.vects)
  sets      = names(datalist)
  for (i in c(1:(no.vects-1) ) ){
   for (j in c((i+1):(no.vects) ) ){
       if(usewilcox){
          ijpval      = wilcoxtest4Hugedata(datalist[[i]], datalist[[j]])
       }else{
          tt    = t.test(datalist[[i]], datalist[[j]])
          ijpval=tt$p.value
       }
       datout[i,j] = signif(ijpval, 3)
   }
  }

  # add a column for the set names
  fdatoutA = cbind(sets, datout)
  #print(fdatout)
  #print(sets)
  colnames(fdatoutA) <- c("Wilcoxon test", sets)

  fdatout = rbind(c("Wilcoxon.pvalue", sets), fdatoutA)

  fdatout
}


# transform single point based histograms into step-based representation
#  
transformHisto2Steps=function(mindex, mhist){
  npoints = length(mindex)
  mstep   = mindex[2] - mindex[1]
  hstep   = mstep/2
  newindex=NULL
  newhist =NULL
  for (i in c(1:npoints) ){
       idx = mindex[i] # current value
       iy  = mhist[i]
       newindex = c(newindex, idx-hstep)     
       newindex = c(newindex, idx+hstep)
       newhist  = c(newhist, iy)
       newhist  = c(newhist, iy)
  }
  newpoints = length(newindex)

  #process start and end which are now out of the original limit
  # because of operations + and - hstep
  newindex[1]=mindex[1]
  newindex[newpoints]=mindex[npoints]

  retMatrix = cbind(newindex, newhist)  

  return (retMatrix )
}

# plot multi-histograms in one figure based on the step-style plot
# mindex:     the mean values of intervals
# histograms: the histograms, all were computed in the same intervals
#
plotMultiHistsograms= function(histograms,mxindex,xlabel="",ylabel="",maintitle="",filename, stepshape=T)
{

 if(stepshape){
  stepHists = NULL
  no.hists = dim(histograms)[2]
  for (i in c(1:no.hists) ){
    istephist=transformHisto2Steps(mxindex, histograms[,i])
    stepHists = cbind(stepHists, istephist[,2])
  }
  stepindex = istephist[,1]
  stepHists <- as.matrix(stepHists)
  colnames(stepHists) <- colnames(histograms)

  plotMultiProfiles(eQTL_profiles=stepHists, xindex=stepindex, 
            xlabel=xlabel,ylabel=ylabel,mainlabel=maintitle, 
            plottype="l", ltype=rep(1,no.hists), filename=filename)
 }else{
  plotMultiProfiles(eQTL_profiles=histograms, xindex=mindex,
            xlabel=xlabel,ylabel=ylabel, mainlabel=maintitle, 
            plottype="l", ltype=rep(1,no.hists), filename=filename)
 }
 
}

# inputs:
# networkfname ~ network link file where each entry is a pair of genes, separated a separator
# idxMatrix    ~ a matrix where the 1st column includes the gene names and the 2nd column are the 
#                coresponding index in the correlation matrix
#
# The output is a 1x2 matrix, where each row, containing two integer numbers, 
#         corersponds to a link in the network file
#
matchNetworkLinks2Corrmatrixindex = function(networkfname, idxMatrix)
{
  bbMatrixAll <- read.delim(networkfname, sep="\t", header=F)
  bbMatrixAll <- as.matrix(bbMatrixAll)
  #dim(bbMatrixAll)

  bbMatrix <- bbMatrixAll[,c(1,2)]
  indexMatrixInt = matchNetworkLinks2Corrmatrixindex_Subfunc(netlinks=bbMatrix, 
                 idxMatrix=idxMatrix)

  return (indexMatrixInt)
}

matchNetworkLinks2Corrmatrixindex_Subfunc = function(netlinks, idxMatrix) {
  leftMatched = merge(idxMatrix, netlinks, by.x=1, by.y=1,all=F)
  #dim(leftMatched)

  bothMatched = merge(idxMatrix, leftMatched, by.x=1, by.y=3,all=F)
  #dim(bothMatched)

  indexMatrix  = bothMatched[,c(2,4)]
  indexMatrix  <- as.matrix(indexMatrix)

  indexMatrixS = cbind(as.matrix(bothMatched[,4]), as.matrix(bothMatched[,2]) )
  indexMatrixT = data.frame(indexMatrixS)

  indexMatrixInt = cbind(as.integer(indexMatrix[,1]), as.integer(indexMatrix[,2]))

  return (indexMatrixInt)
}


# turn one row of an adj matrix into pairwise representation
#
adj2pair = function(adjrow, ithresh=0.8, genenames, outfile)
{  
   rlen = length(adjrow)
   sidx = adjrow[rlen] #source node
   if(sidx >=(rlen-1) ){
       return
   }   

   #consider only the upper part of the matrix above diagonal
   mask = c(rep(F, sidx), rep(T, rlen-sidx))

   sel = adjrow >= ithresh # destination nodes
   sel = (sel & mask)

   sel[rlen] = F # index should be reset   

   if(sum(sel)>0){
     selgenenames <- genenames[sel]
     pairnet=paste(genenames[sidx], selgenenames, sep="\t")
   
     write.table(as.matrix(pairnet), outfile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

     return ( cbind(genenames[sidx], selgenenames) )
   }
   return (NULL) 
}



adj2pairSimple = function(adjrow, sidx, genenames, outfile, directed =F)
{  
   sel  = adjrow ==1 # destination nodes

   # for undirected network, we consider only the upper diaganal part
   #  the trick is to mask out the lower part of the adjrow
   #
   if( !directed) {
     maskoutIdx = c(1:sidx)
     sel[maskoutIdx ] = F
   }

   if(sum(sel)>0){
     selgenenames <- genenames[sel]
     pairnet=paste(genenames[sidx], selgenenames, sep="\t")
   
     write.table(as.matrix(pairnet), outfile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

     return ( pairnet )
   }
   return (NULL) 
}


adjmatrix2linkpairs = function(adjmatrix, genenames, outfile, directed = F)
{
   #clean up
   yfile <- file(outfile, "w+")
   close(yfile)
  
   nrows = dim(adjmatrix)[1]
   for(i in c(1:nrows) ){
       adj2pairSimple(adjrow=adjmatrix[i,], sidx=i, genenames=genenames, outfile=outfile, directed=directed)
   }
}




# here the saving process is decomposed into many small steps, each of which saves only
# a small matrix into the file so that it doesn't need a lot of memory for the whole operation
# but with a little bit more time
#
saveHugeMatrix =function(hugematrix, titlerow=NULL, outfilename,  use1minus=F){

    #title row
    if (!is.null(titlerow)){
       write.table(t(as.matrix(titlerow)), outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    }

    # write into file in many times for saving memory
    no.rows=dim(hugematrix)[1]
    step   =100
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=hugematrix[irows, ]

       if (use1minus){
         write.table(imatrix, outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
       }else{
         write.table(1-imatrix, outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
       }
    }

}

#******************************************************************************************
#********************* FOR converting pair link list into adjacency matrix ****************
#******************************************************************************************
#

############################## links pairs to adj maTRIX #########################
# merge to get indices of two genes forming a link
# netmatrix, mappingMatrix are bot Nx2 matrix
# netmatrix, each row contains source and dest nodes
# mappingMatrix, each rwo represents the mapping between node name and index
#
getNetworknodeIdxByNames = function(netmatrix, mappingMatrix)
{
  # node on the left side
  mergedleft = merge(netmatrix, mappingMatrix, by.x=1, by.y=1, all=F)

  # node on the right side
  mergedright = merge(mergedleft, mappingMatrix, by.x=2, by.y=1, all=F)
  mergedright = as.matrix(mergedright)
  indices2D = cbind(as.integer(mergedright[,3]), as.integer(mergedright[,4]) )
  indices2D
}

# note that name2idxMatrix is a global variable
# coding == NULL: use the third column in inetmatrix
#
makeAjacencyMatrix = function(inetmatrix, coding, matrixsize, directed=F, myname2idxMatrix)
{
  # initialization, no links between any two nodes
  adjMatrix = matrix(0, ncol=matrixsize, nrow=matrixsize)

  # find nodes of each edge, and replace the corresponding element with the edge code
  edgeIndices=getNetworknodeIdxByNames(netmatrix=inetmatrix[,c(1,2)], mappingMatrix=myname2idxMatrix)
  if( !is.null(coding) ) {
     adjMatrix[edgeIndices] = coding[2]
  }else{
     adjMatrix[edgeIndices] = inetmatrix[,3]
  }

  if(!directed){ # for undirected graph, the ajacency matrix is symmetric
    transposedIndices = edgeIndices[,c(2,1)]

    if( !is.null(coding) ) {
       adjMatrix[transposedIndices] = coding[2]
    }else{
       adjMatrix[transposedIndices] = inetmatrix[,3]
    }
  }

  adjMatrix
}


# return adj lists for all nodes
# if returnUniqueNames==T, return only the node names
#
# notice that a node without a link is marked as -1
#
makeAjacencyLists=function(inetmatrix, returnUniqueNames=F, directed=F)
{
    edgesInNet = dim(inetmatrix)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(inetmatrix[,1]) )
    allnodenames = c(allnodenames, as.character(inetmatrix[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    if (returnUniqueNames){
       return ( list(uniquenames, no.uniquenames) )
    }

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, no.uniquenames) )

    # find nodes of each edge, and replace the corresponding element with the edge code
    edgeIndices = getNetworknodeIdxByNames(netmatrix=inetmatrix[,c(1,2)], mappingMatrix=name2idxMatrix)
    no.edges    = dim(edgeIndices)[1]  

    for (k in c(1:no.edges) ){
       adjlists[[ edgeIndices[k,1] ]] = c(adjlists[[ edgeIndices[k,1] ]], edgeIndices[k,2])
       if (!directed){
          adjlists[[ edgeIndices[k,2] ]] = c(adjlists[[ edgeIndices[k,2] ]], edgeIndices[k,1])
       }
    }

    # remove redundant elements & -1
    #
    for (i in c(1:no.uniquenames)){
       if ( length(adjlists[[i]]) >1){
          adjlists[[i]]=setdiff(adjlists[[i]], -1)
       }
    }

    return (adjlists)
}


 
# input link pairs are node indices, so we construct 
# adjlist for the nodes from 1 to the maximum index 
#
makeAjacencyListsFromLinksIndex = function(linksindex, excludedNodes=NULL)
{
    no.edges = dim(linksindex)[1]
    maxIndex = max( max(linksindex[,1]), max(linksindex[,2]) )

    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, maxIndex) )

    for (k in c(1:no.edges) ){
       adjlists[[ linksindex[k,1] ]] = c(adjlists[[ linksindex[k,1] ]], linksindex[k,2])
       adjlists[[ linksindex[k,2] ]] = c(adjlists[[ linksindex[k,2] ]], linksindex[k,1])
    }

    # remove redundant -1
    #
    for (i in c(1:maxIndex )){
       if ( length(adjlists[[i]]) >1){
          adjlists[[i]]=setdiff(adjlists[[i]], -1)
       }
    }

    return (adjlists)
}

# make adj lists for a subset of nodes
#
makeAjacencyListsFromSubsetnodesindex = function(orgAdjLists, subsetNodes=NULL)
{
    if(is.null(subsetNodes)) {
       return(orgAdjLists)
    }

    no.orgnodes = length(orgAdjLists)
    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, no.orgnodes) )
    selNodes = rep(F, no.orgnodes)
    selNodes = T

    for (i in  subsetNodes){
       ioverlap = intersect(orgAdjLists[[i]],subsetNodes)
       if( length(ioverlap)>0) {
          adjlists[[ i ]] = ioverlap
       }
    }

    return (adjlists)
}

# convert dichotmize matrix into link pairs
#

dichoMatrix_to_linkpairs=function(dichotCor, genenames, pairExprnet)
{
  no.genes = dim(dichotCor)[1]
  
  for (i in 1:(no.genes) ){
     noleft     = no.genes - i
     iselpool   = c(rep(F, i), rep(T, noleft) )
     iconnected = (dichotCor[i,] & iselpool)

     if(sum(iconnected )==0){
        next
     }

     ipairs     = paste(genenames[i], "\t", genenames[iconnected],sep="") 
  
     if(i==1){
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
     }else{
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
     }
  }
}




#********************************** INPUT **********************************
# each gene pair in "netpairMatrix" represents a link with two nodes
# for directed networks, the first node (element) is the source and 
#    the second is the destination
#
# codepair=c(0,1):  #[1] for no connection, [2] for connection; ==NULL, use confidence value from the third column of netpairMatrix
# directed= T: #TRUE    #FALSE # directed network or not

convert_links2adjmatrix= function(netpairMatrix,directed= T, codepair=c(0,1), keywords="", fadjname=NULL) {

    #keywords = getFileName(fname)
    #fadjname = paste(outputDir, keywords, ".adj", sep="")

    edgesInNet = dim(netpairMatrix)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(netpairMatrix[,1]) )
    allnodenames = c(allnodenames, as.character(netpairMatrix[,2]) )
     
    #nametable = table(allnodenames)
    #length(nametable)
    #uniquenames     = names(nametable)
    
    uniquenames     = union(allnodenames , NULL)
    uniquenames     = sort(uniquenames)

    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    adjmatrix = makeAjacencyMatrix(inetmatrix=netpairMatrix, coding=codepair,
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix = name2idxMatrix)

    #edgeIndices=getNetworknodeIdxByNames(netmatrix=netpairMatrix[,c(1,2)], mappingMatrix=name2idxMatrix)

    if(is.null(fadjname)){
      return (adjmatrix)
    }


    #---------------------------- output adj matrix ------------------------------------------
    #title row
    write.table(t(as.matrix(uniquenames)), fadjname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, fadjname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    #---------------------------- output statistics --------------------------------------------
    # matched percentage
    #
    linkspernode = as.integer(100*edgesInNet/no.uniquenames)/100
    matchMatrix  = cbind(keywords, edgesInNet, no.uniquenames, linkspernode)
    fmatrix      = rbind( c("network", "# of edges", "# of nodes", "edges per node"), matchMatrix)

    #write.table(fmatrix, flogname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    fmatrix
}

# sna format: 
#*vertices 760		
#nodename	color	shape
#YAL005C	white	50
#YAL034C	white	50
#......
#YAL005C	YAL034C	YJL141C	YAL054C
#0	0	0	0
#0	0	1	0


# normColor/normShape for default nodes
# highColor/highShape for highlighted nodes (signature)
# koColor/koShape for knockout
#
#
makeSNA_from_netpairs=function(netpairs, highlightNodes=NULL, 
           knockouts=NULL,
           koColor  ="green",koShape  ="4",
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, snafile="tmp.sna")
{
    edgesInNet = dim(netpairs)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(netpairs[,1]) )
    allnodenames = c(allnodenames, as.character(netpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = sort(names(nametable))
    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    
    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix
    verticesMatrix = NULL
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      #shape
      if (length(highShape)==1){#use the same shape for all nodes
          nodesHighShape = cbind(highlightNodes,
                                 rep(highShape,length(highlightNodes)) )
      }else{
          nodesHighShape = cbind(highlightNodes,highShape)
      }

      # color
      if (length(highColor)==1){#use the same shape for all nodes
          nodesHighColor = cbind(highlightNodes,
                                 rep(highColor,length(highlightNodes)) )
      }else{
          nodesHighColor = cbind(highlightNodes,highColor)
      }

      #align highlighted nodes to the network index
      #
      colored = mergeTwoMatricesByKeepAllPrimary(cbind(uniquenames),nodesHighColor, 
                                            missinglabel=normColor,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
      verticesMatrix = mergeTwoMatricesByKeepAllPrimary(colored,nodesHighShape, 
                                            missinglabel=normShape,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
    }
    verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")

    # change the settings of knockout
    #
    if( !is.null(knockouts) ){
       for(eko in knockouts) {
          if(!is.element(eko, uniquenames)){
             next
          }
          koIdx = c(1:no.uniquenames)[eko == uniquenames]
          verticesMatrix[koIdx,2]= koColor
          verticesMatrix[koIdx,3]= koShape
       }
    }

    # 3. adjacency matrix
    adjmatrix = makeAjacencyMatrix(inetmatrix=netpairs, coding=c(0,1),
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix = name2idxMatrix)
    
    #edgeIndices=getNetworknodeIdxByNames(netmatrix=netpairMatrix[,c(1,2)], mappingMatrix=name2idxMatrix)

   
    #---------------------------- output adj matrix ------------------------------------------
    #
    write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    # title row
    write.table(t(as.matrix(uniquenames)), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNA from matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# nodenames are sort and in the same order as in the col/row names
# adjmatrix define the color of edges or simply connected or not
#
makeSNA_from_Adjmatrix=function(adjmatrix, nodenames, highlightNodes=NULL, 
           edgecolorlevels,
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, snafile="tmp.sna")
{
    uniquenames     = nodenames
    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    
    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix
    verticesMatrix = NULL
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      #shape
      if (length(highShape)==1){#use the same shape for all nodes
          nodesHighShape = cbind(highlightNodes,
                                 rep(highShape,length(highlightNodes)) )
      }else{
          nodesHighShape = cbind(highlightNodes,highShape)
      }

      # color
      if (length(highColor)==1){#use the same shape for all nodes
          nodesHighColor = cbind(highlightNodes,
                                 rep(highColor,length(highlightNodes)) )
      }else{
          nodesHighColor = cbind(highlightNodes,highColor)
      }

      #align highlighted nodes to the network index
      #
      colored = mergeTwoMatricesByKeepAllPrimary(cbind(uniquenames),nodesHighColor, 
                                            missinglabel=normColor,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
      verticesMatrix = mergeTwoMatricesByKeepAllPrimary(colored,nodesHighShape, 
                                            missinglabel=normShape,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
    }
    verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")
    
   
    #---------------------------- output adj matrix ------------------------------------------
    #
    write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # title row
    write.table(t(as.matrix(uniquenames)), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# the third column of netpairsWtype is the link type 
#  
# highlightNodes can be a vector of genes or a list of gene sets with
#    highColor and highShape for each gene set, i.e.,
#
#  highlightNodes[[i]] <-  highColor[i] and highShape[i]
#
makeSNP =function(netpairsWtype, highlightNodes=NULL,
           edgecolorlevels,
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, snafile="tmp.sna")
{

    pcols = dim(netpairsWtype)[2]

    # consider both columns
    uniquenames    = union(as.character(netpairsWtype[,1]), as.character(netpairsWtype[,2]) )
    uniquenames    = sort(uniquenames)
    no.uniquenames = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    

    # 0. make link index: A1 A2 T & A I
    # A1 A2 T & A I ==> A1 A2 T I1
    #
    leftIdx = merge(netpairsWtype, name2idxMatrix, by.x=1, by.y=1, all=F) 
    
    #  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
    #
    allIdx  = merge(leftIdx, name2idxMatrix, by.x=2, by.y=1, all=F)

    no.pairs = dim(allIdx)[1]

    if(pcols==2){
      no.elevels = length(edgecolorlevels)
      greyIdx    = c(1:no.elevels) [edgecolorlevels=="grey"]
      linksIdxMatrix = cbind( allIdx[,c(3,4)], rep(greyIdx,no.pairs) )
    } else{
      linksIdxMatrix = allIdx[,c(4,5,3)]
    }


    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix 
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      verticesMatrix = matrix("", no.uniquenames, 3)
      verticesMatrix[,1] = uniquenames
      verticesMatrix[,2] = rep(normColor,no.uniquenames)
      verticesMatrix[,3] = rep(normShape,no.uniquenames)

      if ( !is.list(highlightNodes) ){
          # get indices of highlight nodes
          #
          merged = merge(name2idxMatrix, highlightNodes, by.x=1, by.y=1,all=F)
          merged = as.matrix(merged)
          highNodesIdx   = as.integer(merged[,2])

          # set color and shape for highlighted nodes
          #
          verticesMatrix[highNodesIdx, 2] = highColor
          verticesMatrix[highNodesIdx, 3] = highShape
      }else{
         no.highnodeSets = length(highlightNodes)
         for(il in c(1:no.highnodeSets) ) {
             # get indices of highlight nodes
             #             
             merged = merge(name2idxMatrix, highlightNodes[[ il ]], by.x=1, by.y=1,all=F)
             merged = as.matrix(merged)
             highNodesIdx   = as.integer(merged[,2])

             # set color and shape for highlighted nodes
             #
             verticesMatrix[highNodesIdx, 2] = highColor[il]
             verticesMatrix[highNodesIdx, 3] = highShape[il]
         }
      }

    }

    #verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")

    #**************************************************************************
    #
    # 3. output indexed netpairs
    #
    write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # link pairs based index   
    #
    write.table(rbind(c("src","dst", "type")), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    write.table(linksIdxMatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
}


if(F){  # change A at positions on (row,col) pairs with value on the third column
  a=matrix(1,3,3)
  a[1,2]=3
  a[3,3]=4
  a[2,1]=2
  a[1,3]=-1
  a[2,3]=-2
  a[3,3]=-3
  b[a[,c(1,2)]] = a[,3]
}



#
#**********************************  END  **********************************




# Here we assume that both input tables don't have the title row
#
# transform each pairof gene ids in each link in a pairlist file into
#  the coressponding gene symbols, if no match, we use the original IDs
#  so node link will be missed in the out put
#
mapPairlinkIDS2genesymbol = function(fpairlink, inpath="",fmapping, myoutdir="")
{

    # --------- 1. read in network (pair list) ---------------------------------
    #
    fullname = paste(inpath, fpairlink, sep="")
    keywords = getFileName(fpairlink)

    ofname = paste(myoutdir, keywords, "-gsymb.pair", sep="")

    aaMatrix <- read.delim(fullname, sep="\t", header=F)
    aaMatrix =  as.matrix(aaMatrix)
    dim(aaMatrix)

    # --------- 2. read in mapping table (MMT-id to Gene Symbol) ---------------
    #
    mapMatrix <- read.delim(fmapping, sep="\t", header=colheader)
    mapMatrix =  as.matrix(mapMatrix)

    # --------- 3. merge the two t?wice (left node and right node) --------------
    #[A1,A2] & [B1,B2] = [A1,A2,B2L]

    mergeleft = merge(aaMatrix, mapMatrix, by.x=1, by.y=1, all.x=1)
    mergeleft <- as.matrix(mergeleft)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    #
    mergeleft[,3] = ifelse(is.na(mergeleft[,3]), mergeleft[,1], mergeleft[,3])
    dim(mergeleft)


    #[A1,A2,B2] & [B1, B2] = [A2,A1,B2L, B2R]
    #
    mergeright = merge(mergeleft, as.matrix(mapMatrix), by.x=2, by.y=1, all.x=1)
    dim(mergeright)

    mergeright <- as.matrix(mergeright)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    mergeright[,4] = ifelse(is.na(mergeright[,4]), mergeright[,1], mergeright[,4])

    dim(mergeright)

    #---------------------------- output adj matrix ------------------------------------------
    #title row
    
    write.table(mergeright[,c(3,4)], ofname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    return (ofname)
}

# here network can have additional information for each link
#
mapNetworkIDS2genesymbol = function(fnetwork, inpath="",colheader=F, fmapping, myoutdir="")
{

    # --------- 1. read in network (pair list) ---------------------------------
    #
    fullname = paste(inpath, fnetwork, sep="")

    keywords = getFileName(fnetwork)
    extension= getFileExtension(fnetwork)

    ofname = paste(myoutdir, keywords, "-gsymbol.", extension, sep="")

    aaMatrix <- read.delim(fullname, sep="\t", header=colheader)
    aaMatrix =  as.matrix(aaMatrix)
    dim(aaMatrix)

    nocols = dim(aaMatrix)[2]

    # --------- 2. read in mapping table (MMT-id to Gene Symbol) ---------------
    #
    mapMatrix <- read.delim(fmapping, sep="\t", header=F)
    mapMatrix =  as.matrix(mapMatrix)

    # --------- 3. merge the two t?wice (left node and right node) --------------

    colnames( mapMatrix) <- c("dst", "dstGeneSymbol")
    mergeright = merge( mapMatrix, aaMatrix, by.x=1, by.y=2, all.y=1)
    mergeright <- as.matrix(mergeright)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    #
    mergeright[,2] = ifelse(is.na(mergeright[,2]), mergeright[,1], mergeright[,2])
    dim(mergeright)


    #[A1,A2,B2] & [B1, B2] = [A2,A1,B2L, B2R]
    #
    colnames( mapMatrix) <- c("src", "srcGeneSymbol")
    mergeleft = merge(mapMatrix, mergeright, by.x=1, by.y=3, all.y=1)
    dim(mergeleft)

    mergeleft <- as.matrix(mergeleft)
    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    mergeleft[,2] = ifelse(is.na(mergeleft[,2]), mergeleft[,1], mergeleft[,2])

    dim(mergeleft)

    #title row
    
    write.table(mergeleft, ofname, sep="\t",quote=FALSE, col.names=colheader, row.names=FALSE, append=F)

    return (ofname)
}



# perform union and overlap networks between two networks in form of 
#  pair lists and encode the union and overlap networks
#
merge2networksBYpair = function(fnetlist, directeds, keywords, outputDir="") 
{
    codepairs= cbind( c(-1, 1), c(-3,3) )

    flogname = paste(outputDir, keywords, "_log.xls", sep="")
    foutcomb = paste(outputDir, keywords, "_union.adj", sep="")
    foutshare= paste(outputDir, keywords, "_overlap.adj", sep="")
    fshareids= paste(outputDir, keywords, "_overlap-nodeIDs.txt", sep="")
    foutsharpajek= paste(outputDir, keywords, "_overlap4pajek.xls", sep="")

    totalfiles = length(fnetlist)

    nodesInNets=rep(0, totalfiles+2) #plus the union and overlap networks
    edgesInNets=rep(0, totalfiles+2)
    namesofNets=rep("", totalfiles+2)

    allnodenames = NULL
    indexrange  = c(0)
    for (each in fnetlist){
      
      aaMatrix <- read.delim(each, sep="\t", header=F)
      aaMatrix = as.matrix(aaMatrix)

      # consider both columns
      allnodenames = c(allnodenames, as.character(aaMatrix[,1]) )
      allnodenames = c(allnodenames, as.character(aaMatrix[,2]) )

      indexrange  = c(indexrange, length(allnodenames) )

    }

    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    #analyze how many genes overlapped in the two networks
    mergedGmatrix=NULL
    for (i in c(1:totalfiles) ) {
       sidx = indexrange[i] + 1
       eidx = indexrange[i+1]
       
       inametable   = table(allnodenames[c(sidx:eidx) ])
       iuniquenames = as.matrix( names(inametable) )
       if(is.null(mergedGmatrix)){
           mergedGmatrix = iuniquenames
       }else{
           mergedGmatrix = merge(mergedGmatrix, iuniquenames, by.x=1, by.y=1, all=F)
       }
    }


    # initialize matrix with the size equal to no.uniquenames
    # for the first network, edge is represented by 1 and non-edge by -1
    # for the 2nd network, edge is represented by 3 and non-edge by -3
    # so element in the ajacency matrix of the combined network:
    #   = 4: link exists in both networks 
    #   =-2: link exists in 1st network but not 2nd network
    #   = 2: link exists in 2nd network but not 1st network
    #   =-4: link doesn't exist in any of the two networks
    #
    cname2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    combinedMatrix = NULL
    for (i in c(1:totalfiles) ){
      each     = fnetlist[i]
      fname    = getFileName(each)
      outfname = paste(outputDir, fname, sep='')

      aaMatrix <- read.delim(each, sep="\t", header=F)
      dim(aaMatrix)

      iadj = makeAjacencyMatrix(inetmatrix=aaMatrix, coding=codepairs[,i], 
                                matrixsize=no.uniquenames, 
                                directed=directeds[i], 
                                name2idxMatrix=cname2idxMatrix)

      if (is.null(combinedMatrix)){
        combinedMatrix = iadj
      }else{
        combinedMatrix = combinedMatrix + iadj
      }

      # we output the summary information about the current network
      aaMatrix = as.matrix(aaMatrix)
      inodenames = c(as.character(aaMatrix[,2]), as.character(aaMatrix[,1])) 
      inametable = table(inodenames)
      iuniquenames     = names(inametable)
      ino.uniquenames  = length(iuniquenames)

      nodesInNets[i] = ino.uniquenames
      edgesInNets[i] = dim(aaMatrix)[1]
      namesofNets[i] = fname

      # output the node Ids for the current input networks
      #  ifids    = paste(outputDir, fname, "-nodeIDs.txt", sep='')
      ifids    = paste(fname, "-nodeIDs.txt", sep='')
      write.table(as.matrix(iuniquenames), ifids, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)

    }


    # statistics about the union network
    nodesInNets[totalfiles + 1] = no.uniquenames
    edgesInNets[totalfiles + 1] = sum( (combinedMatrix != -4) )
    namesofNets[totalfiles + 1] = "union network"

    # change the coding for simplicity
    # element in the final ajacency matrix of the combined network:
    #   = 4: ( 4)link exists in both networks 
    #   = 1: (-2) link exists in 1st network but not 2nd network
    #   = 2: ( 2) link exists in 2nd network but not 1st network
    #   = 0: (-4) link doesn't exist in any of the two networks
    #
    #need a lot of memory, so we put it in the later part
    #
    #combinedMatrix = ifelse(combinedMatrix==-4, 0, 
    #                        ifelse(combinedMatrix==-2, 1, combinedMatrix) )


    #---------------------------- output combined adj matrix ------------------------------------------

    #title row
    write.table(t(as.matrix(uniquenames)), foutcomb, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    # write into file in many times for saving memory
    no.rows=dim(combinedMatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=combinedMatrix[irows, ]
       imatrix=ifelse(imatrix==-4, 0, 
                            ifelse(imatrix==-2, 1, imatrix) )

       write.table(imatrix, foutcomb, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

       #matched edges
       bmatrix= imatrix==4
       
       total_matched = total_matched + sum(bmatrix)
    }


    #---------------------------- output shared adj matrix ----------------------------------------
    bool.matrix = combinedMatrix==4

    rows.sum    = apply(bool.matrix, 1, sum)
    cols.sum    = apply(bool.matrix, 2, sum)
    shared.bool = (rows.sum>0) | (cols.sum>0)

    sharedNames    = uniquenames[shared.bool]
    no.sharednodes = length(sharedNames)
    sharedMatrix   = combinedMatrix[shared.bool, shared.bool]
    #sum(sharedMatrix==4)

    sharedMatrix = ifelse(sharedMatrix==4, 1, 0)
    #sum(sharedMatrix==1)

    #title row
    write.table(t(as.matrix(sharedNames)), foutshare, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(sharedMatrix, foutshare, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # node IDs
    write.table(as.matrix(sharedNames), fshareids, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)


    #-------------------------- shared matrix for Pajekificator ---------------------------------
    if(1==2){
    colors=rep(13, no.sharednodes)
    cc    =rep(0.5,no.sharednodes)
    pajekMatrix = cbind(sharedNames, cc, colors, as.matrix(sharedMatrix) )

    titlerow=c("pajek","cc","colorCode", sharedNames)

    #title row
    write.table(t(as.matrix(titlerow)), foutsharpajek, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(pajekMatrix, foutsharpajek, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    # statistics about the shared network
    nodesInNets[totalfiles + 2] = no.sharednodes
    edgesInNets[totalfiles + 2] = sum( sharedMatrix )
    namesofNets[totalfiles + 2] = "overlap network"

    #---------------------------- output statistics --------------------------------------------
    # matched percentage
    #
    matchedpercents = signif(total_matched/edgesInNets, 4)
    matchedpercents = matchedpercents*100
    matcheds        = rep(total_matched, totalfiles+2)
    matchMatrix     = cbind(namesofNets, nodesInNets, edgesInNets, matcheds, matchedpercents)
    fmatrix = rbind( c("network", "# of nodes", "# of edges", 
                     "# of matched edges", "matched percentage (%)") ,
                     matchMatrix)

    write.table(fmatrix, flogname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    sharedgenes = dim(mergedGmatrix)[1]
    mystring = paste("\nnumber of common nodes in two networks =", sharedgenes)
    appendStringToFile(flogname, mystring)


    #mystring = paste("\nnumber of matched edges=", total_matched)
    #appendStringToFile(flogname, mystring)
    #cat("total matched edges are ", total_matched, "\n")
    fmatrix
}




############################# Causality test #######################################
#

# examine two models 1) G->T 2) T->G; the assumption here is that G and T are both
#   controlled by L, eg, 
# for each model, we use the following rationale to 
# for the 1st model, if lm(residual(G->T) ~ L) is significant, then we REJECT the 
#  model G->T
#
getPv <- function(L, H, G, T) {
            # G --> T
            lm1T <- lm( T ~ G ) ;
            #lm1Ta <- anova(lm1T) ;
            resT <- resid(lm1T) ;
            lm2T <- lm( resT ~ L + H ) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

            # T --> G
            lm1G <- lm( G ~ T ) ;
            #lm1Ga <- anova(lm1G) ;
            resG <- resid(lm1G) ;
            lm2G <- lm( resG ~ L + H ) ;
            lm2Gs <- summary(lm2G) ; 
            #rsqrG <- lm2Gs$r.squared ; #ntmp <- length(lm2Gs$residuals) ;
            p2G <- 1 - pf(lm2Gs$fstatistic[1], lm2Gs$fstatistic[2], lm2Gs$fstatistic[3])
   return(c(p2G, p2T))
}

getRst <- function(p2G, p2T, pvth1, pvth2) {
            rst <- 0
            if (p2T > pvth1 & p2G < pvth2) rst <- 1 ; #accept G->T & reject T->G
            if (p2T < pvth2 & p2G > pvth1) rst <- 2 ; #accept T->G & reject G->T
            if (p2T < pvth2 & p2G < pvth2) rst <- 3 ; #reject both G->T & T->G, G & T are independent, but both dependent on L
   return(rst)
}

bootstrap <- function(L, H, G, T, callrst, id = 1, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

    ntmp <- length(T) ; 
    sst1 <- NULL ; sst2 <- NULL ; sst3 <- NULL ;
    for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- getPv(L[nni], H[nni], G[nni], T[nni]) ; 
      sst1 <- c(sst1, getRst(tmp[1], tmp[2], .05, .05)) 
     }

    relia <- NULL
    relia$eric <- sum(sst1 == callrst$eric)/nkkk
    if (id >= 2) relia$aic <- sum(sst2 == callrst$aic)/nkkk
    if (id >= 3) relia$bic <- sum(sst3 == callrst$bic)/nkkk
    return(relia)
}
strlookup <- function(str1, str2) {
  q <- vector("numeric", length(str1))
  for (i in 1:length(str1)) {
    qq <- which(str1[i] == str2) ;
    if (length(qq) > 0) q[i] <- qq[1] 
  }
  return(q)
}

# return (3,-2) means that too many missing values
#
causalityTest=function(LL, HH, GG, TT, Bootstrap=F, noboots=10){

  # here we have to preprocess the data to amke sure that no missing values there
  sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  L=LL[sel]
  H=HH[sel]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(LL)/2){
    ret = c(3, -2)
    return(ret) 
  }

  # causality test
  pv <- getPv(L, H, G, T) ; 
  callrst <- NULL;
  callrst$eric <- getRst(pv[1], pv[2], .05, .05);

  #bootstrap
  if (Bootstrap){
     bs <- bootstrap(L, H, G, T, callrst, 1, noboots);
     ret=c(callrst$eric, bs$eric)
  }else{
     ret = c(callrst$eric, -1)
  }
  
  return (ret)
}

##########################################################################
# identify QTL for one gene at one chromosome
# lod score should be ranked from small position to large
# th1 = 1, th2 = 2 for default
# rst has 3 fields: LOD score, max position and index for each QTL
#
# th1 is used to define the bell shape boundary
# th2 is used to define the min QTL LOD
#
identifyQTLorg <- function(lod, pos, th1 = 1, th2 = 2) {

  q1 <- lod > th1 ;
  q2 <- lod <= th1 ;

  tmp <- c(0, q2[1:(length(q2)-1)]) ;
  tmp1 <- which(tmp + q1 == 2) ;
  tmp <- c(q2[2:length(q2)], 0) ;
  tmp2 <- which(tmp + q1 == 2) ;

  if (lod[1] > th1) tmp1 <- c(1, tmp1)
  if (lod[length(lod)] > th1) tmp2 <- c(tmp2, length(lod))

  rst <- NULL;
  rst$lod <- NULL ; rst$pos <- NULL ; rst$idx <- NULL ;
  rst$idx1 <- NULL ; rst$idx2 <- NULL ;
  if (length(tmp1) == 0) return(rst)
  for (i in 1:length(tmp1)) {
    tmp <- max(lod[tmp1[i]:tmp2[i]]) ;
    I   <- which.max(lod[tmp1[i]:tmp2[i]]) ;
    if (tmp > th2) {

        j <- tmp1[i] + I - 1;
        j1 <- max(c(1, j-1)) ;
        j2 <- min(c(j+1, length(pos)))

        if (sum(lod[j1:j2] > th1) < length(j1:j2)) next

        rst$lod  <- c(rst$lod, tmp) ;
        rst$pos  <- c(rst$pos, pos[j]) ;
        rst$idx  <- c(rst$idx, j) ;
        rst$idx1 <- c(rst$idx1, tmp1[i]) ;
        rst$idx2 <- c(rst$idx2, tmp2[i]) ;
    }
  }
  return(rst)
}

# ##############################################################################
# 
# returned results are a matrix of QTL info and the columns are in the order of
# 
# QTLlod QTLchrom  QTLpos    QTLmarkeridx   QTLmarkerLeftBoundIdx     QTLmarkerRightBoundIdx
#

identifyQTL <- function(lod, pos, chrom, lodth1 = 1, lodth2 = 2) {

  q1 <- lod > lodth1 ;
  q2 <- lod <= lodth1 ;

  tmp <- c(0, q2[1:(length(q2)-1)]) ;
  tmp1 <- which(tmp + q1 == 2) ;
  tmp <- c(q2[2:length(q2)], 0) ;
  tmp2 <- which(tmp + q1 == 2) ;

  if (lod[1] > lodth1) tmp1 <- c(1, tmp1)
  if (lod[length(lod)] > lodth1) tmp2 <- c(tmp2, length(lod))

  rst <- NULL;

  # return if no LOD bigger than the thresholds
  if (length(tmp1) == 0) return(rst)
  peaks = (lod > lodth2)
  if ( sum(peaks) == 0) return(rst)

  for (i in 1:length(tmp1)) {
    tmp <- max(lod[tmp1[i]:tmp2[i]]) ;
    I   <- which.max(lod[tmp1[i]:tmp2[i]]) ;
    if (tmp > lodth2) {

        j  <- tmp1[i] + I - 1;
        j1 <- max(c(1, j-1)) ;
        j2 <- min(c(j+1, length(pos)))

        if (sum(lod[j1:j2] > lodth1) < length(j1:j2)) next

        # notice that here cbind generate a matrix instead of a vector
        #
        #            lod  chrom  pos    idx   idx1     idx2
        ires = cbind(tmp, chrom, pos[j], j,  tmp1[i],  tmp2[i])
        rst  = rbind(rst, ires)
    }
  }

  return(rst)
}

############################################################################################
# 
# search for genome-wide QTLs of a given gene
#
# returned results are a matrix of QTL info and the columns are in the order of
# 
# QTLlod QTLchrom QTLpos QTLmarkIdx QTLmarkLeftIdx QTLmarkRightIdx geneChrom genePosBeg genePosEnd
#

identifyQTLGenome <- function(LOD, POS, CHROM, geneChrom=NA, genePosBeg=NA, genePosEnd=NA,Llod=1,Hlod=2,maxcisdist=20000) {
  allres <- NULL;
  
  # how many chromosomes here
  clevels = levels(factor(CHROM))
  chroms  = as.integer(clevels)

  #check LOD by chromosome
  #
  index = c( 1:length(LOD) )# marker index
  for (ec in chroms){
     sel     = (CHROM == ec)
     selnona = sel & (!is.na(sel) )

     marksIdx = index[selnona]

     # returned results are a matrix of QTL info and the columns are in the order of
     # 
     # QTLlod QTLchrom  QTLpos    QTLmarkeridx   QTLmarkerLeftBoundIdx     QTLmarkerRightBoundIdx

     ecRes = identifyQTL(lod=LOD[selnona], pos=POS[selnona], chrom=ec, lodth1 = Llod, lodth2 = Hlod)

     # get absolute Index of markers
     for (mi in c(4:6) ){
         ecRes[,mi] = marksIdx [ ecRes[,mi] ]
     }

     allres= rbind (allres, ecRes)
  }

  if (is.null(allres) ){
     return (NULL)
  }

  noQTLs = dim(allres)[1]

  extracols = cbind( rep(geneChrom,noQTLs), rep(genePosBeg,noQTLs),  rep(genePosEnd,noQTLs) )

  finalres = cbind(allres, extracols)

  return (finalres)
}

# compare two intervals to see if they significantly overlap each other
#
qtlOverlap = function(intervalA, intervalB, minoverlap=0.6)
{
  # 1) there are no overlap at all
  if ( (intervalA[2] <= intervalB[1]) | (intervalA[1]>=intervalB[2]) ){
     return (F)
  }

  # 2) one is inside the other
  #
  if ( (intervalA[1]>=intervalB[1]) & (intervalA[2]<=intervalB[2]) ){
     return (T)
  }
  if ( (intervalB[1]>=intervalA[1]) & (intervalB[2]<=intervalA[2]) ){
     return (T)
  }
  
  #3) there are some overlap, then look at the percentage of overlap
  oS = max(intervalA[1], intervalB[1])
  oE = min(intervalA[2], intervalB[2])
  od = oE-oS

  ovlpA = od/(intervalA[2]- intervalA[1])
  ovlpB = od/(intervalB[2]- intervalB[1])

  if ( (ovlpA>minoverlap) | (ovlpB>minoverlap) ){
      return (T)#(c(T, ovlpA, ovlpB ) )
  }else{
      return (F)#(c(F, ovlpA, ovlpB ) )
  }

}

#
#
############################# END of Causality Test #######################################

# plot SNA format
#
# so, the adj matrix in SNA can be either "the link colors" or binary values 
#  representing connection/disconnection, we need differentiate the two situations
#

plotNetwork = function(input, directed, fimg="tmp.png", disphubs=3, plotfigure=T, nodenamecolor="purple", labelpos=0){

    #filename = getFileName(input)
    #fimg     = paste(filename, ".png",sep="")

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
    nonodes = as.integer(splitted[2])
    nonodes

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    vertexShape = as.integer(as.matrix(vertices$shape) )
    
    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep(nodenamecolor, nonodes)

    # 3) get the link information
    #

    firstrow <- read.delim(input, sep="\t", header=F, skip=nonodes+2, nrows=1)
    firstrow <- as.character(as.matrix(firstrow))

    # so, the adj matrix shows index of the link colors instead of binary values for 
    #  connection/disconnection
    #
    if (firstrow[1]== "edge_color_levels"){

        edge_color_levels= firstrow[-1] #remove flag of "edge_color_levels"

        edges <- read.delim(input, sep="\t", header=T, skip=nonodes+3) #**
        dim(edges)

        enum = as.numeric(as.matrix(edges))
        edges = matrix(enum, nonodes, nonodes)
        
        edgeColor = matrix(edge_color_levels[enum+1], nonodes, nonodes)
        #edges = ifelse(edges ==0, 0, 1)

    } else{
        edges <- read.delim(input, sep="\t", header=T, skip=nonodes+2)
        dim(edges)
        
        edges = matrix(as.numeric(as.matrix(edges)), nonodes, nonodes)

        # consider python vector index starts from '0'
        edgeColor = matrix("white", nonodes, nonodes)
        for (i in c(2:length(linkcolors)) ){
          edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
        }
    }

    dim(edgeColor)

    # find hubs
    edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix
    inlinks   = degree(edges,cmode="indegree")
    outlinks  = degree(edges,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(disphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }


      #make title
      #concatenate hubs and no of links
      nohubs = length(top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = top2hubs[k]
              links= top2links[k]
          }else{
              hubs = paste(hubs,  ", ", top2hubs[k], sep="")
              links= paste(links, ", ", top2links[k], sep="")
          }
      }
      hubsNnets = c(hubs, links)

      if (!plotfigure) {
         rm(edges)
         collect_garbage()
         return (hubsNnets)
      }

      if (disphubs>0){
         ititle = paste("Top hubs: ",hubs,sep="")
      }else{
         ititle = ""
      }

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 3600){
      imgsize = 3600
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)


    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)

    openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    if (disphubs>0){
      par(mfrow=c(1,1), mar=c(0,1,2,1) )
    }else{
      par(mfrow=c(1,1), mar=c(0,1,0,1) )
    }
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=4/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 4/vercexscale,     
             object.scale  = objscale, main=ititle)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()

    rm(edges)
    collect_garbage()

    return (hubsNnets)
}


plotMultiNetworks = function(inputs, directeds, fimg="tmp_multinets.png", disphubs=3, plotfigure=T, labelpos=0){

    if (length(inputs)==1){
       plotNetwork(inputs, directeds)
       return
    }

    nnets    = length(inputs)
    imgsize  = 1000
    nets     = as.list(rep(0, nnets) )

    hubsNnets = NULL

    for (i in c(1:nnets) ) {
      nets[[i]] = plotNetworkCore(inputs[i], directeds[i], idisphubs=disphubs, rimagesize = imgsize)

      #concatenate hubs and no of links
      nohubs = length(nets[[i]]$top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = nets[[i]]$top2hubs[k]
              links= nets[[i]]$top2links[k]
          }else{
              hubs = paste(hubs,  ", ", nets[[i]]$top2hubs[k], sep="")
              links= paste(links, ", ", nets[[i]]$top2links[k], sep="")
          }
      }
      hubsNnets = c(hubsNnets, hubs, links)
    }

    if (!plotfigure){
       rm(nets)
       #collect_garbage()
       return (hubsNnets)
    }

    openImgDev(fimg,iwidth = nnets*imgsize, iheight = imgsize, ipointsize = 12)

    if (disphubs>0){
      par(mfrow=c(1,1), mar=c(0,1,2,1) )
    }else{
      par(mfrow=c(1,1), mar=c(0,1,0,1) )
    }

    #plotNetworkCore(inputs[i], directeds[i], imagesize = imgsize)

    for (i in c(1:nnets) ) {
      #make title
      if (disphubs>0){
         ititle = paste("Top hubs: ",hubsNnets[1+(i-1)*2],sep="")
      }else{
         ititle = ""
      }

      gplot( nets[[i]]$edges,                  diag=T,
             gmode=nets[[i]]$gmod,  #mode="target",#"kamadakawai",
             pad=2,                label.pad=0.1,
             vertex.cex=4/nets[[i]]$vercexscale, vertex.sides=nets[[i]]$vertexShape,
             vertex.col=nets[[i]]$vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=nets[[i]]$edgeColor,
             edge.lty  = 0,          edge.lwd     = 1,
             displaylabels=TRUE,     label.bg="white",label.col=nets[[i]]$labColor,
             label=nets[[i]]$vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=0.8/nets[[i]]$labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/nets[[i]]$vercexscale,
             arrowhead.cex = 2/nets[[i]]$vercexscale,     
             object.scale  = nets[[i]]$objscale,main=ititle) 
      #text(0, 150, ititle, cex =nets[[i]]$labcexscale)
    }
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    dev.off()

    rm(nets)
    collect_garbage()

    return (hubsNnets)
}


plotNetworkCore = function(input, directed, idisphubs=3, rimagesize=1000){

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
    nonodes = as.integer(splitted[2])
    nonodes

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    vertexShape = as.integer(as.matrix(vertices$shape) )

    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep("purple", nonodes)

    # 3) get the link information
    #

    edges <- read.delim(input, sep="\t", header=T, skip=nonodes+2)
    dim(edges)

    edges = matrix(as.numeric(as.matrix(edges)), nonodes, nonodes)

    # consider python vector index starts from '0'
    edgeColor = matrix("white", nonodes, nonodes)
    for (i in c(2:length(linkcolors)) ){
      edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
    }
    dim(edgeColor)

    # find hubs
    edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix
    inlinks   = degree(edges,cmode="indegree")
    outlinks  = degree(edges,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(idisphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }


    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 3600){
      imgsize = 3600
    }

    pr = rimagesize/imgsize

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(rimagesize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes*pr-5)

    objcnst     = 0.08
    objinter    = pr*(0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes*pr-5)
    objscale    = objscale/(1+vercnst^0.5)

    netpara=list(edges=edges,
               gmod=gmod,
               vercexscale =vercexscale,
               vertexShape =vertexShape,
               vertexColor=vertexColor,
               edgeColor   =edgeColor,
               labColor    =labColor,
               vertexName  =vertexName,
               labcexscale =labcexscale,
               objscale    =objscale,
               top2hubs    =top2hubs,
               top2links   =top2links)
               
   return (netpara)
}

#------------------------ plot from a matrix with pairs -------------------------------------
#
# Here we provide the gene information matrix and specify the gene symbol column which to be
#  shown in the network
#
# Also output in/out/total links information as "*_inoutlinks.xls"
#

plotNetwork_inPairs_DisplayGeneSymbol = function(linkpairs, directed, geneinfo=NULL, genesymbolCol=2,
                   nodecolor="blue", nodeshape=30, nodenamecolor="purple", labelpos=1,
                   rankhubby="totallinks",
                   disphubs=3, plotfigure=T, fimg="tmp.png", 
                   saveDegrees=T, colorNodes=NULL, shapeNodes=NULL, shownodename=T, placemode = "fruchtermanreingold", default_title=""){

    #filename = getFileName(input)
    #fimg     = paste(filename, ".png",sep="")

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    netrows <- dim(linkpairs)[1]
    netcols <- dim(linkpairs)[2]

    #------------------------------------------------------------------------------------
    #                        make network adjacency matrix
    #

    # find unique node names 
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    # find corrsponding gene symbols
    if (!is.null(geneinfo)) {
       unmatrix <- cbind(uniquenames)
       colnames(unmatrix) <- c("mmtid")

       ordergeneinfo  = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=unmatrix, 
                                         minorMatrix=geneinfo, missinglabel="", 
                                         keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
       ordergeneinfo = as.matrix(ordergeneinfo)

       # use MMTid if there is no symbol found
       ordergeneinfo[,genesymbolCol] = ifelse(is.na(ordergeneinfo[,genesymbolCol]),
                                              ordergeneinfo[,1],
                                              ordergeneinfo[,genesymbolCol])

       uniquesymbols = as.character(ordergeneinfo[,genesymbolCol])
    } else{
       ordergeneinfo = cbind(uniquenames)
       uniquesymbols = uniquenames
    }

    # look at the third column which is the confidence value 
    # if missing, add one more column 
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    if (netcols ==3) {
       mycoding =NULL
    } else{
       mycoding = c(0, 1)
    }    

    adjmatrix = makeAjacencyMatrix(inetmatrix=linkpairs, 
                                   coding=mycoding,
                                   matrixsize=no.uniquenames, 
                                   myname2idxMatrix= name2idxMatrix,
                                   directed=directed)
    
    # -------------- find hubs --------------------------
    #
    inlinks   = degree(adjmatrix,cmode="indegree")
    outlinks  = degree(adjmatrix,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }

    # rank by different ways
    if (rankhubby=="inlinks"){
       hubidx    = order(-inlinks)
    }else if (rankhubby=="outlinks"){
       hubidx    = order(-outlinks)
    }else{
       hubidx    = order(-totallinks)
    }

    #-- make title using top hubs ---------------------
    #
    #concatenate hubs and no of links
    hubnames = uniquesymbols[hubidx]
    if(disphubs>0){
       tophubs  = concatenate(hubnames[c(1:disphubs)],", ")

      if (rankhubby=="outlinks"){
         ititle = paste("Top causal hubs: ", tophubs,sep="")
      } else{
        ititle = paste("Top hubs: ", tophubs,sep="")
      }
    }else{
        ititle = default_title
    }

    # ++++++++++  save the link matrix +++++++++++++++++++++++
    #
    linkmatrix = cbind(ordergeneinfo, inlinks, outlinks, totallinks)
    colnames(linkmatrix) <-c( colnames(ordergeneinfo), "inlinks", "outlinks", "total_links")
    linkmatrix = linkmatrix [hubidx,]
    if (!is.null(fimg) & saveDegrees ) {
       myfn = getFileName(fimg)
       flink= paste(myfn, "_inoutlinks.xls",sep="")
       write.table(linkmatrix, flink, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }


    # 1) get the number of nodes in the network 
    #
    nonodes = no.uniquenames
    nonodes

    # 2) get the node information
    #
    if ( shownodename){
       vertexName  = uniquesymbols
    } else{
       vertexName  = rep("",nonodes)
    }

    if(length(nodecolor)>1){
         vertexColor = nodecolor
    }else{
         vertexColor         = rep(nodecolor, nonodes)
         if (is.null(colorNodes)) {
            if (disphubs>0){
               vertexColor[hubidx[c(1:disphubs)] ] = "red" #highlight hubs
            }#else{
             #  vertexColor[hubidx[1] ] = "red" #highlight hubs
             #}
          }else{
            colormerged = merge(name2idxMatrix, colorNodes, by.x=1,by.y=1,all=F)
            colormerged = as.matrix(colormerged)
            coloridx    = as.integer(colormerged[,2])
            vertexColor[coloridx] = "red"
          }
    }

    if(length(nodeshape)>1){
         vertexShape = nodeshape
    } else {
         vertexShape = rep(nodeshape, nonodes)
         if (!is.null(shapeNodes)) {
            shapemerged = merge(name2idxMatrix, shapeNodes, by.x=1,by.y=1,all=F)
            shapemerged = as.matrix(shapemerged)
            shapeidx    = as.integer(shapemerged[,2])
            vertexShape [shapeidx] = 3
         }
    }   

    if ( length(nodenamecolor)>1){
       labColor    = nodenamecolor
    }else{
       labColor    = rep(nodenamecolor, nonodes)
    }

    # 3) get the link information
    #
    edges = matrix(as.numeric(as.matrix(adjmatrix)), nonodes, nonodes)

    # consider python vector index starts from '0'
    if (no.uniquenames < 4000) {# to save memory
       edgeColor = matrix("white", nonodes, nonodes)
       for (i in c(2:length(linkcolors)) ){
         edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
       }
       dim(edgeColor)
    } else{
       edgeColor = "gray"
    }

    rm(adjmatrix)
    collect_garbage()

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 2000){ #3600
      imgsize = 2000
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)


    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)


    if (!is.null(fimg)) {
       openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    }
    par(mfrow=c(1,1), mar=c(0,1,2,1) )
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=3/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 3/vercexscale,     
             object.scale  = objscale, main=ititle, mode = placemode)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if (!is.null(fimg)) {
       dev.off()
    }

    rm(edges)
    collect_garbage()

    return (linkmatrix )
}


#********************************************************************************************
#------------------------ plot network in form of matrix -------------------------------------
#
# 1) different from network inpairs where we can only show nodes with at least one link
#    here we can show isolated nodes
#
# 2) geneinfo contains the columns in the order of: nodeIndex, node name, node size
#
# 3) here the link color is proportional to the value
#
# 3) the rest of output is the same as  plotNetwork_inPairs_DisplayGeneSymbol 
#

plotNetwork_inMatrix_DisplayGeneSymbol = function(adjmatrix, directed, geneinfo=NULL, genesymbolCol=2,
                   nodecolor="blue", nodeshape=30, nodenamecolor="purple", labelpos=1,
                   rankhubby="totallinks",
                   disphubs=3, plotfigure=T, fimg="tmp.png", 
                   saveDegrees=T, colorNodes=NULL, shownodename=T, scaleVertex=F,
                   placemode = "fruchtermanreingold", default_title=""){

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    no.uniquenames= dim(adjmatrix)[1]

    # scale up the size of vertex
    #
    nodesize      = rep(1, no.uniquenames)
    if (!is.null(geneinfo) ){
      if (scaleVertex & (dim(geneinfo)[2]==3) ) { 
         nodesize      = as.integer(geneinfo[,3])
         nodesize      = as.integer(5*nodesize/max(nodesize))
         nodesize      = ifelse(nodesize <1, 1, nodesize)
      }
      ordergeneinfo = geneinfo
      uniquesymbols = geneinfo[,2]
      uniquenames   = geneinfo[,2]
    } else{      
      uniquesymbols = colnames(adjmatrix)
      uniquenames   = colnames(adjmatrix)
      ordergeneinfo = cbind(uniquesymbols , uniquesymbols)
    }

    # -------------- find hubs --------------------------
    #
    inlinks   = degree(adjmatrix,cmode="indegree")
    outlinks  = degree(adjmatrix,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }

    # rank by different ways
    if (rankhubby=="inlinks"){
       hubidx    = order(-inlinks)
    }else if (rankhubby=="outlinks"){
       hubidx    = order(-outlinks)
    }else{
       hubidx    = order(-totallinks)
    }

    #-- make title using top hubs ---------------------
    #
    #concatenate hubs and no of links
    hubnames = uniquesymbols[hubidx]
    if(disphubs>0){
       tophubs  = concatenate(hubnames[c(1:disphubs)],", ")

      if (rankhubby=="outlinks"){
         ititle = paste("Top causal hubs: ", tophubs,sep="")
      } else{
        ititle = paste("Top hubs: ", tophubs,sep="")
      }
    }else{
        ititle = default_title
    }

    # ++++++++++  save the link matrix +++++++++++++++++++++++
    #
    linkmatrix = cbind(ordergeneinfo, inlinks, outlinks, totallinks)
    colnames(linkmatrix) <-c( colnames(ordergeneinfo), "inlinks", "outlinks", "total_links")
    linkmatrix = linkmatrix [hubidx,]
    if (!is.null(fimg) & saveDegrees ) {
       myfn = getFileName(fimg)
       flink= paste(myfn, "_inoutlinks.xls",sep="")
       write.table(linkmatrix, flink, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }


    # 1) get the number of nodes in the network 
    #
    nonodes = no.uniquenames
    nonodes

    # 2) get the node information
    #
    if ( shownodename){
       vertexName  = uniquesymbols
    } else{
       vertexName  = rep("",nonodes)
    }

    if(length(nodecolor)>1){
         vertexColor = nodecolor
    }else{
         vertexColor         = rep(nodecolor, nonodes)
         if (is.null(colorNodes)) {
            if (disphubs>0){
               vertexColor[hubidx[c(1:disphubs)] ] = "red" #highlight hubs
            }#else{
             #  vertexColor[hubidx[1] ] = "red" #highlight hubs
             #}
          }else{
            colormerged = merge(name2idxMatrix, colorNodes, by.x=1,by.y=1,all=F)
            colormerged = as.matrix(colormerged)
            coloridx    = as.integer(colormerged[,2])
            vertexColor[coloridx] = "red"
          }
    }

    if(length(nodeshape)>1){
         vertexShape = nodeshape
    } else {
         vertexShape = rep(nodeshape, nonodes)
    }   

    if ( length(nodenamecolor)>1){
       labColor    = nodenamecolor
    }else{
       labColor    = rep(nodenamecolor, nonodes)
    }

    # 3) get the link information
    #
    edges = matrix(as.numeric(as.matrix(adjmatrix)), nonodes, nonodes)

    # consider python vector index starts from '0'
    if (no.uniquenames < 4000) {# to save memory
       edgeColor = matrix("white", nonodes, nonodes)
       for (i in c(2:length(linkcolors)) ){
         edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
       }
       dim(edgeColor)
    } else{
       edgeColor = "gray"
    }

    rm(adjmatrix)
    collect_garbage()

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 2000){ #3600
      imgsize = 2000
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)

    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)

    if (!is.null(fimg)) {
       openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    }

    par(mfrow=c(1,1), mar=c(0,1,2,1) )
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=nodesize*3/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 3/vercexscale,     
             object.scale  = objscale, main=ititle, mode = placemode)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if (!is.null(fimg)) {
       dev.off()
    }

    rm(edges)
    collect_garbage()

    return (linkmatrix )
}

#---------------------- plot SNP ----------------------------------------
plotNetworkSNP = function(input, directed, fimg="tmp.png", disphubs=3, 
                          nodenamecolor="purple", plotfigure=T, labelpos=0){

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
    nonodes = as.integer(splitted[2])
    nonodes

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    vertexShape = as.integer(as.matrix(vertices$shape) )
    
    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep(nodenamecolor, nonodes)

    # 3) get the link information
    #
    firstrow <- read.delim(input, sep="\t", header=F, skip=nonodes+2, nrows=1)
    firstrow <- as.character(as.matrix(firstrow))

    # so, the adj matrix shows index of the link colors instead of binary values for 
    #  connection/disconnection
    #
    edge_color_levels= firstrow[-1] #remove flag of "edge_color_levels"

    linksIdx <- read.delim(input, sep="\t", header=T, skip=nonodes+3) #**
    dim(linksIdx)
    li.cols = dim(linksIdx)[2]
    li.rows = dim(linksIdx)[1]
    linksIdx= matrix(as.integer(as.matrix(linksIdx )), li.rows, li.cols)

    # assign 1 for each link
    #
    edges = matrix(0, nonodes, nonodes)
    edges[ linksIdx[,c(1:2)] ] = 1

    edgeColor = matrix("white", nonodes, nonodes)
    #edges = ifelse(edges ==0, 0, 1)
    edgeColor[ linksIdx[,c(1:2)] ] =  edge_color_levels[ linksIdx[,3] ]

    #dim(edgeColor)

    # find hubs
    #edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix
    inlinks   = degree(edges,cmode="indegree")
    outlinks  = degree(edges,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(disphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }

      #make title
      #concatenate hubs and no of links
      nohubs = length(top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = top2hubs[k]
              links= top2links[k]
          }else{
              hubs = paste(hubs,  ", ", top2hubs[k], sep="")
              links= paste(links, ", ", top2links[k], sep="")
          }
      }
      hubsNnets = c(hubs, links)

      if (!plotfigure) {
         rm(edges)
         collect_garbage()
         return (hubsNnets)
      }

      if (disphubs>0){
         ititle = paste("Top hubs: ",hubs,sep="")
      }else{
         ititle = ""
      }

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 3600){
      imgsize = 3600
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)


    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)

    edgewid = 1.5
    if(nonodes <200){
      edgewid = 2
    } else if (nonodes >= 3000){
      edgewid = 0.8
    }

    openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    if (disphubs>0){
      par(mfrow=c(1,1), mar=c(0,1,2,1) )
    }else{
      par(mfrow=c(1,1), mar=c(0,1,0,1) )
    }
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=4/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = edgewid,
             displaylabels=T,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 4/vercexscale,     
             object.scale  = objscale, main=ititle)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()

    rm(edges)
    collect_garbage()

    return (hubsNnets)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end of plot networks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### read in Roy's combined expression, trait and genotype data #########
#
ReadBinary<-function(bin.fn)
{
    #
    # read previously created bin file; process; save it as text
    #
    now<-proc.time()[3]
    bin.names<-load(bin.fn)
    t.readTime<-proc.time()[3]-now
    #
    # verify data structure name
    #
    tName<-c("all.data")
    if (bin.names[1]!=tName)
    {
        stop(paste("bin data struct ", tName , " not found.\n\t(",
                    bin.names[1],")", sep=""))
    }
    return(all.data)
}

# get Batch ID from cross, tissue, and sex
#
getBatchId = function(sqlserver, dbname, icross, itissue, isex, flag="CorrRegress", shortenBatchName=F){

  # make a string as a combination of icross, itissue, isex, flag
  ilabel = paste( icross, " ", itissue, " ", isex, " ", flag, sep="")  
  ilabel = tolower(ilabel)

  TBbatch   = "Batch"
  econdition= NULL
  selectedfields = "Batch_name, Batch_id"
  batchTable <- accessWholeTableFromDB(sqlserver, dbname, TBbatch,
                 condition=econdition, outfields = selectedfields)

  # no such a table
  if( is.null(batchTable )) {
     return (-1)
  }

  batchTable <- as.matrix(batchTable)
  no.batches <- dim(batchTable)[1]

  # no.words in label
  lwords = splitString(ilabel,sep=" ")
  no.words= length(lwords)

  # to lower
  for (i in c(1:no.batches) ) {
     #batchTable[i,1] = tolower(batchTable[i,1])

     # here we use only the first "no.words" in the batch name
     # since for mci_bxa, it is hard to specify the long batch name
     #
     # for instance, ilabel="mci_bxa Liver All W10 CorrRegress"
     #
     # mci_bxa Liver  Male W10 CorrRegress Express Adjusted
     # mci_bxa Liver All W10 CorrRegress Adjusted Data
     if (shortenBatchName) {
        iwords = splitString(batchTable[i,1],sep=" ")
        no.min = min(no.words, length(iwords))
        iNwords= concatenate(iwords[1:no.min], " ")
        itmp = tolower(iNwords)
     } else{
        itmp = tolower(batchTable[i,1])
     }

     batchTable[i,1] = replaceChars(itmp, "  ", " ")
     
  }

  merged = merge(batchTable, ilabel, by.x=1,by.y=1, all=F)
  merged = as.matrix(merged)

  # no match  
  if ( dim(merged)[1]==0) {
      return (-1)
  }

  # we take the largest Batch_ID if there are multiple matches
  #
  no.matches = dim(merged)[1]

  return (as.integer(merged[no.matches,2]))
}


# direct means the table is directly uner the connectDB specified
#
accessWholeTableFromDB = function(servername, dbname, tablename, condition=NULL, outfields = "*"){
   channel <- odbcConnect(servername)
   if (is.null(condition)){
     query   <- paste("select ", outfields, " from ", dbname, ".dbo.[", tablename, "]", sep="")
   }else{
     query   <- paste("select ", outfields, " from ", dbname, ".dbo.[", tablename, "] where ", condition, sep="")
   }
   
   tb      <- sqlQuery(channel, query)
   odbcClose(channel)
   return (tb)
}


findMarkersInQTLBycM = function (eQTLonegene, markerinfoMx, markerindex)
{ 
  # look at markers on the same chromosome as the given QTL gene
  #
  #selchrom = eQTLonegene$chrom == markerinfoMx$marker_chrom
  #dist     = abs(eQTLonegene$qtl_max_pos - markerinfoMx$marker_pos[selchrom])
  selchrom  = (eQTLonegene[1] == markerinfoMx$marker_chrom)
  selchrom  = ifelse(is.na(selchrom), F, selchrom)

  if ( sum(selchrom,na.rm=T) == 0 ) {
     return ( c(-1,-1,-1) )
  }

  dist     = abs(eQTLonegene[2] - markerinfoMx$marker_pos[selchrom])

  selindex = markerindex[selchrom]

  # find the nearest marker
  #
  mIdx     = which.min(dist)
  #selindex[mIdx]

  myret = c(selindex[mIdx], (markerinfoMx$marker_pos[selchrom])[mIdx], 
                            (markerinfoMx$Base_pair[selchrom])[mIdx] )

  return ( myret)

}


############# network properties ###############################################

# examine the total number of links, # of nodes, average links per node, scale free fitting R^2 & slope
#
# the input is the first two columns in the network file
# 
# if removeRedantLinks==T, then, we first remove redundant links and then do subsequent analysis
#
getNetwork_Statistics_Scalefree_GeneLink = function(networkinput, colheader=F, directed=F, removeRedantLinks=F, outputDir=""){

    imagetype   ="png"
    fname       =getFileName(networkinput)
    extname     =getFileExtension(networkinput)

    # use different tag, if from non-redundant network
    #
    if(removeRedantLinks) {
      fname = paste(fname, "_uniq", sep="")
    }

    logFname    =paste(outputDir, fname, "-log.xls",       sep='')
    imgScaleFree=paste(outputDir, fname, "-imgScaleFree.",imagetype,sep='')

    codepair= c(0,1)  #[1] for no connection, [2] for connection

    #0) read in pairwise network

    aaMatrix   <- read.delim(networkinput, sep="\t", header=colheader)
    
    #1) remove redundant links and save into corresponding files
    #
    if(removeRedantLinks) {
      alllinks = paste(aaMatrix[,1], aaMatrix[,2], sep="\t")
      
      alinktable = table(alllinks)
      length(alinktable)
    
      uniquepairs = names(alinktable)

      uniquefnet  = paste(fname, ".pair",sep="")
      write.table(uniquepairs, uniquefnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
      aaMatrix <- read.delim(uniquefnet, sep="\t", header=F)
    }

    aaMatrix   <- as.matrix(aaMatrix)    
    dim(aaMatrix)    

    edgesInNet = dim(aaMatrix)[1]
    
    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(aaMatrix[,1]) )
    allnodenames = c(allnodenames, as.character(aaMatrix[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    #name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    #if (no.uniquenames<= 4000){
    if (1==2){
      adjmatrix = makeAjacencyMatrix(inetmatrix=aaMatrix, coding=codepair,
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix=name2idxMatrix)
      edgeIndices=getNetworknodeIdxByNames(netmatrix=aaMatrix, mappingMatrix=name2idxMatrix)

      # find hubs
      inlinks   = degree(adjmatrix,cmode="indegree")
      outlinks  = degree(adjmatrix,cmode="outdegree")
      if (directed){
        totallinks= inlinks + outlinks
      }else{
        totallinks= (inlinks + outlinks)/2 #because of matrix symmetry for undirectecd adj matrix
      }

    }

      totallinks= as.integer(nametable) # no of links for each node

      totalmatrix = cbind(names(nametable),  totallinks)

      if(directed){
        # outlines
        dnodenames = as.character(aaMatrix[,1])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        outlinks = as.integer(as.matrix(iolinks[,2])) 


        # inlines
        dnodenames = as.character(aaMatrix[,2])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        inlinks = as.integer(as.matrix(iolinks[,2])) 

      }else{
        inlinks   = totallinks
        outlinks  = totallinks
      }

    hubidx    = order(-totallinks)


    # let's check whether there is a scale free topology
    suminfor   =ScaleFreePlot(totallinks,no.breaks=40, mtitle="",truncated1=T, tofile="")
    scaleResult=ScaleFreePlot(totallinks,no.breaks=40, mtitle="",truncated1=T, tofile=imgScaleFree, outputFitness=T)

    # ------------------ output log file -----------
    #appendStringToFile(logFname, suminfor)

    # output network statistics 
    #
    mtitle = c("# of links", "# of genes", "links/gene", "scale R^2", "trunc. R^2", "slope")
    
    linkpernode = 2*edgesInNet/no.uniquenames

    finalMatrix = rbind( c(edgesInNet,no.uniquenames, linkpernode, scaleResult) )
    finalMatrix = round(finalMatrix, 2)
    
    colnames(finalMatrix) <- mtitle
    write.table(finalMatrix, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE,append=F)


    # output in/out links for each gene
    #
    linksMatrix            = cbind(uniquenames, inlinks, outlinks, totallinks)
    colnames(linksMatrix) <- c("gene", "inlinks", "outlinks", "totallinks")

    write.table(linksMatrix[hubidx, ], logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    rm(linksMatrix)
    rm(inlinks)
    rm(outlinks)
    rm(totallinks)
    rm(aaMatrix)
    collect_garbage()

    return (finalMatrix)
}


#------------------------------- centrality by pairs ---------------------------------------------
#
# to perform centrality analysis on large scale network, we avoid adjacency matrix
# the output is ordered based on total links
#
degree_ByLinkPairs = function(linkpairs, directed=F, cleangarbage=F){

    codepair= c(0,1)  #[1] for no connection, [2] for connection

    edgesInNet = dim(linkpairs)[1]
    
    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)
    
    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    totallinks  = as.integer(nametable) # no of links for each node
    totalmatrix = cbind(names(nametable),  totallinks)

    if(directed){
        # outlines
        dnodenames = as.character(linkpairs[,1])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        outlinks = as.integer(as.matrix(iolinks[,2])) 


        # inlines
        dnodenames = as.character(linkpairs[,2])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        inlinks = as.integer(as.matrix(iolinks[,2])) 

    }else{
        inlinks   = totallinks
        outlinks = totallinks
    }

    #hubidx    = order(-totallinks)

    # output in/out links for each gene
    #
    linksMatrix            = cbind(inlinks, outlinks, totallinks)
    colnames(linksMatrix) <- c("inlinks", "outlinks", "totallinks")
    rownames(linksMatrix) <- uniquenames

    rm(inlinks)
    rm(outlinks)
    rm(totallinks)

    if(cleangarbage) {
        collect_garbage()
    }

    return ( data.frame(linksMatrix) )
}


centralities=function(adjmatrix, gmode="graph", weighted=F)
{
   # 1. degree
   outdegree = apply(adjmatrix, 1, sum)    
   indegree  = apply(adjmatrix, 2, sum)
   if(gmode=="digraph") {
      totaldegree = outdegree + indegree
   }else{
      totaldegree = outdegree
   }

   #2. betweenness(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE,
   #         tmaxdev=FALSE, cmode="directed", geodist.precomp=NULL, 
   #         rescale=FALSE)
   betweenness = betweenness(adjmatrix, gmod=gmode) 

   #3. closeness(graph, v=V(graph), mode = "all") #"in", "out", "all"
   #
   closeness = closeness(adjmatrix, gmod=gmode, cmode="undirected") 
   
   #4. Clustering coefficient
   clusteringcoef = computeClusterCoefficient(adjmatrix, weighted=weighted)
   
   centrls = data.frame(cbind(indegree, outdegree, 
                              totaldegree, betweenness, closeness, clusteringcoef) )
   return(centrls)
}


# pair a given integer with each element in a vector
#
pairI_to_neighbors = function(si, ivect, removeMinus1=T){
  ino    = length(ivect)
  ipairs = cbind( rep(si,ino), ivect)
  isel   = (ipairs[,1] !=-1) & (ipairs[,2] !=-1)

  if(removeMinus1){
     if(sum(isel)==0){
        return (NULL)
     }else{
        return (ipairs[isel, ]) 
     }
  } else{
     return (ipairs)     
  }
}

pairs_To_Index = function(linkpairs, nodeIndex){
  mergedleft=mergeTwoMatricesByKeepAllPrimary(linkpairs, minorMatrix=nodeIndex, 
                                 missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
  mergedright=mergeTwoMatricesByKeepAllPrimary(mergedleft[,c(2,1,3)], minorMatrix=nodeIndex, 
                                  missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
  linksIndex = mergedright[,c(3,4)]
  linksIndex = as.matrix(linksIndex)
  ino.links  = dim(linksIndex)[1]
  linksIndex = matrix(as.integer(linksIndex), ino.links,2)
  colnames(linksIndex ) =c("src","dst")
  return (linksIndex)
}

# remove duplications including reversed duplications, and self links
#
removeDuplicatedLinks= function(linkpairs, directed=F){

    links = paste(linkpairs[,1], linkpairs[,2],sep="\t")

    # 1. remove duplications 
    #
    cleanedlinkMatrix = union(links, NULL)
    length(cleanedlinkMatrix)

    ofname ="tmp.txt"
    write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)

    # 2. remove inversed duplications
    #    

    linkMatrix <- read.delim(ofname, sep="\t", header=F)
    dim(linkMatrix)
    linkMatrix  = as.matrix(linkMatrix)

    if(directed){
       return(linkMatrix)
    }


    #  first, remove self-interactions are also removed
    #
    selSelfLinks      = linkMatrix[,1]==linkMatrix[,2]
    linkMatrix        = linkMatrix[!selSelfLinks,]
    cleanedlinkMatrix = cleanedlinkMatrix[!selSelfLinks]

    reversedLinks = paste(linkMatrix[,2], linkMatrix[,1],sep="\t")

    no.links = length(reversedLinks)
    reversedLinksWithIndex = cbind(reversedLinks, c(1:no.links) )

    merged = merge(cleanedlinkMatrix, reversedLinksWithIndex, by.x=1,by.y=1,all=F)
    dim(merged)

    if (dim(merged)[1]>0) {
        merged      = as.matrix(merged)
        removedCols = as.integer(merged[,2])

        # cosntruct non-duplicated interactions
        #
        dupLinks    = cleanedlinkMatrix[removedCols]
        dupLinksRev = reversedLinksWithIndex[removedCols]
        uniques     = NULL
        for(i in c(1:length(dupLinks)) )
        {
           found    = is.element(dupLinks[i], uniques)
           foundRev = is.element(dupLinksRev[i], uniques)
           combined = found | foundRev
           if (!combined){
               uniques=c(uniques,dupLinks[i])
           }
        }
        length(uniques)
        write.table( cleanedlinkMatrix[-removedCols], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        write.table( uniques, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE,append=T)
    }else{
        write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
    }

    linkMatrix <- read.delim(ofname, sep="\t", header=F)
    dim(linkMatrix)
    linkMatrix  = as.matrix(linkMatrix)
    return(linkMatrix)
}


# get the network with the links between all nodes in specified by subnet
# keep the original order
# the third column can be index, so the linkpair index is also returned
getSubnetwork_LinkPairs = function(linkpairs, subnetNodes){
   mergeright = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
   if( dim(mergeright)[1] == 0) {
     return (NULL)
   }

   mergeleft  = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

   #mergeright = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   #mergeleft2 = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

   if( dim(mergeleft)[1] == 0) {
     return (NULL)
   }

   return (as.matrix(mergeleft) )   
}


findNLayerNeighbors_LinkPairs = function(linkpairs, subnetNodes, nlayers=1){

   mergeright = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   mergeleft  = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
   mergeright <- as.matrix(mergeright)
   mergeleft  <- as.matrix(mergeleft)
   mergeleft  <- mergeleft[,c(2,1)] # keep the original link direction

   if (nlayers==1){
      final=removeDuplicatedLinks(rbind(mergeright, mergeleft))
      return (final)
   }
   
   ineighbors =union(as.character(mergeright), as.character(mergeleft))

   ret=findNLayerNeighbors_LinkPairs(linkpairs,ineighbors, nlayers-1)

   return (ret)
}

mapPairlinkIDS2genesymbol2=function(fpairlink, inpath, fmapping, myoutdir)
{
   infullname = paste(inpath, fpairlink, sep="")
   linkMatrix <- read.delim(infullname,sep="\t", header=F)
   dim(linkMatrix)
   
   transMatrix <- read.delim(fmapping,sep="\t", header=F)
   dim(transMatrix)

   # src dst srcgs   
   leftmatrix = merge(linkMatrix, transMatrix, by.x=1, by.y=1, all=F)

   # src dst srcgs   & yid gs = dst src srcgs dstgs
   rightmatrix = merge(leftmatrix, transMatrix, by.x=2, by.y=1, all=F)

   fname     =getFileName(fpairlink)
   ofname    =paste(myoutdir, fname, "_gs.pair", sep='')

   write.table( rightmatrix[,c(3,4)], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
       
   return (rightmatrix[,c(3,4)])
}

# Here we use hierarchical clustering to idenitify Connected Components
#
find_ConnectedComponents_linkpairs = function(linkpairs, minModuleSize=10, cclabelstart=1){

    # find unique node names 
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )

    #*-------------------------------------------------------------------------------------
    #* prepare matrix
    #*   
    adjmatrix = makeAjacencyMatrix(inetmatrix=linkpairs[,c(1,2)], 
                                   coding=c(0,1),
                                   matrixsize = no.uniquenames, 
                                   myname2idxMatrix= name2idxMatrix,
                                   directed=F)

    colnames(adjmatrix) <- uniquenames
    rownames(adjmatrix) <- uniquenames

    #"single" ensures that if a node is connected to at least one memeber in a candidate cluster
    # then this node will have a disctance 0 to the cluster, i.e, the concept of connect component
    #
    h1row <- hclust(as.dist(1-adjmatrix),method="single")

    #collect_garbage()
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
    #plot(h1row, xlab="",ylab="",main="",sub="")

    # ----------------- output Hierarchical Clustering image ------------------------------
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")

    #*-------------------------------------------------------------------------------------
    #* detect and label modules in TOM based on the hierarchical clustering dendrogram
    #*              
    myheightcutoff    = 0.99
    colcode.reduced   = moduleDetectLabel(hiercluster=h1row, myheightcutoff, 
                                          minsize1=minModuleSize, 
                                          useNumberAsLabel=T, startlabel=cclabelstart)
    table(colcode.reduced)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)

    cc= cbind(uniquenames, colcode.reduced)
    #colnames(cc) <- c("node", "ConnectedComponent")

    return(cc)
}


# here we return index
#
find_ConnectedComponents_adjmatrix = function(adjmatrix, minModuleSize=10, cclabelstart=1){

    dims = dim(adjmatrix)[1]

    # single item
    if (is.null(dims)){ # single item
          ipad = patchZeros(1)
          return ( rbind(c(1,ipad) ) )
    }

    # assign unique IDs
    no.uniquenames  = dims
    uniquenames     = c(1:no.uniquenames)

    colnames(adjmatrix) <- as.character(uniquenames)
    rownames(adjmatrix) <- as.character(uniquenames)

    #"single" ensures that if a node is connected to at least one memeber in a candidate cluster
    # then this node will have a disctance 0 to the cluster, i.e, the concept of connect component
    #
    h1row <- hclust(as.dist(1-adjmatrix),method="single")

    #collect_garbage()
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
    #plot(h1row, xlab="",ylab="",main="",sub="")

    # ----------------- output Hierarchical Clustering image ------------------------------
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")

    #*-------------------------------------------------------------------------------------
    #* detect and label modules in TOM based on the hierarchical clustering dendrogram
    #*              
    myheightcutoff    = 0.99
    colcode.reduced   = moduleDetectLabel(hiercluster=h1row, myheightcutoff, 
                                          minsize1=minModuleSize, 
                                          useNumberAsLabel=T, startlabel=cclabelstart)
    table(colcode.reduced)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)

    cc= cbind(uniquenames, colcode.reduced)
    #colnames(cc) <- c("node", "ConnectedComponent")

    return(cc)
}










##################################### K-CORES k-cores #################################################
#
# Recursively remove all vertices of degree less than K, UNTIL all vertices in the remaining# graph have at least  degree K
#
# return a matrix with the 1st column as node names
#                      the 2nd column as k-shell
#                      the rest columns as connected components in each core
#                               missing values (NA) means that these nodes are not in the cores
#                               -1 means that the nodes are isolated in the cores
#
# if fkey is not null, then it saves the node-kcore matrix and the link-kcore matrix into files
#
#
find_kcores = function (linkpairs, min_core=1, minCCsize=1, returnNodeCores=F, fkey=NULL) {
    kcoresCC      = NULL
    kshell        = NULL

    newlinkpairs  = linkpairs
    degrees       = degree_ByLinkPairs (linkpairs, directed=F)

    i = min_core
    while (T) {

       # i-core, recursively remove nodes of degree less than i
       #  until all nodes in the remaining graph have at least degree i
       #
       while ( T) {
           # get node names
           currnodenames = rownames(degrees)

           icoreSel  = degrees$totallinks >= i
           icoreSize = sum(icoreSel)

           # i-shell (actually i-1 shell)
           ishellSel  = !icoreSel
           ishellSize = sum(ishellSel, na.rm=T)
           if(ishellSize >0){
              ishell     = cbind(currnodenames[ishellSel],  rep(i-1,ishellSize) )
              kshell     = rbind(kshell, ishell)
           }

           # no more nodes with i-core
           if ( icoreSize <=1 ){
             break
           }

           # update adj matrix & network properties
           #
           newlinkpairs = getSubnetwork_LinkPairs(linkpairs=newlinkpairs, 
                                                  subnetNodes=currnodenames[icoreSel])
           if (is.null(newlinkpairs)) {
               icoreSize = 0 #no i-core found
               break
           }

           degrees      = degree_ByLinkPairs (newlinkpairs, directed=F)

           if (min(degrees$totallinks) >= i){
              break
           }
       }

       # no more nodes with i-core
       if ( icoreSize ==0){
         break
       }

       mystr=paste("Shell=", as.character(i), " core=", as.character(min(degrees$totallinks)), "\n",sep="")
       print(mystr)
       #if(i>=8){
       #   print(newlinkpairs )
       #   print(degrees)
       #}

       # finally, we got i-core nodes
       #
       # update node names in i-core
       currnodenames = rownames(degrees)

       icoreSize = length(currnodenames)
       icore     = cbind(currnodenames,  rep(i, icoreSize) )

       # get connected components for the current i-core
       #
       icorelinkpairs = getSubnetwork_LinkPairs(newlinkpairs, currnodenames)

       # nodes in the current icore have no links
       # so we assign a "grey" CC label
       #
       if( dim(icorelinkpairs)[1]==0){
         icclabel = rep("-1", icoreSize)
         icoreCC  = cbind(icore, icclabel)
         colnames(icoreCC) <- c("node", paste("core",as.character(i),sep=""),
                                        paste("cc_core",as.character(i),sep="") )

       }else{ # dectect cc, return is the (nodename, cclabel)
         icclabel = find_ConnectedComponents_linkpairs(linkpairs=icorelinkpairs, 
                                          minModuleSize=minCCsize,
                                          cclabelstart=1)
         icoreCC = merge(icore, icclabel, by.x=1, by.y=1, all.x=T)
         icoreCC = as.matrix(icoreCC)
         icoreCC = ifelse(is.na(icoreCC), "-1", icoreCC ) # -1: not in any CC
         icoreCC = ifelse(icoreCC=="grey", -1,  icoreCC )

         colnames(icoreCC) <- c("node", paste("core",as.character(i),sep=""),
                                        paste("cc_core",as.character(i),sep="") )
       }

       # so, CC's in i-core are put as a column, the missing vALUE indicates 
       #  these nodes are not in the current i-core
       #
       if (is.null(kcoresCC)){
          kcoresCC = icoreCC[,c(1,3)]
       } else{
          kcoresCC = merge(kcoresCC, icoreCC[,c(1,3)], by.x=1, by.y=1, all=T)
       }

       #print(paste(as.character(i), " coreness"))
       #print(icoreCC)

       # update adj matrix & network properties
       #
       #newlinkpairs = getSubnetwork_LinkPairs(linkpairs=newlinkpairs, subnetNodes=currnodenames[icoreSel])
       #degrees      = degree_ByLinkPairs (newlinkpairs, directed=F)

       i=i+1
    }

    colnames(kshell) <- c("node", "k_shell")
    kshell_coresCC = merge(kshell, kcoresCC, by.x=1, by.y=1, all=T)

    # coreness of each node
    # node	k_shell	cc_core1	cc_core2	cc_core3
    # 1a	1	1		NA		NA
    #
    if(!is.null(fkey)){
       fnodecores = paste(fkey, ".xls", sep="")
       write.table(kshell_coresCC, fnodecores, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
    }

    # core based link pairs
    #  a link belongs to a CC in a core iff both nodes are in the CC of the core
    #
    ncols      = dim(kshell_coresCC )[2]
    totalcores = ncols -2
    fnodecores = paste(fkey, "-linkpairs.xls", sep="")
    mergeleft  = merge(linkpairs, kshell_coresCC, by.x=2, by.y=1)
    mergeright = merge(mergeleft, kshell_coresCC, by.x=2, by.y=1)
    mergeright = as.matrix(mergeright )

    # src   dst  k_shell  cc_core1 cc_core2  cc_core3 k_shell cc_core1	cc_core2 cc_core3
    # 1     2    3                                     
    # link's shell is the min of the shells of the two nodes
    #
    newCCinCores = NULL
    shellmatrix = cbind( as.integer(mergeright[,3]), as.integer(mergeright[,totalcores+4]) )
    linkshell   = apply(shellmatrix, 1, min)
    newCCinCores = cbind(mergeright[,c(1:2)], linkshell)

    # look at whether two nodes of each link are in the same CC of the same core
    #
    for (i in c(1:totalcores) ){
        icorelinkflag = mergeright[,i+3]==mergeright[,totalcores+4+i]
        newCCinCores  = cbind(newCCinCores, mergeright[icorelinkflag, i+3] )
    }

    colnames(newCCinCores) <- c("src","dst", "kshell", colnames(kshell_coresCC)[3:(ncols)] )

    if(!is.null(fkey)){
        flinkcores = paste(fkey, "-cc-pairs.xls", sep="")
        write.table(newCCinCores, flinkcores, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
     }

    if(returnNodeCores){
        return(kshell_coresCC)
    }else{
        return (newCCinCores)
    }
}


#
#**********************************  END of K-cores  **********************************





##################################### K-clique #################################################
#
# Recursively remove all vertices of degree less than K, UNTIL all vertices in the remaining
# graph have at least  degree K
#
# return a matrix with the 1st column as node names
#                      the 2nd column as k-shell
#                      the rest columns as connected components in each core
#                               missing values (NA) means that these nodes are not in the cores
#                               -1 means that the nodes are isolated in the cores
#
# if fkey is not null, then it saves the node-kclique matrix and the link-kclique matrix into files
#
#

# check if setA is a subset of setB
setsub = function(setA,setB){
    if (length(setA) > length(setB)){
        return (F)
    }
    setAB = union(setA,setB)
    equal = setequal(setAB,setB)
    return (equal)
}

# check if a set is a subset of any set in a set list
#
setInSets = function (setC, setlist){
   nosets = length(setlist)
   for (i in c(1:nosets) ){
       isSub = setsub(setC, setlist[[i]])
       if (isSub){ #setC is a subset of setlist[[i]]
          return (T)
       }
   }
   return (F)
}

# globalCliquesMatrix is a mtrix with cols as gene index
#   and rows as cliques (sets)
# So, we first get a matrix with only the cols given by setC 
#  then, compute the sum of each row
# If one sum value is equal to the number of elements in setC,
#  then, setC is a subset of setInSets_matrix, otherwise, not
#
setInSets_matrix = function (setC) {
   rsum  = apply(globalCliquesMatrix[,setC], 1, sum)
   found = length(setC)==rsum
   asum  = sum(found)
   return (asum>=1) 
}


# ------ global variables ------------

#~~~~~~~ input ~~~~~~~~~~~~~~~~~~~~
# globalAdjLists
# globalLinkpairsIndex
# globalFullLinked
# globalAdjSize

#~~~~~~~ output ~~~~~~~~~~~~~~~~~~~
# globalCliques
# globalCliqueSize

kcore_subnets = function(linksindex, nodesindex=NULL, kcore=3, name2integer=T) {

    if(!is.null(nodesindex)) {
        mright  = merge(linksindex, nodesindex, by.x=2, by.y=1, all=F)
        if( dim(mright)[1]==0){
            return (list(NULL, NULL, F) )
        }

        orgLinks= merge(mright,     nodesindex, by.x=2, by.y=1, all=F)
        if( dim(orgLinks)[1]==0){
            return (list(NULL, NULL, F) )
        }
    }else{
        orgLinks= linksindex
    }

 
    while(T) {
      no.orglinks   = dim(orgLinks)[1]
      degrees       = degree_ByLinkPairs (orgLinks, directed=F)
      dim(degrees)
      selected    = degrees$totallinks >= kcore
      no.selnodes = sum(selected)

      if (no.selnodes==0) {
          return (list(NULL, NULL, F) )
      }

      selectedNodes    = rownames(degrees)[selected]
      selected_nolinks = degrees$totallinks[selected]/no.selnodes

      # =========== important ======================
      if(name2integer){
          selectedNodes = sort( as.integer(selectedNodes) )
      }

      if (no.selnodes < kcore+1){ # no enough nodes, so no enough links for each node
          return (list(NULL, NULL, F))
      }

      selectedLinks = getSubnetwork_LinkPairs(orgLinks, subnetNodes=selectedNodes)
      if (is.null(selectedLinks)) {
          return (list(NULL, NULL, F))
      }

      selectedLinks = as.matrix(selectedLinks)
      no.newlinks   = dim(selectedLinks)[1]

      if (no.newlinks==0) {
          return (list(NULL, NULL, F))
      }

      orgLinks = selectedLinks

      if(no.newlinks==no.orglinks){
          break
      }

    }

    # check if these nodes are fully linked
    #
    fulllinks  = no.selnodes*(no.selnodes -1)/2
    fulllinked = no.newlinks == fulllinks

    rets = list(selectedNodes, selectedLinks, fulllinked)

    return (rets)
}

#v=775
#ret=kcore_subnets(selectedLinksIndex, c(v, globalAdjLists[[v]]), kcore=5)

subnet_is_clique = function(nodes, links){
   maxlinks = nodes*(nodes-1)/2
   ratio    = links/maxlinks
   isclique = ratio==1
   return (isclique)
}

####################################################################################
################ find all cores and return the core nodes ##########################
#
#  nodeMembershipOnly==T
#       return the T/F corenodesMatrix
#  nodeMembershipOnly==F  
#       return( list(allcores, coresnodes, coreslinks, oressizes, iscliques) )
#
#  returnLinkID=T: this means that the third column of linksindex is IDs of links
#
KCores_ALL = function(linksindex, nodesindex=NULL, minkcore=3, name2integer=T, nodeMembershipOnly=T, minkcoreOnly=F, returnLinkID=F) {

    if(!is.null(nodesindex)) {
        mright  = merge(linksindex, nodesindex, by.x=2, by.y=1, all=F)
        if( dim(mright)[1]==0){
            #return (list(NULL, NULL, F) )
            return (list(NULL, NULL, NULL, NULL, NULL) )
        }

        orgLinks= merge(mright,     nodesindex, by.x=2, by.y=1, all=F)
        if( dim(orgLinks)[1]==0){
            #return (list(NULL, NULL, F) )
            return (list(NULL, NULL, NULL, NULL, NULL) )
        }
    }else{
        orgLinks= linksindex
    }

    degrees      = degree_ByLinkPairs (orgLinks, directed=F)

    kcore        = minkcore
    reachMaxcore = F
    maxindex     = max(degrees$totallinks) #max(linksindex)

    #nodeIdxSequence = c(1:maxindex)
    corenodesMatrix = NULL
    allcores        = NULL

    if(!nodeMembershipOnly) {
        kcoreMemberNodes = as.list(rep(0,maxindex))
        kcoreMemberLinks = as.list(rep(0,maxindex))
    }
    kcoreFullLinked  = rep(0,maxindex)
    kcoresSizes      = rep(0,maxindex)

    corecnt= 0
    while(T) { #all cores

        while(T) { #kcore
          no.orglinks   = dim(orgLinks)[1]
          degrees       = degree_ByLinkPairs (orgLinks, directed=F)
          dim(degrees)
          selected    = degrees$totallinks >= kcore
          no.selnodes = sum(selected)

          if (no.selnodes==0) {
              reachMaxcore=T
              break
          }

          selectedNodes    = rownames(degrees)[selected]
          selected_nolinks = degrees$totallinks[selected]/no.selnodes

          # =========== important ======================
          if(name2integer){
              selectedNodes = sort( as.integer(selectedNodes) )
          }

          if (no.selnodes < kcore+1){ # no enough nodes, so no enough links for each node
              reachMaxcore=T
              break
          }

          selectedLinks = getSubnetwork_LinkPairs(orgLinks, subnetNodes=selectedNodes)
          if (is.null(selectedLinks)) {
              reachMaxcore=T
              break
          }

          selectedLinks = as.matrix(selectedLinks)
          no.newlinks   = dim(selectedLinks)[1]

          if (no.newlinks==0) {
              reachMaxcore=T
              break
          }

          orgLinks = selectedLinks

          # record the nodes in the core
          if(no.newlinks==no.orglinks){
              corecnt   = corecnt + 1
              allcores  = c(allcores, kcore)

              if(nodeMembershipOnly) {
                 coreMembership      = rep(F, maxindex)
                 coreMembership[selectedNodes ] = T
                 corenodesMatrix = cbind(corenodesMatrix, coreMembership)
              }else{
                 kcoreMemberNodes[[corecnt]] = selectedNodes
                 if(returnLinkID) {
                     kcoreMemberLinks[[corecnt]] = as.integer(selectedLinks[,3])
                 }else{
                     kcoreMemberLinks[[corecnt]] = selectedLinks
                 }
              }      

              # clique or not
              kcoreFullLinked[corecnt] = subnet_is_clique(length(selectedNodes), dim(selectedLinks)[1])
              kcoresSizes[corecnt]     = length(selectedNodes)

              break
          }

        } #kcore


        if (reachMaxcore){
          break
        }

        # run only on the minkcore
        #
        if( minkcoreOnly){
           break
        }

        kcore = kcore + 1

    }

    if(nodeMembershipOnly) {
       colnames(corenodesMatrix) <- allcores
       return (corenodesMatrix)   
    } else{
       # NULL is OK, list(NULL) ==> [[1]] NULL
       #
       if(!is.null(allcores) ){
           length(kcoreMemberNodes) <- corecnt
           length(kcoreMemberLinks) <- corecnt
           length(kcoresSizes)      <- corecnt
           length(kcoreFullLinked)  <- corecnt
           finalList = list(corenesses= allcores, 
                        coresnodes= kcoreMemberNodes, 
                        coreslinks= kcoreMemberLinks,
                        coressizes= kcoresSizes, 
                        iscliques = kcoreFullLinked)
       } else{
           finalList = list(corenesses= NULL,
                        coresnodes= NULL, 
                        coreslinks= NULL,
                        coressizes= NULL, 
                        iscliques = NULL)
       }

       return (finalList)
    }
}

#####################################################################################
#  kcores with individual nodes's neighbors being cored
#
#  inputs are global variables:
#      selectedLinksIndex
#      globalAdjLists
#      no.nodes
#
# Notice that if we want k-core, then for the neighborhood of a node we need only
#       (k-1)-core. Keep in mind that
#        K-1 is recorded in global_nodeneighbor_cores_sizes
#
# *Notice that to assignm value to individual elements in a global variable, 
#   you have to use <<- instead of "=", otherwise, there is no change, check the
#   following example
#
if(F){
 # global variable
 a=c(1:10)

 gvar=function(){
     b<<-as.list(rep(10, 5))
 }
 cf = function(x){
     b[[3]] <<- x
 } #good

 cf2 = function(x){a[[6]]=x}
 #cf(10), cf2(11), cf(20)
}

KCores_ALL_in_neighbors = function(minkcore=3, returnLinkID=F) {

    # global variables
    #
    # initialization for each node
    # 1) coreness node/link members
    # 2) no of nodes in each coreness
    # 3) min/max coreness
    #
    global_nodenghr_cores       <<- as.list( rep(NA, no.nodes) )
    global_nodenghr_minmaxcore  <<- matrix(0, no.nodes, 2)

    # estimated max coreness
    #    
    avglinks      = no.links/no.nodes
    maxcore_rough = as.integer(4*avglinks)
    
    # for each core of each node
    #
    #tmplist = as.list( rep(NA, maxcore_rough) )
    #global_nodenghr_cores  = ifelse(is.na(global_nodeneighbor_cores), tmplist)
    #for (i in c(1:no.nodes) ){
    #    global_nodenghr_cores[[i]]      = tmplist
    #}
    #source("C:/Zhang/BinBlast/R-tomfunctions.R")
    #start.time <- proc.time()[3]

    print(paste("find all cores for the neighbors of each of ", no.nodes, " nodes", sep="") )
    cnt=0

    # count how many links in neighbor's cores
    #
    global_nolinks_nghbr_cores <<- rep(0, maxcore_rough)
    
    for (v in c(1:no.nodes)){
       print(as.character(v))
       
       if(returnLinkID) {
         vCores = KCores_ALL(linksindex=linkpairsIndexed, 
                   nodesindex=globalAdjLists[[v]], minkcore=minkcore-1, name2integer=T, 
                   nodeMembershipOnly=F, returnLinkID=returnLinkID)
       }else{
         vCores = KCores_ALL(linksindex=selectedLinksIndex, 
                   nodesindex=globalAdjLists[[v]], minkcore=minkcore-1, name2integer=T, 
                   nodeMembershipOnly=F, returnLinkID=returnLinkID)
       }


       if(!is.null(vCores$corenesses) ){
          #print (v)
          cnt= cnt + 1
          global_nodenghr_cores[[v]]     <<- vCores
          global_nodenghr_minmaxcore[v,] <<- c(min(vCores$corenesses), max(vCores$corenesses) )
          
          # for each core, we count how many links there
          #
          vcoresizes = vCores$coressizes
          length(vcoresizes) = maxcore_rough
          vcoresizes = ifelse(is.na(vcoresizes), 0, vcoresizes)
          global_nolinks_nghbr_cores <<- global_nolinks_nghbr_cores + vcoresizes
       }
    
    }

    # count of links is duplicated
    global_nolinks_nghbr_cores <<- global_nolinks_nghbr_cores

    #preprocessing_Time <- proc.time()[3]-start.time
    #preprocessing_Time
    #cnt

    collect_garbage()

    return (cnt)
}

# remove one node in the neighbors
#
#
# global_nodenghr_cores[[]]$corenesses
# global_nodenghr_cores[[]]$coresnodes
# global_nodenghr_cores[[]]$coreslinks
# global_nodenghr_cores[[]]$coressizes
# global_nodenghr_cores[[]]$iscliques
#
# global_nolinks_nghbr_cores
# global_nodenghr_minmaxcore
#

update_nodenghr_cores=function(nodenghr_cores, removegene=NULL, coreidx, use_link_id=F)
{

           curcoreness = nodenghr_cores$corenesses[coreidx]

           #print(curcoreness)
           #print(nodenghr_cores$coresnodes[[coreidx]])

           finalList   = list(corenesses= NULL,
                        coresnodes= NULL, 
                        coreslinks= NULL,
                        coressizes= NULL, 
                        iscliques = NULL)

           if(curcoreness==0){
               return (NULL)
           }

           # the node to be removed is not a neighbor of the current node
           #
           isneighbor= is.element(removegene, nodenghr_cores$coresnodes[[coreidx]])
           if(!isneighbor){
               return (NULL)
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           if(use_link_id){
             # get actual links and link IDs from the node's link IDs
             #
             ilinks     = linkpairsIndexed[ nodenghr_cores$coreslinks[[coreidx]], ]
             selLeft    = ilinks[,1] != removegene
             selRight   = ilinks[,2] != removegene             
           } else{
             selLeft    = nodenghr_cores$coreslinks[[coreidx]][,1] != removegene
             selRight   = nodenghr_cores$coreslinks[[coreidx]][,2] != removegene
           }

           sel        = selLeft & selRight
           no.remains = sum(sel)
           
           if(no.remains==0){
               return (finalList)
           }

          # redo the core analysis on the updated neighbors but only look at the same coreness
          #
          if(use_link_id){
            vCores = KCores_ALL(linksindex=rbind(ilinks[sel, ]), 
                   nodesindex=NULL, minkcore=curcoreness, name2integer=T, 
                   nodeMembershipOnly=F, minkcoreOnly=T, returnLinkID=use_link_id)
          } else{
            vCores = KCores_ALL(linksindex=rbind(nodenghr_cores$coreslinks[[coreidx]][sel,]), 
                   nodesindex=NULL, minkcore=curcoreness, name2integer=T, 
                   nodeMembershipOnly=F, minkcoreOnly=T, returnLinkID=use_link_id)
          }
          return (vCores)
}


find_maxcoreness=function(orgLinks, returnCoreSize=F){

    if(is.null(orgLinks)) {
         return (0)
    }

    newlinks= orgLinks
    degreesOrg = degree_ByLinkPairs (newlinks, directed=F)

    # find max core
    icore = max( degreesOrg$totallinks)
    no.orglinks = dim(newlinks)[1]
    while(T) {
          selected    = degreesOrg$totallinks >= icore
          nosels      = sum(selected)

          # icore is too big, then there is not enough members
          #  an icore requres at least (icore + 1) nodes
          #
          if (nosels < icore+1) {
               icore = icore - 1
               next
          }else{
               break
          }
    }

    if (icore ==1 ){
       return (icore)
    }

    while (T) { # choose k-core

        # from the very start
        newlinks = orgLinks
        found    = F
        firstrun = T

        while(T) {
          if (firstrun) {
              degrees = degreesOrg # no need to compute links again
              firstrun = F
          } else{
              degrees = degree_ByLinkPairs (newlinks, directed=F)
          }

          no.orglinks = dim(newlinks)[1]
          selected    = degrees$totallinks >= icore
          nosels      = sum(selected)

          # if no nodes was selected or only one node was selected, stop
          if (nosels <=1 ){
             break
          }

          selectedNodes = rownames(degrees)[selected]
          selectedLinks = getSubnetwork_LinkPairs(orgLinks[,c(1,2)],subnetNodes=selectedNodes)
          if (is.null(selectedLinks)) {
              break
          }

          selectedLinks = as.matrix(selectedLinks)
          no.newlinks   = dim(selectedLinks)[1]

          newlinks = selectedLinks

          if(no.newlinks==no.orglinks){
             found=T
             break
          }

        } #while
        #dim(selectedLinks)

        if(found){
           break
        }
        #if (icore <=2){
        #    return (icore)
        #}

        icore = icore - 1

   } #while

   if ( !returnCoreSize) {
      return (icore)    
   } else{
      return ( c(icore,length(selectedNodes), no.newlinks) )
   }
}


#********************************************************************************
#
# * notice that union(list(c(1,2),c(1,3)),list(c(1,2)) ) is different from
#               union(list(c(1,2),c(1,3)),list(c(1:2)) ) 
# the first yields right answer
# also          union(list(c(1,2),c(1,3)),list(c(1:2)) ) gives right answer
# so, internally in R, c(1,2) is different from c(1:2)
#
#--------------------------------------------------------------------------------
#
# initially, setA is empty and setB includes all neighbors of curnode
#  Kclique is the number of nodes needed from setB in order to form a k-clique
#   which includes curnode
#
# global_nodenghr_cores[[]]$corenesses
# global_nodenghr_cores[[]]$coresnodes
# global_nodenghr_cores[[]]$coreslinks
# global_nodenghr_cores[[]]$coressizes
# global_nodenghr_cores[[]]$iscliques
#

recursive_Cliques_nodenghr = function(curnode, setA, setB, neighborlinks, Kclique=3, coreidx, use_link_id=F){

       #print (curnode)
       # at the beginning, global_nodenghr_cores[[v]]$corenesses[coreidx] is not NULL, but it may be NULL after
       # several runs when remove the nodes which have been visited and have core(s)
       #
       # NULL means that the current node don't have enough neighbors (>=Kclique)
       #  so it cannot be in a clique of size = Kclique
       #

       if( is.null(global_nodenghr_cores[[curnode]]$corenesses[coreidx]) ){
           return (0)
       }

       #print(paste("curnode ", curnode, " coresness:"))
       #print(global_nodenghr_cores[[curnode]]$corenesses)

       # this node doesn't have a neighborhood including k-core with k < coreidx-th
       #
       if( length(global_nodenghr_cores[[curnode]]$corenesses)<coreidx){
           return (0)
       }

       if( global_nodenghr_cores[[curnode]]$corenesses[coreidx] ==0 ){
           return (0)
       }
       
       # move up v's neighbors to setA
       #
       setAA      = union(setA, curnode)
       size.setAA = length(setAA)

       #print("Set A, , currnode, B: ")
       #print(setA)
       #print(curnode)
       #print(setB)

       if (is.null(setB) ){
          setBB = global_nodenghr_cores[[curnode]]$coresnodes[[coreidx]]
       }else{
          setBB = setdiff(setB, curnode) # remove v from setB
          #update nodes
          setBB = intersect(setB, global_nodenghr_cores[[curnode]]$coresnodes[[coreidx]] )
       }

       size.setBB = length(setBB)

       # check if it is already a clique
       #
       if (size.setBB <= 1 ){#setBB is null or has only one node, setAA+setBB already form a clique
          isclique = TRUE
       }else{
          # update links
          newneighborlinks = getSubnetwork_LinkPairs(neighborlinks, setBB)
          if(is.null(newneighborlinks)){
              isclique = -1
          } else{
              isclique = subnet_is_clique(nodes=length(setBB), links=dim(newneighborlinks)[1])
          }
       }

       if(isclique==T){
           #print("Early: ")
           #print(sort(union(setAA, setBB) ))
           if( size.setBB == Kclique-size.setAA){
               currclique  = sort(union(setAA, setBB) )
               currcliqueL = list(currclique)

               #set global variable
               #
               #found   = setInSets(currclique, globalCliques)
               found   = setInSets_matrix(currclique)

               #print("isclique==T")
               #print(currclique)
               if(!found){
                  globalCliques   <<- c(globalCliques,  currcliqueL)

                  xcliq = globalNodeZeroIndex
                  xcliq[currclique] = 1
                  globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                  return (1)
               }
               return (1)
           }
           return (0)
       }

       # the remaining nodes in setBB don't have any link with each other
       # so we stop the recursive process
       #
       if(isclique==-1){
           if( size.setAA == Kclique){
               currclique = sort(setAA)
               currcliqueL = list(currclique)
               #set global variable
               #
               #found   = setInSets(currclique, globalCliques)
               found   = setInSets_matrix(currclique)

               #print("isclique==-1")
               #print(currclique)

               if(!found){
                  globalCliques   <<- c(globalCliques,  currcliqueL)

                  xcliq = globalNodeZeroIndex
                  xcliq[currclique] = 1
                  globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                  return (1)
               }
               return (1)
           }else if(size.setAA == Kclique-1) {

               for(ek in setBB) {
                   currclique = sort(c(setAA, ek) )
                   currcliqueL = list(currclique)
                   #set global variable
                   #
                   #found   = setInSets(currclique, globalCliques)
                   found   = setInSets_matrix(currclique)

                   if(!found){
                      globalCliques   <<- c(globalCliques,  currcliqueL)

                      xcliq = globalNodeZeroIndex
                      xcliq[currclique] = 1
                      globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)
                   }
               }
               
               return (1)
           }

           return (0)
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )       
       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }
       
       # do they have enough neighbors
       setAB = sort(union(setAA,setBB))
       if (length(setAB)<Kclique ) {
          return (0)
       }

       # check if seA+setB is already identified as a cliuqe
       if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             #found   = setInSets(setAB, globalCliques)
             found   = setInSets_matrix(setAB)

             if (found){
                #print("Already found (first)")
                return(1)
             }
       }
      
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       nlinks = dim(newneighborlinks)[1]
       nnodes = length(setBB)
       ncc    = 2*nlinks/nnodes/(nnodes-1)

       if(F){
       #if(ncc <=tuo_coreCC  & nnodes >= tuo_coreSZ) {
       ret=kcore_subnets(linksindex=newneighborlinks, nodesindex=setBB,
                               kcore=bkclique-1 )
       setBB               = ret[[1]]
       newneighborlinks    = ret[[2]]
       neighborsfulllinked = ret[[3]]

       # if all remaining neighbors are fully linked, and they form a bkclique-clique,
       #   we simply return the union of AA+BB
       #  or we reject this clique because of the bigger/smaller size
       #
       if(neighborsfulllinked) {
                #print( paste("SetBBcore=", concatenate(as.character(sort(setBB)),","), sep="") )

                if ( length(setBB)==bkclique) {
                    sortedAB   = sort(union(setAA, setBB) )
                    currclique = list(sortedAB)
                    #print( paste("CLIQUEfull=", concatenate(as.character(currclique),","), sep="") )

                    #found    = setInSets(sortedAB, globalCliques)
                    found    = setInSets(sortedAB)

                    if(!found) {
                        #currCliques          = union(cliquesFound,  currclique)
                        #options(globalCliques= currCliques)
                        globalCliques       <<- c(globalCliques, currclique)

                        xcliq = globalNodeZeroIndex
                        xcliq[sortedAB] = 1
                        globalCliquesMatrix <<- rbind(globalCliquesMatrix, xcliq)

                        return (1)
                    }
                }
                #currclique = list(sort(union(setAA, setBB) ) )
                #print( paste("CLIQUEfull -- =", concatenate(as.character(currclique),","), sep="") )
                return (1)
       }     

       # the neighbors don't form a (bkclique-2)-core, so curnode and its neighbors cannot
       #  form a Kclique-clique
       if(is.null(setBB)) {
          #print( "no cores exisit to form a cilque")
          return (0)
       }

      }



       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          #endIdx = bbsize -1

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques_nodenghr(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    neighborlinks=newneighborlinks, 
                                    Kclique      =Kclique, coreidx=coreidx)
             total_cliques = total_cliques + count
          }


       return (total_cliques)
}

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


recursive_Cliques_Super=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)
       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA      = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )
       
       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
          setAB = sort(union(c(setAA,setBB), NULL))
          #setAB = union(setAA,setBB) 

          if (length(setAB)<Kclique ) {
            return (0)
          }

          # check if seA+setB is already identified as a cliuqe
          if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             found   = setInSets(setAB, globalCliques)
             if (found){
                #print("Already found (first)")
                return(0)
             }
          }

       if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          if( length(currclique) == Kclique ){
             currcliqueL          = list(currclique)

             #set global variable
             #
             #globalCliques       <<- union(globalCliques,  currcliqueL)
             globalCliques       <<- c(globalCliques,  currcliqueL)

             #options(globalCliques= currCliques)

             return (1)
          }
          return (0)

       }else if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             currclique           = sort(setAA)
             currcliqueL          = list(currclique)

             #set global variable
             #
             #currCliques          = union(globalCliques,  currcliqueL)
             #options(globalCliques= currCliques)
             #globalCliques       <<- union(globalCliques,  currcliqueL)
             globalCliques       <<- c(globalCliques,  currcliqueL)

             return (1)
          }
          return (0)
       }
      
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

       # special process
       #if (bkclique==2 & Kclique<=4){
       if (bkclique==2 & Kclique<= 5){
           blinks = getSubnetwork_LinkPairs(linkpairs=neighborlinks, subnetNodes=setBB)
           if( is.null(blinks) ){
               return (0)
           }

           no.blinks = dim(blinks)[1]
           qcliques  = cbind(blinks, rep(curnode,no.blinks))
           qcliquesO = apply(qcliques, 1, sort) #notice that the array is transposed after sort

           #print(paste("row,col=", dim(qcliques)[1],  dim(qcliques)[2],sep=" " ))
           #print(paste("row,col (O)=", dim(qcliquesO)[1],  dim(qcliquesO)[2],sep=" " ))

           #cliquesFound = NULL
           for (q in no.blinks){
             #cliquesFound = union(cliquesFound, list(qcliquesO[, q]) )
             found         = setInSets(list(qcliquesO[,q]), globalCliques)
             if(!found) {
                #print( paste("CLIQUEnormal2=", concatenate(as.character(currclique),","), sep="") )
                #cliquesFound = c(cliquesFound, list(qcliquesO[,q]) )
                
                globalCliques <<- c(globalCliques,  list(qcliquesO[,q]))
                
             }
           }
           #
           #globalCliques       <<- union(globalCliques,  cliquesFound)
           #options(globalCliques= cliquesFound)

       }else{
          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques_Super(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       }

       return (total_cliques)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

recursive_Cliques=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")
       #adjlists     = getOption("")

       #sytr= paste("A=",  concatenate(as.character(setA),","), 
       #            " B=", concatenate(as.character(setB),",") )
       #print(sytr)

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)

       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

       # bigger clique, not interest
       #
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
       if (!is.null(setA) | !is.null(setB)) {
          setAB = sort(union(c(setAA,setBB), NULL))
          if (length(setAB)<Kclique ) {
            return (0)
          }
          # check if seA+setB is already identified as a cliuqe
          if (!is.null(cliquesFound)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             found   = setInSets(setAB, globalCliques)

             if (found){
                #print("Already found (first)")
                return(0)
             }
          }
       }
       
       if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             sortedAA = sort(setAA)
             found    = setInSets(sortedAA, globalCliques)
             if(!found) {
               currclique = list(sortedAA)
               #print( paste("CLIQUEnormal1=", concatenate(as.character(currclique),","), sep="") )
               #currCliques          = union(cliquesFound,  currclique)
               #options(globalCliques= currCliques)

               globalCliques       <<- union(globalCliques,  currclique )

             }
             return (1)

          }
          return (0)

       } else if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          currcliqueL = list(currclique)

          if( length(currclique) == Kclique ){
             found    = setInSets(currclique, globalCliques)
             if(!found) {
                #print( paste("CLIQUEnormal2=", concatenate(as.character(currclique),","), sep="") )
                #currCliques          = union(cliquesFound,  currcliqueL)
                #options(globalCliques= currCliques)

                globalCliques       <<- union(globalCliques, currcliqueL)

                return (1)
             }
             return (1)
          }
          return (0)

       }

       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA
       ret=kcore_subnets(linksindex=neighborlinks, nodesindex=setBB,
                               kcore=bkclique-1 )

       setBB               = ret[[1]]
       neighborlinks2      = ret[[2]]
       neighborsfulllinked = ret[[3]]

       # if all remaining neighbors are fully linked, and they form a bkclique-clique,
       #   we simply return the union of AA+BB
       #  or we reject this clique because of the bigger/smaller size
       #
       if(neighborsfulllinked) {
                #print( paste("SetBBcore=", concatenate(as.character(sort(setBB)),","), sep="") )

                if ( length(setBB)==bkclique) {
                    sortedAB   = sort(union(setAA, setBB) )
                    currclique = list(sortedAB)
                    #print( paste("CLIQUEfull=", concatenate(as.character(currclique),","), sep="") )

                    found    = setInSets(sortedAB, cliquesFound)
                    if(!found) {
                        #currCliques          = union(cliquesFound,  currclique)
                        #options(globalCliques= currCliques)
                        globalCliques       <<- union(globalCliques, currclique)

                        return (1)
                    }
                }
                #currclique = list(sort(union(setAA, setBB) ) )
                #print( paste("CLIQUEfull -- =", concatenate(as.character(currclique),","), sep="") )
                return (1)
       }              

       # the neighbors don't form a (bkclique-2)-core, so curnode and its neighbors cannot
       #  form a Kclique-clique
       if(is.null(setBB)) {
          #print( "no cores exisit to form a cilque")
          return (0)
       }

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       currCliques   = NULL
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize       |
       #  **********************
       #             | bkclique |
       #
       #             ^
       #             |
       #            endIdx

       if(1==1) {
          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks2, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       } else{
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          # try all possible combinations
          for (each in setBB){
             count = recursive_Cliques(curnode =each, setA=setAA, setB=setBB, 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks2, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       }
       return (total_cliques)
}



###############################################################################################################
############################### Start of normal clique identification #########################################
#

recursive_CliquesNormal=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)
       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )

       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
          setAB = sort(union(c(setAA,setBB), NULL))
          if (length(setAB)<Kclique ) {
            return (0)
          }
          # check if seA+setB is already identified as a cliuqe
          if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )
             found   = setInSets(setAB, globalCliques)
             if (found){
                #print("Already found (first)")
                return(0)
             }
          }

       if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             currclique           = sort(setAA)
             currcliqueL          = list(currclique)
             found   = setInSets(currclique, globalCliques)
             if(!found) {
                 globalCliques        <<- union(globalCliques,  currcliqueL)
                 #options(globalCliques= currCliques)
             }
             return (1)
          }
          return (0)

       } else if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          if( length(currclique) == Kclique ){
             currcliqueL          = list(currclique)
             found   = setInSets(currclique, globalCliques)
             if(!found) {
                globalCliques        <<- union(globalCliques,  currcliqueL)
                #currCliques          = union(cliquesFound,  currcliqueL)
                #options(globalCliques= currCliques)
             }
             return (1)
          }
          return (0)
       }
       
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       currCliques   = NULL
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_CliquesNormal(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       return (total_cliques)
}


# used global variables
#      selectedNodesName
#      selectedLinksIndex
#      degrees
#      no.nodes
#      no.links


############################################################################################################
#
#                                  Super-fast clique identification
# used global variables
#      selectedNodesName
#      selectedLinksIndex
#      degrees
#      no.nodes
#      no.links
#
#  global_nodenghr_cores: each element corresponding to a node is a list of following items 
#     for the node's neighbors:
#           finalList = list(corenesses= allcores, 
#                        coresnodes= kcoreMemberNodes, 
#                        coreslinks= kcoreMemberLinks,
#                        coressizes= kcoresSizes, 
#                        iscliques = kcoreFullLinked)
#    
#  global_nodenghr_minmaxcore: no.nodes x 2 matrix with (min coreness, max coreness)
#
#
# IF MemoyForSpeed= T, then the core identification will return all nodes and links
#    in each core so that there is no need to perform merge operation to get link members
#    in each core
#

find_kcliques = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, MemoyForSpeed=T, cliqueOnly=T, use_link_id=F) {
    DEBUG=F

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",       sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram", sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    linkpairsIndexed <<- cbind(selectedLinksIndex, c(1:no.links) )
    colnames(linkpairsIndexed) <<-c("src","dst","id")

    nulcoresngbr = list(corenesses= NULL,coresnodes= NULL, 
                        coreslinks= NULL, coressizes= NULL, iscliques = NULL)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # find all cores in each node's neighbors
    #
    # k-core is defined as  number of links for each node is >= k
    # so a k+1 clique is contained in (k-1)-core
    #
    nodes_with_cores = KCores_ALL_in_neighbors(minkcore=min_clique-1, returnLinkID=use_link_id)

    allcores.readTime <- proc.time()[3]-now

    # as the coreness is based on node's neighbor, so we need add the node in
    #  
    maxcoreness  = max(global_nodenghr_minmaxcore[,2]) + 1

    coreindex    = c( (min_clique-1):maxcoreness)
    cliqueindex  = coreindex   + 1 
    max_clique   = maxcoreness + 1

    maxcliqueVect= rep(max_clique, no.nodes)
    no.cores     = length(coreindex)

    allSelnodesIndex = c(1:no.nodes)

    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (i in c(no.cores:1) ) {

        iclq= cliqueindex[i]
        iclq.now<-proc.time()[3]

        # here coreness (global_nodenghr_minmaxcore) is based on a node's neighbor, 
        # so including the node itself leads to  (global_nodenghr_minmaxcore + 1)-core
        #
        iclqNodesSel   = global_nodenghr_minmaxcore[,2] + 1 >= iclq-1

        iclqNodesIndex = allSelnodesIndex[iclqNodesSel]
        no.nodesINiclq = length(iclqNodesIndex)

        if(no.nodesINiclq < iclq){
           # recording time
           iclq.readTime <- proc.time()[3]-iclq.now
           iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, 
                            global_nolinks_nghbr_cores[i],iclq.readTime) )
           next
        }

        # iterate each node to identify all iclq-cliques including this node
        #
        jjindex = c(1:no.nodesINiclq)
        for ( j in jjindex ) {

           cliques_count = 0
           v = iclqNodesIndex[j]
           
           mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
           print(mystr)


           # global_nodenghr_cores[[]]$coresnesses
           # global_nodenghr_cores[[]]$coresnodes
           # global_nodenghr_cores[[]]$coreslinks
           # global_nodenghr_cores[[]]$coressizes
           # global_nodenghr_cores[[]]$iscliques
           #neighbors    = global_nodenghr_cores[[v]]$coresnodes[[i]]

          # first call, whether its neighbors form a clique already has been checked
          #
          # when the current node has Kclique-1 neighbors, there are two scenarios:
          #  1) they form a (Kclique-1)-clique
          #  2) they don't form a (Kclique-1)-clique
          #
          if(global_nodenghr_cores[[v]]$coressizes[i] == iclq-1) {
             # 1)
             if(global_nodenghr_cores[[v]]$iscliques[i]){
                 currclique  = sort(c(v, global_nodenghr_cores[[v]]$coresnodes[[i]]) )
                 currcliqueL = list(currclique)

                 #set global variable
                 #
                 #globalCliques  <<- union(globalCliques,  currcliqueL)
                 #found   = setInSets(currclique, globalCliques)
                 found   = setInSets_matrix(currclique)

                 if(!found){
                    globalCliques <<- c(globalCliques,  currcliqueL)

                    xcliq = globalNodeZeroIndex
                    xcliq[currclique] = 1
                    globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                 }
                 cliques_count = 1
             }
             cliques_count = 1

          } else if(global_nodenghr_cores[[v]]$coressizes[i] < iclq-1) {
             # if the currrent node does't have enough neighbors, it should be removed
             #  from the neighbors of each other node
             # Notice that global_nodenghr_cores[[v]]$coressizes is dynamically changed
             #
             cliques_count = 1
          } else{

             if(use_link_id) {
                    ijlinks = selectedLinksIndex[ global_nodenghr_cores[[v]]$coreslinks[[i]], ]
                    cliques_count = recursive_Cliques_nodenghr(curnode=v, setA=NULL, setB=NULL, 
                              neighborlinks= ijlinks,
                              Kclique=iclq, coreidx=i, use_link_id=use_link_id)#-1)
             }else{
                    cliques_count = recursive_Cliques_nodenghr(curnode=v, setA=NULL, setB=NULL, 
                              neighborlinks= global_nodenghr_cores[[v]]$coreslinks[[i]],
                              Kclique=iclq, coreidx=i, use_link_id=use_link_id)#-1)
             }

          }

           #cliques_count

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }


           # remove the current node and its links
           #
           # 1) reset the i-th core of the current node v
           #
           if(DEBUG){
           global_nodenghr_cores[[v]]$corenesses[i]    <- 0
           global_nodenghr_cores[[v]]$coresnodes[[i]]  <- NA
           global_nodenghr_cores[[v]]$coreslinks[[i]]  <- NA
           global_nodenghr_cores[[v]]$coressizes[i]    <- 0
           global_nodenghr_cores[[v]]$iscliques[i]     <- F
           }else{
           global_nodenghr_cores[[v]]$corenesses[i]    <<- 0
           global_nodenghr_cores[[v]]$coresnodes[[i]]  <<- NA
           global_nodenghr_cores[[v]]$coreslinks[[i]]  <<- NA
           global_nodenghr_cores[[v]]$coressizes[i]    <<- 0
           global_nodenghr_cores[[v]]$iscliques[i]     <<- F
           }

           # 2) remove the current node from LinkPairs and AdjList of
           #     the i-th core of each other node (don't include the current node v)
           #    so that the same cliques including v will not be found again
           #
           for ( m in jjindex[-j] ) {
                qq = iclqNodesIndex[m]
                isneighbor = is.element(v, global_nodenghr_cores[[qq]]$coresnodes[[i]])
                if(!isneighbor){
                    next
                }

                icorenew =update_nodenghr_cores(nodenghr_cores=global_nodenghr_cores[[qq]], 
                                                removegene=v, coreidx=i, use_link_id=use_link_id)

                # qq's neighbors don't form a core or v is not in qq's neighbors
                if(is.null(icorenew)){
                   next
                }


                # elements of List can be null, but not for vector
                #
                # IMPORTANT: update list entries with NA instead of NULL,
                #
                if (is.null(icorenew$corenesses)){
                   if(DEBUG){
                   global_nodenghr_cores[[qq]]$corenesses[i]   <- 0
                   global_nodenghr_cores[[qq]]$coressizes[i]   <- 0
                   global_nodenghr_cores[[qq]]$iscliques[i]    <- F
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <- NA
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <- NA
                   }else{
                   global_nodenghr_cores[[qq]]$corenesses[i]   <<- 0
                   global_nodenghr_cores[[qq]]$coressizes[i]   <<- 0
                   global_nodenghr_cores[[qq]]$iscliques[i]    <<- F
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <<- NA
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <<- NA
                   }

                }else{
                   if(DEBUG){
                   global_nodenghr_cores[[qq]]$corenesses[i]   <- icorenew$corenesses[1]
                   global_nodenghr_cores[[qq]]$coressizes[i]   <- icorenew$coressizes[1]
                   global_nodenghr_cores[[qq]]$iscliques[i]    <- icorenew$iscliques[1]
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <- icorenew$coresnodes[[1]]
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <- icorenew$coreslinks[[1]]
                   }else{
                   global_nodenghr_cores[[qq]]$corenesses[i]   <<- icorenew$corenesses[1]
                   global_nodenghr_cores[[qq]]$coressizes[i]   <<- icorenew$coressizes[1]
                   global_nodenghr_cores[[qq]]$iscliques[i]    <<- icorenew$iscliques[1]
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <<- icorenew$coresnodes[[1]]
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <<- icorenew$coreslinks[[1]]
                   }
                }
           }#qq


        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, 
                            global_nolinks_nghbr_cores[i],iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes",  "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(allcores.readTime/60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(allcores.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")

    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)    


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }
    
    # return if clique only
    if(cliqueOnly) {
        return (clique.readTime)
    }


    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_kcliquesNormal = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, cliqueOnly=T) {

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram",sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    linkpairsIndexed = cbind(selectedLinksIndex, c(1:no.links) )

    # find the maximum coreness, which will be used for inferring max(k)-clique
    #
    #maxcoreness  = find_maxcoreness(selectedLinksIndex)
    #max_clique   = maxcoreness + 1

    # max clique
    #
    degrees       = degree_ByLinkPairs (linkpairsIndexed, directed=F)
    maxdegree     = max(degrees$totallinks)

    maxcliqueVect = rep(maxdegree, no.nodes)    
    max_clique    = maxdegree

    #globalAdjLists  = makeAjacencyListsFromLinksIndex(selectedLinksIndex)

    adjlist.readTime <- proc.time()[3]-now


    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (iclq in c(max_clique:min_clique) ) {
        
        iclq.now<-proc.time()[3]

        # for normal clique identification, we use all links & nodes
        #
        iclqNodesIndex     = c(1:no.nodes) #selected nodes (index)
        iclqLinkpairsIndex = linkpairsIndexed       #selected links (indices)

        # the index of the lists correspond to the original node index
        #
        #iclqAdjlists   = makeAjacencyListsFromLinksIndex(linkpairsIndexed)
        iclqAdjlists   = globalAdjLists

        no.linksINiclq     = dim(linkpairsIndexed)[1]
        no.nodesINiclq     = no.nodes        # iterate each node to identify all iclq-cliques including this node

        #
        for ( j in c(1:no.nodesINiclq) ) {
           cliques_count = 0
           v = iclqNodesIndex[j]
           
           if ( length(iclqAdjlists[[v]])  >= iclq-1 ){

               mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
               print(mystr)


               cliques_count = recursive_CliquesNormal(curnode=v, setA=NULL, setB=NULL,
                                         adjlists     =iclqAdjlists,
                                         neighborlinks=iclqLinkpairsIndex, Kclique=iclq)#-1)
                  
              #cliquesFoundAfter = getOption("globalCliques")

           }#if ( length(iclqAdjlists[[v]]) >= iclq-1

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           selLeft    = iclqLinkpairsIndex[,1] != v
           selRight   = iclqLinkpairsIndex[,2] != v
           sel        = selLeft & selRight
           no.remains = sum(sel)

           # no enough links left, stop the current clique finding
           #
           if ( no.remains < iclq-1 ) {
              break
           } else {
              if(no.remains ==1) {
                 iclqLinkpairsIndex = rbind(iclqLinkpairsIndex[sel,])
              } else{
                 iclqLinkpairsIndex = iclqLinkpairsIndex[sel,]
              }
           }

           # ii.2) remove it from AdjList
           #
           iclqAdjlists[[v]] = -1           
           for ( m in c(1:no.nodesINiclq) ) {
                qq = iclqNodesIndex[m]
                mneighbors = setdiff(iclqAdjlists[[qq]], v)
                if (length(mneighbors)==0) {
                     iclqAdjlists[[qq]] = -1
                }else{
                     iclqAdjlists[[qq]] = mneighbors
                }
            }

            #iclqAdjlists
            #iclqLinkpairsIndex 

        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, no.linksINiclq, iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes", "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(adjlist.readTime /60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(adjlist.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")
    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    # return if clique only
    if(cliqueOnly) {
        return (1)
    }
        
    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    return (globalCliques)
}

find_kcliques_use_globalCORES = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, MemoyForSpeed=T, cliqueOnly=T) {

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram",sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    if(1==2) {
    newlinkpairs  = linkpairs
    degrees       = degree_ByLinkPairs (newlinkpairs, directed=F)
    rets          = makeAjacencyLists(inetmatrix=linkpairs, returnUniqueNames=T)
    #adjlists  = makeAjacencyLists(inetmatrix=linkpairs, returnUniqueNames=F)
    nodenames = rets[[1]]
    no.nodes  = rets[[2]]
    # put link index into the 3rd column, for k-community links output
    #
    no.links  = dim(newlinkpairs)[1]
    }

    linkpairsIndexed = cbind(selectedLinksIndex, c(1:no.links) )

    # find the maximum coreness, which will be used for inferring max(k)-clique
    #
    #maxcoreness  = find_maxcoreness(selectedLinksIndex)
    CORES = KCores_ALL(linksindex=selectedLinksIndex, 
                       nodesindex=NULL, minkcore=min_clique-1, name2integer=T, nodeMembershipOnly=F)#!(MemoyForSpeed) )
    #globalAdjLists  = makeAjacencyListsFromLinksIndex(selectedLinksIndex)

    allcores.readTime <- proc.time()[3]-now

    print(paste("time for searching cores: ", allcores.readTime, sep="") )

    if(MemoyForSpeed){
       coreindex = CORES[[1]]
    } else{
       coreindex = as.integer(colnames(CORES))
    }

    cliqueindex  = coreindex  + 1 
    max_coreness = max(coreindex )
    max_clique   = max(cliqueindex) + 1
    maxcliqueVect= rep(max_clique, no.nodes)
    no.cores     = length(coreindex)
    
    allSelnodesIndex=c(1:no.nodes)

    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (i in c(no.cores:1) ) {

        iclq= cliqueindex[i]
        iclq.now<-proc.time()[3]

        # remove the nodes with the total links less than min_clique
        #
        #rets = kcore_subnets(linksindex=selectedLinksIndex, 
        #                     nodesindex=NULL, kcore=iclq-1)
        #
        #if ( !is.null(rets[[1]]) ){
        #     iclqNodesIndex     = sort(rets[[1]]) #selected nodes (index)
        #     iclqLinkpairsIndex = rets[[2]]       #selected links (indices)
        #}else{
        #     next
        #}

        # nodes in (iclq-1)-core
        if (i==1) {
           iclqNodesIndex     = c(1:no.nodes)
           iclqLinkpairsIndex = selectedLinksIndex
        } else{
          if(MemoyForSpeed) {
             iclqNodesIndex     = (CORES[[2]])[[i]]
             iclqLinkpairsIndex = (CORES[[3]])[[i]]
          }else{
             iclqNodesIndex     = allSelnodesIndex[ CORES[,i] ]
             merged1= merge(iclqNodesIndex,selectedLinksIndex,by.x=1,by.y=1,all=F)
             iclqLinkpairsIndex = merge(iclqNodesIndex,merged1,by.x=1,by.y=2,all=F)
          }
        }

        # the index of the lists correspond to the original node index
        #
        #iclqAdjlists   = makeAjacencyListsFromLinksIndex(iclqLinkpairsIndex)
        if (i==1) {
           iclqAdjlists   = globalAdjLists
        } else{
           iclqAdjlists   = makeAjacencyListsFromSubsetnodesindex(
                            orgAdjLists=globalAdjLists, 
                            subsetNodes=iclqNodesIndex)
        }

        no.nodesINiclq = length(iclqNodesIndex)
        no.linksINiclq = dim(linkpairsIndexed)[1]

        # iterate each node to identify all iclq-cliques including this node
        #
        for ( j in c(1:no.nodesINiclq) ) {
           cliques_count = 0
           v = iclqNodesIndex[j]
           
           #if ( length(iclqAdjlists[[v]]) >= iclq-1 & maxcliqueVect[v]>=iclq ){
           if ( length(iclqAdjlists[[v]]) >= iclq-1){

               mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
               print(mystr)

               #cliques_count = recursive_CliquesNormal(curnode=v, setA=NULL, setB=NULL,
               cliques_count = recursive_Cliques_Super(curnode=v, setA=NULL, setB=NULL,
                                         adjlists     = iclqAdjlists,
                                         neighborlinks=iclqLinkpairsIndex, Kclique=iclq)#-1)

           }#if ( length(iclqAdjlists[[v]]) >= iclq-1
           cliques_count

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           selLeft    = iclqLinkpairsIndex[,1] != v
           selRight   = iclqLinkpairsIndex[,2] != v
           sel        = selLeft & selRight
           no.remains = sum(sel)

           # no enough links left, stop the current clique finding
           #
           if ( no.remains < iclq-1 ) {
              break
           } else {
              if(no.remains ==1) {
                 iclqLinkpairsIndex = rbind(iclqLinkpairsIndex[sel,])
              } else{
                 iclqLinkpairsIndex = iclqLinkpairsIndex[sel,]
              }
           }

           # ii.2) remove it from AdjList
           #
           iclqAdjlists[[v]] = -1           
           for ( m in c(1:no.nodesINiclq) ) {
                qq = iclqNodesIndex[m]
                mneighbors = setdiff(iclqAdjlists[[qq]], v)
                if (length(mneighbors)==0) {
                     iclqAdjlists[[qq]] = -1
                }else{
                     iclqAdjlists[[qq]] = mneighbors
                }
            }

            #iclqAdjlists
            #iclqLinkpairsIndex 

        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, length(iclqLinkpairsIndex), iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes", "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(allcores.readTime/60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(allcores.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")

    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)    


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }
    
    # return if clique only
    if(cliqueOnly) {
        return (globalCliques)
    }


    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    return (globalCliques)
}



############################### end of normal clique identification #########################################
#############################################################################################################











binaryVectorSimilarity = function(n1,n2, n12, measure="Correlation"){
   s11=n12/2
   s00=n12/2
   s10=n1-n12
   s01=n2-n12
   if (measure=="Correlation") {
      alpha      = sqrt((s10+s11)*(s01+s00)*(s11+s01)*(s00+s10))
      similarity = (s11*s00-s10*s01)/alpha
      similarity = (similarity +1)/2.0

   }else if (measure=="Tanimoto" | measure=="Russell-Rao"| measure=="Sokal-Michener"|measure=="Jaccard-Needham") {
      similarity = n12/(n1+n2-n12)

   }else if (measure=="Dice") {
      s11=n12
      similarity = 2*s11/(2*s11+s01+s10)

   }else if (measure=="Yule") {
      similarity = (s11*s00-s10*s01)/(s11*s00+s10*s01)
      similarity = (similarity +1)/2.0

   }else if (measure=="Rogers-Tanimoto" | measure=="Kulzinsky") {
      s11=n12
      similarity = s11/(s11+2*s10+2*s01)

   } else if (measure=="Sigmoid") {
      x = (s11*s00-s10*s01)/(s11*s00+s10*s01)
      x = (x +1)/2.0
      similarity = 1/(1+exp(-5*(x-0.5)) )

   } else if (measure=="SigmoidDice") {
      s11=n12
      x = 2*s11/(2*s11+s01+s10)
      similarity = 1/(1+exp(-5*(x-0.5)) )

   }else{ # Tanimoto
      similarity = n12/(n1+n2-n12)
   }
   
   return(similarity)
}

options(a=c(1:10))
test_globalvariable = function()
{
  a=getOption("a")
  print(a)
  a=c(1:3)
  print (a)
  a=getOption("a")
  print(a)

  a=getOption("a")  
  print (a)
}



#
#**********************************  END of k-cliques **********************************




##################################### FUNCTION ##########################
#
# Either by CDF of histgram count
#
# 1) get PDF by either histogram or scdf function
# 2) for each input data point, look at the nearest index in PDF and get 
#    its probability (P)
# 3) for the given P, check qnorm to get the inverse normal transformed value
#
# toy example
#    xx=c(1,1,2,2,3,3,4,4,4,4,4.5)
#    xx2=ecdf(xx)
#    cbind(xx2$x,xx2$y)
#

inverseNormalTransform = function(input, which="cdf", breakpoints=100) {

    x = round(input, 5)
    N = length(x)

    ######################### BY CDF #################################
    #
    if(which =="cdf"){
        #ecdf(x, datadensity='density')
        #ecdf(x, datadensity='hist')
        xecdf     = ecdf(x, xlab="test", datadensity='density', pl=F)

        xecdf.x   = xecdf$x
        xecdf.y   = xecdf$y
        xecdf.np  = length(xecdf.y)

    } else {

        xhist  = hist(x, br = breakpoints)
        #xhist2 = cbind(xhist$mids, xhist$counts, xhist$density)
        #xhist2

        xecdf.y = cumsum(xhist$counts)/N
        xecdf.x = xhist$mids
        xecdf.np  = length(xecdf.y)
    }

    # handle the boundary
    #
    xecdf.y[1]= (xecdf.y[2] + xecdf.y[1])/2
    xecdf.y[xecdf.np]= (xecdf.y[xecdf.np-1] + xecdf.y[xecdf.np])/2

    xmean = mean(input, na.rm=T)
    xsd   = sd(input,   na.rm=T)

    invernorm = qnorm(xecdf.y, mean=xmean, sd=xsd)
    xhist3    = cbind(xecdf.x, xecdf.y, invernorm)
    xhist3 


    # 2. inverse normal distribution
    #
    xp = rep(0, N)
    for(i in c(1:N) ){
      xp[i] = invernorm[ which.min(abs(x[i]-xecdf.x) ) ]
    }

    #xop = cbind(x, xp)
    #xop

    return (xp)
}

#
###################### END of Inverse Normal Transform ##############################




# ************************* Useful Functions ******************************************
#
#
adjusted_cortest<-function(x1, tt1, gender1, age1){
   
   nnaTT = !(is.na(tt1))
   nnaXX = !(is.na(x1))
   nna   = nnaTT & nnaXX
   x=x1[nna]
   tt=tt1[nna]

   gender=gender1[nna]
   age=age1[nna]

   res<-cor.test(resid(lm(tt~gender+age)),resid(lm(as.numeric(as.vector(x))~gender+age)))
   return( c(res$p.value, res$estimate) )
}

adjusted_cor_multitraits <- function(x, ys, gender, age){
    no.ys = dim(ys)[2]
    mres  = NULL
    for (j in c(1:no.ys)){
        jres = adjusted_cortest(x, ys[,j], gender, age)
        mres = c(mres, jres)
    }
    return (mres)
}


