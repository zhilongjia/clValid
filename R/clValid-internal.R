####################################################################
## Internal functions for clValid
####################################################################


####################################################################
## vClusters() function
####################################################################
## Arguments
####################################################################
## mat - data matrix
## clMethod - clustering method
## nClust - minimun number of clusters to evaluate
## nclustMax - maximun number of clusters to evaluate
## ... - arguments to pass to clustering functions
####################################################################
## Value
####################################################################
## List with:
## clusterObj - clustering results
## measures - validation measures
####################################################################


vClusters <- function(mat,clMethod,nClust,nclustMax, validation,
                      Distmat, method, metric, annotation, GOcategory,
                      goTermFreq, neighbSize, dropEvidence, ncore, verbose, ... ) {
  ###########################################################
  #   #only for test
    require(clValid)
    data(mouse)
    express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
    rownames(express) <- mouse$ID[1:25]
    obj=express
    mat=as.matrix(obj)
    nClust = 2:4
    nclustMax = max(nClust)
    ncore=2
    goTermFreq=0.05
    dropEvidence=NULL
    verbose=FALSE
    clMethods=c("hierarchical","kmeans","pam")
    clMethod = clMethods[2]
    validation="biological"
    metric="euclidean"
    method="average"
    GOcategory="all"
    neighbSize=10
    fc <- tapply(rownames(express),mouse$FC[1:25], c)
    fc <- fc[-match( c("EST","Unknown"), names(fc))]
    annotation=annotationListToMatrix(fc,rownames(express))
  ####################################################
  measNames <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                 if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                 if("biological"%in%validation) c("BHI","BSI"))
  measures <- vector(length=length(measNames), mode="numeric")
  names(measures) <- measNames

  switch(clMethod,
         hierarchical = {
           clusterObj <- hclust(Distmat,method)
         },
         diana = {
           clusterObj <- diana(Distmat)
           #clusterObj <- diana(Distmat, ...)
         },
         kmeans = {
           clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust
           clusterObjInit <- hclust(Distmat,method)
         },
         agnes = {
           clusterObj <- agnes(Distmat, method=method)
           #clusterObj <- agnes(Distmat, method=method, ...)
         },
         ## otherwise - sota, fanny, som, model, pam, clara
         { clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust })

  library(doMC)
  registerDoMC(ncore)
  clusterList = list()
  clusterList <- foreach (nc = nClust) %dopar% {
      #ind tracks number clusters
      #ind=which(nClust==nc)
      #ind=1
      switch(clMethod,
           kmeans = {
             initial <- tapply(mat, list(rep(cutree(clusterObjInit,nc),ncol(mat)),col(mat)),
                               function(x) mean(x, na.rm=TRUE))
             ######################
#              aa <- list(matrix(rep(cutree(clusterObjInit,nc),ncol(mat)), ncol=ncol(mat)), col(mat))
#              #aa <- list(rep(cutree(clusterObjInit,nc),ncol(mat)), col(mat))
#              tapply(mat, aa, function(x) mean(x, na.rm=TRUE))
             #####################################
             if(length(dup <- which(duplicated(initial)))>0) {
               for(dupi in dup) 
                 initial[dupi,] <- initial[dupi,] + jitter(initial[dupi,])
             }
             dimnames(initial) <- list(NULL,dimnames(mat)[[2]])
             clusterObj <- kmeans(mat,initial)
             #clusterObj <- kmeans(mat,initial,...)
             cluster <- clusterObj$cluster
           },
           fanny = {
             clusterObj <- fanny(Distmat, nc)
             #clusterObj <- fanny(Distmat, nc, ...)
             cluster <- clusterObj$clustering
           },
           model = {
             clusterObj <- Mclust(mat,nc)
             #clusterObj <- Mclust(mat,nc, ...)
             cluster <- clusterObj$classification
           },
           som = {
             clusterObj <- som(mat, grid=somgrid(1,nc))
             #clusterObj <- som(mat, grid=somgrid(1,nc), ...)
             cluster <- clusterObj$unit.classif
           },
           pam = {
             clusterObj <- pam(Distmat, nc)
             #clusterObj <- pam(Distmat, nc, ...)
             cluster <- clusterObj$clustering
           },
           clara = {
             clusterObj <- clara(mat, nc, metric=ifelse(metric=="correlation","euclidean",metric))
             #clusterObj <- clara(mat, nc, metric=ifelse(metric=="correlation","euclidean",metric), ...)
             cluster <- clusterObj$clustering
           },
           sota = {
             clusterObj <- sota(mat,nc-1)
             cluster <- clusterObj$clust

           },
           ## otherwise - hierarchical, diana, agnes
           {cluster <- cutree(clusterObj,nc)}) 
    if(length(table(cluster))!=nc) {
      warning(paste(clMethod, "unable to find",nc,"clusters, returning NA for these validation measures"))
      measures <- NA
      #ind <- ind+1
      #next()
    }
    
    ## internal validation measures
    if ("internal"%in%validation) {
      measures["Dunn"] <- dunn(Distmat ,cluster)
      measures["Silhouette"] <- mean(silhouette(cluster, dmatrix=as.matrix(Distmat))[,3])
      measures["Connectivity"] <- connectivity(Distmat ,cluster, neighbSize=neighbSize)
      if(verbose) print(paste("Finished internal validation,", clMethod, nc, "clusters"))
    }
    
    if("biological"%in%validation) {
      measures["BHI"] <- BHI(cluster,annotation=annotation, names=rownames(mat),
                                 category=GOcategory, dropEvidence=dropEvidence)
      if(verbose & "biological"%in%validation)
        print(paste("Finished BHI,", clMethod, nc, "clusters"))      
    }

    ## stability validation measures
    if ("stability"%in%validation | "biological"%in%validation) {
      co.del <- 0 ## for use in verbose printing of progress
      for (del in 1:ncol(mat)) {
        matDel <- mat[,-del]               ## matDel <- as.matrix(matDel)
        if(metric=="correlation") {
          DistDel <- as.dist(1-cor(t(matDel), use="pairwise.complete.obs"))
        } else {
          DistDel <- dist(matDel,method=metric)
        }
        switch(clMethod,
               hierarchical = clusterObjDel <- hclust(DistDel,method),
               kmeans = clusterObjInitDel <- hclust(DistDel,method),
               #diana = clusterObjDel <- diana(DistDel, ...),
               #agnes = clusterObjDel <- agnes(DistDel, method=method, ...),
               #clara = clusterObjDel <- clara(matDel,nc,metric=ifelse(metric=="correlation","euclidean",metric), ...)
               diana = clusterObjDel <- diana(DistDel),
               agnes = clusterObjDel <- agnes(DistDel, method=method),
               clara = clusterObjDel <- clara(matDel,nc,metric=ifelse(metric=="correlation","euclidean",metric)))


        switch(clMethod,
               kmeans = {
                 initialDel <- tapply(matDel, list(rep(cutree(clusterObjInitDel,nc),
                                                       ncol(matDel)), col(matDel)),
                                      function(x) mean(x, na.rm=TRUE))
                 if(length(dup <- which(duplicated(initialDel)))>0) {
                   for(dupi in dup) 
                     initialDel[dupi,] <- initialDel[dupi,] + jitter(initialDel[dupi,])
                 }
                 dimnames(initialDel) <- list(NULL,dimnames(matDel)[[2]])
                 kmdel <- kmeans(matDel,initialDel)
                 #kmdel <- kmeans(matDel,initialDel, ...)
                 clusterDel <- kmdel$cluster
               },
               fanny = {
                 hfdel <- fanny(DistDel, nc)
                 #hfdel <- fanny(DistDel, nc, ...)
                 clusterDel <- hfdel$clustering
               },
               model = {
                 clusterDel <- Mclust(matDel,nc)$classification
                 #clusterDel <- Mclust(matDel,nc, ...)$classification
               },
               som = {
                 hsdel <- try(som(matDel, grid=somgrid(1,nc)))
                 #hsdel <- try(som(matDel, grid=somgrid(1,nc), ...))
                 clusterDel <- hsdel$unit.classif
               },
               pam = {
                 clusterDel <- pam(DistDel, nc, cluster.only=TRUE)
                 #clusterDel <- pam(DistDel, nc, cluster.only=TRUE, ...)
               },
               clara = {
                 clusterDel <- clusterObjDel$clustering
               },
               sota = {
                 clusterDel <- sota(matDel,nc-1)$clust
               },
               ## otherwise - hierarchical, diana, agnes
               {clusterDel <- cutree(clusterObjDel,nc)})

        if("stability"%in%validation) {
          stabmeas <- stability(mat, Distmat, del, cluster, clusterDel)
          measures["APN"] <- measures["APN"] + stabmeas["APN"]
          measures["AD"]  <- measures["AD"]  + stabmeas["AD"]
          measures["ADM"] <- measures["ADM"] + stabmeas["ADM"]
          measures["FOM"] <- measures["FOM"] + stabmeas["FOM"]
        }
        if("biological"%in%validation) {
          tmp <- BSI(cluster,clusterDel,annotation=annotation,
                     names=rownames(mat), category=GOcategory, goTermFreq=goTermFreq,
                     dropEvidence=dropEvidence)
          measures["BSI"] <- measures["BSI"] + tmp
        }
        ## VERBOSE printing
        if (del/ncol(mat) > 0.25 & co.del==0) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 25% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 25% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
        else if (del/ncol(mat) > 0.50 & co.del==1) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 50% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 50% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
        else if (del/ncol(mat) > 0.75 & co.del==2) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 75% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 75% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
      } #END OF del LOOP
      if(verbose & "stability"%in%validation)
        print(paste("Finished stability validation,", clMethod, nc, "clusters"))
      if(verbose & "biological"%in%validation)
      print(paste("Finished BSI,", clMethod, nc, "clusters"))
    } #END of STABILITY measures
    #ind <- ind+1  
    ## if(verbose) print(paste("Finished with", nc, "clusters"))
    list(clusterObj=clusterObj, measures=measures)
  } #END OF NC LOOP
  
  
  if ("stability"%in%validation) {
    measures["APN"] <- measures["APN"]/ncol(mat)
    measures["AD"] <-  measures["AD"]/ncol(mat)
    measures["ADM"] <- measures["ADM"]/ncol(mat)
    measures["FOM"] <- measures["FOM"]/ncol(mat)  ## little different from Yeung paper (doesn't do this)
  }
  if ("biological"%in%validation) {
    measures["BSI"] <- measures["BSI"]/ncol(mat)
  }
  ########################################################
  #combine the paralled results
  measuresComb <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  rownames(measuresComb) <- measNames
  colnames(measuresComb) <- nClust
  for (j in seq(length(clusterList))){
    clusterObj[[j]] <- clusterList[[j]]$clusterObj  
    measuresComb[,j] <- clusterList[[j]]$measures} 
  list(clusterObj=clusterObj, measures=measuresComb)
}

