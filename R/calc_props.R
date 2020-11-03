#' @title Calculate network properties
#'
#' @description Calculates network properterties for a given adjacency matrix
#'
#' @param adjaMat adjacency matrix
#' @param dissMat dissimilarity matrix
#' @param assoMat association matrix
#' @param sPathDisconnected indicates how to handle disconnected components
#' @param sPathNorm logical. If TRUE, shortest paths are normalized by 
#'   average dissimilarity.
#' @param sPathAlgo character indicating the algorithm used for shortest path
#'   calculation.
#' @param weighted indicated whether the network is weighted
#' @param isempty indicator whether the network contains any edges.
#' @param clustMethod character indicating the clustering algorithm.
#' @param clustPar list with parameters passed to the clustering functions.
#' @param hubPar character vector with one or more centrality measures used
#'   for identifying hub nodes. Possible values are \code{degree},
#'   \code{betweenness}, \code{closeness}, and \code{eigenvector}.
#' @param hubQuant quantile used for determining hub nodes.
#' @param lnormFit hubs are nodes with a centrality value above the 95\%
#'   quantile of the fitted log-normal distribution (if \code{lnormFit = TRUE})
#'   or of the empirical distribution of centrality values.
#' @param connectivity logical indicating whether edge and vertex connectivity 
#'   should be calculated. Might be disabled to reduce execution time.
#' @param weightDeg if \code{TRUE}, the weighted degree is used.
#' @param normDeg,normBetw,normClose,normEigen if \code{TRUE}, a normalized
#'   version of the respective centrality values is returned. By default, all
#'   centralities are normalized.
#' @param jaccard shall the Jaccard index be calculated?
#' @param jaccQuant quantile for the Jaccard index
#'
#'
#' @importFrom igraph graph_from_adjacency_matrix decompose.graph
#' @importFrom stats hclust as.dist cutree qlnorm quantile


calc_props <- function(adjaMat, dissMat, assoMat, sPathNorm, sPathAlgo,
                       connectivity, normNatConnect, weighted, isempty, 
                       clustMethod, clustPar, hubPar, hubQuant, lnormFit, 
                       weightDeg, normDeg, normBetw, normClose, normEigen,
                       centrLCC, jaccard = FALSE, jaccQuant = NULL){
  
  if(isempty){
    output <- list(clust = NULL, tree = NULL, deg = 0, deg_unnorm = 0,
                   betw = 0, betw_unnorm = 0, close = 0, close_unnorm = 0,
                   eigen = 0, eigen_unnorm = 0, lccSize = 0,
                   avDiss = 0, avDiss_lcc = 0, avPath = 0, avPath_lcc = 0,
                   vertconnect = 0, vertconnect_lcc = 0, edgeconnect = 0, 
                   edgeconnect_lcc = 0, clustCoef = 0, clustCoef_lcc = 0,
                   density = 0, density_lcc = 0, modul = 0, modul_lcc = 0,
                   pnRatio = 0, pnRatio_lcc = 0,  hubs = NULL, topdeg = NULL, 
                   topbetw = NULL, topclose = NULL, topeigen = NULL)
    return(output)
  }

  #== Create igraph objects and decompose graph ================================
  # Create graph from adjacency matrix
  net <- igraph::graph_from_adjacency_matrix(adjaMat, weighted=T,
                                             mode="undirected", diag=F)
  
  dg_net <- igraph::decompose.graph(net, min.vertices=2)
  nComp <- length(dg_net)
  nNodes <- ncol(adjaMat)
  
  # Largest connected component (LCC)
  indLCC <- which.max(unlist(lapply(dg_net, function(x) length(V(x)))))
  
  net_lcc <- dg_net[[indLCC]]
  lccSize <- length(V(net_lcc))
  lccSizeRel <- lccSize / nNodes

  # Adjacency of the LCC
  adjaMat_lcc <- as.matrix(as_adjacency_matrix(net_lcc, attr="weight"))
  diag(adjaMat_lcc) <- 1
  
  # Names of nodes in the LCC
  lccNames <- colnames(adjaMat_lcc)
  
  if(!weighted){
    dissMat[!is.infinite(dissMat)] <- 1
    diag(dissMat) <- 0
  }
  
  # Dissimilarity of the LCC
  dissMat_lcc <- dissMat[lccNames, lccNames]
  
  #== Graph objects of dissimilarity matrices ==================================
  
  dissMatnoInf <- dissMat
  dissMatnoInf[is.infinite(dissMatnoInf)] <- 0
  
  dissMatnoInf_lcc <- dissMat_lcc
  dissMatnoInf_lcc[is.infinite(dissMatnoInf_lcc)] <- 0
  
  # Whole network
  dissnet <- igraph::graph_from_adjacency_matrix(dissMatnoInf, 
                                                 weighted=T,
                                                 mode="undirected", 
                                                 diag=F)
  
  # LCC
  dissnet_lcc <- igraph::graph_from_adjacency_matrix(dissMatnoInf_lcc, 
                                                    weighted=T,
                                                    mode="undirected", 
                                                    diag=F)
  
  #== Clustering ============================================================

  clust <- clust_lcc <- NULL
  tree <- tree_lcc <- NULL
  
  if(clustMethod != "none"){
    if(clustMethod == "hierarchical"){

      if(is.null(clustPar$method)){
        clustPar$method <- "average"
      }

      dissMat.tmp <- dissMat
      dissMat.tmp[is.infinite(dissMat.tmp)] <- 1
      
      dissMat_lcc.tmp <- dissMat_lcc
      dissMat_lcc.tmp[is.infinite(dissMat_lcc.tmp)] <- 1
      
      tree <- hclust(as.dist(dissMat.tmp), method = clustPar$method)
      tree_lcc <- hclust(as.dist(dissMat_lcc.tmp), method = clustPar$method)
      rm(dissMat.tmp, dissMat_lcc.tmp)

      if(is.null(clustPar$k) & is.null(clustPar$h)){
        clustPar$k <- 3
      }
      
      clust <- do.call(cutree, list(tree = tree, k = clustPar$k, 
                                    h = clustPar$h))
      clust_lcc <- do.call(cutree, list(tree = tree_lcc, k = clustPar$k, 
                                        h = clustPar$h))
      
      names(clust) <- rownames(adjaMat)
      names(clust_lcc) <- rownames(adjaMat_lcc)

    } else{

      if(clustMethod == "cluster_edge_betweenness"){
        clustres <- do.call(clustMethod,  c(list(graph = net,
                                                 weights = E(dissnet)$weight),
                                            clustPar))
      } else{
        clustres <- do.call(clustMethod, c(list(net), clustPar))
      }

      clust <- clustres$membership
      names(clust) <- clustres$names
      
      # LCC
      clust_lcc <- clust[lccNames]
    }

    cltab <- table(clust)
    cltab_lcc <- table(clust_lcc)
    
    # cluster 0 assigned to elements in a cluster with only one element
    clust[clust %in% which(cltab == 1)] <- 0
    clust_lcc[clust_lcc %in% which(cltab_lcc == 1)] <- 0
  }
  

  
  #== Shortest paths ===========================================================
  
  # Whole network
  sPath <- distances(dissnet, algorithm = sPathAlgo)
  
  # LCC
  sPath_lcc <- distances(dissnet_lcc, algorithm = sPathAlgo)

  #== Average dissimilarity/distance ===========================================
  
  if(weighted){
    dissVec <- dissMat[lower.tri(dissMat)]
    dissVec[is.infinite(dissVec)] <- NA
    
    dissVec_lcc <- dissMat_lcc[lower.tri(dissMat_lcc)]
    dissVec_lcc[is.infinite(dissVec_lcc)] <- NA
    
    ### average dissimilarity
    avDiss <- mean(dissVec, na.rm = TRUE)
    avDiss_lcc <- mean(dissVec_lcc, na.rm = TRUE)

    # "normalized" shortest paths (units with average dissimilarity)
    if(sPathNorm){
      sPath <- sPath / avDiss
      sPath_lcc <- sPath_lcc / avDiss_lcc
    }
  } else{
    avDiss <- avDiss_lcc <- 1
  }

  #== global network properties ================================================

  ### Average shortest path length 
  
  sPathVec <- sPath[lower.tri(sPath)]
  sPathVec[is.infinite(sPathVec)] <- NA
  avPath <- mean(sPathVec, na.rm = TRUE)
  if(is.na(avPath)) avPath <- 0

  sPathVec_lcc <- sPath_lcc[lower.tri(sPath_lcc)]
  avPath_lcc <- mean(sPathVec_lcc, na.rm = TRUE)
  if(is.na(avPath_lcc)) avPath_lcc <- 0
  
  #-------------------------------
  ### Connectivity 
  
  if(connectivity){
    # vertex connectivity
    vertconnect <- vertex_connectivity(net)
    vertconnect_lcc <- vertex_connectivity(net_lcc)
    
    # edge connectivity
    edgeconnect <- edge_connectivity(net)
    edgeconnect_lcc <- edge_connectivity(net_lcc)
  } else{
    vertconnect <- vertconnect_lcc <- NA
    edgeconnect <- edgeconnect_lcc <- NA
  }
  
  #-------------------------------
  ### Natural connectivity

  natConnect <- pulsar::natural.connectivity(adjaMat, norm = normNatConnect)
  natConnect_lcc <- pulsar::natural.connectivity(adjaMat_lcc,
                                                 norm = normNatConnect)
  
  #-------------------------------
  ### Clustering coefficient
  
  # Complete network
  clustCoef <- transitivity(net, type = "global")
  if(is.na(clustCoef)) clustCoef <- 0
  
  # LCC
  clustCoef_lcc <- transitivity(net_lcc, type = "global")
  if(is.na(clustCoef_lcc)) clustCoef_lcc <- 0
  
  #-------------------------------
  ### Modularity
  
  # Complete network
  if(clustMethod != "none"){
    modul <- modularity(net, (clust+1))
  } else{
    modul <- NA
  }
  
  # LCC
  if(clustMethod != "none"){
    modul_lcc <- modularity(net_lcc, (clust_lcc+1))
  } else{
    modul_lcc <- NA
  }
  
  #-------------------------------
  ### Density (relative number of edges)
  density <- edge_density(net)
  density_lcc <- edge_density(net_lcc)
  
  #-------------------------------
  ### Ratio of positive to negative associations

  # Complete network
  if(is.null(assoMat)){ # no negative dissimilarities possible
    pnRatio <- 1
    pnRatio_lcc <- 1
    
  } else{
    posAsso <- sum(assoMat[lower.tri(assoMat)] > 0)
    negAsso <- sum(assoMat[lower.tri(assoMat)] < 0)
    pnRatio <- posAsso / negAsso
    if(is.infinite(pnRatio) || is.nan(pnRatio)){
      pnRatio <- posAsso
    }
    
    # LCC
    assoMat_lcc <- assoMat[rownames(dissMat_lcc), colnames(dissMat_lcc)]
    posAsso <- sum(assoMat_lcc[lower.tri(assoMat_lcc)] > 0)
    negAsso <- sum(assoMat_lcc[lower.tri(assoMat_lcc)] < 0)
    pnRatio_lcc <- posAsso / negAsso
    if(is.infinite(pnRatio_lcc) || is.nan(pnRatio_lcc)){
      pnRatio_lcc <- posAsso
    }
  }

  #== Centrality measures ===================================================

  ### degree
  if(weightDeg){
    deg <- deg_unnorm <- strength(net)
    deg_lcc <- strength(net_lcc)

  } else{
    deg <- degree(net, normalized = normDeg)
    deg_unnorm <- degree(net)
  }
  
  if(centrLCC){
    deg[!names(deg) %in% lccNames] <- 0
    deg_unnorm[!names(deg_unnorm) %in% lccNames] <- 0
  }

  #-------------------------------
  ### betweenness centrality (based on distances)

  if(centrLCC){
    betw <- betw_unnorm <- rep(0, ncol(adjaMat))
    names(betw) <- names(betw_unnorm) <- colnames(adjaMat)
    
    betw.tmp <- betweenness(dissnet_lcc, normalized = normBetw)
    betw_unnorm.tmp <- betweenness(dissnet_lcc)
    
    betw[names(betw.tmp)] <- betw.tmp
    betw_unnorm[names(betw_unnorm.tmp)] <- betw_unnorm.tmp
    
  } else{
    betw <- betweenness(dissnet, normalized = normBetw)
    betw_unnorm <- betweenness(dissnet)
  }
  
  betw[is.nan(betw)] <- 0
  betw_unnorm[is.nan(betw_unnorm)] <- 0

  #-------------------------------
  ### closeness centrality (based on distances)

  if(centrLCC){
    close_unnorm <- rep(0, ncol(adjaMat))
    names(close_unnorm) <- colnames(adjaMat)
    
    sPath_rev <- 1/sPath_lcc
    diag(sPath_rev) <- NA
    
    close_unnorm.tmp <- sapply(1:nrow(sPath_rev), function(i){
      sum(sPath_rev[i,], na.rm = TRUE)
    })
    names(close_unnorm.tmp) <- lccNames
    close_unnorm[lccNames] <- close_unnorm.tmp
    
    if(normClose){
      # normalize closeness by n-1 
      close <- close_unnorm / (lccSize - 1)
    } else{
      close <- close_unnorm
    }
    
  } else{
    #if(nComp > 1 && sPathNorm){
      #warning("Normalized shortest paths depend on average dissimilarity, which
      #may be biased because infinite paths between unconnected nodes are ignored. 
      #Hence, closeness centrality may also be biased and should either be based 
      #on non-normalized shortest paths or be computed only for the LCC.")
    #}
    
    sPath_rev <- 1/sPath
    diag(sPath_rev) <- NA
    
    close_unnorm <- sapply(1:nrow(sPath_rev), function(i){
      sum(sPath_rev[i,], na.rm = TRUE)
    })
    names(close_unnorm) <- colnames(adjaMat)
    
    if(normClose){
      # normalize closeness by n-1 
      close <- close_unnorm / (nNodes - 1)
    } else{
      close <- close_unnorm
    }
  }

  #-------------------------------
  ### Eigenvector centrality

  if(nComp > 1 && !centrLCC){
    dgcount <- unlist(lapply(dg_net, vcount))
    dgcount[dgcount == 1] <- 0
    vnumb <- sum(unlist(dgcount))

    ev <- numeric(0)
    for(i in seq_along(dg_net)){
      ev <- c(ev, eigen_centrality(dg_net[[i]], scale = FALSE)$vector * (dgcount[i] / vnumb))
    }

    eigen_unnorm <- ev[colnames(adjaMat)]

    if(normEigen){
      eigen <- eigen_unnorm / max(eigen_unnorm)
    } else{
      eigen <- eigen_unnorm
    }

  } else{
    
    eigen.tmp <- eigen_centrality(net_lcc, scale = normEigen)$vector
    eigen_unnorm.tmp <- eigen_centrality(net_lcc, scale = FALSE)$vector
    
    if(nComp > 1){
      eigen <- eigen_unnorm <- numeric(nNodes)
      names(eigen) <- names(eigen_unnorm) <- colnames(adjaMat)
      eigen[names(eigen.tmp)] <- eigen.tmp
      eigen_unnorm[names(eigen_unnorm.tmp)] <- eigen_unnorm.tmp
    } else{
      eigen <- eigen.tmp
      eigen_unnorm <- eigen_unnorm.tmp
    }
  }

  #== hubs and Jaccard index ================================================

  if(lnormFit){
    # identify nodes with highest centrality value
    pdeg <- try(MASS::fitdistr(deg_unnorm[deg_unnorm>0], "lognormal")$estimate, silent = TRUE)
    if(class(pdeg) == "try-error"){
      topdeg <- hubdeg <- NULL
    } else{
      hubdeg <- names(deg_unnorm[deg_unnorm > qlnorm(hubQuant, pdeg[1], pdeg[2])])
      if(jaccard){
        topdeg <- names(deg_unnorm[deg_unnorm > qlnorm(jaccQuant, pdeg[1], pdeg[2])])
      } else{
        topdeg <- NULL
      }
    }

    pbetw <- try(MASS::fitdistr(betw_unnorm[betw_unnorm>0], "lognormal")$estimate, silent = TRUE)
    if(class(pbetw) == "try-error"){
      topbetw <- hubbetw <- NULL
    } else{
      hubbetw <- names(betw_unnorm[betw_unnorm > qlnorm(hubQuant, pbetw[1], pbetw[2])])
      if(jaccard){
        topbetw <- names(betw_unnorm[betw_unnorm > qlnorm(jaccQuant, pbetw[1], pbetw[2])])
      } else{
        topbetw <- NULL
      }
    }

    pclose <- try(MASS::fitdistr(close_unnorm[close_unnorm>0], "lognormal")$estimate, silent = TRUE)
    if(class(pclose) == "try-error"){
      topclose <- hubclose <- NULL
    } else{
      hubclose <- names(close_unnorm[close_unnorm > qlnorm(hubQuant, pclose[1], pclose[2])])
      if(jaccard){
        topclose <- names(close_unnorm[close_unnorm > qlnorm(jaccQuant, pclose[1], pclose[2])])
      } else{
        topclose <- NULL
      }
    }

    peigen <- try(MASS::fitdistr(eigen_unnorm[eigen_unnorm>0], "lognormal")$estimate, silent = TRUE)
    if(class(peigen) == "try-error"){
      topeigen <- hubeigen <- NULL
    } else{
      hubeigen <- names(eigen_unnorm[eigen_unnorm > qlnorm(hubQuant, peigen[1], peigen[2])])
      if(jaccard){
        topeigen <- names(eigen_unnorm[eigen_unnorm > qlnorm(jaccQuant, peigen[1], peigen[2])])
      } else{
        topeigen <- NULL
      }
    }

  } else{
    hubdeg <- names(deg[deg > quantile(deg, hubQuant)])
    if(jaccard){
      topdeg <- names(deg[deg > quantile(deg, jaccQuant)])
    } else{
      topdeg <- NULL
    }

    hubbetw <- names(betw[betw > quantile(betw, hubQuant)])
    if(jaccard){
      topbetw <- names(betw[betw > quantile(betw, jaccQuant)])
    } else{
      topbetw <- NULL
    }

    hubclose <- names(close[close > quantile(close, hubQuant)])
    if(jaccard){
      topclose <- names(close[close > quantile(close, jaccQuant)])
    } else{
      topclose <- NULL
    }

    hubeigen <- names(eigen[eigen > quantile(eigen, hubQuant)])
    if(jaccard){
      topeigen <- names(eigen[eigen > quantile(eigen, jaccQuant)])
    } else{
      topeigen <- NULL
    }
  }

  ###############
  # identify hub nodes
  hubparams <-  list()
  if("degree" %in% hubPar){
    hubparams <- c(hubparams, list(hubdeg))
  }
  if("betweenness" %in% hubPar){
    hubparams <- c(hubparams, list(hubbetw))
  }
  if("closeness" %in% hubPar){
    hubparams <- c(hubparams, list(hubclose))
  }
  if("eigenvector" %in% hubPar){
    hubparams <- c(hubparams, list(hubeigen))
  }

  hubs <- Reduce(intersect, hubparams)

  #========================================================================

  output <- list(nComp = nComp, lccNames = lccNames,
                 clust = clust, tree = tree, 
                 clust_lcc = clust_lcc, tree_lcc = tree_lcc,
                 deg = deg, deg_unnorm = deg_unnorm,
                 betw = betw, betw_unnorm = betw_unnorm,
                 close = close, close_unnorm = close_unnorm,
                 eigen = eigen, eigen_unnorm = eigen_unnorm,
                 lccSize = lccSize, lccSizeRel = lccSizeRel,
                 avDiss = avDiss, avDiss_lcc = avDiss_lcc,
                 avPath = avPath, avPath_lcc = avPath_lcc,
                 vertconnect = vertconnect, vertconnect_lcc = vertconnect_lcc,
                 edgeconnect = edgeconnect, edgeconnect_lcc = edgeconnect_lcc,
                 natConnect = natConnect, natConnect_lcc = natConnect_lcc,
                 clustCoef = clustCoef, clustCoef_lcc = clustCoef_lcc,
                 density = density, density_lcc = density_lcc,
                 modul = modul, modul_lcc = modul_lcc,
                 pnRatio = pnRatio, pnRatio_lcc = pnRatio_lcc,
                 hubs = hubs, 
                 topdeg = topdeg, topbetw = topbetw, 
                 topclose = topclose, topeigen = topeigen)
  
  return(output)
}
