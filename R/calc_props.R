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
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats hclust as.dist cutree qlnorm quantile


calc_props <- function(adjaMat, dissMat, assoMat, sPathDisconnected,
                       sPathNorm, sPathAlgo, weighted, isempty, clustMethod, 
                       clustPar, hubPar, hubQuant, lnormFit, connectivity,
                       weightDeg, normDeg, normBetw, normClose, normEigen,
                       jaccard = FALSE, jaccQuant = NULL){
  
  if(isempty){
    output <- list(clust = NULL, tree = NULL, deg = 0, deg_unnorm = 0,
                   betw = 0, betw_unnorm = 0, close = 0, close_unnorm = 0,
                   eigen = 0, eigen_unnorm = 0, lcSize = 0,
                   avDiss = 0, avDiss_lc = 0, avPath = 0, avPath_lc = 0,
                   vertconnect = 0, vertconnect_lc = 0, edgeconnect = 0, 
                   edgeconnect_lc = 0, clustCoef = 0, clustCoef_lc = 0,
                   density = 0, density_lc = 0, modul = 0, modul_lc = 0,
                   pnRatio = 0, pnRatio_lc = 0,  hubs = NULL, topdeg = NULL, 
                   topbetw = NULL, topclose = NULL, topeigen = NULL)
    return(output)
  }

  #== create igraph objects and decompose graph ================================
  # create graph from adjacency matrix

  net <- igraph::graph_from_adjacency_matrix(adjaMat, weighted=T,
                                             mode="undirected", diag=F)
  
  dg_net <- decompose.graph(net)
  
  # graph with only the largest component
  net_lc <- dg_net[[1]]
  lcSize <- length(V(net_lc))
  
  # adjacency of largest component
  adjaMat_lc <- as.matrix(as_adjacency_matrix(net_lc, attr="weight"))
  
  
  if(!weighted){
    dissMat[!is.infinite(dissMat)] <- 1
    diag(dissMat) <- 0
  }
  
  # dissimilarity of the largest component
  dissMat_lc <- dissMat[rownames(adjaMat_lc), colnames(adjaMat_lc)]
  
  #== graph objects of dissimilarity matrices ==================================
  
  dissMatnoInf <- dissMat
  dissMatnoInf[is.infinite(dissMatnoInf)] <- 0
  
  dissMatnoInf_lc <- dissMat_lc
  dissMatnoInf_lc[is.infinite(dissMatnoInf_lc)] <- 0
  
  # whole network
  dissnet <- igraph::graph_from_adjacency_matrix(dissMatnoInf, 
                                                 weighted=T,
                                                 mode="undirected", 
                                                 diag=F)
  
  # largest component
  #dissMatEdit_lc <- dissMatEdit[rownames(adjaMat_lc), colnames(adjaMat_lc)]
  
  dissnet_lc <- igraph::graph_from_adjacency_matrix(dissMatnoInf_lc, 
                                                    weighted=T,
                                                    mode="undirected", 
                                                    diag=F)
  
  
  #== shortest paths ===========================================================
  
  # whole network
  sPath <- distances(dissnet, algorithm = sPathAlgo)
  # largest component
  sPath_lc <- distances(dissnet_lc, algorithm = sPathAlgo)
  
  nNodes <- ncol(adjaMat)
  
  if(sPathDisconnected == "maxPath"){
    sPath[is.infinite(sPath)] <- nNodes
  }

  #== average dissimilarity/distance ===========================================

  if(weighted){
    dissVec <- dissMat[lower.tri(dissMat)]
    dissVec_lc <- dissMat_lc[lower.tri(dissMat_lc)]
    
    dissVec[is.infinite(dissVec)] <- NA
    dissVec_lc[is.infinite(dissVec_lc)] <- NA
    
    ### average dissimilarity
    avDiss <- mean(dissVec, na.rm = TRUE)
    avDiss_lc <- mean(dissVec_lc, na.rm = TRUE)
    
    # "normalized" shortest paths (steps with average dissimilarity)
    if(sPathNorm){
      sPath <- sPath / avDiss
      sPath_lc <- sPath_lc / avDiss_lc
    }
  } else{
    avDiss <- avDiss_lc <- 1
  }

  #== clustering ============================================================

  clust <- NULL
  tree <- NULL
  if(clustMethod != "none"){
    if(clustMethod == "hierarchical"){

      if(is.null(clustPar$method)){
        clustPar$method <- "average"
      }

      dissMat.tmp <- dissMat
      dissMat.tmp[is.infinite(dissMat.tmp)] <- 1
      tree <- hclust(as.dist(dissMat.tmp), method = clustPar$method)

      if(is.null(clustPar$k) & is.null(clustPar$h)){
        clustPar$k <- 3
      }
      clust <- do.call(cutree, list(tree = tree, k = clustPar$k, h = clustPar$h))
      names(clust) <- rownames(adjaMat)

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

      cltab <- table(clust)

      # cluster 0 assigned to elements in a cluster with only one element
      clust[clust %in% which(cltab == 1)] <- 0
    }

  }
  
  # largest component
  clust_lc <- clust[colnames(adjaMat_lc)]

  #== centrality measures ===================================================

  ### degree
  if(weightDeg){
    deg <- deg_unnorm <- strength(net)

  } else{
    deg <- degree(net, normalized = normDeg)
    deg_unnorm <- degree(net)
  }

  #-------------------------------
  ### betweenness centrality (based on distances)

  betw <- betweenness(dissnet, normalized = normBetw)
  betw_unnorm <- betweenness(dissnet)

  #-------------------------------
  ### closeness centrality (based on distances)

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

  #-------------------------------
  ### Eigenvector centrality

  if(length(dg_net) > 1){
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
    eigen <- eigen_centrality(net, scale = normEigen)$vector
    eigen_unnorm <- eigen_centrality(net, scale = FALSE)$vector
  }

  #== global network properties ================================================

  ### average shortest path length
  # complete network
  sPathVec <- sPath[lower.tri(sPath)]
  sPathVec[is.infinite(sPathVec)] <- NA
  avPath <- mean(sPathVec, na.rm = TRUE)
  if(is.na(avPath)) avPath <- 0
  
  # largest component
  sPathVec_lc <- sPath_lc[lower.tri(sPath_lc)]
  avPath_lc <- mean(sPathVec_lc, na.rm = TRUE)
  if(is.na(avPath_lc)) avPath_lc <- 0

  #-------------------------------
  # connectivity
  
  if(connectivity){
    # vertex connectivity
    vertconnect <- vertex_connectivity(net)
    vertconnect_lc <- vertex_connectivity(net_lc)
    
    # edge connectivity
    edgeconnect <- edge_connectivity(net)
    edgeconnect_lc <- edge_connectivity(net_lc)
    
  } else{
    vertconnect <- NA
    vertconnect_lc <- NA
    edgeconnect <- NA
    edgeconnect_lc <- NA
  }
  
  #-------------------------------
  ### clustering coefficient
  # complete network
  clustCoef <- transitivity(net, type = "global")
  if(is.na(clustCoef)) clustCoef <- 0
  
  # largest component
  clustCoef_lc <- transitivity(net_lc, type = "global")
  if(is.na(clustCoef_lc)) clustCoef_lc <- 0

  #-------------------------------
  ### modularity
  
  # complete network
  if(clustMethod != "none"){
    modul <- modularity(net, (clust+1))
  } else{
    modul <- NA
  }
  
  # largest component
  if(clustMethod != "none"){
    modul_lc <- modularity(net_lc, (clust_lc+1))
  } else{
    modul_lc <- NA
  }
  
  #-------------------------------
  ### density (relative number of edges)
  density <- edge_density(net)
  density_lc <- edge_density(net_lc)

  #-------------------------------
  ### ratio of positive to negative associations
  # complete network
  posAsso <- sum(assoMat[lower.tri(assoMat)] > 0)
  negAsso <- sum(assoMat[lower.tri(assoMat)] < 0)
  pnRatio <- posAsso / negAsso
  
  # complete network
  assoMat_lc <- assoMat[rownames(dissMat_lc), colnames(dissMat_lc)]
  posAsso <- sum(assoMat_lc[lower.tri(assoMat_lc)] > 0)
  negAsso <- sum(assoMat_lc[lower.tri(assoMat_lc)] < 0)
  pnRatio_lc <- posAsso / negAsso

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

  output <- list(clust = clust, tree = tree, 
                 deg = deg, deg_unnorm = deg_unnorm,
                 betw = betw, betw_unnorm = betw_unnorm,
                 close = close, close_unnorm = close_unnorm,
                 eigen = eigen, eigen_unnorm = eigen_unnorm,
                 lcSize = lcSize,
                 avDiss = avDiss, avDiss_lc = avDiss_lc,
                 avPath = avPath, avPath_lc = avPath_lc,
                 vertconnect = vertconnect, vertconnect_lc = vertconnect_lc,
                 edgeconnect = edgeconnect, edgeconnect_lc = edgeconnect_lc,
                 clustCoef = clustCoef, clustCoef_lc = clustCoef_lc,
                 density = density, density_lc = density_lc,
                 modul = modul, modul_lc = modul_lc,
                 pnRatio = pnRatio, pnRatio_lc = pnRatio_lc,
                 hubs = hubs, 
                 topdeg = topdeg, topbetw = topbetw, 
                 topclose = topclose, topeigen = topeigen)
  
  return(output)
}
