calc_diff_props <- function(adja1, adja2, dissMat1, dissMat2, assoMat1, assoMat2,
                            sPathDisconnected, sPathNorm, sPathAlgo,
                            weighted, clustMethod, clustPar, clustPar2, hubPar,
                            hubQuant, jaccQuant, lnormFit, connectivity,
                            weightDeg, normDeg, normBetw, normClose, normEigen,
                            testJacc = TRUE, jaccTestGreater = FALSE,
                            testRand = TRUE, nPermRand = 1000){

  isempty1 <- all(adja1[lower.tri(adja1)] == 0)
  isempty2 <- all(adja2[lower.tri(adja2)] == 0)

  if(isempty1 & isempty2) stop("There are no connected nodes in both networks.")

  # calculate network properties
  props1 <- calc_props(adjaMat = adja1, dissMat = dissMat1, assoMat = assoMat1, 
                       sPathDisconnected = sPathDisconnected,
                       sPathNorm = sPathNorm, sPathAlgo = sPathAlgo,
                       weighted = weighted, isempty = isempty1, 
                       clustMethod = clustMethod, clustPar = clustPar, 
                       hubPar = hubPar, hubQuant = hubQuant, 
                       lnormFit = lnormFit, connectivity = connectivity,
                       weightDeg = weightDeg, normDeg = normDeg,
                       normBetw = normBetw, normClose = normClose,
                       normEigen = normEigen, 
                       jaccard = TRUE, jaccQuant = jaccQuant)

  props2 <- calc_props(adjaMat = adja2, dissMat = dissMat2, assoMat = assoMat2, 
                       sPathDisconnected = sPathDisconnected,
                       sPathNorm = sPathNorm, sPathAlgo = sPathAlgo,
                       weighted = weighted, isempty = isempty2, 
                       clustMethod = clustMethod, clustPar = clustPar2, 
                       hubPar = hubPar, hubQuant = hubQuant, 
                       lnormFit = lnormFit, connectivity = connectivity,
                       weightDeg = weightDeg, normDeg = normDeg,
                       normBetw = normBetw, normClose = normClose,
                       normEigen = normEigen, 
                       jaccard = TRUE, jaccQuant = jaccQuant)


  #== differences in network properties ========================================
  
  ### centralities
  
  diffdeg <- sort(props1$deg - props2$deg[names(props1$deg)])
  diffbetw <- sort(props1$betw - props2$betw[names(props1$betw)])
  diffclose <- sort(props1$close - props2$close[names(props1$close)])
  diffeigen <- sort(props1$eigen - props2$eigen[names(props1$eigen)])

  # absolute differences
  absdiffdeg <- abs(props1$deg - props2$deg[names(props1$deg)])
  absdiffbetw <- abs(props1$betw - props2$betw[names(props1$betw)])
  absdiffclose <- abs(props1$close - props2$close[names(props1$close)])
  absdiffeigen <- abs(props1$eigen - props2$eigen[names(props1$eigen)])

  #--------------------------------------------------------------------------
  ### Average dissimilarity
  
  avDiss1 <- props1$avDiss
  avDiss2 <- props2$avDiss
  
  avDiss1_lc <- props1$avDiss_lc
  avDiss2_lc <- props2$avDiss_lc
  
  if(is.na(avDiss1)){
    avDiss1 <- 0
    avDiss1_lc <- 0
  } 
  
  if(is.na(avDiss2)){
    avDiss2 <- 0
    avDiss2_lc <- 0
  }
  
  diffdiss <- abs(avDiss1 - avDiss2)
  diffdiss_lc <- abs(avDiss1_lc - avDiss2_lc)
  
  #--------------------------------------------------------------------------
  ### Average Path Length (or Mean Distance)

  avPath1 <- props1$avPath
  avPath2 <- props2$avPath
  
  avPath1_lc <- props1$avPath_lc
  avPath2_lc <- props2$avPath_lc

  if(is.na(avPath1)){
    avPath1 <- 0
    avPath1_lc <- 0
  } 
  
  if(is.na(avPath2)){
    avPath2 <- 0
    avPath2_lc <- 0
  }

  diffpath <- abs(avPath1 - avPath2)
  diffpath_lc <- abs(avPath1_lc - avPath2_lc)

  #--------------------------------------------------------------------------
  ### Clustering

  clustCoef1 <- props1$clustCoef
  clustCoef2 <- props2$clustCoef
  
  clustCoef1_lc <- props1$clustCoef_lc
  clustCoef2_lc <- props2$clustCoef_lc

  if(is.na(clustCoef1)){
    clustCoef1 <- 0
    clustCoef1_lc <- 0
  } 
  
  if(is.na(clustCoef2)){
    clustCoef2 <- 0
    clustCoef2_lc <- 0
  } 
  
  diffclustcoef <- abs(clustCoef1 - clustCoef2)
  diffclustcoef_lc <- abs(clustCoef1_lc - clustCoef2_lc)

  #--------------------------------------------------------------------------
  # differential Modularity
  if(clustMethod != "none"){
    modul1 <- props1$modul
    modul2 <- props2$modul
    modul1_lc <- props1$modul_lc
    modul2_lc <- props2$modul_lc
    diffmod <- abs(modul1 - modul2)
    diffmod_lc <- abs(modul1_lc - modul2_lc)
  } else{
    modul1 <- modul2 <- modul1_lc <- modul2_lc <- NULL
    diffmod <- diffmod_lc <- NULL
  }

  #--------------------------------------------------------------------------
  if(connectivity){
    # vertex connectivity
    diffvertconnect <- abs(props1$vertconnect - props2$vertconnect)
    diffvertconnect_lc <- abs(props1$vertconnect_lc - props2$vertconnect_lc)
    
    # edge connectivity
    diffedgconnect <- abs(props1$edgeconnect - props2$edgeconnect)
    diffedgconnect_lc <- abs(props1$edgeconnect_lc - props2$edgeconnect_lc)
  } else{
    diffvertconnect <- diffvertconnect_lc <- NULL
    diffedgconnect <- diffedgconnect_lc <- NULL
  }

  #--------------------------------------------------------------------------
  # relative number of edges(density)
  diffdensity <- abs(props1$density - props2$density)
  diffdensity_lc <- abs(props1$density_lc - props2$density_lc)
  
  #--------------------------------------------------------------------------
  # positive-to-negative ratio
  diffpnratio <- abs(props1$pnRatio - props2$pnRatio)
  diffpnratio_lc <- abs(props1$pnRatio_lc - props2$pnRatio_lc)

  #--------------------------------------------------------------------------
  # size of the largest connected component
  difflcsize <- abs(props1$lcSize - props2$lcSize)

  #--------------------------------------------------------------------------
  # adjusted Rand index for measuring similarity between two clusterings
  clust1 <- props1$clust
  clust2 <- props2$clust
  if(isempty1 || isempty2){
    randInd <- NA
  } else{
    randInd <- c(value = randIndex(table(clust1, clust2), adjust = TRUE), 
                 pval = NA)
  }


  # significance test for Rand index
  if(testRand){
    randPerm <- numeric(nPermRand)
    for(i in 1:nPermRand){
      clust1.tmp <- gtools::permute(clust1)
      clust2.tmp <- gtools::permute(clust2)
      randPerm[i] <- randIndex(table(clust1.tmp, clust2.tmp), adjust = TRUE)
    }
    randMean <- mean(randPerm)
    randSD <- sd(randPerm)
    normRandPerm <- (randPerm - randMean) / randSD
    normRand <- (randInd[1] - randMean) / randSD
    randInd["pval"] <- (sum(normRandPerm >= abs(normRand)) + 
                          sum(normRandPerm <= -abs(normRand))) / nPermRand
  }

  #--------------------------------------------------------------------------
  # Jaccard Index

  jaccDeg <-    calc_jaccard(props1$topdeg, props2$topdeg, sigTest = testJacc)
  jaccBetw <-   calc_jaccard(props1$topbetw, props2$topbetw, sigTest = testJacc)
  jaccClose <-  calc_jaccard(props1$topclose, props2$topclose, sigTest = testJacc)
  jaccEigen <-  calc_jaccard(props1$topeigen, props2$topeigen, sigTest = testJacc)
  jaccHub <-    calc_jaccard(props1$hubs, props2$hubs, sigTest = testJacc)

  output <- list(jaccDeg = jaccDeg,
                 jaccBetw =jaccBetw,
                 jaccClose = jaccClose,
                 jaccEigen = jaccEigen,
                 jaccHub = jaccHub, 
                 randInd = randInd,
                 diffsGlobal = list(diffDiss = diffdiss,
                                    diffPath = diffpath,
                                    diffDensity = diffdensity,
                                    diffVertConnect = diffvertconnect,
                                    diffEdgeConnect = diffedgconnect,
                                    diffpnRatio = diffpnratio,
                                    diffClustCoef = diffclustcoef,
                                    diffModul = diffmod),
                 diffsGlobalLC = list(difflcSize = difflcsize,
                                      diffDiss = diffdiss_lc,
                                      diffPath = diffpath_lc,
                                      diffDensity = diffdensity_lc,
                                      diffVertConnect = diffvertconnect_lc,
                                      diffEdgeConnect = diffedgconnect_lc,
                                      diffpnRatio = diffpnratio_lc,
                                      diffClustCoef = diffclustcoef_lc,
                                      diffModul = diffmod_lc),
                 diffsCentr = list(diffDeg = diffdeg, 
                                   diffBetw = diffbetw,
                                   diffClose = diffclose, 
                                   diffEigen = diffeigen),
                 absDiffsCentr = list(absDiffDeg = absdiffdeg, 
                                      absDiffBetw = absdiffbetw,
                                      absDiffClose = absdiffclose, 
                                      absDiffEigen = absdiffeigen),
                 props = list(deg1 = props1$deg, 
                              deg2 = props2$deg,
                              betw1 = props1$betw, 
                              betw2 = props2$betw,
                              close1 = props1$close, 
                              close2 = props2$close,
                              eigen1 = props1$eigen, 
                              eigen2 = props2$eigen,
                              hubs1 = props1$hubs, 
                              hubs2 = props2$hubs2,
                              avDiss1 = avDiss1,
                              avDiss2 = avDiss2,
                              avPath1 = avPath1, 
                              avPath2 = avPath2,
                              density1 = props1$density, 
                              density2 = props2$density,
                              vertConnect1 = props1$vertconnect, 
                              vertConnect2 = props2$vertconnect,
                              edgeConnect1 = props1$edgeconnect, 
                              edgeConnect2 = props2$edgeconnect,
                              pnRatio1 = props1$pnRatio,
                              pnRatio2 = props2$pnRatio,
                              clustCoef1 = clustCoef1, 
                              clustCoef2 = clustCoef2,
                              modul1 = modul1, 
                              modul2 = modul2,
                              clust1 = clust1, 
                              clust2 = clust2),
                 propsLC = list(lcSize1 = props1$lcSize,
                                lcSize2 = props2$lcSize,
                                avDiss1 = avDiss1_lc,
                                avDiss2 = avDiss2_lc,
                                avPath1 = avPath1_lc, 
                                avPath2 = avPath2_lc,
                                density1 = props1$density_lc, 
                                density2 = props2$density_lc,
                                vertConnect1 = props1$vertconnect_lc, 
                                vertConnect2 = props2$vertconnect_lc,
                                edgeConnect1 = props1$edgeconnect_lc, 
                                edgeConnect2 = props2$edgeconnect_lc,
                                pnRatio1 = props1$pnRatio_lc,
                                pnRatio2 = props2$pnRatio_lc,
                                clustCoef1 = clustCoef1_lc, 
                                clustCoef2 = clustCoef2_lc,
                                modul1 = modul1_lc, 
                                modul2 = modul2_lc)
  )
  
  return(output)
}
