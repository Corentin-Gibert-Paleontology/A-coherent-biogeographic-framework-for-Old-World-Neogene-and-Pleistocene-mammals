
######### Function used to compute HCPC and UPGMA Map and Tree on map

HCPC_UPGMA_Home  <- function(mds.Faune, Dist.Faune, i, id, sortie = TRUE, IllustrationMAP = TRUE)
{
  #HCPC analysis, production of HCPC Tree, on nMDS (argument 1)
  hcpc.Faune <- HCPC(as.data.frame(mds.Faune$points), nb.clust = i)
  if(sortie == TRUE){
    ExtractTitle <- paste("HCPC on nMDS ", id, " Nbre Clust ",max(as.numeric(hcpc.Faune$data.clust$clust)), "_K_", (length(hcpc.Faune$data.clust)-1), ".pdf", sep = "")
    dev.print(device = pdf, file = ExtractTitle)
  }
  if(IllustrationMAP == TRUE)
  {
    #x11()
    Color_HCPC <- IllustrationCluster(hcpc.Faune, etik = TRUE, cexEtik = 0.1)
    axis(side = 1, labels = "HCPC 'initial' clusters", at = 0)
  }
  if(sortie == TRUE){
    ExtractTitle <- paste("MAP HCPC ", id, " Nbre Clust ",max(as.numeric(hcpc.Faune$data.clust$clust)), "_K_", (length(hcpc.Faune$data.clust)-1), ".pdf", sep = "")
    dev.print(device = pdf, file = ExtractTitle, width = 10)
  }
  
  #UPGMA analysis on the imported (argument 2) Distance matrix
  UPGMA.Faune <- hclust(Dist.Faune, method = "ward.D2")
  
  #Label preparation
  FrameUPGMAonNMDS <- data.frame(Label = hcpc.Faune$call$t$tree$labels[hcpc.Faune$call$t$tree$order])
  FrameUPGMAonNMDS <- cbind(FrameUPGMAonNMDS, clust = rep(0, length(FrameUPGMAonNMDS$Label)))
  FrameUPGMA <- data.frame(Label = UPGMA.Faune$labels[UPGMA.Faune$order])
  FrameUPGMA <- cbind(FrameUPGMA, clust = rep(0, length(FrameUPGMA$Label)))
  for(i in 1:length(hcpc.Faune$data.clust$clust))
  {
    for(j in 1:length(FrameUPGMAonNMDS$Label))
    {
      if(dimnames(hcpc.Faune$data.clust)[[1]][i] == FrameUPGMAonNMDS$Label[j])
      {
        FrameUPGMAonNMDS[j,"clust"]  <- hcpc.Faune$data.clust$clust[i]
      }
    }
    for(j in 1:length(FrameUPGMA$Label))
    {
      if(dimnames(hcpc.Faune$data.clust)[[1]][i] == FrameUPGMA$Label[j])
      {
        FrameUPGMA[j,"clust"]  <- hcpc.Faune$data.clust$clust[i]
      }
    }
  }
  
  #HCPC colors on UPGMA tree
  T1 <- fviz_dend(UPGMA.Faune, cex = 0.3, k = max(as.numeric(hcpc.Faune$data.clust$clust)), label_cols = FrameUPGMA$clust)
  #UPGMA colors on UPGMA tree
  T1b <- fviz_dend(UPGMA.Faune, cex = 0.3, k = max(as.numeric(hcpc.Faune$data.clust$clust)), color_labels_by_k = TRUE)
  
  ### Illustration UPGMA on Map
  UPGMA_List <- list()
  data.clust <- data.frame(clust = T1b$plot_env$data$labels$col) #here are retained the cluster ID of each locality in UPGMA tree
  dimnames(data.clust)[[1]] <- UPGMA.Faune$labels[UPGMA.Faune$order]
  UPGMA_List <- list(data.clust = data.clust)
  LabelUPGMA <-   attr(Dist.Faune, "Labels")
  ClusterUPGMA <- c()
  Nclust <- c()
  for(i in 1:length(table(UPGMA_List$data.clust)))
  {
    if(i == 1)
    {
      Nclust <- c(Nclust, rep(i, sum(UPGMA_List$data.clust == UPGMA_List$data.clust[i,])))
    }
    if(i > 1)
    {
      Nclust <- c(Nclust, rep(i, sum(UPGMA_List$data.clust == UPGMA_List$data.clust[length(Nclust)+1,])))
    }
  }
  UPGMA_List$data.clust <- cbind(UPGMA_List$data.clust, Nclust)
  for(i in 1:length(LabelUPGMA))
  {
    for(j in 1:length(rownames(UPGMA_List$data.clust)))
    {
      if(LabelUPGMA[i] == rownames(UPGMA_List$data.clust)[j])
      {
        ClusterUPGMA[i] <- UPGMA_List$data.clust[j,"Nclust"]
      }
    }
  }
  
  if(IllustrationMAP == TRUE){
    #x11()
    ColUPGMA <- Map_UPGMA(UPGMA_List)
    axis(side = 1, labels = "UPGMA 'initial' clusters", at = 0)
    if(sortie == TRUE){
      ExtractTitle <- paste("MAP_UPGMA_TrueCol ", id, " nb.clust_", max(as.numeric(hcpc.Faune$data.clust$clust)), "_K_", (length(hcpc.Faune$data.clust)-1), ".pdf", sep = "")
      dev.print(device = pdf, file = ExtractTitle, width = 10)
    }
  }
  T2 <- fviz_dend(hcpc.Faune$call$t$tree, cex = 0.3, k = max(as.numeric(hcpc.Faune$data.clust$clust)), label_cols = FrameUPGMAonNMDS$clust)
  #HCPC Tree with real color
  if(IllustrationMAP == TRUE)
  {
    results <- list(hcpc.Faune = hcpc.Faune, UPGMA.Faune = UPGMA.Faune, UPGMAwithNMDScolor = T1, 
                    UPGMAtrueColor = T1b, UPGMAofHCPC = T2, UPGMA_List = UPGMA_List, Color_UPGMA = ColUPGMA, Color_HCPC = Color_HCPC,
                    ClusterUPGMA = ClusterUPGMA)
  }
  if(IllustrationMAP == FALSE)
  {
    results <- list(hcpc.Faune = hcpc.Faune, UPGMA.Faune = UPGMA.Faune, UPGMAwithNMDScolor = T1, 
                    UPGMAtrueColor = T1b, UPGMAofHCPC = T2, UPGMA_List = UPGMA_List, Color_UPGMA = ColUPGMA, ClusterUPGMA = ClusterUPGMA)  
  }
  return(results)
}


Map_UPGMA <- function(UPGMA_test)
{
  
  LAT <- strsplit(dimnames(UPGMA_test$data.clust)[[1]], "_")
  
  Lat.Monde <- c()
  Long.Monde <- c()
  for(i in 1:length(LAT))
  {
    Lat.Monde <- c(Lat.Monde, LAT[[i]][1])
    Long.Monde <- c(Long.Monde, LAT[[i]][2])
  }
  Lat.Monde <- as.numeric(Lat.Monde)
  Long.Monde <- as.numeric(Long.Monde)
  
  map("worldHires", col='gray85', fill=T, xlim=c(min(Long.Monde) -15, max(Long.Monde) +15), ylim=c(min(Lat.Monde) -15,max(Lat.Monde) +15))
  points(Long.Monde, Lat.Monde, col = UPGMA_test$data.clust$clust, pch = 16, cex = 0.75)
  CoordandCol_UPGMA <- data.frame(Lat = Lat.Monde, Long = Long.Monde, Color = UPGMA_test$data.clust$clust)
  return(CoordandCol_UPGMA)
}

######### Function used to compute HCPC on map

IllustrationCluster <- function(ListeFaunique, etik = FALSE, cexEtik = 0.5)
{
  library(maps)
  library(mapdata)
  
  Plot.EU1 <- data.frame()
  LAT <- strsplit(dimnames(ListeFaunique$data.clust)[[1]], "_")
  
  Lat.EU1 <- c()
  Long.EU1 <- c()
  for(i in 1:length(LAT))
  {
    Lat.EU1 <- c(Lat.EU1, LAT[[i]][1])
    Long.EU1 <- c(Long.EU1, LAT[[i]][2])
  }
  Lat.EU1 <- as.numeric(Lat.EU1)
  Long.EU1 <- as.numeric(Long.EU1)
  
  Color_C <- c()
  for(i in 1:length(ListeFaunique$data.clust$clust))
  {
    if(ListeFaunique$data.clust$clust[i] == 1){Color_C <- c(Color_C, "black")}
    if(ListeFaunique$data.clust$clust[i] == 2){Color_C <- c(Color_C, "brown1")}
    if(ListeFaunique$data.clust$clust[i] == 3){Color_C <- c(Color_C, "chartreuse3")}
    if(ListeFaunique$data.clust$clust[i] == 4){Color_C <- c(Color_C, "blue")}
    if(ListeFaunique$data.clust$clust[i] == 5){Color_C <- c(Color_C, "cyan")}
    if(ListeFaunique$data.clust$clust[i] == 6){Color_C <- c(Color_C, "pink")}
    if(ListeFaunique$data.clust$clust[i] == 7){Color_C <- c(Color_C, "gold")}
    if(ListeFaunique$data.clust$clust[i] == 8){Color_C <- c(Color_C, "aquamarine2")}
    if(ListeFaunique$data.clust$clust[i] == 9){Color_C <- c(Color_C, "bisque2")}
    if(ListeFaunique$data.clust$clust[i] == 10){Color_C <- c(Color_C, "blue2")}
    if(ListeFaunique$data.clust$clust[i] == 11){Color_C <- c(Color_C, "cornflowerblue")}
    if(ListeFaunique$data.clust$clust[i] == 12){Color_C <- c(Color_C, "darkgoldenrod1")}
    if(ListeFaunique$data.clust$clust[i] == 13){Color_C <- c(Color_C, "brown4")}
    if(ListeFaunique$data.clust$clust[i] == 14){Color_C <- c(Color_C, "azure3")}
    if(ListeFaunique$data.clust$clust[i] == 15){Color_C <- c(Color_C, "chocolate")}
    if(ListeFaunique$data.clust$clust[i] == 16){Color_C <- c(Color_C, "coral")}
    if(ListeFaunique$data.clust$clust[i] == 17){Color_C <- c(Color_C, "cornsilk2")}
    if(ListeFaunique$data.clust$clust[i] == 18){Color_C <- c(Color_C, "white")}
    if(ListeFaunique$data.clust$clust[i] == 19){Color_C <- c(Color_C, "blueviolet")}
    if(ListeFaunique$data.clust$clust[i] == 20){Color_C <- c(Color_C, "darkorchid")}
  }
  PlotEUROPE1 <- data.frame(Lat.EU1, Long.EU1, Color_C)
  PlotEUROPE1$Color_C <- as.character(PlotEUROPE1$Color_C)
  plot(PlotEUROPE1$Long.EU1, PlotEUROPE1$Lat.EU1, col = PlotEUROPE1$Color_C, ylim = c(), xlim = c(), cex = 0.5, pch = 16)
  if(etik == TRUE)
  {
    text(PlotEUROPE1$Long.EU1, PlotEUROPE1$Lat.EU1, labels=dimnames(ListeFaunique$data.clust)[[1]], cex = cexEtik, pos = 3)
  }
  map("worldHires", col='gray85', fill=T, xlim=c(min(Long.EU1)-15,max(Long.EU1)+15), ylim=c(min(Lat.EU1)-15,max(Lat.EU1)+15))
  points(PlotEUROPE1$Long.EU1, PlotEUROPE1$Lat.EU1, col = PlotEUROPE1$Color_C, pch = 16, cex = 0.75)
  Sortie_Illustration <- data.frame(Lat = PlotEUROPE1$Lat.EU1, Long = PlotEUROPE1$Long.EU1, Color = PlotEUROPE1$Color_C)
  return(Sortie_Illustration)
}


######### Function used to transform faunal list to occurrence matrix

Presence_Absence_matrix <- function(ListeFaunique, type = "Genus", singletons = TRUE, min5 = FALSE){
  library(stringr)
  
  if(type == "Species")
  {
    ch <- str_detect(ListeFaunique$SPECIES, "_")
    for(i in 1:length(ch))
    {
      if(ch[i] == TRUE)
      {
        TempSplit <- str_split(ListeFaunique$SPECIES[i], "_")
        if(length(TempSplit[[1]]) < 3){
          ListeFaunique$SPECIES[i] <- TempSplit[[1]][2]}
        if(length(TempSplit[[1]]) > 2){
          ListeFaunique$SPECIES[i] <- TempSplit[[1]][length(TempSplit[[1]])]}
      }
    }
    Species <- subset(ListeFaunique, ListeFaunique$SPECIES != "indet." & ListeFaunique$SPECIES != "_" & ListeFaunique$SPECIES != "sp." & ListeFaunique$SPECIES != "ssp.") #ajouter _sp
    SpeciesConcat <- paste(Species$GENUS, Species$SPECIES, sep = "_") 
    AllGenera <- c(as.character(unique(SpeciesConcat)))
    AllLocality <- c(as.character(unique(Species$NAME))) 
  }
  else if(type == "Genus")
  {
    AllGenera <- c(as.character(unique(ListeFaunique$GENUS))) 
    AllLocality <- c(as.character(unique(ListeFaunique$NAME))) 
  }
  matrix.PreAbs <- matrix(ncol=length(AllGenera), nrow=length(AllLocality), data = 0)
  dimnames(matrix.PreAbs) <- list(AllLocality, AllGenera) 
  
  for(i in 1:length(matrix.PreAbs[1,]))  
  {
    if(type == "Species")
    {
      LocalityOfGenera <- subset(ListeFaunique, AllGenera[i] == paste(ListeFaunique$GENUS, ListeFaunique$SPECIES, sep = "_"))
    }
    else if(type == "Genus")
    {
      LocalityOfGenera <- subset(ListeFaunique, AllGenera[i] == ListeFaunique$GENUS) 
    }
    vecTemp <- as.character(LocalityOfGenera$NAME)
    
    for(j in 1:length(matrix.PreAbs[,1]))                 
    {
      for(k in 1:length(vecTemp))                      
      {
        if(vecTemp[k] == AllLocality[j]){matrix.PreAbs[j,i] <- 1}  
      }
    }
  }
  
  if(singletons == FALSE)
  {
    tempMat <- apply(matrix.PreAbs, 2, sum)
    want <- which(tempMat == 1)
    tempMat <- matrix.PreAbs[,-want]
    matrix.PreAbs <- tempMat
    tempMat <- apply(matrix.PreAbs, 1, sum)
    want <- which(tempMat == 0)
    
    if(length(want) != 0)
    {
      tempMat <- matrix.PreAbs[-want,]
      matrix.PreAbs <- tempMat
    }
    tempMat <- apply(matrix.PreAbs, 2, sum)
    want <- which(tempMat == 0)
    if(length(want) != 0)
    {
      tempMat <- matrix.PreAbs[,-want]
      matrix.PreAbs <- tempMat
    }
  }
  
  if(min5 != FALSE)
  {
    tempMat <- apply(matrix.PreAbs, 1, sum)
    want <- which(tempMat < min5)
    if(length(want) > 0){
      tempMat <- matrix.PreAbs[-want,]
      matrix.PreAbs <- tempMat
      tempMat <- apply(matrix.PreAbs, 2, sum)
      want <- which(tempMat == 0)
      if(length(want) != 0)
      {
        tempMat <- matrix.PreAbs[,-want]
        matrix.PreAbs <- tempMat
      }
    }
  }
  
  
  return(matrix.PreAbs)
}


#Matrix for Polycohort, MN in columns, taxa in rows

Polycohorte_matrix <- function(ListeFaunique, ListeMN, NomMN){
  
  AllGenera <- c(as.character(unique(ListeFaunique$GENUS))) #Taxa name vector
  AllMN <- NomMN #vector of MN/biozones boundaries ages/datum
  matrix.PreAbs <- matrix(ncol=length(NomMN), nrow=length(AllGenera), data = 0)#Empty matrix
  dimnames(matrix.PreAbs) <- list(AllGenera, AllMN) #Copy/paste names in new matrix
  
  for(i in 1:length(matrix.PreAbs[,1])) #Looking for taxa
  {
    LocalityOfGenera <- subset(ListeFaunique, AllGenera[i] == ListeFaunique$GENUS)
    vecTempMAX <- max(LocalityOfGenera$MAX_AGE)
    vecTempMIN <- min(LocalityOfGenera$MIN_AGE)
    cibleMN <- c()
    ## Temporal occurrence are counted to fill empty matrix ##
    for(j in 1:length(matrix.PreAbs[1,]))                 #Loop through MN/biozones
    {
      if(ListeMN[j] > vecTempMIN & ListeMN[j] < vecTempMAX & ListeMN[j+1] >= vecTempMAX
         | ListeMN[j] <= vecTempMIN & ListeMN[j+1] > vecTempMIN & ListeMN[j+1] < vecTempMAX
         | ListeMN[j] > vecTempMIN & ListeMN[j+1] < vecTempMAX
         | ListeMN[j] <= vecTempMIN & ListeMN[j+1] >= vecTempMAX)
      {cibleMN <- c(cibleMN, j)}
    }
    matrix.PreAbs[i,cibleMN] <- 1
  }
  return(matrix.PreAbs)
}

