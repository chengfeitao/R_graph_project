library("igraph")
##################problem1####################################################
edgelistFile <- "Y:/ee232e_project/facebook_combined.txt"

g <- read.graph(edgelistFile , format = "ncol" , directed=FALSE)


edgesNum <- ecount(g)
nodesNum <- vcount(g)

is.connected(g)
d <-diameter(g)
degreesVector <- degree(g)
h <- hist(degreesVector, breaks=seq(-0.5, by=1 , length.out=max(degreesVector)+2))
x=h$mids[1:max(degreeVector)+1]
y=h$counts[1:max(degreeVector)+1]
m<-nls(y~I(exp(1)^(a+b*x)),start=list(a=0,b=0))
summary(m)

#draw curve
s=seq(from=0,to=max(degreesVector)+2,by=1)
lines(s, predict(m, list(x = s)), col = "green") #change color

MSE<-(sum((residuals(m))^2))/max(degreesVector)
averageDegree<-mean(degreesVector)

########################problem2################################################

FirstNeighbor<-neighborhood(g,1,nodes=1)
FirstNeighbor<-FirstNeighbor[[1]]                           
nonFirstNeighbor<-which(!((1:vcount(g))%in%FirstNeighbor))
FirstGraph<-delete.vertices(g,nonFirstNeighbor)
ecount(FirstGraph)
vcount(FirstGraph)
#####################problem3#####################################################
CoreDegree = 0
coreNodes <- which(neighborhood.size(g, 1 , nodes=V(g)) > 200)
CN<-length(coreNodes)
for (i in 1:CN) {
  CoreDegree=CoreDegree+neighborhood.size(g,1,nodes=coreNodes[i])-1
}
AvgDegree=CoreDegree/CN;
core <- coreNodes[1]
subGraphNodes <- neighborhood(g , 1 , nodes=core)


subGraphNodes <- subGraphNodes[[1]]

length(subGraphNodes)

nonSubGraphNodes <- which( !( (1:vcount(g)) %in% subGraphNodes)  )
length(nonSubGraphNodes)
subGraph <- delete.vertices(g , nonSubGraphNodes)

V(g)$number <- 1:vcount(g)
subGraph <- delete.vertices(g, nonSubGraphNodes)
subGraphNodes <- sort(subGraphNodes)
#---------plotting personal network---------------
plot(subGraph)

plot(subGraph , vertex.size=2 , vertex.label=NA)
vertexSizeVector = rep(2,vcount(subGraph))


coreIndex<-which(subGraphNodes==core)
vertexSizeVector[coreIndex]=7

plot(subGraph , vertex.size=vertexSizeVector , vertex.label=NA , asp=9/16)
fg<-fastgreedy.community(subGraph)
color = rainbow(length(fg),alpha=0.5)
vcolor<-rep(0,vcount(subGraph))
for(i in 1:vcount(subGraph)) {
  vcolor[i] = color[fg$mem[i]]
}
plot.igraph(subGraph, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)

bt<-edge.betweenness.community(subGraph,directed = FALSE)
color2 = rainbow(length(bt),alpha=0.5)
vcolor2<-rep(0,vcount(subGraph))
for(i in 1:vcount(subGraph)) {
  vcolor2[i] = color2[bt$mem[i]]
}
plot.igraph(subGraph, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)

ifo<-infomap.community(subGraph)
color = rainbow(length(ifo),alpha=0.5)
#vcolor<-rep(0,vcount(subGraph))
for(i in 1:vcount(subGraph)) {
  vcolor[i] = color[ifo$mem[i]]
}
plot.igraph(subGraph, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)
######################problem4#####################################
subGraph_NoCore <- delete.vertices(subGraph, core)
vertexSizeVector = rep(2,vcount(subGraph_NoCore))
plot(subGraph_NoCore , vertex.size=vertexSizeVector , vertex.label=NA , asp=9/16)
fg1<-fastgreedy.community(subGraph_NoCore)
color = rainbow(length(fg1),alpha=0.5)
vcolor<-rep(0,vcount(subGraph_NoCore))
for(i in 1:vcount(subGraph_NoCore)) {
  vcolor[i] = color[fg1$mem[i]]
}
plot.igraph(subGraph_NoCore, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)


bt1<-edge.betweenness.community(subGraph_NoCore,directed = FALSE)
color = rainbow(length(bt1),alpha=0.5)
vcolor<-rep(0,vcount(subGraph_NoCore))
for(i in 1:vcount(subGraph_NoCore)) {
  vcolor[i] = color[bt1$mem[i]]
}
plot.igraph(subGraph_NoCore, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)


ifo1<-infomap.community(subGraph_NoCore)
color = rainbow(length(ifo1),alpha=0.5)
vcolor<-rep(0,vcount(subGraph_NoCore))
for(i in 1:vcount(subGraph_NoCore)) {
  vcolor[i] = color[ifo1$mem[i]]
}
plot.igraph(subGraph_NoCore, vertex.size =vertexSizeVector, vertex.label=NA, vertex.color=vcolor)
###############################problem5###############################
embeddedness = numeric(0)
dispersionTotal = numeric(0)
for (k in 1:CN){
  core <- coreNodes[k]
  subGraphNodes <- neighborhood(g , 1 , nodes=core)
  
  
  subGraphNodes <- subGraphNodes[[1]]
  
  
  nonSubGraphNodes <- which( !( (1:vcount(g)) %in% subGraphNodes)  )
  subGraph <- delete.vertices(g , nonSubGraphNodes)
  coreNew<- which.max(neighborhood.size(subGraph,1, nodes=V(subGraph)))
  subGraph_NoCore <- delete.vertices(subGraph, coreNew)  
  for (i in 1:vcount(subGraph_NoCore)){
    embeddednessNodes<- neighborhood(subGraph_NoCore , 1 , nodes=V(subGraph_NoCore)[i])
    embeddednessNodes<-embeddednessNodes[[1]]
    subGraph_NoCoreNoGoal<-delete.vertices(subGraph_NoCore,V(subGraph_NoCore)[i])
    embeddedness <- c(embeddedness,  neighborhood.size(subGraph_NoCore , 1 , nodes=V(subGraph_NoCore)[i])-1)
    
    #------------------------------------------------------------
    
    
    dispersionTemp <- 0;
    
    if(length(embeddednessNodes)==1){  #subsub cfij +h
      dispersionTemp <- 0
    } else {
      for(j in 2:length(embeddednessNodes)){ #start from second which is h
        neighbor1<-neighborhood(subGraph_NoCoreNoGoal,2,V(subGraph_NoCore)$name[embeddednessNodes[j]]) #index duiying
        neighbor1<-neighbor1[[1]]
        neighbor2<-which( V(subGraph_NoCoreNoGoal)$name[neighbor1] %in% V(subGraph_NoCore)$name[embeddednessNodes]) #c's neighbor in h neighbor
        dispersionTemp = dispersionTemp + length(embeddednessNodes)-length(neighbor2)-1 #total-connected-h
      }
    }
    dispersion <-c(dispersion, dispersionTemp/2) #every distance is caculated twice
  }
}

h <- hist(dispersion, breaks=seq(from=-0.5, by=100 , length.out=max(dispersion)/100+2),labels=TRUE)

    
    
    #--------------------------------------------------------------

h <- hist(embeddedness, breaks=seq(-0.5, by=1 , length.out=max(embeddedness)+2))
embeddedness



#------------------teng plot----------------------------------------
select<-c(2)

for(k in 1:3){
  #analyze the personal network
  core <- coreNodes[k]
  subGraphNodes=neighborhood(g,1,nodes=core)
  subGraphNodes <- subGraphNodes[[1]]
  nonSubGraphNodes <- which(!(1:vcount(g)) %in% subGraphNodes)
  subGraph <- delete.vertices(g , nonSubGraphNodes)
  
  #delete the core of the subGraph
  coreIndex <- which.max(neighborhood.size(subGraph, 1, nodes=V(subGraph)))
  subGraphNoCore <- delete.vertices(subGraph, coreIndex)
  
  embeddedness <- numeric()
  dispersion <- numeric()
  for(j in 1:vcount(subGraphNoCore)){
    #caculate the embeddedness
    embeddedness<-c(embeddedness,neighborhood.size(subGraphNoCore,1,V(subGraphNoCore)[j])-1)
    
    # caculate the dispersion between core nodes and j
    #set up a subgraph for the node you wanna caculate dispersion
    subGraphNoCoreNoGoal <-delete.vertices(subGraphNoCore,V(subGraphNoCore)[j])
    
    subSubGraphNodes <- neighborhood(subGraphNoCore,1,V(subGraphNoCore)[j])
    subSubGraphNodes <- subSubGraphNodes[[1]]
    
    dispersionTemp <- 0;
    
    if(length(subSubGraphNodes)==1){
      dispersionTemp <- 0
    } else {
      for(i in 2:length(subSubGraphNodes)){
        neighbor1<-neighborhood(subGraphNoCoreNoGoal,2,V(subGraphNoCore)$name[subSubGraphNodes[i]])
        neighbor1<-neighbor1[[1]]
        neighbor2<-which( V(subGraphNoCoreNoGoal)$name[neighbor1] %in% V(subGraphNoCore)$name[subSubGraphNodes])
        dispersionTemp = dispersionTemp + length(subSubGraphNodes)-length(neighbor2)-1
      }
    }
    dispersion <-c(dispersion, dispersionTemp/2) #every distance is caculated twice
  }
  
  #initialize the vertex size vector
  vertexSizeVector_Dispersion = rep(2,vcount(subGraph))
  vertexSizeVector_Embeddedness = rep(2,vcount(subGraph))
  vertexSizeVector_Ratio = rep(2,vcount(subGraph))
  
  #highlight the core
  vertexSizeVector_Dispersion[coreIndex] = 7
  vertexSizeVector_Embeddedness[coreIndex] = 7
  vertexSizeVector_Ratio[coreIndex] = 7
  
  #highlight the maximum dispersion point
  maxDispersionIndex <- which(V(subGraph)$name == V(subGraphNoCore)$name[which.max(dispersion)])
  vertexSizeVector_Dispersion[maxDispersionIndex] = 4
  EdgeDispersionIndex <- get.edge.ids(subGraph, c(coreIndex, maxDispersionIndex))
  
  edgeColorVector_Dispersion <- rep("lightgray", ecount(subGraph))
  edgeColorVector_Dispersion[EdgeDispersionIndex]<-"blue"
  edgeWidthVector_Dispersion <- rep("0.2", ecount(subGraph))
  edgeWidthVector_Dispersion[EdgeDispersionIndex] <- "12"
  
  #highlight the maximum embeddedness point
  maxEmbeddednessIndex <- which(V(subGraph)$name == V(subGraphNoCore)$name[which.max(embeddedness)])
  vertexSizeVector_Embeddedness[maxEmbeddednessIndex] = 4
  EdgeEmbeddednessIndex <- get.edge.ids(subGraph, c(coreIndex, maxEmbeddednessIndex))
  
  edgeColorVector_Embeddedness <- rep("lightgray", ecount(subGraph))
  edgeColorVector_Embeddedness[EdgeEmbeddednessIndex]<-"blue"
  edgeWidthVector_Embeddedness <- rep("0.2", ecount(subGraph))
  edgeWidthVector_Embeddedness[EdgeEmbeddednessIndex] <- "12"
  
  #highlight the maximum ratio point
  maxRatioIndex <- which(V(subGraph)$name == V(subGraphNoCore)$name[which.max(dispersion/embeddedness)])
  vertexSizeVector_Ratio[maxRatioIndex] = 4
  EdgeRatioIndex <- get.edge.ids(subGraph, c(coreIndex, maxRatioIndex))
  
  edgeColorVector_Ratio <- rep("lightgray", ecount(subGraph))
  edgeColorVector_Ratio[EdgeRatioIndex]<-"blue"
  edgeWidthVector_Ratio <- rep("0.2", ecount(subGraph))
  edgeWidthVector_Ratio[EdgeRatioIndex] <- "12"
  
  #fastgreedy algorithm
  communityFg <- fastgreedy.community(subGraph, merges=TRUE, modularity=TRUE,membership=TRUE, weights=NULL)
  colorFg = rainbow(length(communityFg),alpha=0.5)
  vcolorFg <- rep(0,vcount(subGraph))
  for(i in 1:vcount(subGraph)) {
    vcolorFg[i] = colorFg[communityFg$membership[i]]
  }
  
  plot.igraph(subGraph, vertex.size=vertexSizeVector_Dispersion, vertex.label=NA, vertex.color=vcolorFg, 
              edge.color=edgeColorVector_Dispersion, edge.width=edgeWidthVector_Dispersion)
  plot.igraph(subGraph, vertex.size=vertexSizeVector_Embeddedness, vertex.label=NA, vertex.color=vcolorFg, 
              edge.color=edgeColorVector_Embeddedness, edge.width=edgeWidthVector_Embeddedness)
  plot.igraph(subGraph, vertex.size=vertexSizeVector_Ratio, vertex.label=NA, vertex.color=vcolorFg, 
              edge.color=edgeColorVector_Ratio, edge.width=edgeWidthVector_Ratio)
}


###################problem6#########################################
g <- read.graph(edgelistFile , format = "ncol" , directed=FALSE)
coreNodes <- which(neighborhood.size(g, 1 , nodes=V(g)) > 200)


community_l <- numeric(0)
family_ratio <- numeric(0)  #family
family <- numeric(0)
College_ratio <-numeric(0)    #college
College <-numeric(0)
ratio <- numeric(0)

for (core_i in 1:length(coreNodes))
{
  #core_i <-20
  core <- coreNodes[core_i]
  subGraphNodes <- neighborhood(g , 1 , nodes=core)
  subGraphNodes <- subGraphNodes[[1]]
  nonSubGraphNodes <- which(!((1:vcount(g)) %in% subGraphNodes))
  subGraph <- delete.vertices(g , nonSubGraphNodes)
  
  # Use Infomap to find community whose size is large than 10
  ifo <- infomap.community(subGraph)
  community_i <- 1
  for (i in 1:length(ifo))
  {
    community_name <- V(subGraph)$name[which(ifo$membership==i)]
    if (length(community_name) >= 10){
      community_l[community_i] <- i
      community_i <- community_i +1
    }
  }
  
  # Find the High-School community
  for (j in 1:length(community_l))
  {
    non_communityNodes <- V(subGraph)$name[which(ifo$membership!=community_l[j])]
    subSubGraph <- delete.vertices(subGraph , non_communityNodes)
    dgr_mean <- mean(degree(subSubGraph))
    ratio[j] <- dgr_mean / vcount (subSubGraph)
  }
  ratio_max <-max(ratio)
  ratio_max_i <- which.max(ratio)
  ratio_min <- min(ratio)
  ratio_min_i <- which.min(ratio)
  
  family_ratio[core_i] <- ratio_max
  family[core_i] <-ratio_max_i
  College_ratio[core_i] <-ratio_min
  College[core_i] <- ratio_min_i
  
  
  color = rainbow(3,alpha=0.5)
  vcolor<-rep(0,vcount(subGraph))
  for(i in 1:vcount(subGraph)) {
    if(ifo$mem[i]==ratio_max_i){vcolor[i] = color[1]}
    else if(ifo$mem[i]==ratio_min_i){vcolor[i] = color[2]}
    else {vcolor[i] = color[3]}
    
  }
  #plot.igraph(subGraph, vertex.size =2, vertex.label=NA, vertex.color=vcolor)
  
  community_l <- numeric(0)
}
#plot.igraph(subGraph, vertex.size =2, vertex.label=NA, vertex.color=vcolor)



###################problem7#########################################
filesPath <- "Y:/ee232e_project/gplus/"

#egoNodeId <- "115625564993990145546" # 31 circles
egoNodeId <- "111091089527727420853" # 15 circles
egoNodeId <- "104105354262797387583" # 10 circles
egoNodeId <- "106382433884876652170" # 10b circles





edgelistFile <- paste(filesPath , egoNodeId  , ".edges" , sep="")
circlesFile <- paste(filesPath , egoNodeId , ".circles" , sep="")

g2Raw <- read.graph(edgelistFile , format = "ncol" , directed=TRUE)


#--------------add ego Node-----------------
nonEgoNodes = V(g2Raw)

g2 <- add.vertices(g2Raw,1,name=egoNodeId)
egoNodeIndex <- which(V(g2)$name==egoNodeId)

edgeAppendList <- c()
for (nodeIndex in 1:(vcount(g2)-1)) {
  edgeAppendList <- c(edgeAppendList , c(vcount(g2),nodeIndex))
}

g2 <- add.edges(g2,edgeAppendList)

g2U <- as.undirected(g2)

### reading circles

fileConnection <- file(circlesFile , open="r")
lines <- readLines(fileConnection)
circles <- list()

for (i in 1:length(lines)) {
  sp <- strsplit(lines[i],"\t")
  circles[[i]] <- sp[[1]][-1]
}

close(fileConnection)

#----------detect community walktrap-----------
wk<-walktrap.community(g2)
#-----------------infomap-----------------------
info<-infomap.community(g2)

# circle overlap with walktrap community
overlapPercentageDist<-numeric()
for(i in 1:length(circles)){
  overlapPercentage<-numeric()
  for(j in 1:length(wk)){
    communityNodes<-V(g2)$name[wk$mem==j]
    overlapPercentage<-c(overlapPercentage,length(intersect(communityNodes,circles[[i]]))/length(circles[[i]]))
  }
  overlapPercentageDist<-c(overlapPercentageDist,max(overlapPercentage))
}
hist(overlapPercentageDist)

# circle overlap with infomap community
overlapPercentageDist<-numeric()
for(i in 1:length(circles)){
  overlapPercentage<-numeric()
  for(j in 1:length(info)){
    communityNodes<-V(g2)$name[info$mem==j]
    overlapPercentage<-c(overlapPercentage,length(intersect(communityNodes,circles[[i]]))/length(circles[[i]]))
  }
  overlapPercentageDist<-c(overlapPercentageDist,max(overlapPercentage))
}
hist(overlapPercentageDist)

