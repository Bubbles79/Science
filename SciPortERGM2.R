detach("package:igraph", unload=TRUE)
library(ergm)
library(network)
library(statnet)
library (texreg)
library(latticeExtra)
library(RColorBrewer)
library (texreg) 
library(readxl)

#Beliefs Portland: G0101	G0102	G0103	G0104	G0105	G0106	G0107	G0108	G0109	G0110	G0111	G0112	G0113	G0114	G0115	G0117	G0401	G0403	G0404	G0405	G0406	G0407	G0408	G0409	G0410	G0411
seed <- 12345
set.seed (seed)

#import data
SciPort <- read.csv(file="/Users/paulwagner/Desktop/ERGMs/Portugal/SciPortugal.csv",header=TRUE)
coopPort <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Portugal/collabPortugal.csv",header=TRUE)
InfPort <- read.csv("/Users/paulwagner/Desktop/ERGMS/Portugal/InfPortugal.csv", header = TRUE)
OrgTypePort <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Portugal/OrgTypePortugal.csv", header = TRUE, stringsAsFactors=FALSE)

SciPort=network(SciPort,matrix.type="edgelist",directed=TRUE)
coopPort=network(coopPort,matrix.type="edgelist",directed=TRUE) 
InfPort=network(InfPort,matrix.type="edgelist",directed=TRUE) 

SciPort <- as.matrix(SciPort)
coopPort <- as.matrix(coopPort)
InfPort <- as.matrix(InfPort)

SciPort <-SciPort[-c(1,13,16,20,23,24,27,35,38,43,44,45,48,55,56,60,61,63,64,66,69,71,72,74,82,83,84),-c(1,13,16,20,23,24,27,35,38,43,44,45,48,55,56,60,61,63,64,66,69,71,72,74,82,83,84)]
coopPort <-coopPort[-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,68,70,71,73,81,82),-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,68,70,71,73,81,82)]
InfPort <-InfPort[-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,65,66,68,71,73,74,76,84,85),-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,65,66,68,71,73,74,76,84,85)]
OrgTypePort <-OrgTypePort[-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,65,66,68,71,73,74,76,84,85,86),]

PortOrgNames <- as.vector(colnames(coopPort))

#make sure diagonal is all 0s (no self ties)
diag (coopPort) <- 0 
diag (SciPort)  <- 0 
diag (InfPort) <- 0 

rownames(coopPort) <- PortOrgNames
colnames(coopPort) <- PortOrgNames
rownames(SciPort) <- PortOrgNames
colnames(SciPort) <- PortOrgNames
rownames(InfPort) <- PortOrgNames
colnames(InfPort) <- PortOrgNames

#extract org type data from excel sheet, convert to characters and then convert to a vector 
OrgTypePort <- OrgTypePort[,5]
OrgTypePort <- as.character(OrgTypePort)
OrgTypePort <- as.vector(OrgTypePort)

#Policy Beliefs
PBsPort = read.csv("/Users/paulwagner/Desktop/ERGMs/Portugal/PBsPortugal.csv", header = TRUE)
PBsPort <- PBsPort[-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,65,66,68,71,73,74,76,84,85,86),-c(1,2,4,7,8,11,13,14,16,17,19,20,23,26,28)]
PBsPort[PBsPort == 97] <- 3
Pol.dist.Port <- dist(PBsPort, method = "manhattan")
PrefSimMatPort <- max(Pol.dist.Port) - Pol.dist.Port
Pbs.Port <-as.matrix(PrefSimMatPort)

rownames(Pbs.Port) <- PortOrgNames
colnames(Pbs.Port) <- PortOrgNames

#Scientific Beliefs
SciBs.Port = read.csv("/Users/paulwagner/Desktop/SciERGMs/SciQs/PortSciQs.csv", header = TRUE)
SciBs.Port <- SciBs.Port[-c(1,7,9,17,21,24,25,28,33,36,39,44,45,46,49,55,57,61,62,64,65,66,68,71,73,74,76,84,86),-c(1,2,6,7) ]
SciBs.Port[SciBs.Port == 97] <- 3
SciBs.dist.Port <- dist(SciBs.Port, method = "manhattan")
SciBsSimMatPort <- max(SciBs.dist.Port) - SciBs.dist.Port
SciBs.PortMat <-as.matrix(SciBsSimMatPort)

rownames(SciBs.PortMat) <- PortOrgNames
colnames(SciBs.PortMat) <- PortOrgNames

u <- sort(unique(OrgTypePort))
nodecov <- match(OrgTypePort,u)

"BUS CIV GOV NGO SCI"
nw.SciPort <- network (SciPort) # create network object
set.vertex.attribute (nw.SciPort, "OrgTypePort", OrgTypePort)
set.vertex.attribute (nw.SciPort, "influence", degree (InfPort, cmode = "indegree"))

nw.coopPort <- network (coopPort)
set.vertex.attribute (nw.coopPort , "OrgTypePort", OrgTypePort)
set.vertex.attribute (nw.coopPort, "SciID", degree (SciPort, cmode = "indegree"))
set.vertex.attribute (nw.coopPort, "influence", degree (InfPort, cmode = "indegree"))

Port.sci.sci <- matrix (0 , nrow = nrow (coopPort) , ncol = ncol (coopPort))
for (i in 1: nrow ( Port.sci.sci )) {
  for (j in 1: ncol (  Port.sci.sci)) {
    if (( OrgTypePort [i] == "SCI" && OrgTypePort [j] == "SCI") ||
        ( OrgTypePort [i] == "SCI" && OrgTypePort [j] == "SCI")) {
      Port.sci.sci [i, j] <- 1
      Port.sci.sci [j, i] <- 0
    }
  }
} 

##collab network as dependent variable Aug 27th
#Same models as below, but building complexity
SciPort.m1 <- ergm(nw.coopPort ~ edges, 
                    eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciPort.m1)
plot(mcmc.diagnostics(SciPort.m1))
plot(gof(SciPort.m1))

SciPort.m2 <- ergm(nw.coopPort ~ edges + 
                     edgecov(Port.sci.sci) + nodeifactor ("OrgTypePort", base=-c(5)), 
                    eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciPort.m2)
plot(mcmc.diagnostics(SciPort.m2))
plot(gof(SciPort.m2))

SciPort.m3 <- ergm(nw.coopPort ~ edges + 
                     edgecov(Port.sci.sci) + nodeifactor ("OrgTypePort", base=-c(5)) +
                     edgecov(Pbs.Port) + edgecov(SciBs.PortMat) + edgecov (SciPort) + nodeicov("SciID"), 
                    eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciPort.m3)
plot(mcmc.diagnostics(SciPort.m3))
plot(gof(SciPort.m3))

SciPort.m4 <- ergm(nw.coopPort ~ edges + 
                     edgecov(Port.sci.sci) + nodeifactor ("OrgTypePort", base=-c(5)) +
                     edgecov(Pbs.Port) + edgecov(SciBs.PortMat) + edgecov (SciPort) + nodeicov("SciID") +
                     edgecov(InfPort), 
                    eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciPort.m4)
plot(mcmc.diagnostics(SciPort.m4))
plot(gof(SciPort.m4))

SciPort.m5 <- ergm(nw.coopPort ~ edges + 
                      edgecov(Port.sci.sci) + nodeifactor ("OrgTypePort", base=-c(5)) +
                      edgecov(Pbs.Port) + edgecov(SciBs.PortMat) + edgecov (SciPort) + nodeicov("SciID") +
                      edgecov(InfPort) +  
                      mutual +  twopath + gwodegree(0.5, fixed = TRUE ) +  
                      gwesp(0.7 , fixed = TRUE ), 
                    eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 3000 , MCMC.interval = 3000))
summary(SciPort.m5)
plot(mcmc.diagnostics(SciPort.m5))
plot(gof(SciPort.m5))


texreg::screenreg(list(SciPort.m1,SciPort.m2,SciPort.m3,SciPort.m4,SciPort.m5),single.row = T)

texreg::htmlreg(list(SciPort.m1,SciPort.m2,SciPort.m3,SciPort.m4,SciPort.m5),"/Users/paulwagner/Desktop/SciERGMs/Models Aug 27/PortAug29",single.row = T)


#--------------------------------------------------




#--------------------------------------------

####### Plot Network

library(sna)
# install.packages("GGally")
library(GGally)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("RColorBrewer")
library(RColorBrewer)
# install.packages("intergraph")
library(intergraph)
# install.packages("scales")
library(scales)


ggnet2(nw.SciPort, alpha = 0.5, label = TRUE, label.size = 3, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.5),  node.size = 6, 
       color = "OrgTypePort",
       color.legend = "OrgTypePort", edge.color = "black", 
       palette = c("CIV" = "purple", "GOV" = "red", "BUS" = "blue", "SCI" = "yellow", "NGO" = "green"), 
       size = "indegree", size.min = 1, size.cut = 10, size.legend = "Indegree centrality",
       legend.position = "bottom")
ggsave("/Users/paulwagner/Desktop/ERGMs/SciERGMs/NetworkGraphs/Port_graph.png", width = 29, height = 18,  device = NULL, dpi = 300)    

# Vertex Centrality

# Function for centrality stats

detach("package:statnet",unload=TRUE)

#install.packages("igraph")
library(igraph)

myCentrality <- function(net) {
  
  if (!is.igraph(net)) stop ("Input is not an igraph object")
  ideg <- degree (net, mode="in", loop=F, normalized = T)
  odeg <- degree (net, mode="out", loop=F, normalized = T)
  btw <- betweenness (net, normalized=TRUE)
  clo <- closeness (net, normalized=TRUE) # consider removing
  evc <- evcent(net)$vector
  pgr <- page.rank(net)$vector
  id <- V(net)
  
  ret <-data.frame (cbind(ID=id, IDegree=ideg, ODegree=odeg,
                          Betweenness=btw, Closeness=clo,
                          Eigenvector=evc, PageRank=pgr))
  return(ret)
  
}

# make network graph for igraph

nw.SciPort.igraph <- graph.adjacency(SciPort)
Results1Port <- myCentrality(nw.SciPort.igraph)

# Rounding 

Results1Port$Betweenness <- round(Results1Port$Betweenness, digits =3)
Results1Port$Closeness <- round(Results1Port$Closeness, digits =3)
Results1Port$Eigenvector <- round(Results1Port$Eigenvector, digits =3)
Results1Port$PageRank <- round(Results1Port$PageRank, digits =3)
Results1Port$IDegree <- round(Results1Port$IDegree, digits =3)
Results1Port$ODegree <- round(Results1Port$ODegree, digits =3)

Results1Port 
Results1Port <- as.data.frame(Results1Port)

write.csv2(Results1Port, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/CentralityScores/centrality.scores.Port.csv")

Results1Port$Actor.type.Port <- OrgTypePort

## Calculate means per actor type
Results1Port$ID <- NULL
Results1Port$Actor.type.Port <- as.character(Results1Port$Actor.type.Port)

mean_by_actorPort <- aggregate(Results1Port, by=list(Results1Port$Actor.type.Port), FUN=mean, na.rm=TRUE )
mean_by_actorPort

mean_by_actor$Actor.typePort <- NULL

mean_by_actorPort$Betweenness <- round(mean_by_actorPort$Betweenness, digits =3)
mean_by_actorPort$Closeness <- round(mean_by_actorPort$Closeness, digits =3)
mean_by_actorPort$Eigenvector <- round(mean_by_actorPort$Eigenvector, digits =3)
mean_by_actorPort$PageRank <- round(mean_by_actorPort$PageRank, digits =3)
mean_by_actorPort$IDegree <- round(mean_by_actorPort$IDegree, digits =3)
mean_by_actorPort$ODegree <- round(mean_by_actorPort$ODegree, digits =3)

mean_by_actorPort <- as.data.frame(mean_by_actorPort[c(1:5),c(1:7)])
mean_by_actorPort 
write.csv2(mean_by_actor, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/Port.centrality.scores-mean-by-actor.csv")

#Stats---------
detach("package:igraph",unload=TRUE)
library(sna)
plot(cug.test(SciPort,centralization,cmode="dyad",FUN.arg=list(FUN=degree,cmode="indegree")))

den.Port <- gden(SciPort)
SciPort_in  <- degree(SciPort, cmode="indegree")
SciPort_out <- degree(SciPort, cmode="outdegree")
mean.Port <- mean(SciPort_in)
sd.Port <- sd(SciPort_in)

Centra.Port.In <- centralization(SciPort, degree, cmode="indegree")
Centra.Port.out <- centralization(SciPort, degree, cmode="outdegree")
Centra.Port.evcent <- centralization(SciPort, evcent) 
grecip.Port <- grecip(SciPort) 
gtrans.Port <- gtrans(SciPort)

efficiency.Port <- efficiency(SciPort)
efficiency.krackhard.Port <- hierarchy(SciPort, measure= "krackhardt")
hierarchy.Port <-  hierarchy(SciPort, measure= "reciprocity")
lubness.Port <- lubness(SciPort)

#Influence Stats------------------------------------
library(igraph)

myInfCentrality <- function(net) {
  
  if (!is.igraph(net)) stop ("Input is not an igraph object")
  ideg <- degree (net, mode="in", loop=F, normalized = T)
  odeg <- degree (net, mode="out", loop=F, normalized = T)
  btw <- betweenness (net, normalized=TRUE)
  clo <- closeness (net, normalized=TRUE) # consider removing
  evc <- evcent(net)$vector
  pgr <- page.rank(net)$vector
  id <- V(net)
  
  ret <-data.frame (cbind(ID=id, IDegree=ideg, ODegree=odeg,
                          Betweenness=btw, Closeness=clo,
                          Eigenvector=evc, PageRank=pgr))
  return(ret)
  
}

# make network graph for igraph

nw.InfPort.igraph <- graph.adjacency(InfPort)
Results.InfPort <- myInfCentrality(nw.InfPort.igraph)
#Results.InfPort

# Rounding 

Results.InfPort$Betweenness <- round(Results.InfPort$Betweenness, digits =3)
Results.InfPort$Closeness <- round(Results.InfPort$Closeness, digits =3)
Results.InfPort$Eigenvector <- round(Results.InfPort$Eigenvector, digits =3)
Results.InfPort$PageRank <- round(Results.InfPort$PageRank, digits =3)
Results.InfPort$IDegree <- round(Results.InfPort$IDegree, digits =3)
Results.InfPort$ODegree <- round(Results.InfPort$ODegree, digits =3)

Results.InfPort 
Results.InfPort <- as.data.frame(Results.InfPort)

Results.InfPort$Actor.type.Port <- OrgTypePort

## Calculate means per actor type
Results.InfPort$ID <- NULL
Results.InfPort$Actor.type.Port <- as.character(Results.InfPort$Actor.type.Port)

mean_inf_by_actorPort <- aggregate(Results.InfPort, by=list(Results.InfPort$Actor.type.Port), FUN=mean, na.rm=TRUE )
mean_inf_by_actorPort

mean_inf_by_actorPort$Actor.typePort <- NULL

mean_inf_by_actorPort$Betweenness <- round(mean_inf_by_actorPort$Betweenness, digits =3)
mean_inf_by_actorPort$Closeness <- round(mean_inf_by_actorPort$Closeness, digits =3)
mean_inf_by_actorPort$Eigenvector <- round(mean_inf_by_actorPort$Eigenvector, digits =3)
mean_inf_by_actorPort$PageRank <- round(mean_inf_by_actorPort$PageRank, digits =3)
mean_inf_by_actorPort$IDegree <- round(mean_inf_by_actorPort$IDegree, digits =3)
mean_inf_by_actorPort$ODegree <- round(mean_inf_by_actorPort$ODegree, digits =3)

mean_inf_by_actorPort <- as.data.frame(mean_inf_by_actorPort[c(1:5),c(1:7)])
mean_inf_by_actorPort 


#Cooperation Stats------------------------------------
library(igraph)

coopPort <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Portugal/collabPortugal.csv",header=TRUE)
coopPort=network(coopPort,matrix.type="edgelist",directed=TRUE) 
coopPort <- as.matrix(coopPort)
coopPort <-coopPort[-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,68,70,71,73,81,82),-c(1,7,9,17,21,24,25,28,36,39,44,45,46,49,55,57,61,62,64,68,70,71,73,81,82)]
diag (coopPort) <- 0 

mycoopCentrality <- function(net) {
  
  if (!is.igraph(net)) stop ("Input is not an igraph object")
  ideg <- degree (net, mode="in", loop=F, normalized = F)
  odeg <- degree (net, mode="out", loop=F, normalized = F)
  btw <- betweenness (net, normalized=TRUE)
  clo <- closeness (net, normalized=TRUE) # consider removing
  evc <- evcent(net)$vector
  pgr <- page.rank(net)$vector
  id <- V(net)
  
  ret <-data.frame (cbind(ID=id, IDegree=ideg, ODegree=odeg,
                          Betweenness=btw, Closeness=clo,
                          Eigenvector=evc, PageRank=pgr))
  return(ret)
  
}

# make network graph for igraph

nw.coopPort.igraph <- graph.adjacency(coopPort)
Results.coopPort <- mycoopCentrality(nw.coopPort.igraph)
#Results.InfPort

# Rounding 

Results.coopPort$Betweenness <- round(Results.coopPort$Betweenness, digits =3)
Results.coopPort$Closeness <- round(Results.coopPort$Closeness, digits =3)
Results.coopPort$Eigenvector <- round(Results.coopPort$Eigenvector, digits =3)
Results.coopPort$PageRank <- round(Results.coopPort$PageRank, digits =3)
Results.coopPort$IDegree <- round(Results.coopPort$IDegree, digits =3)
Results.coopPort$ODegree <- round(Results.coopPort$ODegree, digits =3)

Results.coopPort 
Results.coopPort <- as.data.frame(Results.coopPort)

Results.coopPort$Actor.type.Port <- OrgTypePort

## Calculate means per actor type
Results.coopPort$ID <- NULL
Results.coopPort$Actor.type.Port <- as.character(Results.coopPort$Actor.type.Port)

mean_coop_by_actorPort <- aggregate(Results.coopPort, by=list(Results.coopPort$Actor.type.Port), FUN=mean, na.rm=TRUE )
mean_coop_by_actorPort

mean_coop_by_actorPort$Actor.typePort <- NULL

mean_coop_by_actorPort$Betweenness <- round(mean_coop_by_actorPort$Betweenness, digits =3)
mean_coop_by_actorPort$Closeness <- round(mean_coop_by_actorPort$Closeness, digits =3)
mean_coop_by_actorPort$Eigenvector <- round(mean_coop_by_actorPort$Eigenvector, digits =3)
mean_coop_by_actorPort$PageRank <- round(mean_coop_by_actorPort$PageRank, digits =3)
mean_coop_by_actorPort$IDegree <- round(mean_coop_by_actorPort$IDegree, digits =3)
mean_coop_by_actorPort$ODegree <- round(mean_coop_by_actorPort$ODegree, digits =3)

mean_coop_by_actorPort <- as.data.frame(mean_coop_by_actorPort[c(1:5),c(1:7)])
mean_coop_by_actorPort 
