library(ergm)
library(network)
library(statnet)
library (texreg)
library(latticeExtra)
library(RColorBrewer)
library (texreg) 
library(sna)
#source()
#stop()

#Beliefs Ireland: G0101	G0102	G0103	G0104	G0105	G0106	G0107	G0108	G0109	G0110	G0111	G0112	G0113	G0114	G0115	G0117	G0401	G0403	G0404	G0405	G0406	G0407	G0408	G0409	G0410	G0411
seed <- 12345
set.seed (seed)

SciIre <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Ireland/SciIre.csv",header=TRUE)
coopIre <- read.csv(file="/Users/paulwagner/Desktop/RegressionData/coopR.csv",header=TRUE)
InfIre <- read.csv("/Users/paulwagner/Desktop/ERGMS/Ireland/influential.csv", header = TRUE, stringsAsFactors=FALSE)
OrgTypeIre <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Ireland/OrgType.csv", header = FALSE, stringsAsFactors=FALSE)

IreOrgNames <- as.vector(colnames(coopIre))

diag (coopIre) <- 0 
diag (InfIre) <- 0 
diag (SciIre) <- 0 

coopIre <- as.matrix(coopIre)
InfIre <- as.matrix(InfIre)
SciIre <- as.matrix(SciIre)
OrgTypeIre  <- as.matrix(OrgTypeIre)

rownames(coopIre) <- IreOrgNames
colnames(coopIre) <- IreOrgNames
rownames(SciIre) <- IreOrgNames
colnames(SciIre) <- IreOrgNames
rownames(InfIre) <- IreOrgNames
colnames(InfIre) <- IreOrgNames

#beliefs G0101	G0103	G0104	G0107	G0108	G0110	G0113	G0401	G0404	G0405	G0407	G0408	G0410
PBfs = read.csv("/Users/paulwagner/Desktop/ERGMS/Ireland/BeliefsIre.csv", header = TRUE, stringsAsFactors=FALSE)
PBfs <- PBfs[-c(10,14,23,37,41),-c(1,2,4,7,8,11,13,14,16,17,18,20,23,26,28)]
PBfs [PBfs  == 97] <- 3
npol.distR <- dist(PBfs, method = "manhattan")
PrefSimMatR<- max(npol.distR) - npol.distR
Pbs.Ire <-as.matrix(PrefSimMatR)

rownames(Pbs.Ire) <- IreOrgNames
colnames(Pbs.Ire) <- IreOrgNames

#Science Beliefs
SciBs.Ire = read.csv("/Users/paulwagner/Desktop/SciERGMs/SciQs/IreSciQs.csv", header = TRUE)
SciBs.Ire <- SciBs.Ire[-c(10,14,23,37,41),-c(1,2,6,7) ]
SciBs.Ire[SciBs.Ire == 97] <- 3
SciBs.dist.Ire <- dist(SciBs.Ire, method = "manhattan")
SciBsSimMatIre <- max(SciBs.dist.Ire) - SciBs.dist.Ire
SciBs.IreMat <-as.matrix(SciBsSimMatIre)

rownames(SciBs.IreMat) <- IreOrgNames
colnames(SciBs.IreMat) <- IreOrgNames

OrgTypeIre <- OrgTypeIre[,3]
OrgTypeIre <- as.character(OrgTypeIre)
OrgTypeIre <- as.vector(OrgTypeIre)

u <- sort(unique(OrgTypeIre))
nodecov <- match(OrgTypeIre,u)

"BUS CIV GOV NGO SCI"
#-------------
nw.SciIre <- network (SciIre) # create network object
set.vertex.attribute (nw.SciIre, "influence", degree (InfIre, cmode = "indegree"))
set.vertex.attribute (nw.SciIre, "OrgTypeIre", OrgTypeIre)

nw.coopIre <- network (coopIre)
set.vertex.attribute (nw.coopIre, "OrgTypeIre", OrgTypeIre)
set.vertex.attribute (nw.coopIre, "SciID", degree (SciIre, cmode = "indegree"))
set.vertex.attribute (nw.coopIre, "influence", degree (InfIre, cmode = "indegree"))


Ire.sci.sci <- matrix (0 , nrow = nrow (coopIre) , ncol = ncol (coopIre))
for (i in 1: nrow (Ire.sci.sci )) {
  for (j in 1: ncol (Ire.sci.sci)) {
    if (( OrgTypeIre [i] == "SCI" && OrgTypeIre [j] == "SCI") ||
        ( OrgTypeIre [i] == "SCI" && OrgTypeIre [j] == "SCI")) {
      Ire.sci.sci [i, j] <- 1
      Ire.sci.sci [j, i] <- 0
    }
  }
}

##collab network as dependent variable Aug 27th
#Same models as below, but building complexity
SciIre.m1 <- ergm(nw.coopIre ~ edges, 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciIre.m1)
plot(mcmc.diagnostics(SciIre.m1))
plot(gof(SciIre.m1))

SciIre.m2 <- ergm(nw.coopIre ~ edges + 
                    nodeifactor ("OrgTypeIre", base=-c(5)), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciIre.m2)
plot(mcmc.diagnostics(SciIre.m2))
plot(gof(SciIre.m2))

SciIre.m3 <- ergm(nw.coopIre ~ edges + 
                    nodeifactor ("OrgTypeIre", base=-c(5)) + 
                    edgecov(Pbs.Ire) + edgecov(SciBs.IreMat) + edgecov (SciIre) + nodeicov("SciID"), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciIre.m3)
plot(mcmc.diagnostics(SciIre.m3))
plot(gof(SciIre.m3))

SciIre.m4 <- ergm(nw.coopIre ~ edges + 
                    nodeifactor ("OrgTypeIre", base=-c(5)) + 
                    edgecov(Pbs.Ire) + edgecov(SciBs.IreMat) + edgecov (SciIre) + nodeicov("SciID") +
                    edgecov(InfIre), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciIre.m4)
plot(mcmc.diagnostics(SciIre.m4))
plot(gof(SciIre.m4))

SciIre.m5 <- ergm(nw.coopIre ~ edges + 
                     nodeifactor ("OrgTypeIre", base=-c(5)) +
                     edgecov(Pbs.Ire) + edgecov(SciBs.IreMat) + edgecov (SciIre) + nodeicov("SciID") + 
                     edgecov(InfIre) + 
                     mutual + twopath + gwodegree(0.3, fixed = TRUE ) +  
                     gwesp(1.0 , fixed = TRUE ), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 3000 , MCMC.interval = 3000))
summary(SciIre.m5)
plot(mcmc.diagnostics(SciIre.m5))
plot(gof(SciIre.m5))

par(mar=c(1,1,1,1))
par(mar=c(3.1,2.1,4.1,1.1))
    
texreg::screenreg(list(SciIre.m1,SciIre.m2,SciIre.m3,SciIre.m4,SciIre.m5),single.row = T)

texreg::htmlreg(list(SciIre.m1,SciIre.m2,SciIre.m3,SciIre.m4,SciIre.m5),"/Users/paulwagner/Desktop/SciERGMs/Models Aug 27/IreAug29",single.row = T)


#------------

texreg::screenreg(list(SciIre.m5,SciPort.m5,SciFin.m5,SciCzech.m5),single.row = T)
texreg::htmlreg(list(SciIre.m5,SciPort.m5,SciFin.m5,SciCzech.m5),"/Users/paulwagner/Desktop/SciERGMs/ModelsAug27",single.row = T)



#-----------

#https://www.sciencedirect.com/science/article/pii/S0264837716304963#bib0475
#science network as dependent variable
SciIre.m2 <- ergm(nw.SciIre ~ edges + mutual + 
                    nodeifactor ("OrgTypeIre", base=-c(5)) + nodeofactor ("OrgTypeIre", base=-c(5)) + #SCI receives
                    nodematch ("OrgTypeIre") + 
                    edgecov(PSMR) + edgecov (coopIre) + edgecov (Ire.gov.sci) +
                    nodeicov("influence") + 
                    gwdsp(0.1 , fixed = TRUE ) + gwesp(0.1 , fixed = TRUE ),
                  eval.loglik = TRUE , check.degeneracy = TRUE , control = control.ergm ( seed = seed , MCMC.samplesize = 2500 , MCMC.interval = 2500))
summary(SciIre.m2)
plot(mcmc.diagnostics(SciIre.m2))
plot(gof(SciIre.m2))

SciIre.m3 <- ergm(nw.PAIre ~ edges + 
                    nodeifactor ("OrgTypeIre") + nodeofactor ("OrgTypeIre") + #SCI receives
                    edgecov(PSMR) + edgecov (coopIre) + edgecov (SciIre) + 
                    mutual + gwdsp(0.1 , fixed = TRUE ) + gwesp(0.1 , fixed = TRUE ),
                  eval.loglik = TRUE , check.degeneracy = TRUE , control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciIre.m3)
plot(mcmc.diagnostics(SciIre.m2))
plot(gof(SciIre.m3))


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


ggnet2(nw.SciIre, alpha = 0.5, label = TRUE, label.size = 3, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.5),  node.size = 6, 
       color = "OrgTypeIre",
       color.legend = "OrgTypeIre", edge.color = "black", 
       palette = c("CIV" = "purple", "GOV" = "red", "BUS" = "blue", "SCI" = "yellow", "NGO" = "green"), 
       size = "indegree", size.min = 1, size.cut = 10, size.legend = "Indegree centrality",
       legend.position = "bottom")
ggsave("/Users/paulwagner/Desktop/ERGMs/SciERGMs/NetworkGraphs/Ire_graph.png", width = 29, height = 18,  device = NULL, dpi = 300)    

# Vertex Centrality

# Function for centrality stats

detach(package:sna)

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

nw.SciIre.igraph <- graph.adjacency(SciIre)
Results1Ire <- myCentrality(nw.SciIre.igraph)

# Rounding 

Results1Ire$Betweenness <- round(Results1Ire$Betweenness, digits =3)
Results1Ire$Closeness <- round(Results1Ire$Closeness, digits =3)
Results1Ire$Eigenvector <- round(Results1Ire$Eigenvector, digits =3)
Results1Ire$PageRank <- round(Results1Ire$PageRank, digits =3)
Results1Ire$IDegree <- round(Results1Ire$IDegree, digits =3)
Results1Ire$ODegree <- round(Results1Ire$ODegree, digits =3)

Results1Ire 
Results1Ire <- as.data.frame(Results1Ire)

write.csv2(Results1Ire, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/CentralityScores/centrality.scores.Ire.csv")

Results1Ire$Actor.type.Ire <- OrgTypeIre

## Calculate means per actor type
Results1Ire$ID <- NULL
Results1Ire$Actor.type.Ire <- as.character(Results1Ire$Actor.type.Ire)

mean_by_actorIre <- aggregate(Results1Ire, by=list(Results1Ire$Actor.type.Ire), FUN=mean, na.rm=TRUE )
mean_by_actorIre

mean_by_actor$Actor.typeIre <- NULL

mean_by_actorIre$Betweenness <- round(mean_by_actorIre$Betweenness, digits =3)
mean_by_actorIre$Closeness <- round(mean_by_actorIre$Closeness, digits =3)
mean_by_actorIre$Eigenvector <- round(mean_by_actorIre$Eigenvector, digits =3)
mean_by_actorIre$PageRank <- round(mean_by_actorIre$PageRank, digits =3)
mean_by_actorIre$IDegree <- round(mean_by_actorIre$IDegree, digits =3)
mean_by_actorIre$ODegree <- round(mean_by_actorIre$ODegree, digits =3)

mean_by_actorIre <- as.data.frame(mean_by_actorIre[c(1:5),c(1:7)])
mean_by_actorIre 
write.csv2(mean_by_actor, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/Ire.centrality.scores-mean-by-actor.csv")

#Stats---------
detach("package:igraph",unload=TRUE)
library(sna)
plot(cug.test(SciIre,centralization,cmode="dyad",FUN.arg=list(FUN=degree,cmode="indegree")))

den.Ire <- gden(SciIre)
SciIre_in  <- degree(SciIre, cmode="indegree")
SciIre_out <- degree(SciIre, cmode="outdegree")
mean.Ire <- mean(SciIre_in)
sd.Ire <- sd(SciIre_in)

Centra.Ire.In <- centralization(SciIre, degree, cmode="indegree")
Centra.Ire.out <- centralization(SciIre, degree, cmode="outdegree")
Centra.Ire.evcent <- centralization(SciIre, evcent) 
grecip.Ire <- grecip(SciIre) 
gtrans.Ire <- gtrans(SciIre)

efficiency.Ire <- efficiency(SciIre)
efficiency.krackhard.Ire <- hierarchy(SciIre, measure= "krackhardt")
hierarchy.Ire <-  hierarchy(SciIre, measure= "reciprocity")
lubness.Ire <- lubness(SciIre)

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

nw.infIre.igraph <- graph.adjacency(infIre)
Results.infIre <- myInfCentrality(nw.infIre.igraph)
#Results.infIre

# Rounding 

Results.infIre$Betweenness <- round(Results.infIre$Betweenness, digits =3)
Results.infIre$Closeness <- round(Results.infIre$Closeness, digits =3)
Results.infIre$Eigenvector <- round(Results.infIre$Eigenvector, digits =3)
Results.infIre$PageRank <- round(Results.infIre$PageRank, digits =3)
Results.infIre$IDegree <- round(Results.infIre$IDegree, digits =3)
Results.infIre$ODegree <- round(Results.infIre$ODegree, digits =3)

Results.infIre 
Results.infIre <- as.data.frame(Results.infIre)

Results.infIre$Actor.type.Ire <- OrgTypeIre

## Calculate means per actor type
Results.infIre$ID <- NULL
Results.infIre$Actor.type.Ire <- as.character(Results.infIre$Actor.type.Ire)

mean_inf_by_actorIre <- aggregate(Results.infIre, by=list(Results.infIre$Actor.type.Ire), FUN=mean, na.rm=TRUE )
mean_inf_by_actorIre

mean_inf_by_actorIre$Actor.typeIre <- NULL

mean_inf_by_actorIre$Betweenness <- round(mean_inf_by_actorIre$Betweenness, digits =3)
mean_inf_by_actorIre$Closeness <- round(mean_inf_by_actorIre$Closeness, digits =3)
mean_inf_by_actorIre$Eigenvector <- round(mean_inf_by_actorIre$Eigenvector, digits =3)
mean_inf_by_actorIre$PageRank <- round(mean_inf_by_actorIre$PageRank, digits =3)
mean_inf_by_actorIre$IDegree <- round(mean_inf_by_actorIre$IDegree, digits =3)
mean_inf_by_actorIre$ODegree <- round(mean_inf_by_actorIre$ODegree, digits =3)

mean_inf_by_actorIre <- as.data.frame(mean_inf_by_actorIre[c(1:5),c(1:7)])
mean_inf_by_actorIre 




#Collaboration  Stats------------------------------------
library(igraph)
coopIre <- read.csv(file="/Users/paulwagner/Desktop/RegressionData/coopR.csv",header=TRUE)
diag (coopIre) <- 0 
coopIre <- as.matrix(coopIre)

coopIre <- as.matrix(coopIre)

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

nw.coopIre.igraph <- graph.adjacency(coopIre)
Results.coopIre <- mycoopCentrality(nw.coopIre.igraph)
#Results.coopIre

# Rounding 

Results.coopIre$Betweenness <- round(Results.coopIre$Betweenness, digits =3)
Results.coopIre$Closeness <- round(Results.coopIre$Closeness, digits =3)
Results.coopIre$Eigenvector <- round(Results.coopIre$Eigenvector, digits =3)
Results.coopIre$PageRank <- round(Results.coopIre$PageRank, digits =3)
Results.coopIre$IDegree <- round(Results.coopIre$IDegree, digits =3)
Results.coopIre$ODegree <- round(Results.coopIre$ODegree, digits =3)

Results.coopIre 
Results.coopIre <- as.data.frame(Results.coopIre)

Results.coopIre$Actor.type.Ire <- OrgTypeIre

## Calculate means per actor type
Results.coopIre$ID <- NULL
Results.coopIre$Actor.type.Ire <- as.character(Results.coopIre$Actor.type.Ire)

mean_coop_by_actorIre <- aggregate(Results.coopIre, by=list(Results.coopIre$Actor.type.Ire), FUN=mean, na.rm=TRUE )
mean_coop_by_actorIre

mean_coop_by_actorIre$Actor.typeIre <- NULL

mean_coop_by_actorIre$Betweenness <- round(mean_coop_by_actorIre$Betweenness, digits =3)
mean_coop_by_actorIre$Closeness <- round(mean_coop_by_actorIre$Closeness, digits =3)
mean_coop_by_actorIre$Eigenvector <- round(mean_coop_by_actorIre$Eigenvector, digits =3)
mean_coop_by_actorIre$PageRank <- round(mean_coop_by_actorIre$PageRank, digits =3)
mean_coop_by_actorIre$IDegree <- round(mean_coop_by_actorIre$IDegree, digits =3)
mean_coop_by_actorIre$ODegree <- round(mean_coop_by_actorIre$ODegree, digits =3)

mean_coop_by_actorIre <- as.data.frame(mean_coop_by_actorIre[c(1:5),c(1:7)])
mean_coop_by_actorIre 

