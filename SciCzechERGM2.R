detach("package:igraph", unload=TRUE)
library(ergm)
library(network)
library(statnet)
library (texreg)
library(latticeExtra)
library(RColorBrewer)
library (texreg) 

#Beliefs Czechland: G0101	G0102	G0103	G0104	G0105	G0106	G0107	G0108	G0109	G0110	G0111	G0112	G0113	G0114	G0115	G0117	G0401	G0403	G0404	G0405	G0406	G0407	G0408	G0409	G0410	G0411
seed <- 12345
set.seed (seed)

#imCzech data
'coopCzech <- read.table(file="/Users/paulwagner/Desktop/ERGMS/Czech Republic/Czech_Support.txt",header=FALSE)
SciCzech  <- read.table(file= "/Users/paulwagner/Desktop/ERGMs/Czech Republic/Czech_ScSc.txt",header=FALSE)
infCzech  <- read.table(file= "/Users/paulwagner/Desktop/ERGMS/Czech Republic/Czech_Inf.txt", header = FALSE)
'
OrgTypeCzech <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Czech Republic/OrgTypeCzech.csv", header = FALSE)

#Import Edgelists from Aasa Dataset
coopCzech <- read.csv(file="/Users/paulwagner/Desktop/ERGMs/Czech Republic/CoopCzechEdgeList.csv",header=TRUE)
SciCzech <- read.csv(file= "/Users/paulwagner/Desktop/ERGMS/Czech Republic/SciCzechEdgeList.csv",header=TRUE)
infCzech <- read.csv(file= "/Users/paulwagner/Desktop/ERGMS/Czech Republic/InfCzechEdgeList.csv", header = TRUE)

coopCzech=network(coopCzech,matrix.type="edgelist",directed=TRUE) 
SciCzech=network(SciCzech,matrix.type="edgelist",directed=TRUE)
infCzech=network(infCzech,matrix.type="edgelist",directed=TRUE) 

coopCzech <- as.matrix(coopCzech)
SciCzech <- as.matrix(SciCzech)
infCzech <- as.matrix(infCzech)

#make sure diagonal is all 0s (no self ties)
diag (coopCzech) <- 0 
diag (SciCzech)  <- 0 
diag (infCzech ) <- 0 

#remove non-respondents in OrgType data
coopCzech <- coopCzech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157),-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157)]
SciCzech <- SciCzech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157),-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157)]
infCzech <- infCzech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157),-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157)]
OrgTypeCzech  <- OrgTypeCzech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132), ]

CzechOrgNames <- as.vector(colnames(coopCzech))

#convert networks to matrices
infCzech <- as.matrix(infCzech)
coopCzech <- as.matrix(coopCzech)
SciCzech  <- as.matrix(SciCzech)

rownames(infCzech) <- CzechOrgNames
colnames(infCzech) <- CzechOrgNames
rownames(coopCzech) <- CzechOrgNames
colnames(coopCzech) <- CzechOrgNames
rownames(SciCzech) <- CzechOrgNames
colnames(SciCzech) <- CzechOrgNames

#extract org type data from excel sheet, convert to characters and then convert to a vector 
OrgTypeCzech <- as.matrix(OrgTypeCzech)
OrgTypeCzech <- OrgTypeCzech[,5]
OrgTypeCzech <- as.character(OrgTypeCzech)
OrgTypeCzech <- as.vector(OrgTypeCzech)

PolBs.Czech = read.csv("/Users/paulwagner/Desktop/ERGMs/Czech Republic/CzechBeliefsJune20.csv", header = TRUE)
PBsCzech <- PolBs.Czech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132),-c(8,9,10,12) ]
PBsCzech[PBsCzech == 97] <- 3
Pol.dist.Czech <- dist(PBsCzech, method = "manhattan")
PrefSimMatCzech <- max(Pol.dist.Czech) - Pol.dist.Czech
Pbs.Czech <-as.matrix(PrefSimMatCzech)

rownames(Pbs.Czech) <- CzechOrgNames
colnames(Pbs.Czech) <- CzechOrgNames

#Science Beliefs
SciBs.Czech = read.csv("/Users/paulwagner/Desktop/SciERGMs/SciQs/CzechSciQs.csv", header = TRUE)
SciBs.Czech <- SciBs.Czech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132),-c(1,2,6,7) ]
SciBs.Czech[SciBs.Czech == 97] <- 3
SciBs.dist.Czech <- dist(SciBs.Czech, method = "manhattan")
SciBsSimMatCzech <- max(SciBs.dist.Czech) - SciBs.dist.Czech
SciBs.CzechMat <-as.matrix(SciBsSimMatCzech)

rownames(SciBs.CzechMat) <- CzechOrgNames
colnames(SciBs.CzechMat) <- CzechOrgNames

#maybe I should remove actor 20 and actor 78. Both have zero beliefs an no outgoing ties. 

u <- sort(unique(OrgTypeCzech))
nodecov <- match(OrgTypeCzech,u)

"BUS CIV GOV NGO SCI"
nw.SciCzech <- network (SciCzech) # create network object
set.vertex.attribute (nw.SciCzech, "OrgTypeCzech", OrgTypeCzech)
set.vertex.attribute (nw.SciCzech, "influence", degree (infCzech, cmode = "indegree"))

nw.coopCzech <- network (coopCzech)
set.vertex.attribute (nw.coopCzech, "OrgTypeCzech", OrgTypeCzech)
set.vertex.attribute (nw.coopCzech, "SciID", degree (SciCzech, cmode = "indegree"))
set.vertex.attribute (nw.coopCzech, "influence", degree (infCzech, cmode = "indegree"))

Cze.sci.sci <- matrix (0 , nrow = nrow (coopCzech) , ncol = ncol (coopCzech))
for (i in 1: nrow (Cze.sci.sci )) {
  for (j in 1: ncol ( Cze.sci.sci)) {
    if (( OrgTypeCzech [i] == "SCI" && OrgTypeCzech [j] == "SCI") ||
        ( OrgTypeCzech [i] == "SCI" && OrgTypeCzech [j] == "SCI")) {
      Cze.sci.sci [i, j] <- 1
      Cze.sci.sci [j, i] <- 0
    }
  }
}

##collab network as dependent variable Aug 27th
#Same models as below, but building complexity
SciCzech.m1 <- ergm(nw.coopCzech ~ edges, 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciCzech.m1)
plot(mcmc.diagnostics(SciCzech.m1))
plot(gof(SciCzech.m1))

SciCzech.m2 <- ergm(nw.coopCzech ~ edges + 
                      edgecov(Cze.sci.sci) + nodeifactor ("OrgTypeCzech", base=-c(5)), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciCzech.m2)
plot(mcmc.diagnostics(SciCzech.m2))
plot(gof(SciCzech.m2))

SciCzech.m3 <- ergm(nw.coopCzech ~ edges + 
                      edgecov(Cze.sci.sci) + nodeifactor ("OrgTypeCzech", base=-c(5)) +
                      edgecov(Pbs.Czech) + edgecov(SciBs.CzechMat) + edgecov (SciCzech) + nodeicov("SciID"), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciCzech.m3)
plot(mcmc.diagnostics(SciCzech.m3))
plot(gof(SciCzech.m3))

SciCzech.m4 <- ergm(nw.coopCzech ~ edges + 
                      edgecov(Cze.sci.sci) + nodeifactor ("OrgTypeCzech", base=-c(5)) +
                      edgecov(Pbs.Czech) + edgecov(SciBs.CzechMat) + edgecov (SciCzech) + nodeicov("SciID") +
                      edgecov(infCzech),
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciCzech.m4)
plot(mcmc.diagnostics(SciCzech.m4))
plot(gof(SciCzech.m4))

SciCzech.m5 <- ergm(nw.coopCzech ~ edges + 
                      edgecov(Cze.sci.sci) + nodeifactor ("OrgTypeCzech", base=-c(5)) +
                      edgecov(Pbs.Czech) + edgecov(SciBs.CzechMat) + edgecov (SciCzech) + nodeicov("SciID") +
                      edgecov(infCzech) +   
                      mutual +  twopath + gwodegree(0.1, fixed = TRUE ) +  
                      gwesp(1.0 , fixed = TRUE ), 
                     eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 3000 , MCMC.interval = 3000))
summary(SciCzech.m5)
plot(mcmc.diagnostics(SciCzech.m5))
plot(gof(SciCzech.m5))


texreg::screenreg(list(SciCzech.m1,SciCzech.m2,SciCzech.m3,SciCzech.m4,SciCzech.m5),single.row = T)

texreg::htmlreg(list(SciCzech.m1,SciCzech.m2,SciCzech.m3,SciCzech.m4,SciCzech.m5),"/Users/paulwagner/Desktop/SciERGMs/Models Aug 27/CzechAug29",single.row = T)
#--------------------------------------------



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


ggnet2(nw.SciCzech, alpha = 0.5, label = TRUE, label.size = 3, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.5),  node.size = 6, 
       color = "OrgTypeCzech",
       color.legend = "OrgTypeCzech", edge.color = "black", 
       palette = c("CIV" = "purple", "GOV" = "red", "BUS" = "blue", "SCI" = "yellow", "NGO" = "green"), 
       size = "indegree", size.min = 1, size.cut = 10, size.legend = "Indegree centrality",
       legend.position = "bottom")
ggsave("/Users/paulwagner/Desktop/ERGMs/SciERGMs/NetworkGraphs/Czech_graph.png", width = 29, height = 18,  device = NULL, dpi = 300)    

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

nw.SciCzech.igraph <- graph.adjacency(SciCzech)
Results1Czech <- myCentrality(nw.SciCzech.igraph)

# Rounding 

Results1Czech$Betweenness <- round(Results1Czech$Betweenness, digits =3)
Results1Czech$Closeness <- round(Results1Czech$Closeness, digits =3)
Results1Czech$Eigenvector <- round(Results1Czech$Eigenvector, digits =3)
Results1Czech$PageRank <- round(Results1Czech$PageRank, digits =3)
Results1Czech$IDegree <- round(Results1Czech$IDegree, digits =3)
Results1Czech$ODegree <- round(Results1Czech$ODegree, digits =3)

Results1Czech 
Results1Czech <- as.data.frame(Results1Czech)

#write.csv2(Results1Czech, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/CentralityScores/centrality.scores.Czech.csv")

Results1Czech$Actor.type.Czech <- OrgTypeCzech

## Calculate means per actor type
Results1Czech$ID <- NULL
Results1Czech$Actor.type.Czech <- as.character(Results1Czech$Actor.type.Czech)

mean_by_actorCzech <- aggregate(Results1Czech, by=list(Results1Czech$Actor.type.Czech), FUN=mean, na.rm=TRUE )
mean_by_actorCzech

mean_by_actor$Actor.typeCzech <- NULL

mean_by_actorCzech$Betweenness <- round(mean_by_actorCzech$Betweenness, digits =3)
mean_by_actorCzech$Closeness <- round(mean_by_actorCzech$Closeness, digits =3)
mean_by_actorCzech$Eigenvector <- round(mean_by_actorCzech$Eigenvector, digits =3)
mean_by_actorCzech$PageRank <- round(mean_by_actorCzech$PageRank, digits =3)
mean_by_actorCzech$IDegree <- round(mean_by_actorCzech$IDegree, digits =3)
mean_by_actorCzech$ODegree <- round(mean_by_actorCzech$ODegree, digits =3)

mean_by_actorCzech <- as.data.frame(mean_by_actorCzech[c(1:5),c(1:7)])
mean_by_actorCzech
#write.csv2(mean_by_actor, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/Czech.centrality.scores-mean-by-actor.csv")

#Stats---------
detach("package:igraph",unload=TRUE)
library(sna)
plot(cug.test(SciCzech,centralization,cmode="dyad",FUN.arg=list(FUN=degree,cmode="indegree")))

den.Czech <- gden(SciCzech)
SciCzech_in  <- degree(SciCzech, cmode="indegree")
SciCzech_out <- degree(SciCzech, cmode="outdegree")
mean.Czech <- mean(SciCzech_in)
sd.Czech <- sd(SciCzech_in)

Centra.Czech.In <- centralization(SciCzech, degree, cmode="indegree")
Centra.Czech.out <- centralization(SciCzech, degree, cmode="outdegree")
Centra.Czech.evcent <- centralization(SciCzech, evcent) 
grecip.Czech <- grecip(SciCzech) 
gtrans.Czech <- gtrans(SciCzech)

efficiency.Czech <- efficiency(SciCzech)
efficiency.krackhard.Czech <- hierarchy(SciCzech, measure= "krackhardt")
hierarchy.Czech <-  hierarchy(SciCzech, measure= "reciprocity")
lubness.Czech <- lubness(SciCzech)


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

nw.infCzech.igraph <- graph.adjacency(infCzech)
Results.infCzech <- myInfCentrality(nw.infCzech.igraph)
#Results.infCzech

# Rounding 

Results.infCzech$Betweenness <- round(Results.infCzech$Betweenness, digits =3)
Results.infCzech$Closeness <- round(Results.infCzech$Closeness, digits =3)
Results.infCzech$Eigenvector <- round(Results.infCzech$Eigenvector, digits =3)
Results.infCzech$PageRank <- round(Results.infCzech$PageRank, digits =3)
Results.infCzech$IDegree <- round(Results.infCzech$IDegree, digits =3)
Results.infCzech$ODegree <- round(Results.infCzech$ODegree, digits =3)

Results.infCzech 
Results.infCzech <- as.data.frame(Results.infCzech)

Results.infCzech$Actor.type.Czech <- OrgTypeCzech

## Calculate means per actor type
Results.infCzech$ID <- NULL
Results.infCzech$Actor.type.Czech <- as.character(Results.infCzech$Actor.type.Czech)

mean_inf_by_actorCzech <- aggregate(Results.infCzech, by=list(Results.infCzech$Actor.type.Czech), FUN=mean, na.rm=TRUE )
mean_inf_by_actorCzech

mean_inf_by_actorCzech$Actor.typeCzech <- NULL

mean_inf_by_actorCzech$Betweenness <- round(mean_inf_by_actorCzech$Betweenness, digits =3)
mean_inf_by_actorCzech$Closeness <- round(mean_inf_by_actorCzech$Closeness, digits =3)
mean_inf_by_actorCzech$Eigenvector <- round(mean_inf_by_actorCzech$Eigenvector, digits =3)
mean_inf_by_actorCzech$PageRank <- round(mean_inf_by_actorCzech$PageRank, digits =3)
mean_inf_by_actorCzech$IDegree <- round(mean_inf_by_actorCzech$IDegree, digits =3)
mean_inf_by_actorCzech$ODegree <- round(mean_inf_by_actorCzech$ODegree, digits =3)

mean_inf_by_actorCzech <- as.data.frame(mean_inf_by_actorCzech[c(1:5),c(1:7)])
mean_inf_by_actorCzech 


#Cooperation Stats------------------------------------
library(igraph)
coopCzech <- read.csv(file="/Users/paulwagner/Desktop/ERGMs/Czech Republic/CoopCzechEdgeList.csv",header=TRUE)
coopCzech=network(coopCzech,matrix.type="edgelist",directed=TRUE) 
coopCzech <- as.matrix(coopCzech)
diag (coopCzech) <- 0 
coopCzech <- coopCzech[-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157),-c(1,3,6,8,16,17,21,32,37,38,40,44,49,50,62,67,70,71,74,81,82,84,86,88,91,95,100,104,105,106,107,109,110,112,115,116,117,118,120,124,125,128,132,133:157)]


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

nw.coopCzech.igraph <- graph.adjacency(coopCzech)
Results.coopCzech <- mycoopCentrality(nw.coopCzech.igraph)
#Results.coopCzech

# Rounding 

Results.coopCzech$Betweenness <- round(Results.coopCzech$Betweenness, digits =3)
Results.coopCzech$Closeness <- round(Results.coopCzech$Closeness, digits =3)
Results.coopCzech$Eigenvector <- round(Results.coopCzech$Eigenvector, digits =3)
Results.coopCzech$PageRank <- round(Results.coopCzech$PageRank, digits =3)
Results.coopCzech$IDegree <- round(Results.coopCzech$IDegree, digits =3)
Results.coopCzech$ODegree <- round(Results.coopCzech$ODegree, digits =3)

Results.coopCzech 
Results.coopCzech <- as.data.frame(Results.coopCzech)

## Calculate means per actor type
Results.coopCzech$Actor.type.Czech <- OrgTypeCzech

Results.coopCzech$ID <- NULL
Results.coopCzech$Actor.type.Czech <- as.character(Results.coopCzech$Actor.type.Czech)

mean_coop_by_actorCzech <- aggregate(Results.coopCzech, by=list(Results.coopCzech$Actor.type.Czech), FUN=mean, na.rm=TRUE )
mean_coop_by_actorCzech

mean_coop_by_actor$Actor.typeCzech <- NULL

mean_coop_by_actorCzech$Betweenness <- round(mean_coop_by_actorCzech$Betweenness, digits =3)
mean_coop_by_actorCzech$Closeness <-   round(mean_coop_by_actorCzech$Closeness, digits =3)
mean_coop_by_actorCzech$Eigenvector <- round(mean_coop_by_actorCzech$Eigenvector, digits =3)
mean_coop_by_actorCzech$PageRank <-    round(mean_coop_by_actorCzech$PageRank, digits =3)
mean_coop_by_actorCzech$IDegree <-     round(mean_coop_by_actorCzech$IDegree, digits =3)
mean_coop_by_actorCzech$ODegree <-     round(mean_coop_by_actorCzech$ODegree, digits =3)

mean_coop_by_actorCzech <- as.data.frame(mean_coop_by_actorCzech[c(1:5),c(1:7)])
mean_coop_by_actorCzech

