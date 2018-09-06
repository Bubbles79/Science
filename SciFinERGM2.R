detach("package:igraph", unload=TRUE)
library(sna)
library(ergm)
seed <- 12345
set.seed (seed)

coopFin <- read.csv(file="/Users/paulwagner/Desktop/ERGMs/Finland/AprilData2018/CollabAprilFin.csv",header=TRUE)
SciFin <- read.csv(file= "/Users/paulwagner/Desktop/ERGMS/Finland/AprilData2018/SciAprilFin.csv",header=TRUE)
infFin <- read.csv(file= "/Users/paulwagner/Desktop/ERGMS/Finland/AprilData2018/InfAprilFin.csv", header = TRUE)
OrgTypeFin <- read.csv(file="/Users/paulwagner/Desktop/ERGMS/Finland/OrgTypeFin.csv", header = TRUE)
coopFin=network(coopFin,matrix.type="edgelist",directed=TRUE) 
SciFin=network(SciFin,matrix.type="edgelist",directed=TRUE)
infFin=network(infFin,matrix.type="edgelist",directed=TRUE) 
coopFin <- as.matrix(coopFin)
SciFin <- as.matrix(SciFin)
infFin <- as.matrix(infFin)
coopFin <- coopFin[c(1:96),c(1:96)]
SciFin <- SciFin[c(1:96),c(1:96)]
infFin <- infFin[c(1:96),c(1:96)]

diag (coopFin) <- 0
diag (SciFin)  <- 0
diag (infFin) <- 0

coopFin <-coopFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85),-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85)]
SciFin  <-SciFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85),-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85)]
infFin  <-infFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85),-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85)]
OrgTypeFin  <-OrgTypeFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85), ]

FinOrgNames <- as.vector(colnames(coopFin))

rownames(coopFin) <- FinOrgNames
colnames(coopFin) <- FinOrgNames
rownames(SciFin) <- FinOrgNames
colnames(SciFin) <- FinOrgNames
rownames(infFin) <- FinOrgNames
colnames(infFin) <- FinOrgNames

OrgTypeFin <- OrgTypeFin[,4]
OrgTypeFin <- as.character(OrgTypeFin)
OrgTypeFin <- as.vector(OrgTypeFin)

#Policy beliefs
PBsFin = read.csv("/Users/paulwagner/Desktop/ERGMs/Finland/BeliefsFin.csv", header = TRUE, stringsAsFactors=FALSE)
PBsFin <- PBsFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85),-c(1,2,10,11,12,13,15) ]
PBsFin[PBsFin == 97] <- 3
Pol.dist.Fin <- dist(PBsFin, method = "manhattan")
PrefSimMatFin <- max(Pol.dist.Fin) - Pol.dist.Fin
Pbs.Fin  <-as.matrix(PrefSimMatFin)

rownames(Pbs.Fin) <- FinOrgNames
colnames(Pbs.Fin) <- FinOrgNames

#Science Beliefs
SciBs.Fin = read.csv("/Users/paulwagner/Desktop/SciERGMs/SciQs/FinSciQs.csv", header = TRUE)
SciBs.Fin <- SciBs.Fin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85,97:121),-c(1,2,6,7) ]
SciBs.Fin[SciBs.Fin == 97] <- 3
SciBs.dist.Fin <- dist(SciBs.Fin, method = "manhattan")
SciBsSimMatFin <- max(SciBs.dist.Fin) - SciBs.dist.Fin
SciBs.FinMat <-as.matrix(SciBsSimMatFin)

rownames(SciBs.FinMat) <- FinOrgNames
colnames(SciBs.FinMat) <- FinOrgNames

u <- sort(unique(OrgTypeFin))
nodecov <- match(OrgTypeFin,u)

#"BUS" "CIV" "GOV" "NGO" "SCI"

nw.coopFin <- network (coopFin)
nw.SciFin <- network (SciFin) # create network object
set.vertex.attribute (nw.SciFin, "OrgTypeFin", OrgTypeFin)
set.vertex.attribute (nw.coopFin, "OrgTypeFin", OrgTypeFin)
set.vertex.attribute (nw.SciFin, "influence", degree (infFin, cmode = "indegree"))
set.vertex.attribute (nw.coopFin, "SciID", degree (SciFin, cmode = "indegree"))
set.vertex.attribute (nw.coopFin, "influence", degree (infFin, cmode = "indegree"))


Fin.sci.sci <- matrix (0 , nrow = nrow (coopFin) , ncol = ncol (coopFin))
for (i in 1: nrow (Fin.sci.sci )) {
  for (j in 1: ncol ( Fin.sci.sci)) {
    if (( OrgTypeFin [i] == "SCI" && OrgTypeFin [j] == "SCI") ||
        ( OrgTypeFin [i] == "SCI" && OrgTypeFin [j] == "SCI")) {
      Fin.sci.sci [i, j] <- 1
      Fin.sci.sci [j, i] <- 0
    }
  }
}


##collab network as dependent variable Aug 27th
#Same models as below, but building complexity
SciFin.m1 <- ergm(nw.coopFin ~ edges, 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciFin.m1)
plot(mcmc.diagnostics(SciFin.m1))
plot(gof(SciFin.m1))

SciFin.m2 <- ergm(nw.coopFin ~ edges + 
                    edgecov(Fin.sci.sci) + nodeifactor ("OrgTypeFin", base=-c(5)), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciFin.m2)
plot(mcmc.diagnostics(SciFin.m2))
plot(gof(SciFin.m2))

SciFin.m3 <- ergm(nw.coopFin ~ edges + 
                    edgecov(Fin.sci.sci) + nodeifactor ("OrgTypeFin", base=-c(5)) +
                    edgecov(Pbs.Fin) + edgecov(SciBs.FinMat) + edgecov (SciFin) + nodeicov("SciID"), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciFin.m3)
plot(mcmc.diagnostics(SciFin.m3))
plot(gof(SciFin.m3))

SciFin.m4 <- ergm(nw.coopFin ~ edges + 
                    edgecov(Fin.sci.sci) + nodeifactor ("OrgTypeFin", base=-c(5)) +
                    edgecov(Pbs.Fin) + edgecov(SciBs.FinMat) + edgecov (SciFin) + nodeicov("SciID") + 
                    edgecov(infFin) , 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 1000 , MCMC.interval = 1000))
summary(SciFin.m4)
plot(mcmc.diagnostics(SciFin.m4))
plot(gof(SciFin.m4))

SciFin.m5 <- ergm(nw.coopFin ~ edges + 
                     edgecov(Fin.sci.sci) + nodeifactor ("OrgTypeFin", base=-c(5)) +
                     edgecov(Pbs.Fin) + edgecov(SciBs.FinMat) + edgecov (SciFin) + nodeicov("SciID") + 
                     edgecov(infFin) +  
                     mutual +  twopath + gwodegree(3, fixed = TRUE ) +  
                     gwesp(0.1 , fixed = TRUE ), 
                   eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 3000 , MCMC.interval = 3000))
summary(SciFin.m5)
plot(mcmc.diagnostics(SciFin.m5))
plot(gof(SciFin.m5))

texreg::screenreg(list(SciFin.m1,SciFin.m2,SciFin.m3,SciFin.m4,SciFin.m5),single.row = T)

texreg::htmlreg(list(SciFin.m1,SciFin.m2,SciFin.m3,SciFin.m4,SciFin.m5),"/Users/paulwagner/Desktop/SciERGMs/Models Aug 27/FinAug29",single.row = T)


#---------------


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

ggnet2(nw.SciFin, alpha = 0.5, label = TRUE, label.size = 3, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.5),  node.size = 6, 
       color = "OrgTypeFin",
       color.legend = "OrgTypeFin", edge.color = "black", 
       palette = c("CIV" = "purple", "GOV" = "red", "BUS" = "blue", "SCI" = "yellow", "NGO" = "green"), 
       size = "indegree", size.min = 1, size.cut = 10, size.legend = "Indegree centrality",
       legend.position = "bottom")
ggsave("/Users/paulwagner/Desktop/ERGMs/SciERGMs/NetworkGraphs/Fin_graph.png", width = 29, height = 18,  device = NULL, dpi = 300)    

# Vertex Centrality

# Function for centrality stats

detach("package:statnet", unload=TRUE)
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

nw.SciFin.igraph <- graph.adjacency(SciFin)
Results1Fin <- myCentrality(nw.SciFin.igraph)

# Rounding 

Results1Fin$Betweenness <- round(Results1Fin$Betweenness, digits =3)
Results1Fin$Closeness <- round(Results1Fin$Closeness, digits =3)
Results1Fin$Eigenvector <- round(Results1Fin$Eigenvector, digits =3)
Results1Fin$PageRank <- round(Results1Fin$PageRank, digits =3)
Results1Fin$IDegree <- round(Results1Fin$IDegree, digits =3)
Results1Fin$ODegree <- round(Results1Fin$ODegree, digits =3)

Results1Fin 
Results1Fin <- as.data.frame(Results1Fin)

write.csv2(Results1Fin, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/CentralityScores/centrality.scores.Fin.csv")

Results1Fin$Actor.type.Fin <- OrgTypeFin

## Calculate means per actor type
Results1Fin$ID <- NULL
Results1Fin$Actor.type.Fin <- as.character(Results1Fin$Actor.type.Fin)

mean_by_actorFin <- aggregate(Results1Fin, by=list(Results1Fin$Actor.type.Fin), FUN=mean, na.rm=TRUE )
mean_by_actorFin

mean_by_actor$Actor.typeFin <- NULL

mean_by_actorFin$Betweenness <- round(mean_by_actorFin$Betweenness, digits =3)
mean_by_actorFin$Closeness <- round(mean_by_actorFin$Closeness, digits =3)
mean_by_actorFin$Eigenvector <- round(mean_by_actorFin$Eigenvector, digits =3)
mean_by_actorFin$PageRank <- round(mean_by_actorFin$PageRank, digits =3)
mean_by_actorFin$IDegree <- round(mean_by_actorFin$IDegree, digits =3)
mean_by_actorFin$ODegree <- round(mean_by_actorFin$ODegree, digits =3)

mean_by_actorFin <- as.data.frame(mean_by_actorFin[c(1:5),c(1:7)])
mean_by_actorFin 
write.csv2(mean_by_actor, file = "/Users/paulwagner/Desktop/ERGMs/SciERGMs/Fin.centrality.scores-mean-by-actor.csv")

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

nw.infFin.igraph <- graph.adjacency(infFin)
Results.infFin <- myInfCentrality(nw.infFin.igraph)
#Results.infFin

# Rounding 

Results.infFin$Betweenness <- round(Results.infFin$Betweenness, digits =3)
Results.infFin$Closeness <- round(Results.infFin$Closeness, digits =3)
Results.infFin$Eigenvector <- round(Results.infFin$Eigenvector, digits =3)
Results.infFin$PageRank <- round(Results.infFin$PageRank, digits =3)
Results.infFin$IDegree <- round(Results.infFin$IDegree, digits =3)
Results.infFin$ODegree <- round(Results.infFin$ODegree, digits =3)

Results.infFin 
Results.infFin <- as.data.frame(Results.infFin)

Results.infFin$Actor.type.Fin <- OrgTypeFin

## Calculate means per actor type
Results.infFin$ID <- NULL
Results.infFin$Actor.type.Fin <- as.character(Results.infFin$Actor.type.Fin)

mean_inf_by_actorFin <- aggregate(Results.infFin, by=list(Results.infFin$Actor.type.Fin), FUN=mean, na.rm=TRUE )
mean_inf_by_actorFin

mean_inf_by_actorFin$Actor.typeFin <- NULL

mean_inf_by_actorFin$Betweenness <- round(mean_inf_by_actorFin$Betweenness, digits =3)
mean_inf_by_actorFin$Closeness <- round(mean_inf_by_actorFin$Closeness, digits =3)
mean_inf_by_actorFin$Eigenvector <- round(mean_inf_by_actorFin$Eigenvector, digits =3)
mean_inf_by_actorFin$PageRank <- round(mean_inf_by_actorFin$PageRank, digits =3)
mean_inf_by_actorFin$IDegree <- round(mean_inf_by_actorFin$IDegree, digits =3)
mean_inf_by_actorFin$ODegree <- round(mean_inf_by_actorFin$ODegree, digits =3)

mean_inf_by_actorFin <- as.data.frame(mean_inf_by_actorFin[c(1:5),c(1:7)])
mean_inf_by_actorFin 


#Cooperation Stats------------------------------------
library(igraph)
coopFin <- read.csv(file="/Users/paulwagner/Desktop/ERGMs/Finland/AprilData2018/CollabAprilFin.csv",header=TRUE)
coopFin=network(coopFin,matrix.type="edgelist",directed=TRUE) 
coopFin <- as.matrix(coopFin)
coopFin <- coopFin[c(1:96),c(1:96)]
coopFin <-coopFin[-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85),-c(4,8,12,25,28,43,44,46,52,53,56,69,73,85)]
diag (coopFin) <- 0

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

nw.coopFin.igraph <- graph.adjacency(coopFin)
Results.coopFin <- mycoopCentrality(nw.coopFin.igraph)
#Results.coopFin

# Rounding 

Results.coopFin$Betweenness <- round(Results.coopFin$Betweenness, digits =3)
Results.coopFin$Closeness <- round(Results.coopFin$Closeness, digits =3)
Results.coopFin$Eigenvector <- round(Results.coopFin$Eigenvector, digits =3)
Results.coopFin$PageRank <- round(Results.coopFin$PageRank, digits =3)
Results.coopFin$IDegree <- round(Results.coopFin$IDegree, digits =3)
Results.coopFin$ODegree <- round(Results.coopFin$ODegree, digits =3)

Results.coopFin 
Results.coopFin <- as.data.frame(Results.coopFin)

Results.coopFin$Actor.type.Fin <- OrgTypeFin














#Stats---------
detach("package:igraph",unload=TRUE)
library(sna)
plot(cug.test(SciFin,centralization,cmode="dyad",FUN.arg=list(FUN=degree,cmode="indegree")))

den.Fin <- gden(SciFin)
SciFin_in  <- degree(SciFin, cmode="indegree")
SciFin_out <- degree(SciFin, cmode="outdegree")
mean.Fin <- mean(SciFin_in)
sd.Fin <- sd(SciFin_in)

Centra.Fin.In <- centralization(SciFin, degree, cmode="indegree")
Centra.Fin.out <- centralization(SciFin, degree, cmode="outdegree")
Centra.Fin.evcent <- centralization(SciFin, evcent) 
grecip.Fin <- mutuality(SciFin) 
gtrans.Fin <- gtrans(SciFin)

efficiency.Fin <- efficiency(SciFin)
efficiency.krackhard.Fin <- hierarchy(SciFin, measure= "krackhardt")
hierarchy.Fin <-  hierarchy(SciFin, measure= "reciprocity")
lubness.Fin <- lubness(SciFin)


