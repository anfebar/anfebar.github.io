# The first section of code fits model 3 (see paper) using R-stan

# Read dataset
# nba.data = read.table("DatReplicationCollapsedIDbyTrade.txt",header = T)
nba.data = read.table(url("https://raw.githubusercontent.com/anfebar/anfebar.github.io/master/Software/BOARS/DatReplicationCollapsedIDbyTrade.txt"),header = T)

# Get game date and home and away team
Home = nba.data$Team1
Away = nba.data$Team2

# Get players
PlayersHnames = as.matrix(nba.data[,16:20])
colnames(PlayersHnames)
PlayersAnames = as.matrix(nba.data[,16:20+5])
colnames(PlayersAnames)

PlayersH = sapply(1:5,function(j) paste(Home,PlayersHnames[,j], sep = "_"))
PlayersA = sapply(1:5,function(j) paste(Away,PlayersAnames[,j], sep = "_"))

FreqPlayers = table(as.vector(cbind(PlayersH,PlayersA)))
NamesPlayers = sort(unique(c(unique(PlayersH),unique(PlayersA))))

# Define design matrix
M = matrix(0,ncol=length(NamesPlayers),nrow=nrow(PlayersA))
colnames(M) = NamesPlayers

for(j in 1:nrow(PlayersH))
{
  M[j,NamesPlayers %in% PlayersH[j,]]=1
  M[j,NamesPlayers %in% PlayersA[j,]]=-1
}
Ind = c(which((colnames(M)%in%names(which.max(colSums(abs(M)))))),
        which((colnames(M)%in%names(which.max(colSums(abs(M)))))==FALSE))
M = M[,Ind]
colnames(M)

# Check rank of M
dim(M)
qr(M)$rank

# Compute difference in points scored
DiffScore = nba.data$TotPoints1 - nba.data$TotPoints2

# "compile model that is subsequently used by stan
library(rstan)
rstan_options(auto_write = TRUE)
#Modelcpp <- stanc(file = "Model3.stan")
S <- "https://raw.githubusercontent.com/anfebar/anfebar.github.io/master/Software/BOARS/Model3.stan"
destfile <- tempfile()
download.file(S, destfile = destfile)
Modelcpp <- stanc(file = destfile)
Model <- stan_model(stanc_ret = Modelcpp, verbose=FALSE)


# Fitting the model
abilities_dat <- list(N = nrow(M),
                      K = ncol(M),
                      y = DiffScore,
                      X = M)

fit = sampling(Model, data = abilities_dat, 
               iter = 100, warmup = 20, thin = 1, cores=1,
               chains = 1, refresh = 10, algorithm="NUTS"	)


# The next section of code produces Tables 1, 2, and 3 of the document
# Install BOARS
install.packages("https://github.com/anfebar/anfebar.github.io/blob/master/Software/BOARS/BOARS_0.1.0.tar.gz?raw=true",
                 repos = NULL, type = "source")

library(BOARS)
library(xtable)



#
# Lineup analysis Table 1
#

# First compare lineups within a team for all playoff teams. 
# This took 7.15 minutes to run with a grid of length 21.


Ind = simplify2array(lapply(
  strsplit(colnames(as.matrix(fit)),"_cen"),length))==2
mcmc.xi = as.matrix(fit)[,Ind]
colnames(mcmc.xi) = colnames(M)
dim(mcmc.xi)
Nplayers = ncol(mcmc.xi)



teams <- unique(sapply(1:ncol(M),function(j)strsplit(colnames(M),"_")[[j]][1]))

playoffs <- c("Yes","Yes","Yes","Yes","Yes","Yes","Yes","No","No","No","No","No","Yes","No","Yes","Yes","No","No","No","No","Yes","Yes","No","Yes","Yes","No","Yes","No","Yes","No")
playoff.teams <- teams[playoffs=="Yes"]

team.lineups <- vector("list", length(playoff.teams))
optimal.specs <- vector("list", length(playoff.teams))
num.lineups <- numeric()
date()
for(ii in 1:length(playoff.teams)){
  cat("ii = ", ii, "\n")
  Indcol = sapply(1:ncol(M),function(j)strsplit(colnames(M),"_")[[j]][1])%in% playoff.teams[ii]
  aux = abs(M[rowSums(M[,Indcol])!=0,Indcol])
  auxlineup = unique(sapply(1:nrow(aux),function(j)paste(names(which(aux[j,]==1)),collapse="@")))

  Mlineup = M[1:length(auxlineup),]-M[1:length(auxlineup),]
  for(j in 1:length(auxlineup)) Mlineup[j,strsplit(auxlineup[j],"@")[[1]]]=1


  mcmc.xi.tmp = t(Mlineup%*%t(mcmc.xi))
  colnames(mcmc.xi.tmp)=auxlineup

  num.lineups[ii] <- length(auxlineup)
  ### Input
  Grid.alpha = sort(seq(0,0.05,len=21),decreasing = TRUE) # grid for alpha
  local.error = 0.1 # grid for t
  Grid.gamma.local = 1-Grid.alpha # grid for gamma
  global.error = local.error # grid for q

  team.lineups[[ii]] <- playoff.teams[ii]

  set.seed(2*ii)
  result.lineup = global_order(mcmc.xi.tmp, Grid.alpha,
                               local.error, Grid.gamma.local,
                               global.error, Pv_threshold=0.0)


  optimal.specs[[ii]] <- result.lineup$Optimal

  if(!is.null(result.lineup$global.statement)){
    tmp <- sapply(1:ncol(result.lineup$global.statement),
                  function(J){c(result.lineup$global.statement[,J]$param ,
                                sum(!is.na(result.lineup$global.statement[,J]$smaller)),
                                sum(!is.na(result.lineup$global.statement[,J]$greater)))})
    team.lineups[[ii]] <- data.frame(lineup=as.character(tmp[1,]),
                                    smaller = as.numeric(tmp[2,]),
                                    greater=as.numeric(tmp[3,]))
  } else {

    team.lineups[[ii]] <- playoff.teams[ii]

  }

}

date()

# Create table
df <- data.frame(team=NULL, lineup=NULL, Below=NULL, Above=NULL, Prob=NULL)
for(tt in 1:(length(team.lineups))){
	tmp <- team.lineups[[tt]]
	if(length(tmp)!=1){
		tmp <- tmp[order(tmp[,2], decreasing=TRUE),]
		lu <- as.character(tmp[,1])
    		lu1 <- strsplit(lu, "@|\\_")
		team <- lu1[[1]][1]
    		lu2 <- paste(lu1[[1]][!lu1[[1]] == team], collapse=", ")

    		df <- rbind(df, data.frame(team=team,
    		                           lineup=lu2, Below=tmp[1,2], Above=tmp[1,3],
    									Prob=as.numeric(optimal.specs[[tt]][5]),
    									stringsAsFactors=FALSE))
	} else {
		df <- rbind(df, data.frame(team=tmp, lineup=NA, Below=NA, Above=NA, Prob=NA, stringsAsFactors=FALSE))
	}
}

print(xtable(df, digits=c(0,0,0,0,0,2)), include.rownames=FALSE, booktabs=TRUE)


#
# Lineup analysis Table 2
#

# Took about 11.6 minutes to fit with the grid of length 21
# all lineups from play off teams compared between teams

Ind = simplify2array(lapply(
  strsplit(colnames(as.matrix(fit)),"_cen"),length))==2
mcmc.xi = as.matrix(fit)[,Ind]
colnames(mcmc.xi) = colnames(M)
dim(mcmc.xi)
Nplayers = ncol(mcmc.xi)

teams <- unique(sapply(1:ncol(M),function(j)strsplit(colnames(M),"_")[[j]][1]))

playoffs <- c("Yes","Yes","Yes","Yes","Yes","Yes","Yes","No","No","No","No","No","Yes","No","Yes","Yes","No","No","No","No","Yes","Yes","No","Yes","Yes","No","Yes","No","Yes","No")
playoff.teams <- teams[playoffs=="Yes"]


Indcol = sapply(1:ncol(M),function(j)strsplit(colnames(M),"_")[[j]][1])%in% playoff.teams
sum(Indcol) # number of players on playoff teams
aux = abs(M[,Indcol][rowSums(M[,Indcol])!=0,])
nrow(aux) # number of encounters among the playoff team lineups
auxlineup = unique(sapply(1:nrow(aux),function(j)paste(names(which(aux[j,]==1)),collapse="@")))
length(auxlineup) # number of unique lineups

Mlineup = M[1:length(auxlineup),]-M[1:length(auxlineup),]
for(j in 1:length(auxlineup)) Mlineup[j,strsplit(auxlineup[j],"@")[[1]]]=1

mcmc.xi.tmp = t(Mlineup%*%t(mcmc.xi))
colnames(mcmc.xi.tmp)=auxlineup
print(dim(mcmc.xi.tmp))


Grid.alpha = sort(seq(0,0.05,len=21),decreasing = TRUE) # grid for alpha
local.error = 0.1 # grid for t
Grid.gamma.local = 1-Grid.alpha # grid for gamma
global.error = local.error # grid for q

set.seed(101)
result.lineup = global_order(mcmc.xi.tmp, Grid.alpha,
                               local.error, Grid.gamma.local,
                               global.error, Pv_threshold=0.0)

result.lineup$Optimal
# Create table
lineup <- apply(result.lineup$global.statement, 2, function(x) x[[1]])
summary.lineup <- t(apply(result.lineup$global.statement, 2, function(x) sapply(x, function(x) sum(!is.na(x)))))
summary.lineup <- summary.lineup[order(lineup),]
team <- sapply(strsplit(lineup, "@|\\_"), function(x) x[[1]])[order(lineup)]
lineup  <- lineup[order(lineup)]

indx <- tapply(summary.lineup[,2], team, which.max)
indx <- c(indx[1], indx[-1] + cumsum(as.numeric(table(team)))[-16])

summary.lineup[indx,]
lineup[indx]
team[indx]

df <- data.frame(team=NULL, lineup=NULL, Below=NULL, Above=NULL)
for(tt in 1:(length(lineup[indx]))){

    lu1 <- strsplit(lineup[indx][tt], "@|\\_")
    lu2 <- paste(lu1[[1]][!lu1[[1]] == team[indx][tt]], collapse=", ")

    df <- rbind(df, data.frame(team=team[indx][tt], lineup=lu2,
    							Below=summary.lineup[indx,][tt,2],
    							Above=summary.lineup[indx,][tt,3],stringsAsFactors=FALSE))
}

print(xtable(df, digits=c(0,0,0,0,0)), include.rownames=FALSE, booktabs=TRUE)






# Individual player analysis

Ind = simplify2array(lapply(strsplit(colnames(as.matrix(fit)),"_cen"),length))==2
mcmc.xi = as.matrix(fit)[,Ind]
colnames(mcmc.xi) = colnames(M)

# Used this for the paper.  
# Total computation time 0.9 minutes using a grid of length 21.
Grid.alpha = sort(seq(0,0.05,len=21),decreasing = TRUE) # grid for alpha
local.error = 0.1 # grid for t
Grid.gamma.local = 1-Grid.alpha # grid for gamma
global.error = local.error # grid for q



set.seed(1)
result.player = global_order(mcmc.xi, Grid.alpha,
                         local.error, Grid.gamma.local,
                         global.error, Pv_threshold=0)


result.player$Optimal
#result$grid.eval
tmp <- sapply(1:ncol(result.player$global.statement),
            function(J){c(result.player$global.statement[,J]$param ,
                          sum(!is.na(result.player$global.statement[,J]$smaller)),
                          sum(!is.na(result.player$global.statement[,J]$greater)))})
play.ord <- data.frame(Player=tmp[1,], smaller = as.numeric(tmp[2,]), greater=as.numeric(tmp[3,]))

play.ord[order(play.ord$smaller,play.ord$greater, decreasing=TRUE),]

ss1 <- play.ord$smaller > play.ord$greater


xtable(cbind(play.ord[ss1,][order(play.ord[ss1,2], decreasing=TRUE),], 
             play.ord[!ss1,][order(play.ord[!ss1,3], decreasing=FALSE),]))
