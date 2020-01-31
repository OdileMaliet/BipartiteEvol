
source("ComputeMetricsEtablissement.R")
source("main.R")
source("beta.R")
library(bipartite)
library(apTreeshape)
library(ade4)


#### define the parameters #####

muP=0.01          # mutation rate clade A
muH=0.01          # mutation rate clade B

alphaP=0.1        # fitness function parameters
alphaH=0.1
rP=10
rH=10

nx=50             # grid size
ny=50

nP=1              # number of updated individuals at each time step (for each clade)
nH=1

NG=100000         # number of generation (you can add some afterwards if not long enough, see at the very end)

dSpace=0.5       # spatial kernel size (put to Inf to remove space, as I did in my simulations with Nicolas)

thin = 100

#### run the model ####
seed=as.integer(1)
set.seed(seed)

mod = sim.BipartiteEvol(nx,ny,NG,dSpace=dSpace,D=3,muP,muH,alphaP=alphaP,alphaH=alphaH,iniP=0,iniH=0,nP=nP,nH=nH,rP=rP,rH=rH,effect=1,
                             verbose=1000,thin=thin,P=NULL,H=NULL)


#### build the genealogy ####

gen = make.gen(mod)

# and plot it
par(mfrow = c(1,1))
plot(gen$H)
# not everything has calesced yet, more iterations should be added -- see at the end

#### buikd the phylogeny ####

# with a species definition threshold s = 1
phy1 = define.species(gen,threshold=1)

# and plot it
par(mfrow = c(1,1))
plot(phy1$Pphylo$tree)


#### a first plot of the result ####

# rmq : each of the pannels can also be ploted independantly
plot.model(gen,phy1, 1)


#### build the network ####

# still s = 1, change the phylogeny to change this parameter
xP = mod$P       # traits of the individuals in clade P
xH = mod$H

net = build.network(xP, xH, gen, phy1, alpha = alphaP, seuil=0, rP=rP, rH=rH)
# you can put a threshold and remove association with to low fitness with seuil, but I used 0 (the default, all associations are kept) in my case

#plot the result
trait.id = 1
par(mfrow = c(1,1))
plot.network(net,phy1,1)

#### a second plot of the model, including the network ####

trait.id = 1
plot.model.network(gen,phy1,trait.id, net,mod)


#### add iterations ####


# you can repeat the following as many times as you which to 
seed=as.integer(1)
set.seed(seed)

mod = sim.BipartiteEvoll(nx,ny,NG,dSpace=dSpace,D=3,muP,muH,alphaP=alphaP,alphaH=alphaH,iniP=0,iniH=0,nP=nP,nH=nH,rP=rP,rH=rH,effect=1,
                    verbose=1000,thin=thin,P=mod$P,H=mod$H)

# update the genealogy
gen = make.gen(out,treeP=gen$P, treeH=gen$H, verbose=T)
par(mfrow = c(1,1))
plot(gen$H)

# update the phylogenies
phy1 = define.species(gen,threshold=1)
par(mfrow = c(1,1))
plot(phy1$Pphylo$tree)

plot.model(gen,phy1, 1)

xP = mod$P       # traits of the individuals in clade P
xH = mod$H
net = build.network(xP, xH, gen, phy1, alpha = alphaP, seuil=0, rP=rP, rH=rH)
plot.model.network(gen,phy1,trait.id, net,mod)