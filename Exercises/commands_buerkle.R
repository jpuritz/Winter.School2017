## R commands for Buerkle lectures and exercises

## section 1.2--#5.e
p<-seq(0,1,0.001)
plot(p, dbeta(p, shape1=160+1, shape2=200-160+1), 
     type ="l", xlab="p", ylab="density")
abline(v=qbeta(p=c(0.025, 0.975), shape1=160+1, shape2=200-160+1))
qbeta(p=c(0.025, 0.975), shape1=160+1, shape2=200-160+1)

## section 1.3
p<-seq(0,1,0.001)
par(mfrow=c(1,3))
plot(p, dbeta(p, shape1=1, shape2=1),
     type="l", main="1,1", xlab="p", ylab="Density")
plot(p, dbeta(p, shape1=100, shape2=100),
     type="l", main="100,100", xlab="p", ylab="Density")
plot(p, dbeta(p, shape1=0.1, shape2=0.1),
     type="l", main="0.1,0.1", xlab="p", ylab="Density")

## section 3.2--#4
p<-seq(0,1,0.001)
plot(p, dbeta(p, shape1 =    0.5  * (-1 + 1/0.4),
                 shape2 = (1-0.5) * (-1 + 1/0.4)),
     col="red", type="l", ylab="density", xlab="p", ylim=c(0,10))
lines(p, dbeta(p,shape1 =    0.5  * (-1 + 1/0.01),
                 shape2 = (1-0.5) * (-1 + 1/0.01)),
      col="blue")
abline(v=0.5, col="red")


## F-model simulation and Bayesian model to recover simulated values
##
## simulate allele frequencies at 100 loci in ancestral population
nloci <- 100
nind <- 25
## generate loci that are likely to be variable with beta(15,15).
## Invariant loci give JAGS trouble. Additionally, drop
## invariant loci in sim.g below. 
sim.pi <- rbeta(nloci, 15, 15)
## simulate allele frequencies in three derived (Fst=0.01) populations
sim.p <- matrix(nrow=nloci, ncol=3)
for(k in 1:3){
  sim.p[,k] <- rbeta(nloci, sim.pi * (-1 + 1/0.01),
                            (1-sim.pi) * (-1 + 1/0.01))
}

plot(sim.pi, sim.p[,1]) ## compare ancestral and derived allele frequencies

## simulate genotypes for nind in each population
sim.g <- array(0, dim=c(nloci,nind,3))
for(k in 1:3){
  sim.g[,,k] <- matrix(rbinom(nloci*nind, 2, prob=sim.p[,k]),
                       nrow=nloci, ncol=nind)
}
## drop invariant loci
todrop<-logical(nloci)
for(i in 1:nloci){
  todrop[i]<-as.logical(sum(apply(sim.g[i,,], 2, function(x){sum(x)/(nind*2)==0 | sum(x)/(nind*2)==1})))
}
nloci<-nloci-sum(todrop)
sim.g<-sim.g[!todrop,,]

## Use R to JAGS interface to estimate Fst from these data
library(rjags)

locusFmodel<-"model {
  for(i in 1:nloci){ 
    for(j in 1:nind){
      ## binomial likelihood for genotype = number allele copies of reference allele
      for(k in 1:npop){
        g[i,j,k] ~ dbinom(p[i,k], 2)  
      }		
    }
  }
  
  for(i in 1:nloci){
    for(k in 1:npop){
      ## population allele frequency in sample populations
      p[i,k] ~ dbeta(0.001+pi[i]*theta[i], 0.001+(1-pi[i])*theta[i])
    }
    theta[i] <- -1 + 1/Fst[i]
    Fst[i] ~  dbeta(0.001+psi*S, 0.001+(1-psi)*S) 
    pi[i] ~  dbeta(1,1)   
  }
  psi ~ dbeta(1, 1)
  S ~ dunif(0.01, 1000) 
}"



mod.jags <- jags.model(textConnection(locusFmodel),
                       data=list(nind=nind, nloci=nloci, npop=3, g=sim.g),
                       n.chains=2)
mod.sam <- jags.samples(model=mod.jags, variable.names=c("Fst", "psi"), n.iter=2000, thin=2)

### plot chains
plot(c(mod.sam$psi[1,,1], mod.sam$psi[1,,2])) ## plot both chains for genome-wide Fst (psi)

plot(mod.sam$Fst[2,,1]) ## for locus 2, chain 1
mean(mod.sam$Fst[2,,1]) ## for locus 2, chain 1
quantile(mod.sam$Fst[2,,1], prob=c(0.025, 0.5, 0.975)) ## for locus 2, chain 1

