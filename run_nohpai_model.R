#############################################################################
# Full posterior run — No-HPAI (pre-outbreak) scenario
#
# Companion to run_full_model.R. This fits the SAME model structure to
# pre-outbreak-only data (2010-2019/2020), giving the "what if HPAIV hadn't
# happened" counterfactual used in the paper's resilience comparison
# (Figure 4: recovery time, inertia). Reproduces Code_NoHPAI.R exactly.
#
# EXPECTED RUNTIME: ~6-7 hours, same order as the HPAI run.
# OUTPUT: pva_nohpai_full.rds (load this alongside pva_hpai_full.rds for the
# resilience/transient comparison plot in the tutorial).
#############################################################################

require(jagsUI)

dUnif <- function(lower, upper) {
  v <- rep(0, upper)
  v[lower:upper] <- 1 / (upper - lower + 1)
  v
}

## ---- 1. Load data --------------------------------------------------------

load("data_nohpai.RData")

jags.data <- list(
  y = data_nohpai$ypool,
  nyears = ncol(data_nohpai$ypool),
  nind = nrow(data_nohpai$ypool),
  f = data_nohpai$funic,
  age = data_nohpai$agesurvpool,
  agerec = data_nohpai$agerecpool,
  agedet = data_nohpai$agedetpool,
  nyscaled = as.numeric(scale(1:ncol(data_nohpai$ypool))),
  FR = data_nohpai$shortvec,
  fr = data_nohpai$shortvec,
  nstates = 6,

  J = data_nohpai$prodc$nfled,
  B = data_nohpai$prodc$nmonitored,

  counts = data_nohpai$cens[18:28],
  pNinit1 = dUnif(1, 70),
  pNinit2 = dUnif(1, 30),
  pNinit3 = dUnif(1, 10),
  pNinit4 = dUnif(1, 3),
  pNinit5 = dUnif(1, 100),
  yearspd = nrow(data_nohpai$prodc),
  yearsproj = 14
)

## ---- 2. JAGS model --------------------------------------------------------
## Identical to Code_NoHPAI.R (note: mean.omega upper bound is 50 here vs.
## 30 in the HPAI model -- matches the original code exactly).

cat(file = "pva2b.jags", "

    model{

    ###########################
    ###Productivity submodel###
    ###########################

  for(t in 1:yearspd){
    J[t] ~ dpois(B[t]*rho[t])
    log.rho[t] ~ dnorm(l.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }

  for(t in (yearspd+1):(yearspd+yearsproj)){
    log.rho[t] ~ dnorm(l.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }

  mean.rho ~ dunif(0,5)
  l.mean.rho <- log(mean.rho)
  sigma.rho ~ dunif(0,5)
  tau.rho <- pow(sigma.rho, -2)

  ##################################
  ###Population dynamics submodel###
  ##################################

  N[1,1] ~ dcat(pNinit1)
  N[2,1] ~ dcat(pNinit2)
  N[3,1] ~ dcat(pNinit3)
  N[4,1] ~ dcat(pNinit4)
  N[5,1] ~ dcat(pNinit3)
  N[6,1] ~ dcat(pNinit4)
  N[7,1] ~ dcat(pNinit4)
  N[8,1] ~ dcat(pNinit5)
  N[9,1] ~ dpois(omega[1])

  for(t in 1:(yearspd-1+yearsproj)){
    N[1,t+1] ~ dpois(rho[t] / 2 * s[1,t+17] * NB[t])
    N[2,t+1] ~ dbin(s[2,t+17]*gam[2], N[1,t])
    N[3,t+1] ~ dbin(s[2,t+17]*(1-gam[2]), N[1,t])
    N[4,t+1] ~ dbin(s[3,t+17]*gam[3], N[3,t])
    N[5,t+1] ~ dbin(s[3,t+17]*(1-gam[3]), N[3,t])
    N[6,t+1] ~ dbin(s[3,t+17]*gam[3], N[5,t])
    N[7,t+1] ~ dbin(s[3,t+17]*(1-gam[3]), N[5,t])
    N[8,t+1] ~ dbin(s[3,t+17], N[2,t]+N[4,t]+N[6,t]+N[7,t]+N[8,t]+N[9,t])
    N[9,t+1] ~ dpois(omega[t+1])

    lambda[t] <- Ntot[t+1]/(Ntot[t]+0.001)
    lambdaNB[t] <- NB[t+1]/NB[t]
  }

  for(t in 1:(yearspd+yearsproj)){
    NB[t] <- N[2,t]+N[4,t]+N[6,t]+N[8,t]+N[9,t]
    Nfloat[t] <- N[3,t]+N[5,t]+N[7,t]
    Nlr[t] <- N[2,t]+N[4,t]+N[6,t]
    Ntot[t] <- sum(N[,t])
    propNB[t] <- NB[t]/Ntot[t]
    propfloat[t] <- Nfloat[t]/Ntot[t]
    propjuv[t] <- N[1,t]/Ntot[t]
    Nfledglings[t] ~ dpois(rho[t] / 2 * NB[t])
  }

  for(t in 1:yearspd){
    counts[t] ~ dpois(NB[t])
  }

  mean.omega ~ dunif(0.001, 50)
  log.mean.omega <- log(mean.omega)
  sigma.omega ~ dunif(0.001,5)
  tau.omega <- pow(sigma.omega, -2)

  for(t in 1:(yearspd+yearsproj)){
    log.omega[t] ~ dnorm(log.mean.omega, tau.omega)T(0, 4.6)
    omega[t] <- exp(log.omega[t])
  }

    #############################
    ###Mark-recapture submodel###
    #############################

  for(a in 1:3){
    mu.s[a] ~ dnorm(0, 0.001)
    mean.s[a] <- ilogit(mu.s[a])
    sigma.s[a] ~ dunif(0.001, 10)
    tau.s[a] <- pow(sigma.s[a],-2)

    mu.p[a] ~ dnorm(0, 0.001)
    sigma.p[a] ~ dunif(0.001, 10)
    mean.p[a] <- ilogit(mu.p[a])
    tau.p[a] <- pow(sigma.p[a],-2)

        for(t in 1:(nyears-1)){
      logit.p[a,t] ~ dnorm(mu.p[a], tau.p[a])
      p[a,t] <- ilogit(logit.p[a,t])
        }
  }

    for(a in 1:3){
      for(t in 1:(nyears-1+yearsproj)){
      logit.s[a,t] ~ dnorm(mu.s[a], tau.s[a])
      s[a,t] <- ilogit(logit.s[a,t])
      }
    }

  gam[1] <- 0
  gam[2] ~ dbeta(1,1)
  gam[3] ~ dbeta(1,1)
  gam[4] <- 1

  for(t in 1:(nyears-1)){
    logit.r[t] <- alpha.r + beta.r * nyscaled[t]
    r[t] <- ilogit(logit.r[t])
  }

  alpha.r ~ dnorm(0, 0.001)
  beta.r ~ dnorm(0, 0.001)

  for(i in 1:nind){
    for (t in f[i]:(nyears-1)){
      psi[1,1,i,t] <- s[age[i,t],t]*(1-gam[agerec[i,t]])
      psi[1,2,i,t] <- 0
      psi[1,3,i,t] <- s[age[i,t],t]*gam[agerec[i,t]]
      psi[1,4,i,t] <- 0
      psi[1,5,i,t] <- 1-s[age[i,t],t]
      psi[1,6,i,t] <- 0

      psi[2,1,i,t] <- 0
      psi[2,2,i,t] <- s[age[i,t],t]*(1-gam[agerec[i,t]])
      psi[2,3,i,t] <- 0
      psi[2,4,i,t] <- s[age[i,t],t]*gam[agerec[i,t]]
      psi[2,5,i,t] <- 1-s[age[i,t],t]
      psi[2,6,i,t] <- 0

      psi[3,1,i,t] <- 0
      psi[3,2,i,t] <- 0
      psi[3,3,i,t] <- s[age[i,t],t]
      psi[3,4,i,t] <- 0
      psi[3,5,i,t] <- 1-s[age[i,t],t]
      psi[3,6,i,t] <- 0

      psi[4,1,i,t] <- 0
      psi[4,2,i,t] <- 0
      psi[4,3,i,t] <- 0
      psi[4,4,i,t] <- s[age[i,t],t]
      psi[4,5,i,t] <- 1-s[age[i,t],t]
      psi[4,6,i,t] <- 0

      psi[5,1,i,t] <- 0
      psi[5,2,i,t] <- 0
      psi[5,3,i,t] <- 0
      psi[5,4,i,t] <- 0
      psi[5,5,i,t] <- 0
      psi[5,6,i,t] <- 1

      psi[6,1,i,t] <- 0
      psi[6,2,i,t] <- 0
      psi[6,3,i,t] <- 0
      psi[6,4,i,t] <- 0
      psi[6,5,i,t] <- 0
      psi[6,6,i,t] <- 1

      po[1,1,i,t] <- p[agedet[i,t], t]
      po[1,2,i,t] <- 0
      po[1,3,i,t] <- 0
      po[1,4,i,t] <- 0
      po[1,5,i,t] <- 0
      po[1,6,i,t] <- 1-p[agedet[i,t], t]

      po[2,1,i,t] <- 0
      po[2,2,i,t] <- 0
      po[2,3,i,t] <- 0
      po[2,4,i,t] <- 0
      po[2,5,i,t] <- 0
      po[2,6,i,t] <- 1

      po[3,1,i,t] <- 0
      po[3,2,i,t] <- 0
      po[3,3,i,t] <- p[agedet[i,t], t]
      po[3,4,i,t] <- 0
      po[3,5,i,t] <- 0
      po[3,6,i,t] <- 1-p[agedet[i,t], t]

      po[4,1,i,t] <- 0
      po[4,2,i,t] <- 0
      po[4,3,i,t] <- 0
      po[4,4,i,t] <- 0
      po[4,5,i,t] <- 0
      po[4,6,i,t] <- 1

      po[5,1,i,t] <- 0
      po[5,2,i,t] <- 0
      po[5,3,i,t] <- 0
      po[5,4,i,t] <- 0
      po[5,5,i,t] <- r[t]
      po[5,6,i,t] <- 1-r[t]

      po[6,1,i,t] <- 0
      po[6,2,i,t] <- 0
      po[6,3,i,t] <- 0
      po[6,4,i,t] <- 0
      po[6,5,i,t] <- 0
      po[6,6,i,t] <- 1
    }
  }

  # Likelihood (Marginalized, Yackulic et al., 2020)

for (i in 1:nind){

  zeta[i,f[i],1] <- equals(1,y[i,f[i]])
  zeta[i,f[i],2] <- equals(2,y[i,f[i]])
  zeta[i,f[i],3] <- equals(3,y[i,f[i]])
  zeta[i,f[i],4] <- equals(4,y[i,f[i]])
  zeta[i,f[i],5] <- 0
  zeta[i,f[i],6] <- 0

  for (t in f[i]:(nyears-1)){
    zeta[i,t+1,1] <- (zeta[i,t,1:6] %*% psi[,1,i,t]) * po[1,y[i,t+1],i,t]
    zeta[i,t+1,2] <- (zeta[i,t,1:6] %*% psi[,2,i,t]) * po[2,y[i,t+1],i,t]
    zeta[i,t+1,3] <- (zeta[i,t,1:6] %*% psi[,3,i,t]) * po[3,y[i,t+1],i,t]
    zeta[i,t+1,4] <- (zeta[i,t,1:6] %*% psi[,4,i,t]) * po[4,y[i,t+1],i,t]
    zeta[i,t+1,5] <- (zeta[i,t,1:6] %*% psi[,5,i,t]) * po[5,y[i,t+1],i,t]
    zeta[i,t+1,6] <- (zeta[i,t,1:6] %*% psi[,6,i,t]) * po[6,y[i,t+1],i,t]
  } #t

  lik[i] <- sum(zeta[i,nyears,])
  fr[i] ~ dbin(lik[i], FR[i])
}

    #Posterior predictive checks for the count and the productivity model

     for(t in 1:yearspd){

    predcounts[t] ~ dpois(NB[t])
    checkcounts[t,1] <- abs((predcounts[t] - NB[t])/(NB[t]+0.0001))
    checkcounts[t,2] <- abs((counts[t] - NB[t])/(NB[t]+0.0001))

    predfledg[t] ~ dpois(B[t]*rho[t])
    checkfledg[t,1] <- pow(predfledg[t] - B[t]*rho[t], 2)/(B[t]*rho[t])
    checkfledg[t,2] <- pow(J[t] - B[t]*rho[t],2)/(B[t]*rho[t])

    }

    chi2[1,1] <- sum(checkcounts[,1])
    chi2[1,2] <- sum(checkcounts[,2])
    chi2[2,1] <- sum(checkfledg[,1])
    chi2[2,2] <- sum(checkfledg[,2])

    }
    ")

## ---- 3. Initial values ----------------------------------------------------

nyears <- ncol(data_nohpai$ypool)
yearspd <- 11
yearsproj <- 14

inits <- function(){list(
  sigma.s = rep(.5,3),
  logit.s = matrix(nrow = 3, ncol = nyears-1+yearsproj, 0),
  mu.s = c(-1,0.2,1),
  mu.p = c(-1,0,1),
  sigma.p = rep(.5,3),
  logit.p = matrix(nrow = 3, ncol = nyears-1, -1),
  alpha.r = -1,
  beta.r = .5,
  gam = c(NA, rep(.5,2), NA),
  mean.rho = 2,
  sigma.rho = 1,
  log.rho = rep(1, yearspd+yearsproj)
)}

## ---- 4. Parameters to monitor --------------------------------------------

parameters <- c("s", "mean.s", "sigma.s", "p", "sigma.p", "mean.p", "r",
                 "alpha.r", "beta.r", "rho", "mean.rho", "sigma.rho",
                 "mean.omega", "sigma.omega", "NB", "Nfloat", "Nlr", "N",
                 "gam", "chi2", "lambda", "lambdaNB", "propNB", "propfloat",
                 "propjuv", "Ntot", "Nfledglings")

## ---- 5. Run (paper's original MCMC settings) ------------------------------

ni <- 80000; nb <- 40000; nc <- 4; nt <- 40

cat("Starting MCMC run (No-HPAI scenario) at", format(Sys.time()), "\n")
cat("Settings: n.chains =", nc, " n.iter =", ni, " n.burnin =", nb, " n.thin =", nt, "\n")
cat("Expected runtime: several hours.\n\n")

t_start <- Sys.time()

out <- jags(jags.data, inits, parameters, "pva2b.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            parallel = TRUE)

t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "mins")
cat("\nMCMC finished at", format(t_end), "\n")
cat("Elapsed time:", round(as.numeric(elapsed), 1), "minutes\n")

## ---- 6. Save results -------------------------------------------------------

saveRDS(out, "pva_nohpai_full.rds")
write.csv(out$summary, "pva_nohpai_full_summary.csv")

cat("\nSaved: pva_nohpai_full.rds (load alongside pva_hpai_full.rds for the\n")
cat("resilience/transient comparison plot in the tutorial)\n")
cat("Max Rhat across monitored parameters:", round(max(out$summary[, "Rhat"], na.rm = TRUE), 3), "\n")
