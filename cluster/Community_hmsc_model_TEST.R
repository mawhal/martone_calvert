# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space


# load libraries
library( Hmsc )



# read models
models <- list.files( path = getwd(), pattern = "*specified.Rdata", recursive = T)
# pick a model
mselect <- models[1]
# load it - will bring "x" into the environment, which has everything we need to run a model
load( mselect ) 


## Run MCMC and save the model
# thin = 100
# samples = 1000
# nChains = 4
thin = 1
samples = 100
transient = 10
nChains = 1
set.seed(1)
# ptm = tic("model run")
m = sampleMcmc(m, samples = samples, thin = thin,
               # adaptNf = rep(ceiling(0.4*samples*thin),m$nr),
               transient = transient,
               nChains = nChains, nParallel = nChains,
               initPar = "fixed effects")
warnings()# computational.time <- toc()
model = "elevxtemp"

filename = file.path(getwd(), paste("model_",as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))
save(m,file=filename)#,computational.time

