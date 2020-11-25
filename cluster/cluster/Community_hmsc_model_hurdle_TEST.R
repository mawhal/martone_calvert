# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space


# load libraries
library( Hmsc )



# read models
model_path <- list.files( path = getwd(), pattern = "*specified.Rdata", recursive = T)
# pick models
mselect <- model_path[1]
# load it - will bring "x" into the environment, which has everything we need to run a model
load( mselect ) 


## Run MCMC and save the model
# thin = 100
# samples = 250
# nChains = 4
thin = 1
samples = 5
transient = ceiling(thin*samples*0.5)
nChains = 1
set.seed(1)
nm = length(models)
for( model in 1:nm ){
  print(paste0("model = ",modelnames[model]))
  m = models[[model]]
  m = sampleMcmc(m, samples = samples, thin = thin,
               # adaptNf = rep(ceiling(0.4*samples*thin),m$nr),
               transient = transient,
               nChains = nChains, #nParallel = nChains,
               )
  models[[model]] = m
}
warnings()# computational.time <- toc()
model = "elevxyear_hurdle_test"
filename = file.path(getwd(), paste("model_",as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))
save(models,file=filename)#,computational.time