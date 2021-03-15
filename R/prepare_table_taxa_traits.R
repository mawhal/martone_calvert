
# need to use comm data from hmsc prep
hmsc <- data.frame( hmsc = T, taxon = colnames(Y))

# need to use comm data from rda script
rank_occurrence <- rank(-apply(dwide[,-c(1:4)],2,function(z) length(z[z>0])), ties.method = "min")
commtest <- gsub(" ",".",colnames(comm))

community <- data.frame( source = 'all', taxon = commtest, rank=rank_occurrence)

# join
taxa_modeled <- full_join(community, hmsc)


# add functional groups, mean shifts?
funs <- read_csv( "Data/taxa/Algae_functional_groups.csv")
lumps <- read_csv( "Data/taxa/TaxonList_corrected_lumped_unique.csv") %>% 
  filter( motile_sessile == "sessile") %>% 
  select( taxon_lumped2, non.alga.flag ) %>% 
  distinct()

funlumps <- left_join(lumps, funs, by = c("taxon_lumped2" = "taxon"))
funlumps$taxon <- gsub(" ",".",funlumps$taxon_lumped2 )

# now merge back with taxa_modeled
taxa_modeled_fun <- left_join( taxa_modeled, funlumps, by = "taxon" ) %>% 
  arrange(rank)

write.csv( taxa_modeled_fun, "R/output/taxa_modeled_funs.csv")
