library(did)
library(did2s)
library(tidyverse)
library(sf)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting
# buffer rings
source("functions.R")

dir.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
cov.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/model.covariates/"
file.dir <- data.frame(file = list.files(path = dir.path))
file.dir <- file.dir %>% mutate(country = sub(".forest.*", "", file),
                                buffer = sub(".buffer.forest.loss.df.RData", "", file),
                                buffer = sub(".forest.mines.", "", buffer),
                                buffer = str_remove(buffer, country))

file.1km <- file.dir %>% filter(buffer == "1km")
file.5km <- file.dir %>% filter(buffer == "5km")
file.10km <- file.dir %>% filter(buffer == "10km")
file.20km <- file.dir %>% filter(buffer == "20km")

ssa.1km <- read.bind.list(file.1km)
ssa.5km <- read.bind.list(file.5km)
ssa.10km <- read.bind.list(file.10km)
ssa.20km <- read.bind.list(file.20km)

ssa.ls <- list("1km" = ssa.1km, 
               "5km" = ssa.5km,
               "10km" = ssa.10km, 
               "20km" = ssa.20km)

## SSA wide --------------------------------------------------------------------

ssa.ls <- lapply(ssa.ls, function(x) {did.prep(x,
                                     lead.time = -23, post.time = 23,
                                     type = "loss", covariates = NULL) %>%
                          group_by(country, CLUSTER_ID) %>%
                          mutate(cluster.country.id = cur_group_id())})



ssa.names <- names(ssa.ls)
ssa.did <- data.frame()

for (i in 1:4) {
  buff.i <- ssa.ls[[i]]
  name.i <- ssa.names[[i]]
  
  gardner.did <- did2s(data = buff.i, yname = "cumulative.forest.loss.perc", 
                   treatment = "treatment",
                   first_stage =  ~ 0 + country | cluster.country.id + year, 
                   second_stage = ~ i(rel.year.first, ref = c(-1)),
                   cluster_var = "cluster.country.id", verbose = TRUE)
  
  gardner.tidy.1km <- did2s.tidy(gardner.did, buff = name.i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)))
  
  csa.did <- att_gt(data = buff.i, yname = "cumulative.forest.loss.perc", tname="year",
                    idname= "cluster.country.id", gname = "first.mine.year",
                    control_group = "notyettreated", base_period = "varying",
                    xformla = ~ country,
                    clustervars = "cluster.country.id", 
                    bstrap=T, cband=T)
  
  csa.tidy.1km <- csa.tidy(csa.did, buff = name.i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)))
  
  did.i <- rbind(gardner.tidy.1km, csa.tidy.1km)
  ssa.did <- rbind(ssa.did, did.i)
  
  
}

save(ssa.did, file = "Outputs/DiD.tables/SSA.tidy.1km.RData")



