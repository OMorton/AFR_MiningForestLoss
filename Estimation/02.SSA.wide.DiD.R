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
cov.dir <- data.frame(file = list.files(path = cov.path, full.names = TRUE))


file.1km <- file.dir %>% filter(buffer == "1km")
file.5km <- file.dir %>% filter(buffer == "5km")
file.10km <- file.dir %>% filter(buffer == "10km")
file.20km <- file.dir %>% filter(buffer == "20km")
file.5km.master <- file.dir %>% filter(buffer == "5km.master")

ssa.1km <- read.bind.list(file.1km)
ssa.5km <- read.bind.list(file.5km)
ssa.10km <- read.bind.list(file.10km)
ssa.20km <- read.bind.list(file.20km)
ssa.5km.master <- read.bind.list(file.5km.master)

ssa.ls <- list("1km" = ssa.1km, 
               "5km" = ssa.5km,
               "10km" = ssa.10km, 
               "20km" = ssa.20km,
               "5km.master" = ssa.5km.master)

## 5km specific
covs.all <- read.covs.list(cov.dir) %>%
  mutate(country = ifelse(country =="Côte d'Ivoire", "Côte d_Ivoire", country))



## SSA wide --------------------------------------------------------------------

ssa.ls <- lapply(ssa.ls, function(x) {did.prep(x,
                                     lead.time = -23, post.time = 23,
                                     type = "loss", covariates = NULL) %>%
                          group_by(country, CLUSTER_ID) %>%
                          mutate(cluster.country.id = cur_group_id())})



ssa.names <- names(ssa.ls)
ssa.did <- data.frame()
ssa.did.covar <- data.frame()

i <- 1
for (i in 1:5) {
  buff.i <- ssa.ls[[i]]
  name.i <- ssa.names[[i]]
 # main 5
  gardner.did <- did2s(data = filter(buff.i, rel.year.first >= -5),
                       yname = "cumulative.forest.loss.perc", 
                   treatment = "treatment",
                   first_stage =  ~ 0 + i(country, year) | cluster.country.id + year, 
                   second_stage = ~ i(rel.year.first, ref = c(-1)),
                   cluster_var = "cluster.country.id", verbose = TRUE)
  
  gardner.tidy <- did2s.tidy(gardner.did, buff = name.i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)),
           pre.period = -5)
  
  # main 10
  gardner.did.10 <- did2s(data = filter(buff.i, rel.year.first >= -10),
                       yname = "cumulative.forest.loss.perc", 
                       treatment = "treatment",
                       first_stage =  ~ 0 + i(country, year) | cluster.country.id + year, 
                       second_stage = ~ i(rel.year.first, ref = c(-1)),
                       cluster_var = "cluster.country.id", verbose = TRUE)
  
  gardner.tidy.10 <- did2s.tidy(gardner.did.10, buff = name.i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)),
           pre.period = -10)
  
  csa.did <- att_gt(data = buff.i, yname = "cumulative.forest.loss.perc", tname="year",
                    idname= "cluster.country.id", gname = "first.mine.year",
                    control_group = "notyettreated", base_period = "varying",
                    xformla = ~ country,
                    clustervars = "cluster.country.id", 
                    bstrap=T, cband=T)
  
  csa.tidy <- csa.tidy(csa.did, buff = name.i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)),
           pre.period = NA)
  
  # main + covariates
  covs.buff.i <- buff.i %>%
    left_join(covs.all) %>%
    group_by(country, CLUSTER_ID) %>%
    mutate(cluster.country.id = cur_group_id()) %>% ungroup()
  covs.buff.i <- scale.tidy(covs.buff.i)
  
  gardner.did.covar <- did2s(data = filter(covs.buff.i, rel.year.first >= -5), 
                             yname = "cumulative.forest.loss.perc", 
                       treatment = "treatment",
                       first_stage =  ~ 0 + slope.z + elevation.z + pop.density.z + travel.time.z +
                         i(country, year) | cluster.country.id + year, 
                       second_stage = ~ i(rel.year.first, ref = c(-1)),
                       cluster_var = "cluster.country.id", verbose = TRUE)
  
  gardner.tidy.covar <- did2s.tidy(gardner.did.covar, buff = name.i) %>%
    mutate(cluster.n = length(unique(covs.buff.i$cluster.country.id)),
           pre.period = -5)
  
  ## 10 
  gardner.did.covar.10 <- did2s(data = filter(covs.buff.i, rel.year.first >= -10), 
                             yname = "cumulative.forest.loss.perc", 
                             treatment = "treatment",
                             first_stage =  ~ 0 + slope.z + elevation.z + pop.density.z + travel.time.z +
                               i(country, year) | cluster.country.id + year, 
                             second_stage = ~ i(rel.year.first, ref = c(-1)),
                             cluster_var = "cluster.country.id", verbose = TRUE)
  
  gardner.tidy.covar.10 <- did2s.tidy(gardner.did.covar.10, buff = name.i) %>%
    mutate(cluster.n = length(unique(covs.buff.i$cluster.country.id)),
           pre.period = -10)
  
  
  did.i <- rbind(gardner.tidy, gardner.tidy.10, csa.tidy)
  ssa.did <- rbind(ssa.did, did.i)
  ssa.did.covar <- rbind(ssa.did.covar, gardner.tidy.covar, gardner.tidy.covar10)
  
  
}

save(ssa.did, file = "Outputs/DiD.tables/SSA.tidy.1km.RData")
save(ssa.did.covar, file = "Outputs/DiD.tables/SSA.tidy.1km.covar.RData")
