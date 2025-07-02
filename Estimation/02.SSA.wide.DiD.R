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
# loop accross pre-treatment periods, the Gardner method can be sensitive to this.
for (t in c(-5, -10)) {
  for (i in 1:5) {
    buff.i <- ssa.ls[[i]]
    name.i <- ssa.names[[i]]
    
    ## Gardner 2022
    gardner.did <- did2s(data = filter(buff.i, rel.year.first >= t),
                         yname = "cumulative.forest.loss.perc", 
                     treatment = "treatment",
                     first_stage =  ~ 0 | cluster.country.id + country^year, 
                     second_stage = ~ i(rel.year.first, ref = c(-1)),
                     cluster_var = "cluster.country.id", verbose = TRUE)
    
    gardner.tidy <- did2s.tidy(gardner.did, buff = name.i) %>%
      mutate(cluster.n = length(unique(buff.i$cluster.country.id)),
             pre.period = t)
    
    ## Gardner 2022 + covariates
    covs.buff.i <- buff.i %>%
      left_join(covs.all) %>%
      group_by(country, CLUSTER_ID) %>%
      mutate(cluster.country.id = cur_group_id()) %>% ungroup()
    covs.buff.i <- scale.tidy(covs.buff.i)
    
    gardner.did.covar <- did2s(data = filter(covs.buff.i, rel.year.first >= t), 
                               yname = "cumulative.forest.loss.perc", 
                               treatment = "treatment",
                               first_stage =  ~ 0 + slope.z + elevation.z + pop.density.z + travel.time.z +
                                 i(country, year) | cluster.country.id + year, 
                               second_stage = ~ i(rel.year.first, ref = c(-1)),
                               cluster_var = "cluster.country.id", verbose = TRUE)
    
    
    gardner.tidy.covar <- did2s.tidy(gardner.did.covar, buff = name.i) %>%
      mutate(cluster.n = length(unique(covs.buff.i$cluster.country.id)),
             pre.period = t)
  
    ## CSA runs optimally when provided with the whole panel.
    ## Thus we dont run at -5 and -10 periods.
    if(t == -10) {
    ## Callaway and Sant'Anna 2021
    csa.did <- att_gt(data = buff.i, yname = "cumulative.forest.loss.perc", tname="year",
                      idname= "cluster.country.id", gname = "first.mine.year",
                      control_group = "notyettreated", base_period = "varying",
                      xformla = ~ country,
                      clustervars = "cluster.country.id", 
                      bstrap=T, cband=T)
    
    csa.tidy.df <- csa.tidy(csa.did, buff = name.i) %>%
      mutate(cluster.n = length(unique(buff.i$cluster.country.id)),
             pre.period = NA)
    

    csa.did.covar <- att_gt(data = covs.buff.i, 
                      yname = "cumulative.forest.loss.perc", tname="year",
                      idname= "cluster.country.id", gname = "first.mine.year",
                      control_group = "notyettreated", base_period = "varying",
                      xformla = ~ country + slope.z + elevation.z + pop.density.z + travel.time.z,
                      clustervars = "cluster.country.id", 
                      bstrap=T, cband=T)
    
    csa.tidy.covar <- csa.tidy(csa.did.covar, buff = name.i) %>%
      mutate(cluster.n = length(unique(covs.buff.i$cluster.country.id)),
             pre.period = NA)
    }
    ## dont write out the CSA results twice for each of the Gardner t periods
    if(t == -5) {
      ssa.did <- rbind(ssa.did, gardner.tidy)
      ssa.did.covar <- rbind(ssa.did.covar, gardner.tidy.covar)
    } else {
      ssa.did <- rbind(ssa.did, gardner.tidy, csa.tidy.df)
      ssa.did.covar <- rbind(ssa.did.covar, gardner.tidy.covar, csa.tidy.covar)  
    }
    
  }
}

save(ssa.did, file = "Outputs/DiD.tables/SSA.tidy.1km.RData")
save(ssa.did.covar, file = "Outputs/DiD.tables/SSA.tidy.1km.covar.RData")
