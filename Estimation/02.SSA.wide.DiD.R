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
covs.all <- read.covs.list(cov.dir)

ssa.covs <- did.prep(ssa.5km.master,
         lead.time = -23, post.time = 23,
         type = "loss", covariates = NULL) %>%
  left_join(covs.all) %>%
  group_by(country, CLUSTER_ID) %>%
  mutate(cluster.country.id = cur_group_id()) %>%
  filter(is.na(first.spillover))

length(unique(ssa.covs$cluster.country.id))
length(unique(ssa.covs$CLUSTER_ID))

gardner.did <- did2s(data = ssa.covs, yname = "cumulative.forest.loss.perc", 
                     treatment = "treatment",
                     first_stage =  ~ 0 + country | cluster.country.id + year, 
                     second_stage = ~ i(rel.year.first, ref = c(-1)),
                     cluster_var = "cluster.country.id", verbose = TRUE)

gardner.tidy.1km <- did2s.tidy(gardner.did, buff = name.i) %>%
  mutate(cluster.n = length(unique(buff.i$cluster.country.id)))

## SSA wide --------------------------------------------------------------------

ssa.ls <- lapply(ssa.ls, function(x) {did.prep(x,
                                     lead.time = -23, post.time = 23,
                                     type = "loss", covariates = NULL) %>%
                          group_by(country, CLUSTER_ID) %>%
                          mutate(cluster.country.id = cur_group_id())})



ssa.names <- names(ssa.ls)
ssa.did <- data.frame()
ssa.grp.did <- data.frame()
ssa.grp.time.did <- data.frame()
i <- 5
for (i in 1:5) {
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
  
  ## temporal
  grp.est <- aggte(csa.did, type = "group", na.rm = TRUE)
  grp.out <- data.frame(t = grp.est$egt, 
                        estimate = grp.est$att.egt, 
                        lci = grp.est$att.egt - grp.est$crit.val.egt*grp.est$se.egt,
                        uci = grp.est$att.egt + grp.est$crit.val.egt*grp.est$se.egt,
                        p.value = NA, method = "Callaway & Sant'Anna 2021", buffer.size = i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)))
  
  grp.time.out <- data.frame(grp = csa.did$group, t = csa.did$t, 
                             estimate = csa.did$att, 
                             lci = csa.did$att - csa.did$c*csa.did$se,
                             uci = csa.did$att + csa.did$c*csa.did$se,
                             p.value = NA, method = "Callaway & Sant'Anna 2021", buffer.size = i) %>%
    mutate(cluster.n = length(unique(buff.i$cluster.country.id)))
  
  ssa.grp.did <- rbind(ssa.grp.did, grp.out)
  ssa.grp.time.did <- rbind(ssa.grp.time.did, grp.time.out)
  
  
}

save(ssa.did, file = "Outputs/DiD.tables/SSA.tidy.1km.RData")
save(ssa.grp.did, file = "Outputs/DiD.tables/SSA.tidy.grp.RData")
save(ssa.grp.time.did, file = "Outputs/DiD.tables/SSA.tidy.grp.time.RData")



