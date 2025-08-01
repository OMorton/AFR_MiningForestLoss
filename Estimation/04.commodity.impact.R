library(did)
library(did2s)
library(tidyverse)
library(sf)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting
# buffer rings
source("functions.R")

## Set up buffer forest data ---------------------------------------------------
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
               "20km" = ssa.20km)
ssa.df <- bind_rows(ssa.ls)

# All commodities likely mined per site ----------------------------------------
comm.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/country.mining.commodity.type/"
comm.dir <- data.frame(file = list.files(path = comm.path, pattern = "5km.buffer"))
comm.df <- data.frame()

# collate commodity data
for (i in 1:nrow(comm.dir)) {
  load(paste0(comm.path,comm.dir$file[[i]]))
  comm.df <- rbind(comm.df, country.mining.commodity.data.df)
}

comm.df <- comm.df %>% separate_longer_delim(material, ",") %>% distinct()

comm.cluster <- comm.df %>%
  group_by(country, CLUSTER_ID) %>%
  mutate(cluster.country.id = cur_group_id()) 
length(unique(comm.cluster$cluster.country.id))
comm.cluster %>% ungroup() %>%
  filter(!material %in% c("No mine in Maus", "Unknown")) %>% summarise(n_distinct(cluster.country.id))

top.comm <- comm.df %>% group_by(material) %>% tally() %>%
  filter(n > 29, material != "No mine in Maus")

comm.ls <- top.comm$material
buff.ls <- c(1, 5, 10, 20)
comm.did.out <- data.frame()

## NOTE small group n prevent inclusion of a country varying effect here.
for (t in c(-5, -10)) {
  for (i in 1:4) {
    for (c in 1:length(comm.ls)) {
      
      buff.i <- buff.ls[[i]]
      comm.c <- comm.ls[[c]]
      
      dat.i.c <- filter(ssa.df, buffer.size == buff.i)
      dat.i.c <- did.prep(dat.i.c, lead.time = -23, post.time = 23,
                          type = "loss", covariates = NULL) %>%
        left_join(comm.df, by = c("CLUSTER_ID", "country"),
                  relationship = "many-to-many") %>%
        filter(material == comm.c) %>%
        group_by(country, CLUSTER_ID) %>%
        mutate(cluster.country.id = cur_group_id(),
               country = as.factor(country))
      
      gardner.did <- did2s(data = filter(dat.i.c, rel.year.first >= t),
                           yname = "cumulative.forest.loss.perc", 
                           treatment = "treatment",
                           first_stage =  ~ 0 | cluster.country.id + year, 
                           second_stage = ~ i(rel.year.first, ref = c(-1)),
                           cluster_var = "cluster.country.id", verbose = TRUE)
      
      gardner.tidy.1km <- did2s.tidy(gardner.did, buff = buff.i) %>%
        mutate(cluster.n = length(unique(dat.i.c$cluster.country.id)),
               pre.period = t)
      
      if (t == -10) {
      csa.did <- att_gt(data = dat.i.c, yname = "cumulative.forest.loss.perc", tname="year",
                        idname= "cluster.country.id", gname = "first.mine.year",
                        control_group = "notyettreated", base_period = "varying",
                        xformla = NULL,
                        clustervars = "cluster.country.id", 
                        bstrap=T, cband=T)
      
      csa.tidy.1km <- csa.tidy(csa.did, buff = buff.i) %>%
        mutate(cluster.n = length(unique(dat.i.c$cluster.country.id)),
               pre.period = NA)
      }
      ## dont write out the CSA results twice for each of the Gardner t periods
      if(t == -5) {
      did.i.c <- gardner.tidy.1km %>% 
        mutate(buffer.size = paste0(buffer.size, "km"),
               commodity = comm.c)
      comm.did.out <- rbind(comm.did.out, did.i.c)
      } else {
        did.i.c <- rbind(gardner.tidy.1km, csa.tidy.1km) %>% 
          mutate(buffer.size = paste0(buffer.size, "km"),
                 commodity = comm.c)
        comm.did.out <- rbind(comm.did.out, did.i.c) 
      }
    }
  }
}

save(comm.did.out, file = "Outputs/DiD.tables/SSA.ALL.commodities.tidy.Jul25.RData")

# Primary commodities likely mined per site ------------------------------------
comm.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/country.mining.commodity.type/"
comm.dir <- data.frame(file = list.files(path = comm.path, pattern = "commodity.5km.buffer.RData"))
comm.df <- data.frame()

# collate commodity data
for (i in 1:nrow(comm.dir)) {
  load(paste0(comm.path,comm.dir$file[[i]]))
  comm.df <- rbind(comm.df, country.mining.commodity.data.df)
}

comm.df <- comm.df %>% separate_longer_delim(material, ",") %>% distinct()

comm.cluster <- comm.df %>%
  group_by(country, CLUSTER_ID) %>%
  mutate(cluster.country.id = cur_group_id()) 
length(unique(comm.cluster$cluster.country.id))
comm.cluster %>% ungroup() %>%
  filter(!material %in% c("No mine in Maus", "Unknown")) %>% summarise(n_distinct(cluster.country.id))

top.comm <- comm.df %>% group_by(material) %>% tally() %>%
  filter(n > 29, material != "No mine in Maus")

comm.ls <- top.comm$material
buff.ls <- c(1, 5, 10, 20)
comm.did.out <- data.frame()

## NOTE small group n prevent inclusion of a country varying effect here.
for (t in c(-5, -10)) {
  for (i in 1:4) {
    for (c in 1:length(comm.ls)) {
      
      buff.i <- buff.ls[[i]]
      comm.c <- comm.ls[[c]]
      
      dat.i.c <- filter(ssa.df, buffer.size == buff.i)
      dat.i.c <- did.prep(dat.i.c, lead.time = -23, post.time = 23,
                          type = "loss", covariates = NULL) %>%
        left_join(comm.df, by = c("CLUSTER_ID", "country"),
                  relationship = "many-to-many") %>%
        filter(material == comm.c) %>%
        group_by(country, CLUSTER_ID) %>%
        mutate(cluster.country.id = cur_group_id(),
               country = as.factor(country))
      
      gardner.did <- did2s(data = filter(dat.i.c, rel.year.first >= t),
                           yname = "cumulative.forest.loss.perc", 
                           treatment = "treatment",
                           first_stage =  ~ 0 | cluster.country.id + year, 
                           second_stage = ~ i(rel.year.first, ref = c(-1)),
                           cluster_var = "cluster.country.id", verbose = TRUE)
      
      gardner.tidy.1km <- did2s.tidy(gardner.did, buff = buff.i) %>%
        mutate(cluster.n = length(unique(dat.i.c$cluster.country.id)),
               pre.period = t)
      
      if (t == -10) {
        csa.did <- att_gt(data = dat.i.c, yname = "cumulative.forest.loss.perc", tname="year",
                          idname= "cluster.country.id", gname = "first.mine.year",
                          control_group = "notyettreated", base_period = "varying",
                          xformla = NULL,
                          clustervars = "cluster.country.id", 
                          bstrap=T, cband=T)
        
        csa.tidy.1km <- csa.tidy(csa.did, buff = buff.i) %>%
          mutate(cluster.n = length(unique(dat.i.c$cluster.country.id)),
                 pre.period = NA)
      }
      ## dont write out the CSA results twice for each of the Gardner t periods
      if(t == -5) {
        did.i.c <- gardner.tidy.1km %>% 
          mutate(buffer.size = paste0(buffer.size, "km"),
                 commodity = comm.c)
        comm.did.out <- rbind(comm.did.out, did.i.c)
      } else {
        did.i.c <- rbind(gardner.tidy.1km, csa.tidy.1km) %>% 
          mutate(buffer.size = paste0(buffer.size, "km"),
                 commodity = comm.c)
        comm.did.out <- rbind(comm.did.out, did.i.c) 
      }
    }
  }
}

save(comm.did.out, file = "Outputs/DiD.tables/SSA.PRIMARY.commodities.tidy.Jul25.RData")
