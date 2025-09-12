library(did)
library(did2s)
library(tidyverse)
library(ggpubr)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting buffer rings
source("functions.R")



## Direct and indirect 0-5 km buffer -------------------------------------------

# extract all the data from the file naming
# data stored on the X: drive
all.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
min.only.path <-  "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers.miningonly/"
cov.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/model.covariates/"
file.dir <- data.frame(file = list.files(path = all.path, pattern = "5km.master")) %>% 
  mutate(country = sub(".forest.*", "", file),
         buffer = "5km.master")%>%
  mutate(mining.only.file = list.files(path = min.only.path, pattern = "5km.master"))

forest.clusters <- read.csv("Outputs/third.forest.clusters.csv")

country.ls <- unique(file.dir$country)
did.all <- list()
j <- 6
t <- -10

for (t in c(-5, -10)) {
  for (j in 1:length(country.ls)) {
    cat(j, "out of", length(country.ls), "\n")
    
    country.j <- country.ls[j]
    nat.buff <- file.dir %>% filter(country == country.j)
    
    # get files for all forest loss
    load(paste0(all.path, nat.buff$file))
    mine.all.forest.loss.df <- mine.buffer.forest.loss.df %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             forest.area.2000, cumulative.forest.loss.area)
    
    # get files for mining only forest loss
    load(paste0(min.only.path, nat.buff$mining.only.file))
    mine.direct.forest.loss.df <- mine.buffer.forest.loss.df %>%
      rename("direct.forest.loss" = "cumulative.forest.loss.area") %>%
      mutate(buffer.size = "5km master") %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             direct.forest.loss)
    
    # # get files for mining only forest loss
    # load(paste0(min.only.path, nat.buff$mining.only.file.1km))
    # direct.1km <- mine.buffer.forest.loss.df %>%
    #   select(CLUSTER_ID, loss.year, cumulative.forest.loss, forest.cells.2000) %>%
    #   rename("cumulative.1km" = "cumulative.forest.loss", 
    #          "forest.cells.2000.1km" = "forest.cells.2000")
    # load(paste0(min.only.path, nat.buff$mining.only.file.5km))
    # direct.5km <- mine.buffer.forest.loss.df
    # direct.5km.master <- left_join(direct.5km, direct.1km) %>%
    #   mutate(cumulative.forest.loss = cumulative.forest.loss + cumulative.1km,
    #          buffer.size = "5km Master")
    # mine.direct.forest.loss.df <- direct.5km.master %>%
    #   rename("direct.forest.loss" = "cumulative.forest.loss") %>%
    #   select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
    #          direct.forest.loss)
    forest.clusters.j <- forest.clusters %>% filter(country == country.j)
    
    
    dir.indir.loss <- left_join(mine.all.forest.loss.df, mine.direct.forest.loss.df,
                                by = c("CLUSTER_ID", "first.mine.year", "buffer.size", "loss.year")) %>%
      # drop years after 2020 as beyond Masosele et al.
      filter(!is.na(direct.forest.loss)) %>%
      mutate(direct.forest.loss = ifelse(direct.forest.loss>cumulative.forest.loss.area,
                                         cumulative.forest.loss.area, direct.forest.loss),
             indirect.forest.loss = cumulative.forest.loss.area - direct.forest.loss,
             indirect.perc = indirect.forest.loss/forest.area.2000 *100,
             direct.perc = direct.forest.loss/forest.area.2000 *100) %>%
      filter(CLUSTER_ID %in% forest.clusters.j$forest.clusters)
    
    load(paste0(cov.path, nat.buff$country, ".50p.forest.final.COVARIATES.RData"))
    covs.i <- 
      select(mine.cluster.points.df, -first.mine.year, -buffer.forest.prop)
    
    # check for any zero cover
    zeroes <- "No"
    if (any(dir.indir.loss$forest.area.2000 ==0)) {
      cat("Warning - 0 forest cover in 2000 detected for", nat.buff$country, "\n")
      dir.indir.loss <- 
        dir.indir.loss %>% filter(forest.area.2000>0)
      zeroes <- "Yes"
    }
    # prep for DiD
    # prep for DiD - using first 10% as start of mining
    i.did.dat <- dir.indir.loss %>% 
      left_join(covs.i, by = "CLUSTER_ID") %>%
      select(-first.mine.year) %>%
      mutate(first.mine.year = year.10p) %>%
      filter(loss.year <= 2020) %>%
      mutate(year = loss.year - 2000,
             treatment = ifelse(year >= first.mine.year, TRUE, FALSE),
             treatment.n = ifelse(year >= first.mine.year, 1, 0),
             rel.year.first = year - first.mine.year) %>%
      filter(rel.year.first >= -23 & rel.year.first <= 23) %>% 
      left_join(covs.i, by = "CLUSTER_ID")
    
    # try statement as some countries cannot isolate a direct effect of mining as the 
    # response ends up constant
    i.did <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                        p.value  = NA, method = NA, y.var = NA)
    try(i.did <- fit.dynamic.DiD(i.did.dat,
                                 yname = c("direct.forest.loss","indirect.forest.loss",
                                           "direct.perc", "indirect.perc"),
                                 xformula = NULL, 
                                 method = c("gardner"), GAR.pre.period = t), silent = TRUE)
    
    
    if (t == -10) {
    # try statement as some countries cannot isolate a direct effect of mining as the 
    # response ends up constant
    i.did.csa <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                        p.value  = NA, method = NA, y.var = NA)
    try(i.did.csa <- fit.dynamic.DiD(i.did.dat,
                                 yname = c("direct.forest.loss","indirect.forest.loss",
                                           "direct.perc", "indirect.perc"),
                                 xformula = NULL, 
                                 method = c("csa"), GAR.pre.period = NA), silent = TRUE)
    }
    
    
    ## dont write out the CSA results twice for each of the Gardner t periods
    if(t == -5) {
      i.did$buffer.size <- "5km.master"
      i.did$country <- country.j
      i.did$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$zeroes <- zeroes
      i.did$pre.length <- t
      did.all <- rbind(did.all, i.did)
    } else {
      i.did.csa$buffer.size <- "5km.master"
      i.did.csa$country <- country.j
      i.did.csa$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did.csa$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did.csa$zeroes <- zeroes
      i.did.csa$pre.length <- NA
      
      i.did$buffer.size <- "5km.master"
      i.did$country <- country.j
      i.did$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$zeroes <- zeroes
      i.did$pre.length <- t
      did.all <- rbind(did.all, i.did, i.did.csa)  
    }
  }
}

save(did.all, file = "Outputs/DiD.tables/all.direct.indirect.AREA.DiD.5km.Sept25.10p.PlantationsRem.RData")

## SSA-Wide: Direct and indirect 0-5 km buffer ---------------------------------

# extract all the data from the file naming
# data stored on the X: drive
all.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
min.only.path <-  "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers.miningonly/"
cov.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/model.covariates/"
file.indir <- data.frame(file = list.files(path = all.path, pattern = "5km.master")) %>% 
  mutate(country = sub(".forest.*", "", file),
         buffer = "5km.master")
file.dir <- data.frame(file = list.files(path = min.only.path, pattern = "5km.master")) %>% 
  mutate(country = sub(".forest.*", "", file),
         buffer = "5km.master")

forest.clusters <- read.csv("Outputs/third.forest.clusters.csv")
forest.clusters <- forest.clusters %>% mutate(forest.id = paste0(country, ".", forest.clusters))

cov.dir <- data.frame(file = list.files(path = cov.path, full.names = TRUE))
covs.all <- read.covs.list(cov.dir) %>%
  mutate(country = ifelse(country =="Côte d'Ivoire", "Côte d_Ivoire", country))

ssa.all.raw  <- data.frame()

for (i in 1:29) {
  file <- file.indir$file[i]
  # Read the text file as a data frame
  load(paste0(all.path,file))  # Change the arguments based on your file format
  mine.buffer.forest.loss.df$country <- file.indir$country[[i]]
  # Bind the data to the combined data frame
  ssa.all.raw  <- rbind(ssa.all.raw , mine.buffer.forest.loss.df)
}
ssa.all.raw  <- ssa.all.raw  %>% filter(loss.year <2021)

ssa.dir.raw  <- data.frame()
for (i in 1:29) {
  file <- file.dir$file[i]
  # Read the text file as a data frame
  load(paste0(min.only.path,file))  # Change the arguments based on your file format
  mine.buffer.forest.loss.df$country <- file.dir$country[[i]]
  # Bind the data to the combined data frame
  ssa.dir.raw  <- rbind(ssa.dir.raw , mine.buffer.forest.loss.df)
}



SSA.did.all <- list()
t <- -10

for (t in c(-5, -10)) {

    # get files for all forest loss
    ssa.all <- ssa.all.raw %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             forest.area.2000, cumulative.forest.loss.area, country)
    
    # get files for mining only forest loss
    ssa.dir <- ssa.dir.raw  %>%
      rename("direct.forest.loss" = "cumulative.forest.loss.area") %>%
      mutate(buffer.size = "5km master") %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             direct.forest.loss, country)
    
    # # get files for mining only forest loss
    # load(paste0(min.only.path, nat.buff$mining.only.file.1km))
    # direct.1km <- mine.buffer.forest.loss.df %>%
    #   select(CLUSTER_ID, loss.year, cumulative.forest.loss, forest.cells.2000) %>%
    #   rename("cumulative.1km" = "cumulative.forest.loss", 
    #          "forest.cells.2000.1km" = "forest.cells.2000")
    # load(paste0(min.only.path, nat.buff$mining.only.file.5km))
    # direct.5km <- mine.buffer.forest.loss.df
    # direct.5km.master <- left_join(direct.5km, direct.1km) %>%
    #   mutate(cumulative.forest.loss = cumulative.forest.loss + cumulative.1km,
    #          buffer.size = "5km Master")
    # mine.direct.forest.loss.df <- direct.5km.master %>%
    #   rename("direct.forest.loss" = "cumulative.forest.loss") %>%
    #   select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
    #          direct.forest.loss)
    
    dir.indir.loss <- left_join(ssa.all, ssa.dir,
                                by = c("CLUSTER_ID", "first.mine.year", "buffer.size", "loss.year", "country")) %>%
      mutate(direct.forest.loss = ifelse(direct.forest.loss>cumulative.forest.loss.area,
                                         cumulative.forest.loss.area, direct.forest.loss),
             indirect.forest.loss = cumulative.forest.loss.area - direct.forest.loss,
             indirect.perc = indirect.forest.loss/forest.area.2000 *100,
             direct.perc = direct.forest.loss/forest.area.2000 *100,
             forest.id = paste0(country, ".", CLUSTER_ID)) %>%
      filter(forest.id %in% forest.clusters$forest.id)
    
    # load(paste0(cov.path, nat.buff$country, ".50p.forest.final.COVARIATES.RData"))
    # covs.i <- 
    #   select(mine.cluster.points.df, -first.mine.year, -buffer.forest.prop)
    # 
    # # check for any zero cover
    # zeroes <- "No"
    # if (any(dir.indir.loss$forest.cells.2000 ==0)) {
    #   cat("Warning - 0 forest cover in 2000 detected for", nat.buff$country, "\n")
    #   dir.indir.loss <- 
    #     dir.indir.loss %>% filter(forest.cells.2000>0)
    #   zeroes <- "Yes"
    # }
    # prep for DiD
    i.did.dat <- dir.indir.loss %>%    
      select(-first.mine.year) %>%
      left_join(covs.all, by = c("country", "CLUSTER_ID")) %>%
      mutate(first.mine.year = year.10p) %>%
      #mining data only runs to 2020
      filter(loss.year <= 2020) %>%
      mutate(year = loss.year - 2000,
             treatment = ifelse(year >= first.mine.year, TRUE, FALSE),
             treatment.n = ifelse(year >= first.mine.year, 1, 0),
             rel.year.first = year - first.mine.year) %>%
      filter(rel.year.first >= -23 & rel.year.first <= 23) %>% 
      group_by(country, CLUSTER_ID) %>%
      mutate(cluster.country.id = cur_group_id()) %>% ungroup()
    
    # try statement as some countries cannot isolate a direct effect of mining as the 
    # response ends up constant
    i.did <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                        p.value  = NA, method = NA, y.var = NA)
    try(i.did <- fit.dynamic.DiD.SSA(i.did.dat,
                                 yname = c("direct.forest.loss","indirect.forest.loss",
                                           "direct.perc", "indirect.perc"),
                                 xformula = NULL, 
                                 method = c("gardner"), GAR.pre.period = t), silent = TRUE)
    
    
    if (t == -10) {
      # try statement as some countries cannot isolate a direct effect of mining as the 
      # response ends up constant
      i.did.csa <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                              p.value  = NA, method = NA, y.var = NA)
      try(i.did.csa <- fit.dynamic.DiD.SSA(i.did.dat,
                                       yname = c("direct.forest.loss","indirect.forest.loss",
                                                 "direct.perc", "indirect.perc"),
                                       xformula = ~ country, 
                                       method = c("csa"), GAR.pre.period = NA), silent = TRUE)
    }
    
    
    ## dont write out the CSA results twice for each of the Gardner t periods
    if(t == -5) {
      i.did$buffer.size <- "5km.master"
      i.did$cluster.n <- length(unique(i.did.dat$cluster.country.id))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$pre.length <- t
      SSA.did.all <- rbind(SSA.did.all, i.did)
    } else {
      i.did.csa$buffer.size <- "5km.master"
      i.did.csa$cluster.n <- length(unique(i.did.dat$cluster.country.id))
      i.did.csa$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did.csa$pre.length <- NA
      
      i.did$buffer.size <- "5km.master"
      i.did$cluster.n <- length(unique(i.did.dat$cluster.country.id))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$pre.length <- t
      SSA.did.all <- rbind(SSA.did.all, i.did, i.did.csa)  
    }
  }


save(SSA.did.all, file = "Outputs/DiD.tables/SSA.direct.indirect.AREA.DiD.5km.Sept25.10p.PlantationsRem.RData")


## Direct and indirect 0-10 km buffer - Not Used -------------------------------------------
# extract all the data from the file naming
# data stored on the X: drive
all.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
min.only.path <-  "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers.miningonly/"
cov.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/model.covariates/"
file.dir <- data.frame(file = list.files(path = all.path, pattern = "10km.master")) %>% 
  mutate(country = sub(".forest.*", "", file),
         buffer = "10km.master")%>%
  mutate(mining.only.file = list.files(path = min.only.path, pattern = "10km.master"))

country.ls <- unique(file.dir$country)
did.all <- list()
j <- 6
t <- -10


for (t in c(-5, -10)) {
  for (j in 1:length(country.ls)) {
    cat(j, "out of", length(country.ls), "\n")
    
    country.j <- country.ls[j]
    nat.buff <- file.dir %>% filter(country == country.j)
    
    # get files for all forest loss
    load(paste0(all.path, nat.buff$file))
    mine.all.forest.loss.df <- mine.buffer.forest.loss.df %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             forest.cells.2000, cumulative.forest.loss, cumulative.forest.loss.prop)
    # get files for mining only forest loss
    load(paste0(min.only.path, nat.buff$mining.only.file))
    mine.direct.forest.loss.df <- mine.buffer.forest.loss.df %>%
      rename("direct.forest.loss" = "cumulative.forest.loss") %>%
      select(CLUSTER_ID, first.mine.year, buffer.size, loss.year, 
             direct.forest.loss)
    
    dir.indir.loss <- left_join(mine.all.forest.loss.df, mine.direct.forest.loss.df,
                                by = c("CLUSTER_ID", "first.mine.year", "buffer.size", "loss.year")) %>%
      # drop years after 2020 as beyond Masosele et al.
      filter(!is.na(direct.forest.loss)) %>%
      mutate(direct.forest.loss = ifelse(direct.forest.loss>cumulative.forest.loss,
                                         cumulative.forest.loss, direct.forest.loss),
             indirect.forest.loss = cumulative.forest.loss - direct.forest.loss,
             indirect.perc = indirect.forest.loss/forest.cells.2000 *100,
             direct.perc = direct.forest.loss/forest.cells.2000 *100)
    
    load(paste0(cov.path, nat.buff$country, ".50p.forest.final.COVARIATES.RData"))
    covs.i <- 
      select(mine.cluster.points.df, -first.mine.year, -buffer.forest.prop)
    
    # check for any zero cover
    zeroes <- "No"
    if (any(dir.indir.loss$forest.cells.2000 ==0)) {
      cat("Warning - 0 forest cover in 2000 detected for", nat.buff$country, "\n")
      dir.indir.loss <- 
        dir.indir.loss %>% filter(forest.cells.2000>0)
      zeroes <- "Yes"
    }
    # prep for DiD - using first 10% as start of mining
    dir.indir.loss <- dir.indir.loss %>% left_join(covs.i, by = "CLUSTER_ID") %>%
      select(-first.mine.year) %>%
      mutate(first.mine.year = year.10p)
    i.did.dat <- did.prep(dir.indir.loss,
                          lead.time = -23, post.time = 23,
                          type = "loss", covariates = NULL)
    
    # try statement as some countries cannot isolate a direct effect of mining as the 
    # response ends up constant
    i.did <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                        p.value  = NA, method = NA, y.var = NA)
    try(i.did <- fit.dynamic.DiD(i.did.dat,
                                 yname = c("direct.forest.loss","indirect.forest.loss",
                                           "direct.perc", "indirect.perc"),
                                 xformula = NULL, 
                                 method = c("gardner"), GAR.pre.period = t), silent = TRUE)
    
    ## CSA runs optimally when provided with the whole panel.
    ## Thus we dont run at -5 and -10 periods.
    if (t == -10) {
      i.did.csa <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                              p.value  = NA, method = NA, y.var = NA)
      try(i.did.csa <- fit.dynamic.DiD(i.did.dat,
                                       yname = c("direct.forest.loss","indirect.forest.loss",
                                                 "direct.perc", "indirect.perc"),
                                       xformula = NULL, 
                                       method = c("csa"), GAR.pre.period = NA), silent = TRUE)
    }
    
    ## dont write out the CSA results twice for each of the Gardner t periods
    if(t == -5) {
      i.did$buffer.size <- "10km.master"
      i.did$country <- country.j
      i.did$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$zeroes <- zeroes
      i.did$pre.length <- t
      did.all <- rbind(did.all, i.did)
    } else {
      i.did.csa$buffer.size <- "10km.master"
      i.did.csa$country <- country.j
      i.did.csa$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did.csa$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did.csa$zeroes <- zeroes
      i.did.csa$pre.length <- NA
      
      i.did$buffer.size <- "10km.master"
      i.did$country <- country.j
      i.did$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
      i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
      i.did$zeroes <- zeroes
      i.did$pre.length <- t
      did.all <- rbind(did.all, i.did, i.did.csa)  
    }
  }
}

save(did.all, file = "Outputs/DiD.tables/all.direct.indirect.DiD.10km.Aug25.10p.RData")
