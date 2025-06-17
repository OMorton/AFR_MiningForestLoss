library(did)
library(did2s)
library(tidyverse)
library(ggpubr)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting buffer rings
source("functions.R")

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
j <- 3
t <- -10
## Estimate DiD ----------------------------------------------------------------

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
    # prep for DiD
    i.did.dat <- did.prep(dir.indir.loss,
                          lead.time = -23, post.time = 23,
                          type = "loss", covariates = covs.i)
    
    # try statement as some countries cannot isolate a direct effect of mining as the 
    # response ends up constant
    i.did <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                        p.value  = NA, method = NA, y.var = NA)
    try(i.did <- fit.dynamic.DiD(i.did.dat,
                                 yname = c("direct.forest.loss","indirect.forest.loss",
                                           "direct.perc", "indirect.perc"),
                                 xformula = NULL, 
                                 method = c("gardner"), GAR.pre.period = t), silent = TRUE)
    
    
    i.did$buffer.size <- "10km.master"
    i.did$country <- country.j
    i.did$cluster.n <- length(i.did.dat$CLUSTER_ID)
    i.did$unique.trt.yrs <- length(unique(i.did.dat$first.mine.year))
    i.did$zeroes <- zeroes
    i.did$pre.length <- t
    did.all <- rbind(did.all, i.did)
  }
}
# gardner.ls <- did.all
# csa.ls <- did.all
# dir.indir <- rbind(gardner.ls, csa.ls)

save(did.all, file = "Outputs/DiD.tables/all.GAR.direct.indirect.DiD.RData")






## Estimate DiD + Spillover ----------------------------------------------------------------

did.spill.ratio <- data.frame()


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
  # prep for DiD
  i.did.dat <- did.prep(dir.indir.loss,
                        lead.time = -23, post.time = 23,
                        type = "loss", covariates = covs.i)
  
  i.spill.dat <- spill.tidy(i.did.dat, buff = i, concentric.ring = "fixed")
    # try statement as some countries cannot isolate a direct effect of mining as the 
  # response ends up constant
  i.spill.DiD.direct <- data.frame(effect = NA, year.since = NA, estimate = NA, lci = NA, uci  = NA,
                            p.value  = NA, method = NA, buffer.size = NA, y.var = NA)
  i.spill.DiD.indirect <- data.frame(effect = NA, year.since = NA, estimate = NA, lci = NA, uci  = NA,
                                   p.value  = NA, method = NA, buffer.size = NA, y.var = NA)
  try(i.spill.DiD.direct <- spillover.dynamic.DiD(i.spill.dat,
                                           yname = "direct.forest.loss"),
      silent = T)
  try(i.spill.DiD.indirect <- spillover.dynamic.DiD(i.spill.dat,
                                                  yname = "indirect.forest.loss"),
      silent = T)
  i.spill.DiD.direct$y.var <- "direct.forest.loss"
  i.spill.DiD.indirect$y.var <- "indirect.forest.loss"
  
  try(i.spill.DiD <- rbind(i.spill.DiD.direct, i.spill.DiD.indirect))
  
  i.spill.DiD$buffer.size <- i
  i.spill.DiD$country <- country.j
  i.spill.DiD$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
  i.spill.DiD$spill.cluster.n <- length(unique(filter(i.spill.dat, 
                                                      !is.na(first.spillover.i))$CLUSTER_ID))
  i.spill.DiD$unique.trt.yrs <- length(unique(i.spill.dat$first.mine.year))
  i.spill.DiD$zeroes <- zeroes
  
  did.spill.ratio <- rbind(did.spill.ratio, i.spill.DiD)
}


save(did.spill.ratio, file = "Outputs/DiD.tables/all.direct.indirect.withspillover.DiD.RData")

indir.sum <- did.spill.ratio %>% 
  filter(y.var %in% c("indirect.forest.loss"), year.since == 5, 
         method == "Gardner 2022", effect == "Mine effect") %>%
  rename("indirect.estimate" = "estimate", "indirect.uci" = "uci",
         "indirect.lci" = "lci") %>% 
  select(country, indirect.estimate, indirect.uci, indirect.lci) %>% add.iso()

dir.sum <- did.spill.ratio %>% 
  filter(y.var %in% c("direct.forest.loss"),
         year.since == 5, method == "Gardner 2022", effect == "Mine effect") %>%
  rename("direct.estimate" = "estimate", "direct.uci" = "uci",
         "direct.lci" = "lci") %>% 
  select(country, direct.estimate, direct.uci, direct.lci)

ratio.all.sum <- left_join(indir.sum, dir.sum, by = "country") %>% 
  mutate(ratio.est = indirect.estimate/direct.estimate) %>%
  group_by(country) %>%
  mutate( # focus only where there are average increases
         ratio.est = ifelse(indirect.estimate < 0 | direct.estimate < 0, NA,
                            ratio.est))

bar.ratio <- ggplot(ratio.all.sum, aes(reorder(iso3, -ratio.est), ratio.est)) +
  geom_col() +
  xlab("Country") +
  ylab("Ratio indirect to direct forest loss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

## Estimate DiD + No Overlap ---------------------------------------------------

did.no.overlap.ratio <- data.frame()

j <- 6
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
  # prep for DiD
  i.did.dat <- did.prep(dir.indir.loss,
                        lead.time = -23, post.time = 23,
                        type = "loss", covariates = covs.i)
  
  i.nospill.dat <- i.did.dat %>% 
    filter(is.na(spillover.10km.master))  # try statement as some countries cannot isolate a direct effect of mining as the 
  # response ends up constant
  i.spill.DiD.direct <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                                   p.value  = NA, method = NA, buffer.size = NA, y.var = NA)
  i.spill.DiD.indirect <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                                     p.value  = NA, method = NA, buffer.size = NA, y.var = NA)
  try(i.spill.DiD.direct <- fit.dynamic.DiD(i.nospill.dat,
                                            yname = "direct.forest.loss",
                                            xformula = NULL, 
                                            method = c("gardner")),
      silent = T)
  try(i.spill.DiD.indirect <- fit.dynamic.DiD(i.nospill.dat,
                                              yname = "indirect.forest.loss",
                                              xformula = NULL, 
                                              method = c("gardner")),
      silent = T)
  i.spill.DiD.direct$y.var <- "direct.forest.loss"
  i.spill.DiD.indirect$y.var <- "indirect.forest.loss"
  
  try(i.spill.DiD <- rbind(i.spill.DiD.direct, i.spill.DiD.indirect))
  
  i.spill.DiD$buffer.size <- i
  i.spill.DiD$country <- country.j
  i.spill.DiD$cluster.n <- length(unique(i.did.dat$CLUSTER_ID))
  i.spill.DiD$no.overlap.cluster.n <- length(unique(i.nospill.dat$CLUSTER_ID))
  i.spill.DiD$unique.trt.yrs <- length(unique(i.spill.dat$first.mine.year))
  i.spill.DiD$zeroes <- zeroes
  
  did.no.overlap.ratio <- rbind(did.no.overlap.ratio, i.spill.DiD)
}


save(did.no.overlap.ratio, file = "Outputs/DiD.tables/all.direct.indirect.no.overlap.DiD.RData")

indir.sum <- did.no.overlap.ratio %>% 
  filter(y.var %in% c("indirect.forest.loss"), year.since == 5, 
         method == "Gardner 2022") %>%
  rename("indirect.estimate" = "estimate", "indirect.uci" = "uci",
         "indirect.lci" = "lci") %>% 
  select(country, indirect.estimate, indirect.uci, indirect.lci) %>% add.iso()

dir.sum <- did.no.overlap.ratio %>% 
  filter(y.var %in% c("direct.forest.loss"),
         year.since == 5, method == "Gardner 2022") %>%
  rename("direct.estimate" = "estimate", "direct.uci" = "uci",
         "direct.lci" = "lci") %>% 
  select(country, direct.estimate, direct.uci, direct.lci)

ratio.all.sum <- left_join(indir.sum, dir.sum, by = "country") %>% 
  mutate(ratio.est = indirect.estimate/direct.estimate) %>%
  group_by(country) %>%
  mutate( # focus only where there are average increases
    ratio.est = ifelse(indirect.estimate < 0 | direct.estimate < 0, NA,
                       ratio.est))

bar.ratio <- ggplot(ratio.all.sum, aes(reorder(iso3, -ratio.est), ratio.est)) +
  geom_col() +
  xlab("Country") +
  ylab("Ratio indirect to direct forest loss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))