library(did)
library(did2s)
library(tidyverse)
library(ggpubr)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting buffer rings
source("functions.R")

# extract all the data from the file naming
# data stored on the X: drive
dir.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
cov.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/model.covariates/"
file.dir <- data.frame(file = list.files(path = dir.path))
file.dir <- file.dir %>% mutate(country = sub(".forest.*", "", file),
                                buffer = sub(".buffer.forest.loss.df.RData", "", file),
                                buffer = sub(".forest.mines.", "", buffer),
                                buffer = str_remove(buffer, country)) %>%
  ## remove low n
  filter(!country %in% c("South Sudan", "Comoros", "Burundi", "Malawi", "Rwanda"))

country.ls <- unique(file.dir$country)
did.ls <- list()
did.covar.ls <- list()

j <- 10
i <- "1km"
t <--10

## Estimate DiD ----------------------------------------------------------------

for (t in c(-5, -10)) {
  for (j in 1:length(country.ls)) {
    cat(j, "out of", length(country.ls), "\n")
    for (i in c("1km", "2.5km", "5km", "10km", "20km", "5km.master")) {
      
      # get ixj combo
      country.j <- country.ls[j]
      nat.buff <- file.dir %>% filter(country == country.j, buffer == i)
      # get files
      load(paste0(dir.path, nat.buff$file))
      load(paste0(cov.path, nat.buff$country, ".50p.forest.final.COVARIATES.RData"))
      covs.i <- 
        select(mine.cluster.points.df, -first.mine.year, -buffer.forest.prop)
    
      # check for any zero cover
      zeroes <- "No"
      if (any(mine.buffer.forest.loss.df$forest.cells.2000 ==0)) {
        cat("Warning - 0 forest cover in 2000 detected for", nat.buff$country, "\n")
        mine.buffer.forest.loss.df <- 
          mine.buffer.forest.loss.df %>% filter(forest.cells.2000>0)
        zeroes <- "Yes"
      }
      # prep for DiD
      i.did.dat <- did.prep(mine.buffer.forest.loss.df,
                            lead.time = -23, post.time = 23,
                              type = "loss", covariates = covs.i)
      
      i.did.fit <- fit.dynamic.DiD(i.did.dat,
                                   yname = "cumulative.forest.loss.perc",
                                   xformula = NULL, 
                                   method = c("csa", "gardner"), GAR.pre.period = t)
      
      # only run the csa covar once per loop (t is irrelevant)
      if(t == -10) {
      i.did.dat <- scale.tidy(i.did.dat)
      
      i.did.covar.fit <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                          p.value  = NA, method = NA, y.var = NA)
      try(i.did.covar.fit <- fit.dynamic.DiD(i.did.dat,
                      yname = "cumulative.forest.loss.perc",
                      xformula = ~ slope.z + elevation.z + pop.density.z + travel.time.z, 
                      method = c("csa"), GAR.pre.period = t),
          silent = TRUE)
      } else {
        i.did.dat <- scale.tidy(i.did.dat)
        
        i.did.covar.fit <- data.frame(year.since = NA, estimate = NA, lci = NA, uci  = NA,
                                      p.value  = NA, method = NA, y.var = NA)
        try(i.did.covar.fit <- fit.dynamic.DiD(i.did.dat,
                                               yname = "cumulative.forest.loss.perc",
                                               xformula = ~ slope.z + elevation.z + pop.density.z + travel.time.z, 
                                               method = c("gardner"), GAR.pre.period = t),
            silent = TRUE)
      }
      i.did.fit$buffer.size <- i
      i.did.fit$country <- country.j
      i.did.fit$cluster.n <- length(mine.cluster.points.df$CLUSTER_ID)
      i.did.fit$unique.trt.yrs <- length(unique(mine.cluster.points.df$first.mine.year))
      i.did.fit$zeroes <- zeroes
      i.did.fit$pre.period <- t
      if (!is.null(i.did.covar.fit)) {
        i.did.covar.fit$buffer.size <- i
        i.did.covar.fit$country <- country.j
        i.did.covar.fit$cluster.n <- length(mine.cluster.points.df$CLUSTER_ID)
        i.did.covar.fit$unique.trt.yrs <- length(unique(mine.cluster.points.df$first.mine.year))
        i.did.covar.fit$zeroes <- zeroes
        i.did.covar.fit$pre.period <- t
      }
      did.ls[[paste0(country.j, ".", i, ".pre", t)]] <- i.did.fit 
      did.covar.ls[[paste0(country.j, ".", i, ".pre", t)]] <- i.did.covar.fit 
      
    }
  }
}

save(did.ls, file = "Outputs/DiD.tables/all.DiD.all.t.Jul25.RData")
save(did.covar.ls, file = "Outputs/DiD.tables/all.covars.DiD.all.t.Jul25.RData")

## Per country plot checking ---------------------------------------------------------------
load("Outputs/DiD.tables/all.not20km.DiD.RData")

plot.ls <- lapply(did.ls, plot.DiD2S.f, pre.yrs = -5, post.yrs = 10, legend = "none",
       method = c("csa", "gardner"),
       pre.mine.col = c("#313695", "#abd9e9"),
       post.mine.col = c("#fdae61", "#d73027"),
       y.label = "Forest loss (% points)")

for (j in 1:length(country.ls)) {
  country.j <- country.ls[j]
  j.ls <- plot.ls[grep(country.j, names(plot.ls))]
  j.ls.names <- names(j.ls)
  
  j.arr <-  ggarrange(j.ls[[1]]$`Gardner 2022`, j.ls[[1]]$`Callaway & Sant'Anna 2021`,
            j.ls[[2]]$`Gardner 2022`, j.ls[[2]]$`Callaway & Sant'Anna 2021`,
            j.ls[[3]]$`Gardner 2022`, j.ls[[3]]$`Callaway & Sant'Anna 2021`,
            j.ls[[4]]$`Gardner 2022`, j.ls[[4]]$`Callaway & Sant'Anna 2021`,
            j.ls[[5]]$`Gardner 2022`, j.ls[[5]]$`Callaway & Sant'Anna 2021`,
            j.ls[[6]]$`Gardner 2022`, j.ls[[6]]$`Callaway & Sant'Anna 2021`,
            nrow = 6, ncol = 2, align = "hv", hjust = 0,
            labels = c(paste0(country.j,".1km.GAR"), paste0(country.j,".1km.CSA"),
                       paste0(country.j,".2.5km.GAR"), paste0(country.j,".2.5km.CSA"),
                       paste0(country.j,".5km.GAR"), paste0(country.j,".5km.CSA"),
                       paste0(country.j,".10km.GAR"), paste0(country.j,".10km.CSA"),
                       paste0(country.j,".20km.GAR"), paste0(country.j,".20km.CSA"),
                       paste0(country.j,".5kmMASTER.GAR"), paste0(country.j,".5kmMASTER.CSA")))
  
  ggsave(path = "Outputs/Figures/Raw/Basic", filename = paste0(country.j, ".all.basic.png"),
         j.arr, bg = "white",
         device = "png", width = 20, height = 35, units = "cm") 
}


