
## Read in lists ---------------------------------------------------------------

 read.bind.list <- function(x) {
   combined.data <- data.frame()
   
   for (i in 1:29) {
  file <- x$file[i]
  # Read the text file as a data frame
  load(paste0(dir.path,file))  # Change the arguments based on your file format
  mine.buffer.forest.loss.df$country <- x$country[[i]]
  # Bind the data to the combined data frame
  combined.data <- rbind(combined.data, mine.buffer.forest.loss.df)
   }
   return(combined.data)
 }
 
## Read covs lists ---------------------------------------------------------------
 
 read.covs.list <- function(x) {
   combined.data <- data.frame()
   
   for (i in 1:29) {
     file <- x$file[i]
     # Read the text file as a data frame
     load(paste0(file))  # Change the arguments based on your file format
     # Bind the data to the combined data frame
     combined.data <- rbind(combined.data, mine.cluster.points.df)
   }
   return(combined.data)
 }

## Add ISO codes ---------------------------------------------------------------
 add.iso <- function(x) {
   x %>% mutate(iso3 = case_when(country == "Angola" ~ "AGO",
                                 country == "Burundi" ~ "BDI",
                                 country == "Cameroon" ~ "CMR",
                                 country == "Central African Republic" ~ "CAF",
                                 country == "Comoros" ~ "COM",
                                 country == "Côte d_Ivoire" ~ "CIV",
                                 country == "Côte d'Ivoire" ~ "CIV",
                                 country == "Democratic Republic of the Congo" ~ "COD",
                                 country == "Equatorial Guinea" ~ "GNQ",
                                 country == "Ethiopia" ~ "ETH",
                                 country == "Gabon" ~ "GAB",
                                 country == "Ghana" ~ "GHA",
                                 country == "Guinea" ~ "GIN",
                                 country == "Guinea-Bissau" ~ "GNB",
                                 country == "Kenya" ~ "KEN",
                                 country == "Liberia" ~ "LBR",
                                 country == "Madagascar" ~ "MDG",
                                 country == "Malawi" ~ "MWI",
                                 country == "Mozambique" ~ "MOZ",
                                 country == "Nigeria" ~ "NGA",
                                 country == "Republic of the Congo" ~ "COG",
                                 country == "Rwanda" ~ "RWA",
                                 country == "Sierra Leone" ~ "SLE",
                                 country == "South Africa" ~ "ZAF",
                                 country == "South Sudan" ~ "SSD",
                                 country == "Swaziland" ~ "SWZ",
                                 country == "Tanzania" ~ "TZA",
                                 country == "Uganda" ~ "UGA",
                                 country == "Zambia" ~ "ZMB",
                                 country == "Zimbabwe" ~ "ZWE",
                                 ))
 }
 
## Tidy raw data and add covariate information----------------------------------
did.prep <- function(x, lead.time = -10, post.time = 15, 
                     type = "loss", covariates = NULL){
  if (type == "loss") {
    out <- x %>% 
      #mining data only runs to 2020
      filter(loss.year <= 2020) %>%
      mutate(year = loss.year - 2000,
             treatment = ifelse(year >= first.mine.year, TRUE, FALSE),
             treatment.n = ifelse(year >= first.mine.year, 1, 0),
             rel.year.first = year - first.mine.year,
             cumulative.forest.loss.perc = cumulative.forest.loss.prop.area*100) %>%
      filter(rel.year.first >= lead.time & rel.year.first <= post.time)
    # left_join(mine.cluster.points.df, by = "CLUSTER")
  }else{
    out <- x %>% 
      #mining data only runs to 2020
      filter(degradation.year <= 2020) %>%
      mutate(year = degradation.year - 2000,
             treatment.n = ifelse(year >= first.mine.year, 1, 0),
             treatment = ifelse(year >= first.mine.year, TRUE, FALSE),
             rel.year.first = year - first.mine.year,
             cumulative.degradation.perc = cumulative.degradation.prop*100) %>%
      filter(rel.year.first >= lead.time & rel.year.first <= post.time)
  }
  if (!is.null(covariates)) { out <- left_join(out, covariates, by = "CLUSTER_ID")}
  return(out)
}

## scale covars ----------------------------------------------------------------
scale.tidy <- function(x) {
  x %>% filter(!is.na(travel.time), !is.na(slope),
               !is.na(pop.density), !is.na(elevation)) %>%
    ungroup() %>%
    mutate(travel.time.z =  (travel.time  - mean(travel.time))/
                 sd(travel.time),
               pop.density.z =  (pop.density  - mean(pop.density))/
                 sd(pop.density),
               slope.z =  (slope  - mean(slope))/
                 sd(slope),
               elevation.z =  (elevation  - mean(elevation))/
                 sd(elevation),
               )
}
## Estimate dynamic DiD --------------------------------------------------------
fit.dynamic.DiD <- function(data = x, yname = "cumulative.forest.loss.perc",
                            xformula = NULL, method = c("csa", "gardner"),
                            GAR.pre.period = -5){
  all.coefs.vars <- NULL
  for (y.var in yname) {
    csa.coefs <- NULL
    gardner.coefs <- NULL
    all.coefs <- NULL
    
    buff <- unique(data$buffer.size)
    cdata <- data
    gdata <- data
    
    if ("gardner" %in% method) {
      if (is.null(xformula)) {
        first_stage = ~ 0 | CLUSTER_ID + year 
        gdata <- gdata %>% filter(rel.year.first >= GAR.pre.period)
      } else {
        first_stage <- stats::as.formula(glue::glue("~ 0 + {xformula} | CLUSTER_ID + year")[[2]])
        gdata <- gdata %>% filter(rel.year.first >= GAR.pre.period)
      }
      
      gardner.did <- did2s(data = gdata, yname = y.var, 
                           treatment = "treatment",
                           first_stage = first_stage, 
                           second_stage = ~ i(rel.year.first, ref = c(-1)),
                           cluster_var = "CLUSTER_ID", verbose = TRUE)
      gardner.coefs <- did2s.tidy(gardner.did, buff = buff)
      gardner.coefs$y.var <- y.var
    }
    
    if ("csa" %in% method) {
      csa.did <- att_gt(data = cdata, yname = y.var, tname="year",
                        idname="CLUSTER_ID", gname = "first.mine.year",
                        control_group = "notyettreated", base_period = "varying",
                        xformla = xformula, clustervars = "CLUSTER_ID", 
                        bstrap=T, cband=T)
      
      csa.coefs <- csa.tidy(csa.did, buff = buff)
      csa.coefs$y.var <- y.var
      
    }
    if (!is.null(csa.coefs) && !is.null(gardner.coefs)) {
      all.coefs <- rbind(csa.coefs, gardner.coefs)
    } else if (!is.null(csa.coefs)) {
      all.coefs <- csa.coefs
    } else if (!is.null(gardner.coefs)) {
      all.coefs <- gardner.coefs
    } else {
      next
    }
    all.coefs.vars <- rbind(all.coefs, all.coefs.vars)
  }
  return(all.coefs.vars)
  
}

## Estimate dynamic DiD - SSA --------------------------------------------------------
fit.dynamic.DiD.SSA <- function(data = x, yname = "cumulative.forest.loss.perc",
                            xformula = NULL, method = c("csa", "gardner"),
                            GAR.pre.period = -5){
  all.coefs.vars <- NULL
  for (y.var in yname) {
    buff <- unique(data$buffer.size)
    cdata <- data
    gdata <- data
    
    if ("gardner" %in% method) {
      if (is.null(xformula)) {
        first_stage = ~ 0 | cluster.country.id + year 
        gdata <- gdata %>% filter(rel.year.first >= GAR.pre.period)
      } else {
        first_stage <- stats::as.formula(glue::glue("~ 0 + {xformula} | cluster.country.id + year")[[2]])
        gdata <- gdata %>% filter(rel.year.first >= GAR.pre.period)
      }
      
      gardner.did <- did2s(data = gdata, yname = y.var, 
                           treatment = "treatment",
                           first_stage = first_stage, 
                           second_stage = ~ i(rel.year.first, ref = c(-1)),
                           cluster_var = "cluster.country.id", verbose = TRUE)
      gardner.coefs <- did2s.tidy(gardner.did, buff = buff)
      gardner.coefs$y.var <- y.var
    }
    
    if ("csa" %in% method) {
      csa.did <- att_gt(data = cdata, yname = y.var, tname="year",
                        idname="cluster.country.id", gname = "first.mine.year",
                        control_group = "notyettreated", base_period = "varying",
                        xformla = xformula, clustervars = "cluster.country.id", 
                        bstrap=T, cband=T)
      
      csa.coefs <- csa.tidy(csa.did, buff = buff)
      csa.coefs$y.var <- y.var
      
    }
    if (exists("csa.coefs") & exists("gardner.coefs")) {
      all.coefs <- rbind(csa.coefs, gardner.coefs)
    } else { if (exists("csa.coefs") & !exists("gardner.coefs")) {
      all.coefs <- csa.coefs }else{all.coefs <- gardner.coefs}}
    all.coefs.vars <- rbind(all.coefs, all.coefs.vars)
  }
  return(all.coefs.vars)
  
}  
## extract DiD2s estimates from a DiD2S fixest object---------------------------
did2s.tidy <- function(x, buff = buff) {
  broom::tidy(x) %>%
    mutate(lci = estimate - std.error*1.96,
           uci = estimate + std.error*1.96,
           year.since = as.numeric(gsub("rel.year.first::", "", term))) %>%
    select(estimate, p.value, uci, lci, year.since) %>%
    rbind(tibble(estimate = 0, p.value= NA, uci= 0, lci= 0, year.since = -1)) %>%
    mutate(method = "Gardner 2022", buffer.size = buff)
}

## extract coefs from CSA methods objects --------------------------------------
csa.tidy <- function(x, buff = buff) {
  x1 <-  aggte(x, type = "dynamic",na.rm=T)
  data.frame(year.since = x1$egt, 
             estimate = x1$att.egt, 
             lci = x1$att.egt - x1$crit.val.egt*x1$se.egt,
             uci = x1$att.egt + x1$crit.val.egt*x1$se.egt,
             p.value = NA, method = "Callaway & Sant'Anna 2021", buffer.size = buff)
}

## plot the ATT as a factor through time----------------------------------------
plot.DiD2S.f <- function(x, 
                         pre.mine.col = c("#C6507EFF", "#0D0887FF"),
                         post.mine.col = c("#EF7F4FFF", "#F0F921FF"),
                         pre.yrs = NULL, post.yrs = NULL, legend = "none",
                         y.label = "Forest loss (% points)", method = c("csa", "gardner")){
  
  precol_magma <- colorRampPalette(pre.mine.col)
  postcol_magma <- colorRampPalette(post.mine.col)
  
  storage.ls <- list()
  for (i in method) {
    if (i == "csa") {i <- "Callaway & Sant'Anna 2021"}
    if (i == "gardner") {i <- "Gardner 2022"}
    
    if (is.null(pre.yrs) & is.null(post.yrs)){
      post.yrs <- max(x$year.since)
      pre.yrs <- min(x$year.since)
    } else {
      x.method <- x[[1]] %>% filter(year.since >= pre.yrs & year.since <= post.yrs,
                               method == i)
    }
    
    
    plt.i <- ggplot(x.method, 
                    aes(year.since, estimate, colour = year.since)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, size = .8) +
      scale_color_gradientn(colors=c(precol_magma(abs(pre.yrs)),postcol_magma(post.yrs + 1)),
                            name="Years since mining",breaks = c(pre.yrs,0,post.yrs)) +
      geom_vline(xintercept = -1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Years since mining") + ylab(paste0(y.label)) +
      theme_minimal(base_size = 8) +
      theme(legend.position = legend, legend.key.width = unit(01, "cm"))
    
    storage.ls[[paste0(i)]] <- plt.i 
  }
  return(storage.ls)
}

## plot the ATT as a ribbon through time----------------------------------------
plot.DiD2S.s <- function(x, pre.yrs = NULL, post.yrs = NULL,
                         pre.col = "#1b7837", post.col = "#762a83", 
                         y.label = "Forest loss (% points)", method = c("csa", "gardner")){
  storage.ls <- list()
  for (i in method) {
    if (i == "csa") {i <- "Callaway & Sant'Anna 2021"}
    if (i == "gardner") {i <- "Gardner 2022"}
    
    if (is.null(pre.yrs) & is.null(post.yrs)){
      post.yrs <- max(x$year.since)
      pre.yrs <- min(x$year.since)
    }
    
    plt.dat.pre <- x %>%
      filter(year.since >= pre.yrs & year.since <= post.yrs, year.since<=-1,
             method == i) %>%
      mutate(treatment = "pre")
    plt.dat.post <- x %>%
      filter(year.since >= pre.yrs & year.since <= post.yrs, year.since>=-1,
             method == i) %>%
      mutate(treatment = "post")
    
    plt.i <-  ggplot(plt.dat.post, aes(year.since, estimate, colour = treatment, fill = treatment)) +
      geom_line() +
      geom_ribbon(aes(ymin = lci, ymax = uci), alpha = .5, linetype = "dotted") +
      geom_line(data = plt.dat.pre) +
      geom_ribbon(data = plt.dat.pre, aes(ymin = lci, ymax = uci), alpha = .5, linetype = "dotted") +
      scale_fill_manual(values = c(post.col, pre.col)) +
      scale_colour_manual(values = c(post.col, pre.col)) +
      geom_vline(xintercept = -1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Years since mining") + ylab(paste0(y.label)) +
      theme_minimal(base_size = 8) +
      theme(legend.position = "none")
    
    storage.ls[[paste0(i)]] <- plt.i 
  }
  return(storage.ls)
}

## Scaling ---------------------------------------------------------------------

scale.did <- function(x) {
  x %>% ungroup() %>%
    mutate(travel.time.z = (travel.time - mean(travel.time))/sd(travel.time),
           pop.density.z = (pop.density - mean(pop.density))/sd(pop.density),
           road.distance.z = (road.distance - mean(road.distance))/sd(road.distance),
           river.distance.z = (river.distance - mean(river.distance))/sd(river.distance),
           elevation.z = (elevation - mean(elevation))/sd(elevation),
           slope.z = (slope - mean(slope))/sd(slope))
}

## Spillover tidying -----------------------------------------------------------

spill.tidy <- function(x, buff, concentric.ring) {
  if (concentric.ring == "fixed") {
    x %>% 
      rename("first.spillover.i" = "spillover.10km.master") %>%
      mutate(rel.spill.year = year - first.spillover.i,
             rel.spill.year = ifelse(is.na(rel.spill.year), -Inf, rel.spill.year),
             spill.treat = ifelse(rel.spill.year >= 0, 1, 0),
             spill.or.treat = ifelse(spill.treat + treatment.n > 0, 1, 0))
  }else{
  x %>% 
    rename("first.spillover.i" = paste0("spillover.", buff)) %>%
    mutate(rel.spill.year = year - first.spillover.i,
               rel.spill.year = ifelse(is.na(rel.spill.year), -Inf, rel.spill.year),
               spill.treat = ifelse(rel.spill.year >= 0, 1, 0),
               spill.or.treat = ifelse(spill.treat + treatment.n > 0, 1, 0))
  }
}

## spillover dynamic DiD -------------------------------------------------------

spillover.dynamic.DiD <- function(x, yname = "cumulative.forest.loss.perc"){
  
  second_stage_spill = ~ i(rel.year.first, ref = c(-1)) + i(rel.spill.year, ref = c(-Inf, -1))
  second_stage_nospill = ~ i(rel.year.first, ref = c(-1))
  
  spill.mod <- did2s::did2s(
    data = x,
    yname = yname,
    first_stage = ~ 0 | CLUSTER_ID + year,
    second_stage = second_stage_spill,
    treatment = "spill.or.treat",
    cluster_var = "CLUSTER_ID")
  
  no.spill.mod <- did2s::did2s(
    data = x,
    yname = yname,
    first_stage = ~ 0 | CLUSTER_ID + year,
    second_stage = second_stage_nospill,
    treatment = "treatment",
    cluster_var = "CLUSTER_ID")
  
  spill.coefs <- broom::tidy(spill.mod) %>%
    mutate(lci = estimate - std.error*1.96,
           uci = estimate + std.error*1.96,
           year.since = as.numeric(gsub("rel.year.first::|rel.spill.year::", "", term)),
           effect = ifelse(grepl("spill", term),"Spillover effect", "Mine effect"))%>%
    select(effect, estimate, p.value, uci, lci, year.since) %>%
    rbind(tibble(estimate = 0, p.value= NA, uci= 0, lci= 0, year.since = -1, 
                 effect = c("Mine effect", "Spillover effect"))) %>%
    mutate(method = "Gardner 2022", buffer.size = "master5000")
  
  no.spill.coefs <- broom::tidy(no.spill.mod) %>%
    mutate(lci = estimate - std.error*1.96,
           uci = estimate + std.error*1.96,
           year.since = as.numeric(gsub("rel.year.first::|rel.spill.year::", "", term)),
           effect = "Total effect")%>%
    select(effect, estimate, p.value, uci, lci, year.since) %>%
    rbind(tibble(estimate = 0, p.value= NA, uci= 0, lci= 0, year.since = -1, 
                 effect = "Total effect")) %>%
    mutate(method = "Gardner 2022", buffer.size = "master5000")
  
  rbind(spill.coefs, no.spill.coefs)
}

## Spillover plot summary ------------------------------------------------------

spill.plot <- function(x, 
                         mine.col = "red", spill.col = "blue",
                         total.col = "black",
                         pre.yrs = NULL, post.yrs = NULL, legend = "none",
                         y.label = "Forest loss (% points)") {
  
    if (is.null(pre.yrs) & is.null(post.yrs)){
      post.yrs <- max(x$year.since)
      pre.yrs <- min(x$year.since)
    } 
  
      x.i <- x %>% filter(year.since >= pre.yrs & year.since <= post.yrs)
    
    plt.i <- ggplot(x.i, 
                    aes(year.since, estimate, colour = effect)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, size = .8) +
      scale_color_manual(values = c(mine.col, spill.col, total.col)) +
      geom_vline(xintercept = -1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Years since mining") + ylab(paste0(y.label)) +
      theme_minimal(base_size = 8) +
      theme(legend.position = legend, legend.key.width = unit(01, "cm"))
}

## spillover basic plot --------------------------------------------------------
plot.spill.DiD2S <- function(x, 
                         pre.mine.col = c("#C6507EFF", "#0D0887FF"),
                         post.mine.col = c("#EF7F4FFF", "#F0F921FF"),
                         pre.yrs = NULL, post.yrs = NULL, legend = "none",
                         y.label = "Forest loss (% points)"){
  
  precol_magma <- colorRampPalette(pre.mine.col)
  postcol_magma <- colorRampPalette(post.mine.col)
  
    if (is.null(pre.yrs) & is.null(post.yrs)){
      post.yrs <- max(x$year.since)
      pre.yrs <- min(x$year.since)
    } else {
      x.method <- x %>% filter(year.since >= pre.yrs & year.since <= post.yrs)
    }
    
    
    ggplot(x.method, 
                    aes(year.since, estimate, colour = year.since)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, size = .8) +
      scale_color_gradientn(colors=c(precol_magma(abs(pre.yrs)),postcol_magma(post.yrs + 1)),
                            name="Years since mining",breaks = c(pre.yrs,0,post.yrs)) +
      geom_vline(xintercept = -1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Years since mining") + ylab(paste0(y.label)) +
      theme_minimal(base_size = 8) +
      theme(legend.position = legend, legend.key.width = unit(01, "cm"))
    
}

## create_sub_exp() ------------------------------------------------------------

## credit to:
# https://www.nber.org/system/files/working_papers/w32054/w32054.pdf
# https://github.com/hollina/stacked-did-weights/blob/main/stacked-example-r-and-stata.qmd

# create_sub_exp() constructs a sub-experimental data set for a specified 
# adoption year. The code that defines the function is below in both Stata and R. 
# 
# The inputs and outputs of the `create_sub_exp()` are:
#   
#   -   inputs:
#   -   `dataset`: underlying panel dataset
# -   `timeID`: name of the variable measuring calendar time
# -   `groupID`: name of the variable indicating group
# -   `adoptionTime`: name of the variable containing the treatment adoption time for each unit; the variable is set to `NA` for those that never adopt
# -   `focalAdoptionTime`: user specified number indicating the focal adoption time for this sub-experiment
# -   `kappa_pre`: user specified number of pre-treatment periods
# -   `kappa_post`: user specified number of post-treatment periods
# -   outputs: A data.table that consists of the treated units that adopt in `timeID` `adoptTime`, the clean controls that adopt no earlier than `adoptTime` + `kappa_post`, and only the calendar time units that fall within $adoptTime - kappa_{pre}$ and $atime + kappa_{post}$. The output table consists of one row for each state-year observation in the sub-experimental data set. There is also a `feasible` variable that indicates whether or not a given state-year observation is "allowed" given our trimming restrictions. 


create_sub_exp <- function(dataset, timeID, groupID, adoptionTime, 
                           focalAdoptionTime, kappa_pre, kappa_post){
  
  # Copy dataset 
  dt_temp = copy(dataset)
  
  # Determine earliest and latest time in the data. 
  # Used for feasibility check later
  minTime = dt_temp[, fmin(get(timeID))]
  maxTime = dt_temp[, fmax(get(timeID))]
  
  # Include only the treated groups and the clean controls that adopt at least kappa_post periods after the focal atime.
  dt_temp = dt_temp[get(adoptionTime) == focalAdoptionTime | get(adoptionTime) > focalAdoptionTime + kappa_post | get(adoptionTime) == TRUE | is.na(get(adoptionTime))]
  
  # Limit to time periods inside the event window defined by the kappas
  dt_temp = dt_temp[get(timeID) %in% (focalAdoptionTime - kappa_pre):(focalAdoptionTime + kappa_post)]
  
  # Make treatment group dummy
  dt_temp[, treat := 0]
  dt_temp[get(adoptionTime) == focalAdoptionTime, treat := 1] 
  
  # Make a post variable
  dt_temp[, post := fifelse(get(timeID) >= focalAdoptionTime, 1, 0)]
  
  # Make event time variable
  dt_temp[, event_time := get(timeID) - focalAdoptionTime]
  
  # Create a feasible variable
  dt_temp[, feasible := fifelse(focalAdoptionTime - kappa_pre >= minTime & focalAdoptionTime + kappa_post <= maxTime, 1, 0)]
  
  # Make a sub experiment ID.
  dt_temp[, sub_exp := focalAdoptionTime]
  
  return(dt_temp)
} 

## compute_weights() -----------------------------------------------------------

## credit to:
# https://www.nber.org/system/files/working_papers/w32054/w32054.pdf
# https://github.com/hollina/stacked-did-weights/blob/main/stacked-example-r-and-stata.qmd

# `compute_weights()` constructs the corrective weights used in the weighted stacked 
# DID regressions. `compute_weights()` takes a stacked data set as input and 
# computes the corrective sample weights for the treated and control groups
# 
# -   inputs: 
#   -   `stack_data`: a stacked dataset created using `create_sub_exp()` and the appending procedure above
# -   `treatedVar`: the name of the variable that indicates whether a unit serves as a treated unit in a given sub-experiment
# -   `eventTimeVar`: the name of the variable that indicates the event time for each sub-experiment
# -   `subexpVar`: the name of the variable that indicates the sub-experiment for each unit
# -   outputs: 
#   - the original `stack_data`, but with a new column of corrective sample weights `stack_weights`

compute_weights <- function(dataset, treatedVar, eventTimeVar, subexpVar) {
  
  # Create a copy of the underlying dataset
  stack_dt_temp = copy(dataset)
  
  # Step 1: Compute stack - time counts for treated and control
  stack_dt_temp[, `:=` (stack_n = .N,
                        stack_treat_n = sum(get(treatedVar)),
                        stack_control_n = sum(1 - get(treatedVar))), 
                by = get(eventTimeVar)
  ]  
  # Step 2: Compute sub_exp-level counts
  stack_dt_temp[, `:=` (sub_n = .N,
                        sub_treat_n = sum(get(treatedVar)),
                        sub_control_n = sum(1 - get(treatedVar))
  ), 
  by = list(get(subexpVar), get(eventTimeVar))
  ]
  
  # Step 3: Compute sub-experiment share of totals
  stack_dt_temp[, sub_share := sub_n / stack_n]
  
  stack_dt_temp[, `:=` (sub_treat_share = sub_treat_n / stack_treat_n,
                        sub_control_share = sub_control_n / stack_control_n
  )
  ]
  
  # Step 4: Compute weights for treated and control groups
  stack_dt_temp[get(treatedVar) == 1, stack_weight := 1]
  stack_dt_temp[get(treatedVar) == 0, stack_weight := sub_treat_share/sub_control_share]
  
  return(stack_dt_temp)
}  
