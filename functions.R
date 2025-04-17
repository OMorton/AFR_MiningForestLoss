
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
             cumulative.forest.loss.perc = cumulative.forest.loss.prop*100) %>%
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

## Estimate dynamic DiD --------------------------------------------------------
fit.dynamic.DiD <- function(data = x, yname = "cumulative.forest.loss.perc",
                            xformula = NULL, method = c("csa", "gardner")){
  
  buff <- unique(data$buffer.size)
  cdata <- data
  gdata <- data
  
  if ("gardner" %in% method) {
    if (is.null(xformula)) {
      first_stage = ~ 0 | CLUSTER_ID + year 
      #gdata <- gdata %>% filter(rel.year.first >-11)
    } else {
      first_stage <- stats::as.formula(glue::glue("~ 0 + {xformula} | CLUSTER_ID + year")[[2]])
      #gdata <- gdata %>% filter(rel.year.first >-11)
    }
    
    gardner.did <- did2s(data = gdata, yname = yname, 
                         treatment = "treatment",
                         first_stage = first_stage, 
                         second_stage = ~ i(rel.year.first, ref = c(-1)),
                         cluster_var = "CLUSTER_ID", verbose = TRUE)
    gardner.coefs <- did2s.tidy(gardner.did, buff = buff)
  }
  
  if ("csa" %in% method) {
    csa.did <- att_gt(data = cdata, yname = yname, tname="year",
                      idname="CLUSTER_ID", gname = "first.mine.year",
                      control_group = "notyettreated", base_period = "varying",
                      xformla = xformula, clustervars = "CLUSTER_ID", 
                      bstrap=T, cband=T)
    
    csa.coefs <- csa.tidy(csa.did, buff = buff)
  }
  if (exists("csa.coefs") & exists("gardner.coefs")) {
    all.coefs <- rbind(csa.coefs, gardner.coefs)
  } else { if (exists("csa.coefs") & !exists("gardner.coefs")) {
    all.coefs <- csa.coefs }else{all.coefs <- gardner.coefs}
    return(all.coefs)
  }
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
      x.method <- x %>% filter(year.since >= pre.yrs & year.since <= post.yrs,
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

spill.tidy <- function(x) {
  x %>% mutate(rel.spill.year = year - first.spillover,
               rel.spill.year = ifelse(is.na(rel.spill.year), -Inf, rel.spill.year),
               spill.treat = ifelse(rel.spill.year >= 0, 1, 0),
               spill.or.treat = ifelse(spill.treat + treatment.n > 0, 1, 0))
}

## spillover dynamic DiD -------------------------------------------------------

spillover.dynamic.DiD <- function(x, yname = "cumulative.forest.loss.perc"){
  
  second_stage_spill = ~ i(rel.year.first, ref = c(-1)) + i(rel.spill.year, ref = c(-Inf, -1))
  second_stage_nospill = ~ i(rel.year.first, ref = c(-1))
  
  spill.mod <- did2s::did2s(
    data = x,
    yname = yname,
    first_stage = ~ 0 | CLUSTER + year,
    second_stage = second_stage_spill,
    treatment = "spill.or.treat",
    cluster_var = "CLUSTER")
  
  no.spill.mod <- did2s::did2s(
    data = x,
    yname = yname,
    first_stage = ~ 0 | CLUSTER + year,
    second_stage = second_stage_nospill,
    treatment = "treatment",
    cluster_var = "CLUSTER")
  
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
