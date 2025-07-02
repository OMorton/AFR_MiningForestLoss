library(did)
library(did2s)
library(tidyverse)
library(ggpubr)
library(rnaturalearth)

# this reads in a suite of functions this script uses predominately to ease the
# the repetitive nature of certain tasks e.g. processing, modelling and plotting buffer rings
source("functions.R")


## Figure 2 - National and SSA-wide temp and coefs -----------------------------
load("Outputs/DiD.tables/all.DiD.all.t.Jul25.RData")
load("Outputs/DiD.tables/SSA.tidy.Jul25.RData")

did.df <- bind_rows(did.ls) %>% add.iso()
did.df %>% distinct(country, cluster.n)

# Remove countries with less than 30 mines
did.df <- did.df %>% 
  filter(cluster.n > 29, pre.period %in% c(-10, -5)) %>%
  # just remove the csa -10 period as its identical to -5 as it was fit on the full panel
  filter(method == "Callaway & Sant'Anna 2021"& pre.period==-10 |
         method == "Gardner 2022"& pre.period==-10 |
         method == "Gardner 2022"& pre.period==-5) %>%
  mutate(sig = sign(lci) == sign(uci), 
         pre.period = ifelse(method == "Callaway & Sant'Anna 2021", "csa", pre.period))
write.csv(did.df, "Outputs/DiD.tables/national.did.all.csv")

did.df %>% group_by(method, buffer.size) %>%
  filter(year.since %in% c(5,10)) %>%
  summarise(sig = sum(sig))
ssa.did <-ssa.did %>% filter(pre.period %in% c(-10,-5, NA)) %>%
  mutate(pre.period = ifelse(is.na(pre.period), "csa", pre.period))

## Summaries
national.summary <- did.df %>% filter(year.since %in% c(5,10),
                  buffer.size %in% c("1km", "5km", "10km", "20km"))
ssa.summary <-ssa.did %>% filter(year.since %in% c(5,10),
                             buffer.size %in% c("1km", "5km", "10km", "20km"))
write.csv(national.summary, "Outputs/DiD.tables/national.did.summary.csv")
write.csv(ssa.summary, "Outputs/DiD.tables/SSA.did.summary.csv")

## plotting set up
method.ls <- unique(did.df$method)
buffer.ls <- c("1km", "5km", "10km", "20km")

did.plt.ls <- list()
yr.since <- 10

t <- -5
i <- 1
j <- 1

GAR.grid <- expand.grid(method.i = "Gardner 2022", 
                        buffer.i = buffer.ls, pre.period = c(-5, -10))
CSA.grid <- expand.grid(method.i = "Callaway & Sant'Anna 2021", 
                        buffer.i = buffer.ls, pre.period = "csa")
plt.grid <- rbind(GAR.grid, CSA.grid)

for(i in 1:nrow(plt.grid)) {

      method.i <- plt.grid[i,]$method.i
      buffer.i <- plt.grid[i,]$buffer.i
      t <- plt.grid[i,]$pre.period
      
      nat.i <- did.df %>% filter(year.since >= -5 & year.since <= yr.since,
                                 pre.period == t,
                                buffer.size == buffer.i, method == method.i)%>%
        select(estimate, p.value, sig, uci, lci, year.since, method, buffer.size, country, iso3)  %>%
        mutate(sig.yr = ifelse(year.since == yr.since, sig, NA)) %>%
        group_by(iso3) %>%
        fill(sig.yr, .direction = "updown") %>% ungroup()
      
      ssa.i <- ssa.did %>% filter(year.since >= -5 & year.since <= yr.since,
                                  pre.period == t,
                                 buffer.size == buffer.i, method == method.i) %>%
        mutate(country = "All", iso3 = "ALL", order = 25,
               sig = sign(lci) == sign(uci))%>%
        select(estimate, p.value, sig, uci, lci, year.since, method, buffer.size, country, iso3, order) %>%
        mutate(sig.yr= sig)
      
      did.all.i <- nat.i %>% filter(year.since == yr.since) %>%
        arrange(estimate) %>%
        mutate(order = 1:n()) %>%
        rbind(filter(ssa.i, year.since == yr.since))
      
  
      # time plot
      temp.plt.i <- ggplot(nat.i, aes(year.since, estimate, group = country)) +
        # geom_ribbon(aes(ymin = lci , ymax = uci, group = country),
        #             linetype = "dashed", alpha = .3) +
       # geom_line(aes(colour = sig.yr), alpha = .4) +
       #  scale_color_manual(values = c("red", "black")) +
        geom_line(colour = "grey25", alpha = .4) +
        geom_ribbon(data = ssa.i, 
                     aes(year.since, estimate, group = NA, 
                         ymin = lci, ymax = uci, colour = NA), fill = "red",
                     alpha = .2) +
        geom_line(data = ssa.i, aes(year.since, estimate, group = NA, colour = NA), 
                   colour = "red", size = 1) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = -1, linetype = "dashed") +
        xlab("Years since mining") +
        ylab("Additional \n deforestation (pp)") +
        theme_minimal() + 
        theme(legend.position  = "none")
      
      coef.horiz.plt.i <- ggplot(did.all.i, 
                           aes(reorder(iso3, -order), estimate, , colour = country, 
                               shape = sig)) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = lci, ymax = uci), width = 0) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("darkred", rep("black", 24))) +
        scale_shape_manual(values = c(1, 19)) +
        ylab("Additional \n deforestation (pp)") +
        theme_minimal() +
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_text(face = c("bold", rep("plain", 25)), angle = 45, vjust = 0),
              legend.position = "none")
  
      if (method.i == "Callaway & Sant'Anna 2021") {
        did.plt.ls[[paste0("temporal.", buffer.i, ".", method.i)]] <- temp.plt.i 
        did.plt.ls[[paste0("coefs.yr.horiz.", buffer.i, ".", method.i)]] <- coef.horiz.plt.i}else{
          did.plt.ls[[paste0("temporal.", buffer.i, ".", method.i, ".pre",t)]] <- temp.plt.i 
          did.plt.ls[[paste0("coefs.yr.horiz.", buffer.i, ".", method.i, ".pre",t)]] <- coef.horiz.plt.i
        } 
}


# GAR5
GAR5.did.horiz.plt.arr <- ggarrange(did.plt.ls$`temporal.1km.Gardner 2022.pre-5`,
                             did.plt.ls$`coefs.yr.horiz.1km.Gardner 2022.pre-5`,
                             did.plt.ls$`temporal.5km.Gardner 2022.pre-5`,
                             did.plt.ls$`coefs.yr.horiz.5km.Gardner 2022.pre-5`,
                             did.plt.ls$`temporal.10km.Gardner 2022.pre-5`,
                             did.plt.ls$`coefs.yr.horiz.10km.Gardner 2022.pre-5`,
                             did.plt.ls$`temporal.20km.Gardner 2022.pre-5`,
                             did.plt.ls$`coefs.yr.horiz.20km.Gardner 2022.pre-5`,
                             ncol = 2, nrow = 4,
                             labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                             widths = c(1, 1.7),
                             align = "hv", hjust = 0)

GAR5.did.horiz.plt.arr2 <- GAR5.did.horiz.plt.arr +
  annotation_custom(text_grob("0-1 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.95, ymax = 1.0) +
  annotation_custom(text_grob("1-5 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.7, ymax = 0.75) +
  annotation_custom(text_grob("5-10 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.45, ymax = 0.5) +
  annotation_custom(text_grob("10-20 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.2, ymax = 0.25)
  
ggsave(path = "Outputs/Figures/DiD", filename = "GAR5.DiD.Fig2.horiz.10yr.png",
       GAR5.did.horiz.plt.arr2, bg = "white",
       device = "png", width = 25, height = 25, units = "cm") 

# GAR10
GAR10.did.horiz.plt.arr <- ggarrange(did.plt.ls$`temporal.1km.Gardner 2022.pre-10`,
                                    did.plt.ls$`coefs.yr.horiz.1km.Gardner 2022.pre-10`,
                                    did.plt.ls$`temporal.5km.Gardner 2022.pre-10`,
                                    did.plt.ls$`coefs.yr.horiz.5km.Gardner 2022.pre-10`,
                                    did.plt.ls$`temporal.10km.Gardner 2022.pre-10`,
                                    did.plt.ls$`coefs.yr.horiz.10km.Gardner 2022.pre-10`,
                                    did.plt.ls$`temporal.20km.Gardner 2022.pre-10`,
                                    did.plt.ls$`coefs.yr.horiz.20km.Gardner 2022.pre-10`,
                                    ncol = 2, nrow = 4,
                                    labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                                    widths = c(1, 1.7),
                                    align = "hv", hjust = 0)

GAR10.did.horiz.plt.arr2 <- GAR10.did.horiz.plt.arr +
  annotation_custom(text_grob("0-1 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.95, ymax = 1.0) +
  annotation_custom(text_grob("1-5 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.7, ymax = 0.75) +
  annotation_custom(text_grob("5-10 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.45, ymax = 0.5) +
  annotation_custom(text_grob("10-20 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.2, ymax = 0.25)

ggsave(path = "Outputs/Figures/DiD", filename = "GAR10.DiD.Fig2.horiz.10yr.png",
       GAR10.did.horiz.plt.arr2, bg = "white",
       device = "png", width = 25, height = 25, units = "cm") 

# CSA
CSA.did.horiz.plt.arr <- ggarrange(did.plt.ls$`temporal.1km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`coefs.yr.horiz.1km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`temporal.5km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`coefs.yr.horiz.5km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`temporal.10km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`coefs.yr.horiz.10km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`temporal.20km.Callaway & Sant'Anna 2021`,
                             did.plt.ls$`coefs.yr.horiz.20km.Callaway & Sant'Anna 2021`,
                             ncol = 2, nrow = 4, widths =  c(1, 1.7),
                             labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                             align = "hv", hjust = 0)

CSA.did.horiz.plt.arr2 <- CSA.did.horiz.plt.arr +
  annotation_custom(text_grob("0-1 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.95, ymax = 1.0) +
  annotation_custom(text_grob("1-5 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.7, ymax = 0.75) +
  annotation_custom(text_grob("5-10 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.45, ymax = 0.5) +
  annotation_custom(text_grob("10-20 km", face = "bold", size= 12), xmin = 0.90, xmax = 1, ymin = 0.2, ymax = 0.25)


ggsave(path = "Outputs/Figures/SM", filename = "CSA.DiD.Fig2.10yr.png",
       CSA.did.horiz.plt.arr2, bg = "white",
       device = "png", width = 25, height = 25, units = "cm")


## Figure 3 direct and indirect effects ----------------------------------------
load("Outputs/DiD.tables/all.GAR.direct.indirect.DiD.RData")


# note the clsuter size of some points
gardner.ratio.sum <- did.all %>% 
  filter(y.var %in% c("direct.forest.loss","indirect.forest.loss"),
         year.since == 5, pre.length == -10, method == "Gardner 2022") %>%
  mutate(ci95.cert.dir = ifelse(uci > 0 & lci > 0, "yes", "no"),
         estimate = ifelse(estimate > 0, estimate, NA)) %>%
  group_by(country) %>%
  mutate(ci95.cert.dir = ifelse(any(ci95.cert.dir == "no"), "no", "yes")) %>%
  pivot_wider(id_cols = c("country", "buffer.size", "cluster.n", 
                          "pre.length" , "ci95.cert.dir"), 
              values_from = "estimate", names_from = "y.var") %>%
  mutate(ratio = indirect.forest.loss/direct.forest.loss,
         small.cluster = ifelse(cluster.n < 30, "yes", "no"),
         ratio = ifelse(ci95.cert.dir == "no", NA, ratio)) %>% 
  add.iso() %>%
  filter(small.cluster == "no")
write.csv(gardner.ratio.sum, "Outputs/DiD.tables/national.ratios.all.csv")

## maps
library(rnaturalearth)
afr.map <- ne_countries(type = "countries", continent = "Africa", returnclass = "sf") %>%
  select(iso_a3,iso_a3_eh, geometry)
ratio.all.iso <- gardner.ratio.sum %>% 
  left_join(afr.map, by = c("iso3" = "iso_a3"))

afr.inc <- afr.map %>% filter(iso_a3 %in% unique(gardner.ratio.sum$iso3))

map.ratio <- ggplot(ratio.all.iso, aes(fill = ratio)) +
  geom_sf(data = afr.map, aes(geometry = geometry), fill = "grey90") +
  geom_sf(aes(geometry = geometry)) +
  scale_fill_gradient(low ="#ffeda0", high = "#800026", 
                      "Ratio of indirect forest loss") +
  theme_void() +
  theme(legend.position = "bottom", legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(10,10,10,10))

## xy plot
load("X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/loss.summaries/all.countries.mining.loss.50p.df.RData")
country.mining.total.loss.df <- country.mining.total.loss.df %>% 
  mutate(country = ifelse(country == "Côte d'Ivoire", "Côte d_Ivoire", country))
country.ratio.area <- left_join(ratio.all.iso, country.mining.total.loss.df, by = "country") %>%
  mutate(missing = ifelse(is.na(ratio), "yes", "no"),
         ratio = ifelse(is.na(ratio), 0, ratio))

(country.xy.plt <- ggplot(filter(country.ratio.area, missing == "no"), 
                          aes(total.mined.area, ratio, shape = small.cluster,
                              size = total.mined.area, colour = ratio)) +
    geom_point(data = filter(country.ratio.area, missing == "yes"), colour = "grey50") +
    geom_point() +
    scale_color_gradient(low = "#ffeda0", high = "#800026") +
    scale_shape_manual(values = c(19, 21)) +
    geom_point() +
    #geom_label(aes(label = iso3)) +
    ggrepel::geom_text_repel(aes(label = iso3), size = 3, colour = "black" ) +
    ggrepel::geom_text_repel(data = filter(country.ratio.area, missing == "yes"),
                             aes(label = iso3), size = 3, colour = "black" ) +
    
    #scale_x_log10() +
    #scale_y_log10() +
    xlab("Mining forest loss (ha)") +
    ylab("Ratio of indirect forest loss") +
    coord_cartesian(ylim = c(-15, 260)) +
    theme_minimal() +
    theme(legend.position = "none"))

map.xy.arr <- ggarrange(map.ratio, country.xy.plt,
                        labels = c("a", "b"), common.legend = TRUE,
                        legend = "bottom", widths = c(1, 1.4))

ggsave(path = "Outputs/Figures/DiD", filename = "map.ratio.xy.png",
       map.xy.arr, bg = "white",
       device = "png", width = 7, height = 4, units = "in") 

## Figure 4 - commodity impacts ------------------------------------------------
dir.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/forest.loss.in.buffers/"
comm.path <- "X:/morton_research/User/bi1om/Research/Mining/AfricaWideMining_ForestLoss/Analysis/Data/NationalMining/country.mining.commodity.type/"

# load commodity data
comm.dir <- data.frame(file = list.files(path = comm.path, pattern = "5km.buffer"))
comm.df <- data.frame()
for (i in 1:nrow(comm.dir)) {
  load(paste0(comm.path,comm.dir$file[[i]]))
  comm.df <- rbind(comm.df, country.mining.commodity.data.df)
}

comm.df <- comm.df %>% separate_longer_delim(material, ",") %>% 
  distinct() %>% 
  filter(material != "No mine in Maus")

top.comm <- comm.df %>% group_by(material) %>% tally() %>%
  filter(n > 29)

nat.comm.df <- comm.df %>%
  group_by(country, material) %>% tally() %>%
  filter(material %in% top.comm$material)

afr.map <- ne_countries(type = "countries", continent = "Africa", returnclass = "sf") %>%
  select(iso_a3,iso_a3_eh, geometry)
nat.comm.iso <- add.iso(nat.comm.df) %>% 
  left_join(afr.map, by = c("iso3" = "iso_a3"))
afr.inc <- afr.map %>% filter(iso_a3 %in% unique(nat.comm.iso$iso3))

# get n
load("Outputs/DiD.tables/SSA.commodities.tidy.new.RData")
comm.n <- comm.did.out %>% group_by(commodity) %>% summarise(n = unique(cluster.n))

col.ls <- c("#3690c0", "#cc4c02", "#bdbdbd", "#fec44f", "#800026", "#9e9ac8",
            "#525252", "black", "#d9f0a3")
comm.map.ls <- list()
i <- 1

for (i in 1:length(top.comm$material)) {
  comm.i <- top.comm$material[i]
  col.i <- col.ls[i]
  nat.comm.i <- nat.comm.iso %>% filter(material == comm.i)
  comm.n.i <- comm.n %>% filter(commodity == comm.i)
  
  comm.map.i <- ggplot() +
    geom_sf(data = afr.map, colour = "grey20", fill = "white") +
    geom_sf(data = afr.inc, fill = "grey85") +
    geom_sf(data = nat.comm.i, aes(geometry = geometry), fill = col.i) +
    annotate("text", x = -5, y = -20, label = paste0("n = ", comm.n.i$n), size = 4.5) +
    theme_void()
  
  comm.map.ls[[paste0(comm.i , ".map")]] <- comm.map.i 
  
}  



## plot DiD
method.ls <- unique(comm.did.out$method)
comm.ls <- unique(comm.did.out$commodity)
buff.ls <- c(1, 5, 10, 20)
col.ls <- c("#3690c0", "#cc4c02", "#bdbdbd", "#fec44f", "#800026", "#9e9ac8",
            "#525252", "black", "#d9f0a3")
did.plt.ls <- list()

i <- 1
c <- 2
j <- 1
comm.plt.ls <- list()
for (m in 1:2) {
  for (i in 1:length(buff.ls)) {
    for (c in 1:length(comm.ls)){
      method.i <- method.ls[[m]]
      comm.c <- comm.ls[[c]]
      col.c <- col.ls[[c]]
      buff.j <- buff.ls[[i]]
      
      comm.i.c <- comm.did.out %>% 
        filter(commodity == comm.c, method == method.i, 
               buffer.size == paste0(buff.j, "km"),
               year.since >= -5 & year.since <= 10)
   
      comm.plt <- ggplot(comm.i.c, aes(year.since, estimate)) +
        geom_ribbon(aes(ymin = lci, ymax = uci), linetype = "dashed",
                    colour = col.c, fill = col.c, alpha = .3) +
        geom_line(size = 1, colour =col.c ) + 
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = -1, linetype = "dashed") +
       # coord_cartesian(ylim = c(-5, 15)) +
        xlab("Years since mining") +
        ylab("Additional \n deforestation (pp)") +
        theme_minimal(base_size = 12) 
      
      comm.plt.ls[[paste0("temporal.",comm.c , ".", buff.j, "km.", method.i)]] <- comm.plt 
    }
  }
}

empty <- ggplot() + theme_void()

comm.map.gardner.arr <- ggarrange(empty,
                                  ggarrange(comm.plt.ls$`temporal.Cobalt.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Cobalt.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Cobalt.10km.Gardner 2022`,
                              comm.map.ls$Cobalt.map,
                              comm.plt.ls$`temporal.Copper.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Copper.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Copper.10km.Gardner 2022`,
                              comm.map.ls$Copper.map,
                              comm.plt.ls$`temporal.Diamonds.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Diamonds.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Diamonds.10km.Gardner 2022`,
                              comm.map.ls$Diamonds.map,
                              comm.plt.ls$`temporal.Gold.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Gold.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Gold.10km.Gardner 2022`,
                              comm.map.ls$Gold.map,
                              comm.plt.ls$`temporal.Iron.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Iron.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Iron.10km.Gardner 2022`,
                              comm.map.ls$Iron.map,
                              comm.plt.ls$`temporal.Manganese.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Manganese.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Manganese.10km.Gardner 2022`,
                              comm.map.ls$Manganese.map,
                              comm.plt.ls$`temporal.Silver.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Silver.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Silver.10km.Gardner 2022`,
                              comm.map.ls$Silver.map,
                              comm.plt.ls$`temporal.Uranium.1km.Gardner 2022`,
                              comm.plt.ls$`temporal.Uranium.5km.Gardner 2022`,
                              comm.plt.ls$`temporal.Uranium.10km.Gardner 2022`,
                              comm.map.ls$Uranium.map,
                              ncol =4, nrow = 8, widths = c(1,1,1,.5),
                              labels = c("a - Cobalt", "", "", "",
                                         "b - Copper", "", "", "",
                                         "c - Diamonds", "", "", "",
                                         "d - Gold", "", "", "",
                                         "e - Iron", "", "", "",
                                         "f - Manganese", "", "", "",
                                         "g - Silver", "", "", "",
                                         "h - Uranium", "", "", ""),
                              hjust = 0, vjust = 0),
                              ncol = 1, heights = c(1, 40))

comm.map.gardner.arr2 <- comm.map.gardner.arr +
  annotation_custom(text_grob("0-1 km", face = "bold", size= 16), xmin = 0.6/3.5, xmax = 0.6/3.5, ymin = 0.97, ymax = 1.0) +
  annotation_custom(text_grob("1-5 km", face = "bold", size= 16), xmin = 1.6/3.5, xmax = 1.6/3.5, ymin = 0.97, ymax = 1.0) +
  annotation_custom(text_grob("5-10 km", face = "bold", size= 16), xmin = 2.6/3.5, xmax = 2.6/3.5, ymin = 0.97, ymax = 1.0)


# ggsave(path = "Outputs/Figures/DiD", filename = "GAR.commodities.map.fixedscale.png",
#        comm.map.gardner.arr2, bg = "white",
#        device = "png", width = 35, height = 40, units = "cm")  

ggsave(path = "Outputs/Figures/DiD", filename = "GAR.commodities.map.png",
       comm.map.gardner.arr2, bg = "white",
       device = "png", width = 35, height = 40, units = "cm")  

