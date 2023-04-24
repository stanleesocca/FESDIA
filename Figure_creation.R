
## =============================================================================
##    Creation of figures in the manuscripts 
## Version: For submission to Journal GMD
## FESDIA: Exploring temporal variations of sediment biogeochemistry under the influence of flood events using numerical modelling
### Stanley: 16/02/2022 
## =============================================================================
library(FESDIA)
library(ggplot2)
library(magrittr)
library(cowplot)
library(ocedata)
library(oce)
library(sf)

## load coastline data 
data("coastlineWorldFine")

## load local scripts and data
source("scripts/helper_func.R")
# data for mapping 
load("data/maps-data.RData")


## ==========================================================
##          Load sensitivity analysis result 
## ==========================================================
load("data/relax_SO4_conc_thick.RData")
load("data/relax_DIC_conc_thick.RData")
load("data/relax_O2_conc_thick.RData")


## ==========================================================
##          Load Bootstrapping result 
## ==========================================================
load("data/bootstrap_relax_est.RData")
load("data/EM_Scenario.RData")
load("data/Model_Obs_Dec_2008.RData")

station_A <- readr::read_delim("data/station_A.txt", delim = "\t")
## special treatment because oxygen data arrived late 
load("data/O2_Data_2008.RData")
## special treatment because oxygen data arrived late
sobsDIC <- readr::read_csv("data/DIC_Dec_2008.csv", col_names = FALSE)

## set the figure directory if does not exist in user's directory
if(!dir.exists("./figures")) dir.create("./figures")

## =============================================================================
##                                Figure 1
## =============================================================================

RESO <- 300
PS <- 10

#Overall units in inches
WIDTHS <- c(2,2,2,2)
HEIGHTS <- c(2,2)
OMA <- c(0.5,0.5,2,0.5)
HEIGHT <- sum(HEIGHTS) + OMA[1]*PS*1/72 + OMA[3]*PS*1/72
WIDTH <- sum(WIDTHS) + OMA[2]*PS*1/72 + OMA[4]*PS*1/72
LOFRANCE <- matrix(c(1,1,1,2), nrow=2, ncol=2, byrow=TRUE)

png("figures/fig01.png", width = 10, height = 8, units = "in", res = RESO)

layout(LOFRANCE, heights=HEIGHTS, widths=WIDTHS)

bathyLon = as.numeric(rownames(b))
bathyLat = as.numeric(colnames(b))
bathyZ = as.numeric(b)
dim(bathyZ) = dim(b)

plot(coastlineWorldFine, clon = mean(c(3, 8)), clat = mean(c(40, 46)),  
     span = 300,  yaxt = "n", xaxt = "n", axes = FALSE)
box()
axis(1, at = axTicks(1), labels = parse(text = paste(axTicks(1), "*degree ~ N", sep = "")), cex.axis = 1.5)
axis(2, at = axTicks(2), labels = parse(text = paste(axTicks(2), "*degree ~ E", sep = "")), cex.axis = 1.5)
points(4.8, 43.3, pch = 24, cex=2, col=3, bg="yellow", lwd=3)
text(4.8, 43.25, label = "Station A", cex = 1)
text(4.8, 43.52, label = "Rhone River", cex = 1, srt = -70)
contour(bathyLon,bathyLat,bathyZ,
        levels = c(-50, -100, -150, -200, -250),
        lwd = c(1, 1, 2, 2, 3),
        lty = c(3, 1, 3, 1, 3),
        drawlabels = T, labcex = 0.8, add = TRUE, col = 'darkgray')
## add bathymetry legend 
segments(x0 = 3.8, x1 = 4, y0 = 42.3, y1 = 42.3, col = "darkgray", lwd = 2) 
text(x = 4.5, y = 42.3, label = "Depth (m)", cex = 1.5)
prettymapr::addnortharrow()
prettymapr::addscalebar()
par(mar=c(5.1, 4.1, 2.1, 2.1))
mapPlot(coastlineCut(coastlineWorldFine),
        longitudelim=c(1,5), latitudelim=c(42, 50),
        projection="+proj=leac +lat_0=44 +lat_1=49 +lon_0=3", col='gray', axes = FALSE, grid = FALSE)

mapText(3, 46, labels = "FRANCE", cex = 1)
mapPolygon(longitude = c(4.143516, 4.143516, 5.70489, 5.70489), latitude = c(44.09636, 42.95620, 42.95620, 44.09636), lwd = 2)

dev.off()


## =============================================================================
##                                Figure 2
## =============================================================================

## created with Inkscape 

## =============================================================================
##                                Figure 3
## =============================================================================

## create with Inkscape

## =============================================================================
##                                Figure 4
## =============================================================================
png("GMD_Pub_2022/figures/fig04.png", width = 15, height = 10, units = "in", res = RESO)
## ==================================================
##                process DIC data 
## ==================================================
names(sobsDIC) <- c("DIC", "Depth")  ## rename the data to match convention
sobsDIC$DIC <- sobsDIC$DIC * 1e3     ## convert DIC to model unit (mmol/m3)

## ==================================================
##               process O2 data 
## ==================================================
## convert depth in mm to cm, and invert the depth scale, and mean of the replicate data
sobsO2 <- data.frame(Depth = -1*Dec08[, 1]/10, O2 = rowMeans(Dec08[, 2:5])) 
## rename colnames 
colnames(sobsO2) <- c("Depth", "O2")  
## take only depth within model domain
sobsO2 = sobsO2[sobsO2$Depth >= 0, ] 

## ==================================================
##        process station observation data 
## ==================================================
## Subset only data in December flood event 
sobsDEC <- data.frame( Depth = station_A$depths_dec_1, 
                       Fe = station_A$fe_dec*1E-3, 
                       NO3 = station_A$no3_dec, 
                       NH3 = station_A$nh4_dec, 
                       SO4 = station_A$so4_dec*1e3, 
                       TOC = station_A$toc_dec
)

## ==================================================
##        Integrate all data into one dataframe 
## ==================================================
## merge all data into one
sobs = merge(sobsDEC, sobsDIC,  all = T)
sobs = merge(sobs, sobsO2, all = T)

## =================================================
##          Data - Model visualization
## ==================================================

## get the longitude and labitude of the bathymetric data 
bathyLon = as.numeric(rownames(b))
bathyLat = as.numeric(colnames(b))
bathyZ = as.numeric(b)
dim(bathyZ) = dim(b)

## create layout grid for plotting the data 
mat <- layout(rbind(c(1, 1, 1, 2, 2), 
                    c(0, 0, 0, 0, 0),
                    c(3, 4, 5, 6, 7)), 
              heights = c(1, lcm(0.25), 2))

## draw map
overlay_map()
## draw descriptive text for data
overlay_context()
## draw data and model 
overlay_data(EVENT_Dec, c("TOC", "SO4", "DIC", "NH3", "O2"), c(314,324, 339), mainObs = sobs, ylim = c(-25, 0))


dev.off()

## =============================================================================
##                                Figure 5
## =============================================================================
png("./GMD_Pub_2022/figures/fig05.png", width = 10, height = 8, units = "in", res = RESO)
## plot depth profile using ggplot2
p1 = plot_relax_depth_profile(EM1, "SO4", c(161, 171))
p2 = plot_relax_depth_profile(EM1, "SO4", c(200, 210))
p3 = plot_relax_depth_profile(EM1, "SO4", c(240, 250))
p4 = plot_relax_depth_profile(EM1, "SO4", c(300, 320))

## plot time series of distance function 
elapsed_EM1 <- ggplot() +
  theme_bw() +
  geom_line(data = data.frame(elapsed_time = seq(161, 700), deviate = ss_so4_EM1),
            mapping = aes(elapsed_time - 161, deviate), col = "black", size = 1) +
  labs(x = "Elapsed time (days)", y = expression(varphi(t))) +
  geom_vline(xintercept = quick.elbow(ss_so4_EM1), lty = "dashed") +
  geom_hline(yintercept = median(ss_so4_EM1), lty = "dashed") +
  scale_x_continuous(breaks = seq(0, 600, 50))

## create inset for bootstrap sim result
inset <- bootstrap_histogram(boot_so4_EM1, remove_outlier = T, n = 10) +
  geom_vline(xintercept = compute_bootCI(boot_so4_EM1, remove_outlier = T, method = "percentile"), col = "red", linetype = "dashed", size = 0.5)

## draw the all plot into one
p = ggdraw(elapsed_EM1 + theme_half_open(12)) +
  draw_plot(inset, .65, .25, .3, 0.8) +
  draw_plot_label(
    label = c("", ""),
    x = c(0, 0.65),
    y = c(1, 0.95),
    hjust = 0.5,
    size = 12
  )


plt_obj <- plot_grid(plot_grid(p1, p2, p3, p4, ncol = 4), p, nrow = 2)

ggsave(filename = "./GMD_Pub_2022/figures/fig05.png", 
       plot = plt_obj,
       width = 10, height = 8, units = "in", dpi = RESO)

## =============================================================================
##                                Figure 6
## =============================================================================

p1 = plot_relax_depth_profile(EM1, "DIC", c(161, 171))
p2 = plot_relax_depth_profile(EM1, "DIC", c(200, 210))
p3 = plot_relax_depth_profile(EM1, "DIC", c(240, 250))
p4 = plot_relax_depth_profile(EM1, "DIC", c(320, 350))

## plot time series of distance function
elapsed_EM1 <- ggplot() +
  theme_bw() +
  geom_line(data = data.frame(elapsed_time = seq(161, 700) , deviate = ss_dic_EM1),
            mapping = aes(elapsed_time - 161, deviate), size = 1) +
  labs(x = "Elapsed time (days)", y = expression(varphi(t))) +
  geom_vline(xintercept = quick.elbow(ss_dic_EM1), lty = "dashed") +
  geom_hline(yintercept = median(ss_dic_EM1), lty = "dashed") +
  scale_x_continuous(breaks = seq(0, 600, 50))

## create inset for bootstrap sim result
inset <- bootstrap_histogram(boot_dic_EM1, remove_outlier = T, n = 10) +
  geom_vline(xintercept = compute_bootCI(boot_dic_EM1, remove_outlier = T, method = "percentile"), col = "red", linetype = "dashed", size = 0.5)

## draw the all plot into one
p = ggdraw(elapsed_EM1 + theme_half_open(12)) +
  draw_plot(inset, .65, .25, .3, 0.8) +
  draw_plot_label(
    label = c("", ""),
    x = c(0, 0.65),
    y = c(1, 0.95),
    hjust = 0.5,
    size = 12
  )

plt_obj <- plot_grid(plot_grid(p1, p2, p3, p4, ncol = 4), p, nrow = 2)

ggsave(filename = "./GMD_Pub_2022/figures/fig06.png", 
       plot = plt_obj,
       width = 10, height = 8, units = "in", dpi = RESO)

## =============================================================================
##                                Figure 7
## =============================================================================

## plot depth profile using ggplot2

  p1 = plot_relax_depth_profile(EM2, "SO4", c(314, 324))
  p2 = plot_relax_depth_profile(EM2, "SO4", c(340, 350))
  p3 = plot_relax_depth_profile(EM2, "SO4", c(400, 420))
  p4 = plot_relax_depth_profile(EM2, "SO4", c(425, 430))
  offset <- 314 



## plot time series of distance function depending on use_equal_time switch

elapsed_EM2 <- ggplot() +
    theme_bw() +
    geom_line(data = data.frame(elapsed_time = seq(314, 700) , deviate = ss_so4_EM2),
              mapping = aes(elapsed_time - offset, deviate), size = 1) +
    labs(x = "Elapsed time (days)", y = expression(varphi(t))) +
    geom_vline(xintercept = quick.elbow(ss_so4_EM2, 25), lty = "dashed") +
    geom_hline(yintercept = median(ss_so4_EM2), lty = "dashed") +
    scale_x_continuous(breaks = seq(0, 600, 50))

## create inset for bootstrap sim result
inset <- bootstrap_histogram(boot_so4_EM2, remove_outlier = T, n = 10, bias = 25) +
  geom_vline(xintercept = compute_bootCI(boot_so4_EM2, remove_outlier = T, method = "percentile", bias = 25), col = "red", linetype = "dashed", size = 0.5)

## draw the all plot into one
p = ggdraw(elapsed_EM2 + theme_half_open(12)) +
  draw_plot(inset, .65, .25, .3, 0.8) +
  draw_plot_label(
    c("", ""),
    c(0, 0.65),
    c(1, 0.95),
    size = 12
  )

plt_obj <- plot_grid(plot_grid(p1, p2, p3, p4, ncol = 4), p, nrow = 2)

ggsave(filename = "./GMD_Pub_2022/figures/fig07.png", 
       plot = plt_obj,
       width = 10, height = 8, units = "in", dpi = RESO)
## =============================================================================
##                                Figure 8
## =============================================================================


  p1 = plot_relax_depth_profile(EM2, "DIC", c(314, 324))
  p2 = plot_relax_depth_profile(EM2, "DIC", c(340, 350))
  p3 = plot_relax_depth_profile(EM2, "DIC", c(400, 420))
  p4 = plot_relax_depth_profile(EM2, "DIC", c(425, 430))
  offset <- 314


## plot time series of distance function depending on use_equal_time switch

  elapsed_EM2 <- ggplot() +
    theme_bw() +
    geom_line(data = data.frame(elapsed_time = seq(314, 700) , deviate = ss_dic_EM2),
              mapping = aes(elapsed_time - offset, deviate), size = 1) +
    labs(x = "Elapsed time (days)", y = expression(varphi(t))) +
    geom_vline(xintercept = quick.elbow(ss_dic_EM2, 25), lty = "dashed") +
    geom_hline(yintercept = median(ss_dic_EM2), lty = "dashed") +
    scale_x_continuous(breaks = seq(0, 600, 50))


## create inset for bootstrap sim result
inset <- bootstrap_histogram(boot_dic_EM2, remove_outlier = T, n = 10, bias = 25) +
  geom_vline(xintercept = compute_bootCI(boot_dic_EM2, remove_outlier = T, method = "percentile", bias = 25), col = "red", linetype = "dashed", size = 0.5)

## draw the all plot into one
p = ggdraw(elapsed_EM2 + theme_half_open(12)) +
  draw_plot(inset, .65, .25, .3, 0.8) +
  draw_plot_label(
    c("", ""),
    c(0, 0.65),
    c(1, 0.95),
    size = 12
  )


plt_obj <- plot_grid(plot_grid(p1, p2, p3, p4, ncol = 4), p, nrow = 2)

ggsave(filename = "./GMD_Pub_2022/figures/fig08.png", 
       plot = plt_obj,
       width = 10, height = 8, units = "in", dpi = RESO)

## =============================================================================
##                                Figure 9
## =============================================================================

out <- sapply(1:5, function(x){
  sapply(1:5, function(y){
    quick.elbow(relax_SO4_conc_thickness[[x]][[y]])
  })
})

## calculate relaxation time for sensitivity analysis of DIC
out_DIC <- sapply(1:5, function(x){
  sapply(1:5, function(y){
    quick.elbow(relax_DIC_conc_thickness[[x]][[y]])
  })
})

## calculate relaxation time for sensitivity analysis of O2
out_O2 <- sapply(1:5, function(x){
  sapply(1:5, function(y){
    quick.elbow(relax_O2_conc_thickness[[x]][[y]])
  })
})

## alpha and thickness grid line 
concfac = 10^seq(-0.5, 2,by = 0.5)
thickness <- seq(2, 30, by = 5)

## create matrix spanning the thickness and alpha resolution to draw the contour map
gridmat <- expand.grid(thickness = thickness[1:5], alpha = concfac[1:5])
## add individual vector in the matrix
gridmat$O2 <- as.vector(out_O2[1:5, 1:5])
gridmat$SO4 <- as.vector(out[1:5, 1:5])
gridmat$DIC <- as.vector(out_DIC[1:5, 1:5])

## draw graph in 1 x 3 plot
plt_obj <- cowplot::plot_grid(draw_contour_map(O2, ylab = TRUE), draw_contour_map(SO4), draw_contour_map(DIC), ncol = 3)

ggsave(filename = "./GMD_Pub_2022/figures/fig09.png", 
       plot = plt_obj,
       width = 15, height = 8, units = "in", dpi = RESO)
