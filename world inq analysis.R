#Libraries----
library(easypackages)
libraries("Hmisc", "rworldmap", "velox", "plm", "fields", "olsrr", "countrycode", "sp", "rgdal", "openxlsx", "ggplot2", "raster", "ineq", "ncdf4", "plyr", "gstat", "spatial.tools", "GGally", "mgcv", "rasterVis", "tictoc", "spatialEco", "usdm", "boot", "mixtools")
setwd("D:/R/Light Inequality/Illuminating global economic inequality using night light data")
load("D:/R/Light Inequality/Illuminating global economic inequality using night light data/world inq analysis.RData")

#Colors
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#Common files----
world <- readOGR(dsn = "world boundaries", layer = "TM_WORLD_BORDERS-0.3")
world_2010 <- readOGR(dsn = "world boundaries", layer = "world_2010")
world_2015 <- readOGR(dsn = "world boundaries", layer = "world_2015")
world_2006 <- readOGR(dsn = "world boundaries", layer = "world_2006")
pop_density <- read.csv("Pop/Pop Density/Pop density country wise/API_EN.POP.DNST_DS2_en_csv_v2.csv")
lights_1992 <- raster("Light/Lights 1992.tif")
lights_1995 <- raster("Light/Lights 1995.tif")
lights_2000 <- raster("Light/Lights 2000.tif")
lights_2005 <- raster("Light/Lights 2005.tif")
lights_2010 <- raster("Light/Lights 2010.tif") 
lights_2015 <- raster("Light/lights 2015.tif") 
pop_count_1990 <- raster("Pop/Pop Count/Pop 1990.asc", crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
pop_count_1995 <- raster("Pop/Pop Count/Pop 1995.asc", crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
pop_count_2000 <- raster("Pop/Pop Count/Pop 2000.tif")
pop_count_2005 <- raster("Pop/Pop Count/Pop 2005.tif")
pop_count_2010 <- raster("Pop/Pop Count/Pop 2010.tif")
pop_count_2015 <- raster("Pop/Pop Count/Pop 2015.tif")
pop_density_2010 <- raster("Pop/Pop Density/2010/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals_2010.tif")
pop_density_2015 <- raster("Pop/Pop Density/2015/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals_2015.tif")
df_inequality_2010 <- read.csv("Wealth Gini/Inequality Data 2010.csv")
df_inequality_2015 <- read.csv("Wealth Gini/Inequality Data 2015.csv")
pop_density_2010_gt0 <- raster("Pop/Pop Density/2010/pop_density_2010_gt0.tif")
pop_count_2010_gt0 <- raster("Pop/Pop Count/2010/pop_count_2010_gt0.tif")
pop_density_2015_gt0 <- raster("Pop/Pop Density/2015/pop_density_2015_gt0.tif")
pop_count_2015_gt0 <- raster("Pop/Pop Count/2015/pop_count_2015_gt0.tif")
pop_density_2010_gt10 <- raster("Pop/Pop Density/2010/pop_density_2010_gt10.tif")
pop_count_2010_gt10 <- raster("Pop/Pop Count/2010/pop_count_2010_gt10.tif")
pop_density_2015_gt10 <- raster("Pop/Pop Density/2015/pop_density_2015_gt10.tif")
pop_count_2015_gt10 <- raster("Pop/Pop Count/2015/pop_count_2015_gt10.tif")
gw_use_0610 <- raster("Water use/global/nag_cru_0610.tif")
sw_use_0610 <- raster("Water use/global/aus_cru_0610.tif")
gw_use_1115 <- raster("Water use/global/nag_cru_1115.tif")
sw_use_1115 <- raster("Water use/global/aus_cru_1115.tif")
lights_1990_calib <- raster("Light/Lights 1990 calib.tif")
lights_1995_calib <- raster("Light/Lights 1995 calib.tif")
lights_2000_calib <- raster("Light/Lights 2000 calib.tif")
lights_2005_calib <- raster("Light/Lights 2005 calib.tif")
lights_2010_calib <- raster("Light/Lights 2010 calib.tif") 
lights_2015_calib <- raster("Light/lights 2015 calib.tif")
pop_calib_1990 <- raster("Pop/Pop Count/Pop 1990 calib.tif")
pop_calib_1995 <- raster("Pop/Pop Count/Pop 1995 calib.tif")
pop_calib_2000 <- raster("Pop/Pop Count/Pop 2000 calib.tif")
pop_calib_2005 <- raster("Pop/Pop Count/Pop 2005 calib.tif")
pop_calib_2010 <- raster("Pop/Pop Count/Pop 2010 calib.tif")
pop_calib_2015 <- raster("Pop/Pop Count/Pop 2015 calib.tif")
lpp_calib_1990 <- raster("LPP/LPP 1990.tif")
lpp_calib_1995 <- raster("LPP/LPP 1995.tif")
lpp_calib_2000 <- raster("LPP/LPP 2000.tif")
lpp_calib_2005 <- raster("LPP/LPP 2005.tif")
lpp_calib_2010 <- raster("LPP/LPP 2010.tif")
lpp_calib_2015 <- raster("LPP/LPP 2015.tif")
lpp_calib_1990_gt10 <- raster("LPP/LPP 1990 gt10.tif")
lpp_calib_1995_gt10 <- raster("LPP/LPP 1995 gt10.tif")
lpp_calib_2000_gt10 <- raster("LPP/LPP 2000 gt10.tif")
lpp_calib_2005_gt10 <- raster("LPP/LPP 2005 gt10.tif")
lpp_calib_2010_gt10 <- raster("LPP/LPP 2010 gt10.tif")
lpp_calib_2015_gt10 <- raster("LPP/LPP 2015 gt10.tif")
agg_map_90 <- raster("Smoothed Gini/agg_gini_1990.tif")
agg_map_95 <- raster("Smoothed Gini/agg_gini_1995.tif")
agg_map_00 <- raster("Smoothed Gini/agg_gini_2000.tif")
agg_map_05 <- raster("Smoothed Gini/agg_gini_2005.tif")
agg_map_10 <- raster("Smoothed Gini/agg_gini_2010.tif")
agg_map_15 <- raster("Smoothed Gini/agg_gini_2015.tif")
bl_gini_map_90 <- raster("Smoothed Gini/bl_gini_1990.tif")
bl_gini_map_95 <- raster("Smoothed Gini/bl_gini_1995.tif")
bl_gini_map_00 <- raster("Smoothed Gini/bl_gini_2000.tif")
bl_gini_map_05 <- raster("Smoothed Gini/bl_gini_2005.tif")
bl_gini_map_10 <- raster("Smoothed Gini/bl_gini_2010.tif")
bl_gini_map_15 <- raster("Smoothed Gini/bl_gini_2015.tif")
bl_cmap_9095 <- raster("Smoothed Gini/c_bl_gini_1990_1995.tif")
bl_cmap_9500 <- raster("Smoothed Gini/c_bl_gini_1995_2000.tif")
bl_cmap_0005 <- raster("Smoothed Gini/c_bl_gini_2000_2005.tif")
bl_cmap_0510 <- raster("Smoothed Gini/c_bl_gini_2005_2010.tif")
bl_cmap_1015 <- raster("Smoothed Gini/c_bl_gini_2010_2015.tif")
bl_cmap_9010 <- raster("Smoothed Gini/c_bl_gini_1990_2010.tif")
map_90 <- raster("Smoothed Gini/focal/gini_1990.tif")
map_95 <- raster("Smoothed Gini/focal/gini_1995.tif")
map_00 <- raster("Smoothed Gini/focal/gini_2000.tif")
map_05 <- raster("Smoothed Gini/focal/gini_2005.tif")
map_10 <- raster("Smoothed Gini/focal/gini_2010.tif")
map_15 <- raster("Smoothed Gini/focal/gini_2015.tif")
cmap_9095 <- raster("Smoothed Gini/focal/c_gini_1990_1995.tif")
cmap_9500 <- raster("Smoothed Gini/focal/c_gini_1995_2000.tif")
cmap_0005 <- raster("Smoothed Gini/focal/c_gini_2000_2005.tif")
cmap_0510 <- raster("Smoothed Gini/focal/c_gini_2005_2010.tif")
cmap_1015 <- raster("Smoothed Gini/focal/c_gini_2010_2015.tif")
cmap_9010 <- raster("Smoothed Gini/focal/c_gini_1990_2010.tif")

#WB data
wb_data <- read.csv("World bank data/WDIData.csv")
wb_data <- melt(wb_data)
wb_data1 <- dcast(wb_data, ?..Country.Name + Country.Code + variable ~ Indicator.Name)
yearlights_1990_calib <- calc(lights_1992, na_fun_light, filename = "Light/Lights 1990 calib.tif", format = "GTiff", overwrite = TRUE)
pop_calib_1990 <- calc(pop_count_1990, na_fun_pop, filename = "Pop/Pop Count/Pop 1990 calib.tif", format = "GTiff", overwrite = TRUE)
pop_calib_1990 <- spatial_sync_raster(pop_calib_1990, lights_1990_calib)
writeRaster(pop_calib_1990, filename = "Pop/Pop Count/Pop 1990 calib.tif", format = "GTiff", overwrite = TRUE)
lpp_calib_1990 <- overlay(lights_1990_calib, pop_calib_1990, fun = div, filename = "LPP/LPP 1990.tif", format = "GTiff", overwrite = TRUE)

#Functions----
dif <- function(x, y) {x - y}
div <- function (x,y) {x/y}
na_fun_pop <- function (x) {
  x[x <= 1] <- NA
  return(x)
} 
na_fun_light <- function (x) {
  x[x <= 0 | x == 255] <- NA
  return(x)
} 
na_fun_gt0 <- function (x) {x[x <= 0] <- NA; return(x)} 
na_fun_gt10 <- function (x) {x[x <= 10] <- NA; return(x)}
topone <- function(x) {
  s <- sort(x, decreasing = T)
  n <- round(length(x) * 0.01)
  return(sum(s[1:n])/sum(s))
} 
standardize <- function(x) {
  return ((x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T)))
}
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
calib <- function(x, c0 = -0.6201, c1 = 1.3504, c2 = -0.0049) {
  x[x <= 0 | x == 255] <- NA
  x <- c0 + c1*x + c2*x^2
  ifelse (x >= 20, 10^(log10(x) + 0.0684*(-1/log10(x/66) - 1.92477)^0.74), x)
}
sat_corr <- function(x) {
  ifelse (x >= 20, 10^(log10(x) + 0.0684*(-1/log10(x/66) - 1.92477)^0.74), x)
}
#Raupach, M. R., Rayner, P. J., & Paget, M. (2010). Regional variations in spatial structure of nightlights, population density and fossil-fuel CO 2 emissions. Energy Policy, 38(9), 4756-4764.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#SWIID data----
swiid_summary$ISO3 <- countrycode(swiid_summary$country, origin = "country.name", destination = "iso3c")
swiid_summary <- swiid_summary[, c(1, 13, 2:12)]
swiid_15 <- swiid_summary[swiid_summary$year == 2015,]
swiid_15 <- merge(swiid_15, df_inequality_2015_c, by = "ISO3", all.x = T)
swiid_10 <- swiid_summary[swiid_summary$year == 2010,]
swiid_10 <- merge(swiid_10, WIID_10, by = "ISO3", all.x = T); swiid_10 <- swiid_10[!is.na(swiid_10$ISO3),]
swiid_05 <- swiid_summary[swiid_summary$year == 2005,]
swiid_06 <- swiid_summary[swiid_summary$year == 2006,]
swiid_06 <- swiid_06[!is.na(swiid_06$ISO3),] 
df <- merge(world_2006, swiid_06, by = "ISO3")
swiid_df <- swiid_summary[swiid_summary$year %in% c(1990, 1995, 2000, 2005, 2010, 2015),]
swiid_df <- merge(swiid_df, swiid_10[,c("ISO3", "year", "light_gini_lpp")], by = c("ISO3", "year"), all.x = T)

#WIID data ----
WIID <- read.xlsx("United Nations University's World Income Inequality Database/WIID3.4_19JAN2017New.xlsx")
WIID_10_exp <- WIID[WIID$Year == 2010, c(1, 2, 4, 5, 6, 7, 8, 9, 29)]
WIID_10 <- unique(WIID_10_exp[,c(1:5)])
WIID_10$mean_income <- as.vector(by(WIID_10_exp$Mean, WIID_10_exp$ISO3, function(x) mean(x, na.rm = T)))
WIID_10$median_income <- as.vector(by(WIID_10_exp$Median, WIID_10_exp$ISO3, function(x) mean(x, na.rm = T)))
WIID_10$Gini <- as.vector(by(WIID_10_exp$Gini, WIID_10_exp$ISO3, function(x) mean(x, na.rm = T)))
WIID_10$Quality <- as.vector(by(WIID_10_exp$Quality, WIID_10_exp$ISO3, Mode))

#Country rasters----
country_rasters <- function(year, fname, address) {
n <- nrow(wb_data[wb_data$Year == year,])
for (i in 1:n) {
  cat("iteration =", i, "\n")
  id <- as.character(wb_data[wb_data$Year == year, ]$ISO3[i])
  ras <- crop(fname, world[world$ISO3 == id,])
  r <- rasterize(world[world$ISO3 == id,], ras)
  l <- mask(ras, r, filename = paste(address, id, ".tif", sep = ""), format = "GTiff", overwrite = T)
}
}

sol_calc <- function(year, address) {
  n <- nrow(swiid_df[swiid_df$year == year,])  
  for (i in 1:n) {
    cat("iteration =", i, "\n")
    id <- swiid_df[swiid_df$year == year,]$ISO3[i]
    l <- raster(paste(address, id, ".tif", sep = ""))
    swiid_df[swiid_df$year == year,]$sol[i] <<- cellStats(l, stat = 'sum')
  }
}

gini_calc <- function(year, address) {
  n <- nrow(wb_data[wb_data$Year == year,])  
  for (i in 1:n) {
    cat("iteration =", i, "\n")
    id <- as.character(wb_data[wb_data$Year == year, ]$ISO3[i])
    l <- raster(paste(address, id, ".tif", sep = ""))
    v <- getValues(l)
    wb_data[wb_data$Year == year,]$light_gini[i] <<- Gini(v)*100
  }
}

#Global Gini Maps by Country----
world_2010 <- merge(world, swiid_df[swiid_df$year == 2010,], by = "ISO3", all.x = T)
spplot(world_2010, c("gini_disp"), col.regions = colorRampPalette(brewer.pal(9, 'YlOrRd'))(20), 
       cuts = 19, col = "grey")
spplot(world_2010, c("light_gini_lpp"), col.regions = colorRampPalette(brewer.pal(9, 'YlOrRd'))(20), 
       cuts = 19, col = "grey")
spplot(world_2010, c("gini_disp_se"), col.regions = colorRampPalette(brewer.pal(9, 'Blues'))(20), 
       cuts = 19, col = "grey")
spplot(world_2010, c("sol"), col.regions = colorRampPalette(brewer.pal(9, 'Blues'))(20), 
       cuts = 19, col = "grey")


#World Wealth Gini distribution
ggplot() + geom_histogram(data = ww_gini, aes(Gini), bins = 165)
ggplot() + geom_histogram(data = ww_gini, aes(log(T_Wealth)), bins = 50)
ggplot() + geom_histogram(data = ww_gini, aes(Wealth_C), bins = 50)
ggplot() + geom_histogram(data = ww_gini, aes(Wealth_A), bins = 50)

#Wealth-Gini
ggplot(data = ww_gini, aes(log(T_Wealth), Gini)) + 
  geom_text(aes(label = Country, color = D_status), size = 2) 
ggplot(data = ww_gini[which(ww_gini$D_status == "T"), ], aes(log(Wealth_C), Gini)) + 
  geom_text(aes(label = Country), size = 2) + geom_smooth(method = "glm")
ggplot(data = ww_gini[which(ww_gini$D_status == "T"), ], aes(log(T_Wealth), Gini)) + 
  geom_text(aes(label = Country), size = 2) + geom_smooth(method = "glm")
ggplot() + 
  geom_histogram(data = ww_gini[which(ww_gini$D_status == "T"), ], aes(Gini), bins = 100)

#Relationship between Gini and Light Gini---- 
corr_gini <- function(year, se = 10, sol = 0) {
  s <- seq(0, 1000000, 100000)
  c1 <- vector(length = length(s))
  for (i in 1:length(s)) {
    c1[i] <- cor.test(~ gini_disp + light_gini_lpp, data = swiid_df[swiid_df$year == year & swiid_df$sol > s[i],])$estimate
  }
  df_cor <- as.data.frame(cbind(s, c1))
  g1 <- ggplot(data = df_cor, aes(s, c1)) + geom_point() + 
    labs(x = "Sum of Light (SOL)", y = "Correlation", title = "Correlations for countries with greater SOL")
  s <- seq(1, 4, 0.1)
  c1 <- vector(length = length(s))
  for (i in 1:length(s)) {
    c1[i] <- cor.test(~ gini_disp + light_gini_lpp, data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < s[i],])$estimate
  }
  df_cor <- as.data.frame(cbind(s, c1))
  g2 <- ggplot(data = df_cor, aes(s, c1)) + geom_point() + 
    labs(x = "SE Income Gini", y = "Correlation", title = "Correlations for countries with lesser SE for Light Gini")
  c <- cor.test(~ gini_disp + light_gini_lpp, data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol > sol,])
  g <- ggplot(data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol > sol,], aes(light_gini_lpp, gini_disp)) + 
    geom_text(aes(label = country), size = 3) + stat_smooth(method = lm, se = F) + 
    labs(x = "Light Gini", y = "Income Gini", title = paste("By Country Inequality", year), 
    subtitle = paste(nrow(swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol > sol,]), " Countries", "; Correlation = ", round(c$estimate, 2), ", p-value < 0.001", sep = ""))
  print(g2)
  print(g1)
  print(g)
  c
}

corr_gini_gt10 <- function(year, se = 10, sol_gt10 = 0) {
  s <- seq(0, 1000000, 100000)
  c1 <- vector(length = length(s))
  for (i in 1:length(s)) {
    c1[i] <- cor.test(~ gini_disp + light_gini_lpp_gt10, data = swiid_df[swiid_df$year == year & swiid_df$sol_gt10 > s[i],])$estimate
  }
  df_cor <- as.data.frame(cbind(s, c1))
  g1 <- ggplot(data = df_cor, aes(s, c1)) + geom_point() + 
    labs(x = "Sum of Light (sol_gt10)", y = "Correlation", title = "Correlations for countries with greater sol_gt10")
  s <- seq(1, 4, 0.1)
  c1 <- vector(length = length(s))
  for (i in 1:length(s)) {
    c1[i] <- cor.test(~ gini_disp + light_gini_lpp_gt10, data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < s[i],])$estimate
  }
  df_cor <- as.data.frame(cbind(s, c1))
  g2 <- ggplot(data = df_cor, aes(s, c1)) + geom_point() + 
    labs(x = "SE Income Gini", y = "Correlation", title = "Correlations for countries with lesser SE for Light Gini")
  c <- cor.test(~ gini_disp + light_gini_lpp_gt10, data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol_gt10 > sol_gt10,])
  g <- ggplot(data = swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol_gt10 > sol_gt10,], aes(light_gini_lpp_gt10, gini_disp)) + 
    geom_text(aes(label = country), size = 3) + stat_smooth(method = lm, se = F) + 
    labs(x = "Light Gini", y = "Income Gini", title = paste("By Country Inequality", year), 
         subtitle = paste(nrow(swiid_df[swiid_df$year == year & swiid_df$gini_disp_se < se & swiid_df$sol_gt10 > sol_gt10,]), " Countries", "; Correlation = ", round(c$estimate, 2), ", p-value < 0.001", sep = ""))
  print(g2)
  print(g1)
  print(g)
  c
}

gini_ts_country <- function(country) {
  ggplot(data = swiid_df[swiid_df$ISO3 == country,]) + 
    geom_point(aes(x = year, y = gini_disp, color = "Gini"), size = 3) + 
    geom_line(aes(x = year, y = gini_disp, color = "Gini")) +
    geom_point(aes(x = year, y = light_gini_lpp, color = "Light Gini"), size = 3) + 
    geom_line(aes(x = year, y = light_gini_lpp, color = "Light Gini")) +
    labs(colour = "Income Inequality") + xlim(1990, 2010)
}

gini_ts_country_sep <- function(country) {
  g1 <- ggplot(data = swiid_df[swiid_df$ISO3 == country,]) + 
    geom_point(aes(x = year, y = gini_disp)) +
    labs(x = "Years", y = "Income Gini", title = paste("Inequality Timeseries", country))
  g2 <- ggplot(data = swiid_df[swiid_df$ISO3 == country,]) +
    geom_point(aes(x = year, y = light_gini_lpp)) +
    labs(x = "Years", y = "Light Gini", title = paste("Inequality Timeseries", country))
  multiplot(g1, g2, cols = 2)
}

ggplot(data = swiid_df[swiid_df$region == "Europe",]) +
  geom_line(aes(x = year, y = light_gini_lpp, group = country, color = country))

model <- lm(gini_disp ~ light_gini_lpp, data = swiid_10[swiid_10$sol > 700000,]); summary(model)
model <- lm(gini_disp ~ light_gini_lpp + as.numeric(factor(ISO3))-1, data = swiid_10[swiid_10$sol > 700000,]); summary(model)

c <- NA
s <- seq(0, 1000000, 100000)
for (i in 1:length(s)) {
  c[i] <- cor(swiid_10[swiid_10$sol > s[i], c("light_gini_lpp", "gini_disp")])[2]
}
df_cor <- data.frame(sol = s, cor = c)
ggplot(data = df_cor, aes(sol, cor)) + geom_point()

c <- NA
s <- seq(0, 3.7, 0.1)
for (i in 1:length(s)) {
  c[i] <- cor(swiid_10[swiid_10$gini_disp_se < s[i], c("light_gini_lpp", "gini_disp")])[2]
}
df_cor <- data.frame(se = s, cor = c)
ggplot(data = df_cor, aes(se, cor)) + geom_point()

Gini_i_10 <- merge(world, swiid_10[swiid_10$sol > 300000, c(1, 4, 17, 18)], by = "ISO3", all.x = T)
spplot(Gini_i_10, c("light_gini_lpp"), main = "Light Gini 2010", sub = "80 Countries", col = "transparent", col.regions = terrain.colors(100), cuts = 99)
spplot(Gini_i_10, c("gini_disp"), main = "Income Gini 2010", sub = "80 Countries", col = "transparent", col.regions = terrain.colors(100), cuts = 99)


#Panel regression with fixed effects
df_fe <- pdata.frame(df_fe, index = c("ISO3", "year"))
fe.mod <- plm(gini_disp ~ light_gini_lpp, data = swiid_df, index = c("ISO3", "year"), model = "within"); summary(fe.mod)
re.mod <- plm(gini_disp ~ light_gini_lpp, data = swiid_df, index = c("ISO3", "year"), model = "random"); summary(re.mod)
phtest(fe.mod, re.mod)

#US Data----
US <- raster("Country rasters/Light 2015/USA_l.tif")
US_b <- readOGR(dsn = "US boundaries", layer = "cb_2016_us_nation_5m")
US_states_10 <- readOGR(dsn = "US boundaries", layer = "gz_2010_us_040_00_500k")
US_states_10 <- spTransform(US_states_10, lights_2010@crs)
US_states_pop <- read.csv("Country rasters/US states/PEP_2016_PEPANNRES_with_ann.csv", skip = 1)
US_ineq_10 <- merge(US_ineq_10, US_states_pop[,3:4], by = "NAME", all.x = T)
US_ineq_10$lpp <- US_ineq_10$sol/US_ineq_10$April.1..2010...Census
US_ineq_10[10, 7] <- NA
US_ineq_10[53, 7] <- NA
cor(US_ineq_10[2:53, c("light_gini_lpp", "Estimate..Gini.Index")], use = "complete.obs")
US_ineq_10$diff <- abs(US_ineq_10$Estimate..Gini.Index - US_ineq_10$light_gini_lpp)
cor(US_ineq_10[order(US_ineq_10$diff),][1:50, c("light_gini_lpp", "Estimate..Gini.Index")], use = "complete.obs")
cor(US_ineq_10[US_ineq_10$sol > 1e6, c("light_gini_lpp", "Estimate..Gini.Index")], use = "complete.obs")
ggplot(data = US_ineq_10[US_ineq_10$sol > 1.1e6,], aes(light_gini_lpp, Estimate..Gini.Index)) + 
  geom_text(aes(label = NAME), size = 3) + stat_smooth(method = lm, se = F) + 
  annotate("text", label = "Correlation = 0.3728564, p-value < 0.05", x = 0.43, y = 0.53, size = 3, colour = "red") + 
  labs(x = "Light Gini", y = "Income Gini", title = "US Inequality by States 2010", subtitle = "49 States")
summary(lm(Estimate..Gini.Index ~ light_gini_lpp, data = US_ineq_10[US_ineq_10$sol > 300000,]))
cor.test(~ Estimate..Gini.Index + light_gini_lpp, data = US_ineq_10[US_ineq_10$sol > 300000,])
colnames(US_ineq_10)[2] <- "GEO_ID"
US_states_10 <- merge(US_states_10, US_ineq_10, by = "GEO_ID")
US_county <- readOGR(dsn = "US boundaries/County boundaries", layer = "gz_2010_us_050_00_500k")
US_gini_county <- read.csv("American Community Survey (ACS)/County data/ACS_10_5YR_B19083_with_ann.csv", skip = 1)
USA <- raster("Country rasters/LPP 2010/USA.tif")
colnames(US_gini_county)[1] <- "GEO_ID"; crs(US_county) <- crs(USA);
US_county <- merge(US_county, US_gini_county, by = "GEO_ID")
n <- length(US_county$GEO_ID)
for (i in 1:n) {
  cat("iteration =", i, "\n")
  id <- US_county$GEO_ID[i]
  ras <- raster(paste("Country rasters/County lpp 10/", id, ".tif", sep = ""))
  p <- rasterToPoints(ras)
  US_county$sol[i] <- sum(p[,3])
}
s <- seq(1000, 1500000, 1000)
c <- NA
for (i in 1:length(s)) {
c[i] <- cor(US_ineq_10[US_ineq_10$sol > s[i], c(4, 6)], use = "p")[2]
}
plot(s, c)


s <- seq(0.3, 0.01, length.out = 30)
c <- NA
for (i in 1:length(s)) {
  c[i] <- cor(x = US_county@data[US_county@data$Margin.of.Error..Gini.Index < s[i], c(11, 9)], use = 'p')[2]
}
plot(s, c)

plot(US_county@data[US_county@data$Margin.of.Error..Gini.Index < 0.02, c(11, 9)])
summary(lm(Estimate..Gini.Index ~ light_gini, data = US_county@data[US_county@data$Margin.of.Error..Gini.Index < 0.02,]))

#Global Gini maps----
pak <- raster(paste("Country rasters/LPP 2010/", "PAK", ".tif", sep = ""))
mat <- matrix(1, nrow = 99, ncol = 99)
map_90 <- focal(lpp_calib_1990, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_1990.tif", format = "GTiff", overwrite = T)
map_95 <- focal(lpp_calib_1995, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_1995.tif", format = "GTiff", overwrite = T)
map_00 <- focal(lpp_calib_2000, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_2000.tif", format = "GTiff", overwrite = T)
map_05 <- focal(lpp_calib_2005, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_2005.tif", format = "GTiff", overwrite = T)
map_10 <- focal(lpp_calib_2010, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_2010.tif", format = "GTiff", overwrite = T)
map_15 <- focal(lpp_calib_2015, mat, fun = Gini, na.rm = T, filename = "Smoothed Gini/gini_2015.tif", format = "GTiff", overwrite = T)
cmap_9095 <- overlay(map_95, map_90, fun = dif, filename = "Smoothed Gini/c_gini_1990_1995.tif", format = "GTiff", overwrite = TRUE)
cmap_9500 <- overlay(map_00, map_95, fun = dif, filename = "Smoothed Gini/c_gini_1995_2000.tif", format = "GTiff", overwrite = TRUE)
cmap_0005 <- overlay(map_05, map_00, fun = dif, filename = "Smoothed Gini/c_gini_2000_2005.tif", format = "GTiff", overwrite = TRUE)
cmap_0510 <- overlay(map_10, map_05, fun = dif, filename = "Smoothed Gini/c_gini_2005_2010.tif", format = "GTiff", overwrite = TRUE)
cmap_1015 <- overlay(map_15, map_10, fun = dif, filename = "Smoothed Gini/c_gini_2010_2015.tif", format = "GTiff", overwrite = TRUE)
cmap_9010 <- overlay(map_10, map_90, fun = dif, filename = "Smoothed Gini/c_gini_1990_2010.tif", format = "GTiff", overwrite = TRUE)

#Gini maps
agg_map_90 <- aggregate(lpp_calib_1990, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_1990.tif", format = "GTiff", overwrite = T)
agg_map_95 <- aggregate(lpp_calib_1995, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_1995.tif", format = "GTiff", overwrite = T)
agg_map_00 <- aggregate(lpp_calib_2000, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_2000.tif", format = "GTiff", overwrite = T)
agg_map_05 <- aggregate(lpp_calib_2005, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_2005.tif", format = "GTiff", overwrite = T)
agg_map_10 <- aggregate(lpp_calib_2010, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_2010.tif", format = "GTiff", overwrite = T)
agg_map_15 <- aggregate(lpp_calib_2015, fact = 99, fun = Gini, filename = "Smoothed Gini/agg_gini_2015.tif", format = "GTiff", overwrite = T)
gini_map_90 <- disaggregate(agg_map_90, fact = 99, filename = "Smoothed Gini/gini_1990.tif", format = "GTiff", overwrite = T)
gini_map_95 <- disaggregate(agg_map_95, fact = 99, filename = "Smoothed Gini/gini_1995.tif", format = "GTiff", overwrite = T)
gini_map_00 <- disaggregate(agg_map_00, fact = 99, filename = "Smoothed Gini/gini_2000.tif", format = "GTiff", overwrite = T)
gini_map_05 <- disaggregate(agg_map_05, fact = 99, filename = "Smoothed Gini/gini_2005.tif", format = "GTiff", overwrite = T)
gini_map_10 <- disaggregate(agg_map_10, fact = 99, filename = "Smoothed Gini/gini_2010.tif", format = "GTiff", overwrite = T)
gini_map_15 <- disaggregate(agg_map_15, fact = 99, filename = "Smoothed Gini/gini_2015.tif", format = "GTiff", overwrite = T)
bl_gini_map_90 <- disaggregate(agg_map_90, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_1990.tif", format = "GTiff", overwrite = T)
bl_gini_map_95 <- disaggregate(agg_map_95, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_1995.tif", format = "GTiff", overwrite = T)
bl_gini_map_00 <- disaggregate(agg_map_00, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_2000.tif", format = "GTiff", overwrite = T)
bl_gini_map_05 <- disaggregate(agg_map_05, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_2005.tif", format = "GTiff", overwrite = T)
bl_gini_map_10 <- disaggregate(agg_map_10, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_2010.tif", format = "GTiff", overwrite = T)
bl_gini_map_15 <- disaggregate(agg_map_15, fact = 99, method = 'bilinear', filename = "Smoothed Gini/bl_gini_2015.tif", format = "GTiff", overwrite = T)

#Change maps
agg_cmap_9095 <- overlay(agg_map_95, agg_map_90, fun = dif, filename = "Smoothed Gini/c_agg_gini_1990_1995.tif", format = "GTiff", overwrite = TRUE)
agg_cmap_9500 <- overlay(agg_map_00, agg_map_95, fun = dif, filename = "Smoothed Gini/c_agg_gini_1995_2000.tif", format = "GTiff", overwrite = TRUE)
agg_cmap_0005 <- overlay(agg_map_05, agg_map_00, fun = dif, filename = "Smoothed Gini/c_agg_gini_2000_2005.tif", format = "GTiff", overwrite = TRUE)
agg_cmap_0510 <- overlay(agg_map_10, agg_map_05, fun = dif, filename = "Smoothed Gini/c_agg_gini_2005_2010.tif", format = "GTiff", overwrite = TRUE)
agg_cmap_1015 <- overlay(agg_map_15, agg_map_10, fun = dif, filename = "Smoothed Gini/c_agg_gini_2010_2015.tif", format = "GTiff", overwrite = TRUE)
agg_cmap_9010 <- overlay(agg_map_10, agg_map_90, fun = dif, filename = "Smoothed Gini/c_agg_gini_1990_2010.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_9095 <- overlay(bl_gini_map_95, bl_gini_map_90, fun = dif, filename = "Smoothed Gini/c_bl_gini_1990_1995.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_9500 <- overlay(bl_gini_map_00, bl_gini_map_95, fun = dif, filename = "Smoothed Gini/c_bl_gini_1995_2000.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_0005 <- overlay(bl_gini_map_05, bl_gini_map_00, fun = dif, filename = "Smoothed Gini/c_bl_gini_2000_2005.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_0510 <- overlay(bl_gini_map_10, bl_gini_map_05, fun = dif, filename = "Smoothed Gini/c_bl_gini_2005_2010.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_1015 <- overlay(bl_gini_map_15, bl_gini_map_10, fun = dif, filename = "Smoothed Gini/c_bl_gini_2010_2015.tif", format = "GTiff", overwrite = TRUE)
bl_cmap_9010 <- overlay(bl_gini_map_10, bl_gini_map_90, fun = dif, filename = "Smoothed Gini/c_bl_gini_1990_2010.tif", format = "GTiff", overwrite = TRUE)

#Average map
gini_stack <- stack(agg_map_90, agg_map_95, agg_map_00, agg_map_05, agg_map_10)
bl_gini_stack <- stack(bl_gini_map_90, bl_gini_map_95, bl_gini_map_00, bl_gini_map_05, bl_gini_map_10)
gini_avg <- overlay(gini_stack, fun  = mean, na.rm = T)
bl_gini_avg <- overlay(bl_gini_stack, fun  = mean, na.rm = T)

country_map <- function(fname, id, des_fname) {
  n <- as.character(world[world$ISO3 == id, ]$NAME)
  ras <- crop(fname, world[world$ISO3 == id,])
  r <- rasterize(world[world$ISO3 == id,], ras)
  l <- mask(ras, r, filename = des_fname, format = "GTiff", overwrite = TRUE)
  levelplot(l, margin = F, par.settings = myTheme, at = brk_g)
}


myTheme <- rasterTheme(jet.colors(10))
myTheme <- rasterTheme(bpy.colors(10))
my.colors = colorRampPalette(c("blue", "lightblue" , "red"))
myTheme <- rasterTheme(my.colors(10))
brk_g <- seq(0, 1, length.out = 100)
brk <- c(-1, seq(-0.5, 0.5, 0.02), 1)
mycolorkey <- list(labels = list(labels = seq(-1, 1, 0.25), at = seq(-1, 1, 0.25)))

gini_plot <- function(name, year) {
  levelplot(name, margin = F, main = paste("Global Gini Map", year), par.settings = myTheme, at = brk_g) + layer(sp.lines(world, lwd = 0.8, col = 'darkgrey'))
}

cmap_plot <- function(name, from, to) {
  levelplot(name, margin = F, main = paste("Global Gini Change Map", from, "-", to), par.settings = myTheme, at = brk, colorkey = mycolorkey) + layer(sp.lines(world, lwd = 0.8, col = 'darkgrey'))
}

#Gini by region
swiid_df[,c("ISO3", "country", "region")] <- lapply(swiid_df[,c("ISO3", "country", "region")], as.factor) 

ggplot(aes(country, light_gini_lpp), data = swiid_df[swiid_df$region == "Europe",]) + geom_point(alpha = 0.5) + 
  stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, color="red")

ggplot(aes(country, light_gini_lpp), data = swiid_df[swiid_df$region == "Asia",]) + geom_point(alpha = 0.5) + 
  stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, color="red")

ggplot(aes(region, light_gini_lpp), data = swiid_df) + geom_point(alpha = 0.25) +
  stat_summary(fun.data = "mean_sdl", colour = "red", geom = "crossbar", fatten = 3, width = 0.2) +
  labs(x = "Regions", y = "Light Gini", title = "Region-wise Inequality Comparison")

ggplot(aes(region, gini_disp), data = swiid_df) + geom_point(alpha = 0.25) +
  stat_summary(fun.data = "mean_sdl", colour = "red", geom = "crossbar", fatten = 3, width = 0.2)

wilcox.test(swiid_df[swiid_df$region == "Europe", "light_gini_lpp"], swiid_df[swiid_df$region == "Asia", "light_gini_lpp"], alternative = "two.sided") #compare two independent groups of samples. It’s used when your data are not normally distributed.
kruskal.test(light_gini_lpp ~ region, data = swiid_df) #Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups. It’s recommended when the assumptions of one-way ANOVA test are not met
pairwise.wilcox.test(swiid_df$light_gini_lpp, swiid_df$region, p.adjust.method = "BH") #pairwise comparisons between group levels with no assumption of equal variances or normality

#WB Data

inequality_plot <- function(var, year) {
ggplot(aes_string(x = var, y = "light_gini"), data = wb_data[wb_data$Year == year, ]) + geom_text(aes(label = ISO3), size = 3) + 
  stat_smooth(method = lm, se = F) + labs(y = "Light Gini", title = paste("Inequality Associations", year))
}

"`GDP growth (annual %)`"
"log(`GDP per capita (constant 2010 US$)`)"
"log(`Population density (people per sq. km of land area)`)"
"`Urban population growth (annual %)`"
"`Tax revenue (% of GDP)`"
wb_data$`Employment in industry (% of total employment)`


c <- cor.test(~ `GINI index (World Bank estimate)` + light_gini, data = wb_data[wb_data$Year == 2010, ])
ggplot(aes(`GINI index (World Bank estimate)`, light_gini), data = wb_data[wb_data$Year == 2010, ]) + geom_text(aes(label = ISO3), size = 3) + 
  stat_smooth(method = lm, se = F) + labs(y = "Light Gini", title = "Economic Growth - Income Inequality 2010")

#Variogram 
lpp_10 <- rasterToPoints(lpp_calib_2010, spatial = T) 
vgram_lpp_10 <- variogram(log(LPP_2010)~1, data = lpp_10)

pak <- raster("Country rasters/LPP 2010/PAK.tif")
v_pak <- rasterToPoints(pak, spatial = T)
vgram_pak <- variogram(log(PAK)~1, data = v_pak)
plot(vgram_pak, pch = 20, cex = 2)

china <- raster("Country rasters/LPP 2010/CHN.tif")
v_chn <- rasterToPoints(china, spatial = T)
vgram_chn <- variogram(log(CHN)~1, data = v_chn)
plot(vgram_pak, pch = 20, cex = 2)

#Distribution of inequality
gini_10 <- rasterToPoints(map_10)
mean(gini_10[,3]); sd(gini_10[,3]) 
gini_10_s <- sample(gini_10[,3], 1e7)
c_gini_9010 <- rasterToPoints(cmap_9010)
mean(c_gini_9010[,3]); sd(c_gini_9010[,3]) 
c_gini_9010_s <- sample(c_gini_9010[,3], 1e7)
mean(c_gini_9010_s); sd(c_gini_9010_s)
gini_10_nm <- normalmixEM(gini_10_s, k = 2)
plot(gini_10_nm, density = T, which = 2)
sdnorm = function(x, mean, sd, lambda){lambda*dnorm(x, mean = mean, sd = sd)}

#Bootstrapping
pak_df <- rasterToPoints(pak)
Gini(pak_df[,3])
gini <- function(d, i) {Gini(d[i])}
b <- boot(pak_df[,3], gini, 100)
plot(b)
ci <- boot.ci(b, type = "basic")
myboot <- function(d) {
  gini <- function(d, i) {Gini(d[i])}
  b <- boot(d, gini, 100)
  return(sd(b$t))
}

xx <- raster(matrix(1:10000, 100, 100))
plot(xx); text(xx)
f_xx <- focal(xx, matrix(1,3,3), Gini)
b_xx <- focal(xx, matrix(1,3,3), myboot)
plot(f_xx); text(f_xx)

#light map
lcolors = colorRampPalette(c("white", "yellow" , "orange"))
lTheme <- rasterTheme(lcolors(10))
brk_g <- seq(0, 1, length.out = 100)
brk <- c(-1, seq(-0.5, 0.5, 0.02), 1)
mycolorkey <- list(labels = list(labels = seq(-1, 1, 0.25), at = seq(-1, 1, 0.25)))
levelplot(lights_2010, margin = F, par.settings = lTheme) + 
  layer(sp.lines(world, lwd = 0.8, col = 'darkgrey'))
country <- raster('Country rasters/Light 2010/ESP.tif')
#country <- setExtent(country, ext = extent(-179, -50, 19, 71), keepres = T) #for US
levelplot(country, margin = F, par.settings = lTheme) + 
  layer(sp.lines(world[world$ISO3 == 'ESP',], lwd = 0.8, col = 'darkgrey'))

##Paper Figure
#Fig 1
spplot(world_2010, c("gini_disp"), col.regions = colorRampPalette(brewer.pal(9, 'YlOrRd'))(20), 
       cuts = 19, col = "grey", par.settings=list(fontsize=list(text=20)))
spplot(world_2010, c("light_gini_lpp"), col.regions = colorRampPalette(brewer.pal(9, 'YlOrRd'))(20), 
       cuts = 19, col = "grey", par.settings=list(fontsize=list(text=20)))
ggplot(data = swiid_df[swiid_df$year == 2010 & swiid_df$gini_disp_se < 10 & swiid_df$sol > 300000,], aes(light_gini_lpp, gini_disp)) + 
  geom_text(aes(label = country), size = 3) + stat_smooth(method = lm, se = F) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
lm_gini <- lm(gini_disp ~ light_gini_lpp, data = swiid_df[swiid_df$year == 2010 & swiid_df$gini_disp_se < 10 & swiid_df$sol > 300000,])
summary(lm_gini)

#Fig 2
hi <- WIID[WIID$income_cat == 'High Income', 'ISO3']
mi <- WIID[WIID$income_cat == 'Middle Income', 'ISO3']
li <- WIID[WIID$income_cat == 'Low Income', 'ISO3']
swiid_df[swiid_df$ISO3 %in% hi, 'income_cat'] <- 'High Income'
swiid_df[swiid_df$ISO3 %in% mi, 'income_cat'] <- 'Middle Income'
swiid_df[swiid_df$ISO3 %in% li, 'income_cat'] <- 'Low Income'
ggplot(data = swiid_df, aes(income_cat, light_gini_lpp)) + 
  geom_boxplot(colour = 'darkblue', fill = 'lightblue', size = 1) +
  labs(x = '', y = '') +
  theme_minimal(base_size = 20) + xlim(c('High Income', 'Low Income', 'Middle Income'))
ggplot(data = swiid_df, aes(income_cat, gini_disp)) + 
  geom_boxplot(colour = 'darkblue', fill = 'lightblue', size = 1) +
  labs(x = '', y = '') +
  theme_minimal(base_size = 20) + xlim(c('High Income', 'Low Income', 'Middle Income'))
ggplot(data = swiid_df, aes(region, light_gini_lpp)) + 
  geom_boxplot(colour = 'red4', fill = 'indianred1', size = 1) +
  labs(x = '', y = '') +
  theme_minimal(base_size = 20) 
ggplot(data = swiid_df, aes(region, gini_disp)) + 
  geom_boxplot(colour = 'red4', fill = 'indianred1', size = 1) +
  labs(x = '', y = '') +
  theme_minimal(base_size = 20) 
wilcox.test(swiid_df[swiid_df$income_cat == "High Income", "light_gini_lpp"], 
            swiid_df[swiid_df$income_cat == "Low Income", "light_gini_lpp"], 
            alternative = "two.sided")
wilcox.test(swiid_df[swiid_df$income_cat == "High Income", "light_gini_lpp"], 
            swiid_df[swiid_df$income_cat == "Middle Income", "light_gini_lpp"], 
            alternative = "two.sided")
by(swiid_df$light_gini_lpp, swiid_df$income_cat, mean)
by(swiid_df$light_gini_lpp, swiid_df$region, mean)

#Fig 3
ggplot(data = US_ineq_10[US_ineq_10$sol > 1.1e6,], aes(light_gini_lpp*100, Estimate..Gini.Index*100)) + 
  geom_text(aes(label = NAME), size = 3) + stat_smooth(method = lm, se = F) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
lm_usgini <- lm(Estimate..Gini.Index ~ light_gini_lpp, data = US_ineq_10[US_ineq_10$sol > 1.1e6,])
summary(lm_usgini)

US_states_10@bbox <- as.matrix(extent(usa_10))
spplot(US_states_10, c("Estimate..Gini.Index"), par.settings=list(fontsize=list(text=20)), 
       col.regions = colorRampPalette(brewer.pal(9, 'Reds'))(50), cuts = 49, col = "transparent")
spplot(US_states_10, c("light_gini"), par.settings=list(fontsize=list(text=20)),
       col.regions = colorRampPalette(brewer.pal(9, 'Reds'))(50), cuts = 49, col = "transparent")
spplot(US_county, c("Estimate..Gini.Index"), par.settings=list(fontsize=list(text=20)),
       col.regions = colorRampPalette(brewer.pal(9, 'Reds'))(50), cuts = 49, col = "transparent")
spplot(US_county, c("Estimate..Gini.Index"), par.settings=list(fontsize=list(text=20)),
       col.regions = colorRampPalette(brewer.pal(9, 'Reds'))(50), cuts = 49, col = "transparent")
spplot(US_county, c("light_gini"), par.settings=list(fontsize=list(text=20)),
       col.regions = colorRampPalette(brewer.pal(9, 'Reds'))(50), cuts = 49, col = "transparent")
reds <- rasterTheme(colorRampPalette(brewer.pal(9, 'Reds'))(50))
usa_10 <- raster("Smoothed Gini/focal/USA/usa_gini_2010.tif")
usa_10 <- setExtent(usa_10, ext = extent(-179, -50, 19, 71), keepres = T)
levelplot(usa_10, margin = F, par.settings = reds, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))

#Fig 4
levelplot(map_10, margin = F, par.settings = myTheme, at = brk_g) + 
  layer(sp.lines(world, lwd = 0.8, col = 'darkgrey'))
levelplot(cmap_9010, margin = F, par.settings = myTheme, at = brk, colorkey = mycolorkey) + 
  layer(sp.lines(world, lwd = 0.8, col = 'darkgrey'))
hist(gini_10[, 3], probability = T, border = 'blue', col = 'cornflowerblue',
     ylab = NULL, xlab = NULL, main = NULL); 
lines(density(gini_10[, 3], adjust = 5), lwd = 2, col = 'darkblue')
hist(c_gini_9010[, 3], probability = T, border = 'blue', col = 'cornflowerblue',
     ylab = "", xlab = "", main = NULL); 
lines(density(c_gini_9010[, 3], adjust = 5), lwd = 2, col = 'darkblue', lty = 2)
grid()
ggplot(data = data.frame(gini_10_s), aes(x = gini_10_s)) + 
  geom_histogram(aes(y = ..density..), fill = 'cornflowerblue', color = 'black') + 
  geom_density(adjust = 5, col = 'darkblue', size = 1, linetype = 'dashed') + 
  stat_function(fun = sdnorm, args = list(mean = gini_10_nm$mu[1], sd = gini_10_nm$sigma[1], lambda = gini_10_nm$lambda[1]), geom = "area", fill = 'red', alpha = 0.5) + 
  stat_function(fun = sdnorm, args = list(mean = gini_10_nm$mu[2], sd = gini_10_nm$sigma[2], lambda = gini_10_nm$lambda[2]), geom = "area", fill = 'green', alpha = 0.5) + 
  labs(x = "", y = "") + theme_minimal(base_size = 20)
ggplot(data = data.frame(c_gini_9010_s), aes(x = c_gini_9010_s)) + 
  geom_histogram(aes(y = ..density..), fill = 'cornflowerblue', color = 'black') + 
  geom_density(adjust = 5, col = 'darkblue', size = 1, linetype = 'dashed') + 
  labs(x = "", y = "") + theme_minimal(base_size = 20)


#Extra figs
usa_90 <- raster("Smoothed Gini/focal/USA/usa_gini_1990.tif")
usa_95 <- raster("Smoothed Gini/focal/USA/usa_gini_1995.tif")
usa_00 <- raster("Smoothed Gini/focal/USA/usa_gini_2000.tif")
usa_05 <- raster("Smoothed Gini/focal/USA/usa_gini_2005.tif")
usa_10 <- raster("Smoothed Gini/focal/USA/usa_gini_2010.tif")
usa_90 <- setExtent(usa_90, ext = extent(-179, -50, 19, 71), keepres = T)
usa_95 <- setExtent(usa_95, ext = extent(-179, -50, 19, 71), keepres = T)
usa_00 <- setExtent(usa_00, ext = extent(-179, -50, 19, 71), keepres = T)
usa_05 <- setExtent(usa_05, ext = extent(-179, -50, 19, 71), keepres = T)
usa_10 <- setExtent(usa_10, ext = extent(-179, -50, 19, 71), keepres = T)
levelplot(usa_90, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(usa_95, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(usa_00, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(usa_05, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(usa_10, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
usa_cmap_10 <- raster("Smoothed Gini/focal/USA/usa_cmap_9010.tif")
usa_cmap_10 <- setExtent(usa_cmap_10, ext = extent(-179, -50, 19, 71), keepres = T)
levelplot(usa_cmap_10, margin = F, par.settings = myTheme, at = brk_g)
chn_90 <- raster("Smoothed Gini/focal/China/chn_gini_1990.tif")
chn_95 <- raster("Smoothed Gini/focal/China/chn_gini_1995.tif")
chn_00 <- raster("Smoothed Gini/focal/China/chn_gini_2000.tif")
chn_05 <- raster("Smoothed Gini/focal/China/chn_gini_2005.tif")
chn_10 <- raster("Smoothed Gini/focal/China/chn_gini_2010.tif")
levelplot(chn_90, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(chn_95, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(chn_00, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(chn_05, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(chn_10, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
bra_90 <- raster("Smoothed Gini/focal/Brazil/bra_gini_1990.tif")
bra_95 <- raster("Smoothed Gini/focal/Brazil/bra_gini_1995.tif")
bra_00 <- raster("Smoothed Gini/focal/Brazil/bra_gini_2000.tif")
bra_05 <- raster("Smoothed Gini/focal/Brazil/bra_gini_2005.tif")
bra_10 <- raster("Smoothed Gini/focal/Brazil/bra_gini_2010.tif")
levelplot(bra_90, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(bra_95, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(bra_00, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(bra_05, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(bra_10, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
ind_90 <- raster("Smoothed Gini/focal/India/ind_gini_1990.tif")
ind_95 <- raster("Smoothed Gini/focal/India/ind_gini_1995.tif")
ind_00 <- raster("Smoothed Gini/focal/India/ind_gini_2000.tif")
ind_05 <- raster("Smoothed Gini/focal/India/ind_gini_2005.tif")
ind_10 <- raster("Smoothed Gini/focal/India/ind_gini_2010.tif")
levelplot(ind_90, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(ind_95, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(ind_00, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(ind_05, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(ind_10, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
pak_90 <- raster("Smoothed Gini/focal/PAK/pak_gini_1990.tif")
pak_95 <- raster("Smoothed Gini/focal/PAK/pak_gini_1995.tif")
pak_00 <- raster("Smoothed Gini/focal/PAK/pak_gini_2000.tif")
pak_05 <- raster("Smoothed Gini/focal/PAK/pak_gini_2005.tif")
pak_10 <- raster("Smoothed Gini/focal/PAK/pak_gini_2010.tif")
levelplot(pak_90, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(pak_95, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(pak_00, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(pak_05, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))
levelplot(pak_10, margin = F, par.settings = myTheme, at = brk_g, xlab = NULL, ylab = NULL, 
          scales = list(draw = F))

ggplot(data = swiid_df[swiid_df$ISO3 %in% c("IND", "BRA", "CHN"),], aes(year, gini_disp, color = ISO3)) + 
  geom_point(size = 4) +
  geom_line(size = 1) + 
  labs(x = 'Year', y = 'Income Gini') + xlim(1990, 2010) +
  theme_minimal() + scale_color_discrete(name = "Country")
ggplot(data = swiid_df[swiid_df$ISO3 %in% c("IND", "BRA", "CHN"),], aes(year, light_gini_lpp, color = ISO3)) + 
  geom_point(size = 4) +
  geom_line(size = 1) + 
  labs(x = 'Year', y = 'Income Gini') + xlim(1990, 2010) +
  theme_minimal() + scale_color_discrete(name = "Country")

ggplot(data = swiid_df[swiid_df$region == "Europe", ], aes(year, light_gini_lpp)) + 
  geom_boxplot() +
  labs(x = 'Year', y = 'Income Gini') +
  theme_minimal() + xlim(c("1990", "1995", "2000", "2005", "2010"))

ggplot(data = swiid_df[swiid_df$region == "Asia", ], aes(year, light_gini_lpp)) + 
  geom_boxplot() +
  labs(x = 'Year', y = 'Income Gini') +
  theme_minimal() + xlim(c("1990", "1995", "2000", "2005", "2010"))


boxplot(light_gini_lpp ~ region, data = swiid_df)
