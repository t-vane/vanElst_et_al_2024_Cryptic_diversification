### Script to generate Supplementary Figure S22 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ggplot2")
library("ggstar")
library("ggsn")
library("ggthemes")
library("mapdata")
library("mapplots")
library("maps")
library("maptools")
library("patchwork")
library("raster")
library("readxl")
library("rgdal")
library("rgeos")
library("scales")
library("sp")

## Read necessary geography data
# World map
worldmap <- map_data('worldHires')
# Shape files (dry and humid forest)
dryforest_df <- fortify(readOGR(shape_folder, "Mada_dry_forest_shape"))
humidforest_df <- fortify(readOGR(shape_folder, "Mada_humid_forest_shape"))

## Read coordinates and plot occurrences for each data set
sets <- c("geno", "morpho", "ENM", "repro", "acoustic")
plots <- list()

for s in sets {
    # Prepare data
    dat <- read_excel(paste0("Microcebus_", s, "_localities.xlsx"), col_names=TRUE)
    dat <- na.omit(dat)
    dat$Latitude <- as.numeric(dat$Latitude)
    dat$Longitude <- as.numeric(dat$Longitude)

    # Set shapes and fill colors
    master <- as.data.frame(read.csv("master.csv", header = TRUE))

    # Construct shape vector
    shapes <- c()
    for (elem in unique(dat$Taxon)[]) {
        for (i in seq(from=1, to=nrow(dat[dat$Taxon==elem,]))) {
            shapes <- c(shapes,master$Symbol[master$Taxon==elem])
        }
    }
    # Construct color vector
    colors <- c()
    for (elem in unique(dat$Taxon)[]) {
        for (i in seq(from=1, to=nrow(dat[dat$Taxon==elem,]))) {
            colors <- c(colors,master$Code[master$Taxon==elem])
        }
    }
    # Construct dot vector
    dot_shapes <- c()
    for (elem in unique(dat$Taxon)[]) {
        for (i in seq(from=1, to=nrow(dat[dat$Taxon==elem,]))) {
            dot_shapes <- c(dot_shapes,master$Inlet[master$Taxon==elem])
        }
    }

    # Plot map
    p <- ggplot() + 
    coord_fixed(xlim = c(43.15,50.75), ylim = c(-25.6075,-12.25)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme_map() +
    geom_polygon(data=dryforest_df, aes(long, lat, group = group), fill="#FFD27F", alpha = 1) +
    geom_polygon(data=humidforest_df, aes(long, lat, group = group), fill = "#7FB17E", alpha = 1) +
    geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill=NA, colour = 'black') +
    scale_y_continuous(breaks = seq(-25, -13, by = 2)) +
    xlab("\nLongitude") + 
    ylab("Latitude\n") +
    ggsn::scalebar(data=worldmap, transform=TRUE, dist=100, dist_unit="km", model="WGS84", location ="topleft",
        anchor = c(x = 48, y = -25), height = 0.001, st.bottom=FALSE, st.dist = 0.0025, st.size =3,
        border.size=0.5)
    # Plot samples
    for (i in seq(from=1, to=nrow(dat))) {
        print(i)
        p <- p + geom_star(data=dat[i,], aes(x=Longitude, y=Latitude), starshape=shapes[i], fill=colors[i], size=4, colour="black", alpha=.65)
        if (!is.na(dot_shapes[i])) {
            p <- p + geom_star(data=dat[i,], aes(x=Longitude, y=Latitude), starshape=dot_shapes[i], fill="white", colour="black", size=1.5, alpha=.65)
        }
    }
    plots[[s]] <- p
}

pdf("Fig_S22.pdf", width = 18, height = 18)
(plots$geno | plots$morpho | plots$ENM) /
(plots$repro | plots$acoustic)
dev.off()
