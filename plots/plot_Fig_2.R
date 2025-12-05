### Script to generate Figure 2 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ape")
library("deeptime")
library("ggplot2")
library("ggsn")
library("ggstar")
library("ggthemes")
library("ggtree")
library("mapdata")
library("mapplots")
library("maps")
library("maptools")
library("phytools")
library("readxl")
library("raster")
library("rgeos")
library("rgdal")
library("scales")
library("sp")
library("treeio")

#################
## A) Plot map ##
#################

## Read necessary geography data
# World map
worldmap <- map_data('worldHires')
# Shape files (dry and humid forest)
dryforest_df <- fortify(readOGR(shape_folder, "Mada_dry_forest_shape"))
humidforest_df <- fortify(readOGR(shape_folder, "Mada_humid_forest_shape"))

## Prepare data
# Read coordinate data
dat <- read_excel(paste0("Microcebus_geno_localities.xlsx"), col_names=TRUE)
dat <- na.omit(dat)
dat$Latitude <- as.numeric(dat$Latitude)
dat$Longitude <- as.numeric(dat$Longitude)

## Set shapes and fill colors
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

## Plot map
svg("Fig_1_A.svg", width = 6, height = 9)
pA <- ggplot() + 
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
        pA <- pA + geom_star(data=dat[i,], aes(x=Longitude, y=Latitude), starshape=shapes[i], fill=colors[i], size=4, colour="black", alpha=.65)
        if (!is.na(dot_shapes[i])) {
            pA <- pA + geom_star(data=dat[i,], aes(x=Longitude, y=Latitude), starshape=dot_shapes[i], fill="white", colour="black", size=1.5, alpha=.65)
        }
    }

# Add north arrow
north2(pA, symbol= 3, x=.7, y=.22, scale=0.08)
dev.off()


#######################
## B) Plot phylogeny ##
#######################

## Set renaming map
list_gsub <- read.csv(sep = ";", text="
old; new
Mirzazaza; italic('Mirza zaza')
Marnholdi; italic('M. arnholdi')                 
Msp1; paste(italic('M.'), ' sp. 1')      
Mmamiratra; italic('M. mamiratra')                  
Mmargotmarshae; italic('M. margotmarshae')                 
Msambiranensis; italic('M. sambiranensis')
Mtavaratra; italic('M. tavaratra')               
Mberthae; italic('M. berthae')
Mrufus; italic('M. rufus')
Mmyoxinus; italic('M. myoxinus')
MlehilahytsaraN; italic('M. lehilahytsara')
Mtanosi; italic('M. tanosi')
Mboraha; italic('M. boraha')
Msimmonsi; italic('M. simmonsi')
Mgerpi; italic('M. gerpi')
Mjollyae; italic('M. jollyae')
Mmarohita; italic('M. marohita')
Mjonahi; italic('M. jonahi')
Mmacarthurii; italic('M. macarthurii')
Mbongolavensis; italic('M. bongolavensis')
Mravelobensis; italic('M. ravelobensis')
Mdanfossi; italic('M. danfossi')
Mmittermeieri; italic('M. mittermeieri')
Mganzhorni; italic('M. ganzhorni')
Mmanitatra; italic('M. manitatra')
MmurinusC; paste(italic('M. murinus'), ' (central)')
MmurinusN; paste(italic('M. murinus'), ' (north)')
Mgriseorufus; italic('M.griseorufus')
")

## Read and clean tree
tree <- read.beast("BPP_mean_FigTree1.tre")
tree <- treeio::drop.tip(tree, tip = c("Mirzazaza", "MlehilahytsaraS"))

## Rename tips
for(x in 1:nrow(list_gsub)) {
  tree@phylo$tip.label <- gsub(list_gsub[x,"old"],list_gsub[x,"new"], tree@phylo$tip.label)
}

## Set node ages and letters
q <- ggtree(tree)
node_ages <- round(max(q$data$x) - q$data$x)
node_letters <- c(rep(NA, 28),letters[1:25])

## Plot
pdf("Fig_2_B.pdf", height = 10, width = 15)
pB <- ggtree(tree) + geom_tree() + geom_tiplab(align=TRUE, parse = TRUE, linetype='dashed', size = 5) + xlim_tree(200) +
    geom_vline(xintercept = c(-1500,-1250,-1000,-750,-500,-250), color="gray",linetype = "dashed") + 
    geom_range("height_95_HPD", color='red', size=3, alpha=.5) +  
    geom_text2(aes(subset=(!isTip), label=node_ages), color = "black", hjust = -0.22, size=5) +
    geom_text2(aes(subset=(node > 27), label=node_letters), color = "dimgray", hjust = 2, vjust = -0.5, size = 5) + 
    theme_tree2() + 
    xlab(label = "Ka ago") + 
    theme(axis.text.x = element_text(size = 15),axis.title.x = element_text(size = 15)) +
    scale_x_continuous(breaks=seq(-1500,0,250), labels=abs(seq(-1500,0,250)))
revts(pB)
dev.off()


######################################################################
## C+D) Plot species delimiation and conservation status rectangles ##
######################################################################

## Set tip limits four current (cc) and revised (rc) classification
tips_cc <- setdiff(1:27, 6)
tips_cc_b <- setdiff(1:27, 5)
tips_cc_m<-vector()
for(i in 1:length(tips_cc)){tips_cc_m<-append(tips_cc_m,mean(c(tips_cc[i],tips_cc_b[i])))}

tips_rc <- setdiff(1:27, c(3, 6, 7, 8, 13, 15, 18, 24))
tips_rc_b <- setdiff(1:27, c(2, 5, 6, 7, 12, 14, 17, 23))
tips_rc_m<-vector()
for(i in 1:length(tips_rc)){tips_rc_m<-append(tips_rc_m,mean(c(tips_rc[i],tips_rc_b[i])))}

## Set conservation status (reverse order as in figure)
cons_status_cc<-c("VU","EN","VU","LC","LC","EN","CR","EN","EN","CR","CR","EN","EN","DD","VU","NE","VU","EN","EN","EN","EN","EN","NT","VU","CR","VU")
cons_status_rc<-c("VU","NT","LC","LC","EN","NT","EN","EN","VU","VU","VU","EN","EN","EN","NT","NT","VU","CR","NT")
cons_status_cc<-factor(cons_status_cc, levels=c("CR","EN","VU","NT","LC","DD"))
cons_status_rc<-factor(cons_status_rc, levels=c("CR","EN","VU","NT","LC","DD"))

## Set required node color and shape data
node_colors <- c("#88CCEE","#882255","#882255","#117733","#CC6677","#CC6677","#CC6677","#CC6677","#999933","#DDDDDD","#AA4499","#DDDDDD","#DDDDDD","#CC6677","#CC6677",
    "#332288","#88CCEE","#88CCEE","#DDDDDD","#AA4499","#CC6677","#999933","#882255","#882255","#332288","#117733","#88CCEE")
cons_status_colors<-c("#c52412","#f28533","#ffc90e","#97c115","#51bc1d","#c3c3c3")
nodes_syn <- c(2,6,7,8,12,15,17,23)
nodes_syn_shape <- c(5,15,13,28,11,5,28,13)

## Plot
lsz=0.4

svg("Fig_2_CD.svg", height = 10, width = 2.5)
pCD <- ggplot() + 
    geom_rect(aes(xmin = c(1.6), xmax =c(2.4), ymin = c(tips_cc)-0.45, ymax = c(tips_cc_b)+0.45), fill = node_colors[tips_cc], alpha = 1, color = "black", size=lsz) + 
    geom_star(aes(2,nodes_syn), starshape=nodes_syn_shape, size=1.75, fill='white') +
    
    geom_rect(aes(xmin = c(2.6), xmax =c(3.4), ymin = c(tips_rc)-0.45, ymax = c(tips_rc_b)+0.45), fill = node_colors[tips_rc], alpha = 1, color = "black", size=lsz) + 
    geom_star(aes(3,nodes_syn), starshape=nodes_syn_shape, size=1.75, fill='white') +
    
    geom_rect(aes(xmin = c(4.6), xmax =c(5.4), ymin = c(tips_cc)-0.45, ymax = c(tips_cc_b)+0.45), fill = cons_status_colors[as.numeric(cons_status_cc)], alpha = 1, color = "black", size=lsz) + 
    annotate("text", x=8, y=tips_cc_m, label= cons_status_cc, size=3.5) +
    
    geom_rect(aes(xmin = c(5.6), xmax =c(6.4), ymin = c(tips_rc)-0.45, ymax = c(tips_rc_b)+0.45), fill = cons_status_colors[as.numeric(cons_status_rc)], alpha = 1, color = "black", size=lsz) + 
    annotate("text", x=9, y=tips_rc_m, label= cons_status_rc, size=3.5) +
    
    xlim(1.55,9.45) + 
    ylim(0,27.6) + 
    theme_void() +
    theme(plot.margin = margin(t = .1, r = .1, b = .1, l = .1, unit = "in")) + 
    annotate("text", x=c(2:9), y=0., label= c("CT","RT","CT","RT"), fontface =c(1,2,1,2))
pCD
dev.off()
