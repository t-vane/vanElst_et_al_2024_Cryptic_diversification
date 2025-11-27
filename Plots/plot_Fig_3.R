### Script to generate Figure 3 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ape")
library("cowplot")
library("deeptime")
library("dplyr")
library("ggplot2")
library("ggridges")
library("ggstar")
library("ggtree")
library("gridExtra")
library("phytools")
library("readxl")
library("svglite")
library("tidytree")
library("treeio")

# Read shapes and fill colors data
master <- as.data.frame(read.csv("master.csv", header = TRUE))


#######################
## A) Plot phylogeny ##
#######################

# Read tree
tree <- read.tree("allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.25.noinv.treefile")

# Keep only specific tips and rename
keep_inds <- c("Msp1_MBB052", "Marn_RMR156", "Mmam_m22y16_lok", "Mmarg_TAF06.5", "Msam_15y02_mahi_S12", "Mtav_D72_2015",
               "Mmit_MBB012", "Mleh_ODY6.9", "Mruf_TRA8.85", "Mber_KIR06.34", "Mmyo_RMR70", "Mtan_RMR209",
               "Mbor_RMR116", "Msim_BET87", "Mger_f51y18_sako", "Mjol_LAKI5.24", "Mmaro_RMR131", "Mjon_A24",
               "Mmac_08y08_hely_S12", "Mgan_0007289B19", "Mman_RMR215", "Mmur_f142y17_jabe", "Mmur_KIR06.31",
               "Mgri_JMR009", "Mrav_f325y18_mari", "Mbon_f09y13_zana", "Mdan_90y04_anji_S12")
tree <- keep.tip(tree, tip = keep_inds)
list_gsub <- read.csv(sep = ";", text="
old; new
Marn_RMR156; italic('M. arnholdi')
Msp1_MBB052; paste(italic('M.'), ' sp. 1')
Mmam_m22y16_lok; italic('M. mamiratra')
Mmarg_TAF06.5; italic('M. margotmarshae')
Msam_15y02_mahi_S12; italic('M. sambiranensis')
Mtav_D72_2015; italic('M. tavaratra')     
Mber_KIR06.34; italic('M. berthae')
Mruf_TRA8.85; italic('M. rufus')
Mmyo_RMR70; italic('M. myoxinus')
Mleh_ODY6.9; italic('M. lehilahytsara')
Mtan_RMR209; italic('M. tanosi')
Mbor_RMR116; italic('M. boraha')
Msim_BET87; italic('M. simmonsi')
Mger_f51y18_sako; italic('M. gerpi')
Mjol_LAKI5.24; italic('M. jollyae')
Mmaro_RMR131; italic('M. marohita')
Mjon_A24; italic('M. jonahi')
Mmac_08y08_hely_S12; italic('M. macarthurii')
Mbon_f09y13_zana; italic('M. bongolavensis')
Mrav_f325y18_mari; italic('M. ravelobensis')
Mdan_90y04_anji_S12; italic('M. danfossi')
Mmit_MBB012; italic('M. mittermeieri')
Mgan_0007289B19; italic('M. ganzhorni')
Mman_RMR215; italic('M. manitatra')
Mmur_KIR06.31; paste(italic('M. murinus'), ' (central)')
Mmur_f142y17_jabe; paste(italic('M. murinus'), ' (north)')
Mgri_JMR009; italic('M.griseorufus')
")
for(x in 1:nrow(list_gsub)) {
  tree$tip.label <- gsub(list_gsub[x,"old"],list_gsub[x,"new"], tree$tip.label)
}

# Assign white edge colors for those behind collapsed triangles
edge_cols=as.list(rep("black", length(tree$edge[,2])+1))
edge_cols[c(10, 11, 13, 14, 23:26, 52, 53)]="white"

# Plot basic tree         
pA <- ggtree(tree, aes(color=edge_cols)) + geom_tree() + xlim_tree(0.3) +
  geom_tiplab(align=TRUE, linetype = NA, parse = TRUE, size = 6, offset = 0.02, color ="black") +
  geom_treescale(x=0.01, fontsize = 5)

# Flip or collapse selected nodes
pA <- pA %>% flip(10,11) %>% flip(3,4) %>% flip(1,2) %>% flip(13,14) %>% flip(16,17) 
pA <- pA %>%
    collapse(42, mode="max") %>% 
    collapse(43, mode="max") %>% 
    collapse(51, mode="max") %>% 
    ggtree::expand(42) %>% 
    ggtree::expand(43) %>% 
    ggtree::expand(51)

# Set tip symbols and colors
taxa <- c(
    "M. arnholdi", "M. sp1", "M. mamiratra", "M. margotmarshae", "M. sambiranensis", "M. tavaratra", "M. berthae", "M. rufus",
    "M. myoxinus", "M. lehilahytsara", "M. mittermeieri", "M. tanosi", "M. boraha", "M. simmonsi", "M. gerpi", "M. jollyae",
    "M. marohita", "M. jonahi", "M. macarthurii", "M. bongolavensis", "M. ravelobensis", "M. danfossi", "M. ganzhorni", "M. manitatra",
    "M. murinus (central)", "M. murinus (north)", "M. griseorufus"
    )
node_colors  <- sapply(taxa, function(t) master$Code[master$Taxon == t])
node_shapes  <- sapply(taxa, function(t) master$Symbol[master$Taxon == t])
node_inlets  <- sapply(taxa, function(t) master$Inlet[master$Taxon == t])

# Plot
nodes_nonsyn <- c(1:9, 11:23, 27)
nodes_syn <- c(2,11,13,17,20,23)
pA <- pA + 
    geom_star(aes(subset=(node %in% nodes_nonsyn)), starshape=node_shapes[nodes_nonsyn], size=4.5, fill=node_colors[nodes_nonsyn], color="black") +
    geom_star(aes(subset=(node %in% nodes_syn)), starshape=node_inlets[nodes_syn], size=1.75, fill="white", color="black") +
    # M. lehilahytsara
    geom_star(aes(x=p1$data$x[match(" italic('M. mittermeieri')", tree$tip.label)],y=p1$data$y[match(" italic('M. lehilahytsara')", tree$tip.label)]), starshape=master$Symbol[master$Taxon == "M. lehilahytsara"], fill=master$Code[master$Taxon == "M. lehilahytsara"],size=4.5,color="black") +
    # M. manitatra
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" italic('M. manitatra')", tree$tip.label)]), starshape=master$Symbol[master$Taxon == "M. manitatra"], fill=master$Code[master$Taxon == "M. manitatra"],size=4.5,color="black") +
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" italic('M. manitatra')", tree$tip.label)]), starshape=master$Inlet[master$Taxon == "M. manitatra"], fill="white", size=1.75,color="black") +
    # M. murinus (north)
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" paste(italic('M. murinus'), ' (north)')", tree$tip.label)]), starshape=master$Symbol[master$Taxon == "M. murinus (north)"], fill=master$Code[master$Taxon == "M. murinus (north)"],size=4.5,color="black") +
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" paste(italic('M. murinus'), ' (north)')", tree$tip.label)]), starshape=master$Inlet[master$Taxon == "M. murinus (north)"], fill="white", size=1.75,color="black") +
    # Murinus (central)
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" paste(italic('M. murinus'), ' (central)')", tree$tip.label)]), starshape=master$Symbol[master$Taxon == "M. murinus (central)"], fill=master$Code[master$Taxon == "M. murinus (central)"],size=4.5,color="black") +
    geom_star(aes(x=p1$data$x[match(" italic('M. ganzhorni')", tree$tip.label)],y=p1$data$y[match(" paste(italic('M. murinus'), ' (central)')", tree$tip.label)]), starshape=master$Inlet[master$Taxon == "M. murinus (central)"], fill="white", size=1.75,color="black")
  
svg("Fig_2_A.svg",height=10,width=12)
pA
dev.off()


########################
## B) Plot clustering ##
########################

## Set taxon subsets, number of clusters, best seed and bind into data frame
subsets <- c("lehmit", "arnsp1", "borsim", "jonmac", "bermyoruf", "sammarmam", "germarjol", "ganmanmur", "bondanrav")
K <- c(2, 2, 2, 2, 3, 3, 3, 3, 3)
seed <-c(8, 8, 6, 9, 10, 2, 10, 7, 6)
overview <- data.frame(subsets, taxa, best)

## Loop over subsets and plot
plots <- list()
subset_colors <- list(
    lehmit = c("M. lehilahytsara", "#AC6488"),
    arnsp1 = c("M. arnholdi", "#ACDBF3"),
    borsim = c("M. simmonsi", "#DB94A0"),
    jonmac = c("M. jonahi", "M. macarthurii"),
    bermyoruf = c("M. rufus","M. berthae","M. myoxinus"),
    sammarmam = c("M. sambiranensis","M. mamiratra","M. margotmarshae"),
    germarjol = c("M. jollyae","M. gerpi","#C47CB8"),
    ganmanmur = c("#DB94A0","M. murinus (north)","#FAF0F1","#EBC2C9"),
    bondanrav = c("M. ravelobensis","#AC6488","M. danfossi")
    )

for i in subsets {
    # Read data and create sorted data fram
    admix <- as.data.frame(read.table(paste0(i, ".ngsAdmix.K", overview[overview$subsets == i,]$K, ".seed", overview[overview$subsets == i,]$best, ".qopt")))
    inds <- read.table(paste0(i, ".popfile.txt"))
    admix$inds <- inds[,1]
    # Format to long
    admix_long <- reshape2::melt(admix, id.vars="inds")
    # Get colors for this subset
    taxa_or_colors <- subset_colors[[i]]
    colors_vec <- sapply(subset_colors[[i]], function(x) {
        if (x %in% master$Taxon) {
            master$Code[master$Taxon == x]
        } else {
            x  # already a hex color
        }
    })
    # Plot
    pB <- ggplot(admix_long, aes(x=inds, y=value, fill=variable)) +
        scale_fill_manual(values=colors_ved) +
        geom_bar(stat = "identity",width=1) + 
        theme_minimal() + 
        scale_x_discrete(limits = admix[order(admix[,1], decreasing=TRUE),]$inds) +
        geom_vline(xintercept=c(0.5,10.5,nrow(inds)+0.5), size=1) +
        geom_hline(yintercept=c(-Inf,Inf), size=1) + 
        theme(panel.grid.major.x = element_blank(),
            axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_line(color = "black", linewidth=1),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            aspect.ratio= 1/3,
            panel.spacing = unit(2, "cm")) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.5), expand = c(0, 0))
    plots[[i]] <- pB
}

svg("Fig_3_B.svg", height=25, width=24/3)
do.call(grid.arrange, c(plots, nrow=9))
dev.off()



#########################################
## C) Plot isolation-by-distance (IBD) ##
#########################################

## Load IBD data
ibd_all <- read_excel("ibd.xlsx")

## Create custom ordering
ibd_all$Pair <- factor(ibd_all$Pair, levels=ibd_all$Pair[order(ibd_all$Order, decreasing = TRUE)])

## Assign colors and shapes and set rectangles for background
taxa_group1 <- c("M. berthae", "M. myoxinus", "M. myoxinus", "M. lehilahytsara",
                 "M. mamiratra", "M. mamiratra", "M. sambiranensis", "M. sambiranensis",
                 "M. sp1", "M. simmonsi", "M. marohita", "M. gerpi",
                 "M. jonahi", "M. ganzhorni", "M. murinus (central)", "M. murinus (north)",
                 "M. murinus (central)", "M. murinus (north)", "M. murinus (north)",
                 "M. bongolavensis", "M. danfossi", "M. danfossi")
color1 <- sapply(taxa_group1, function(t) master$Code[master$Taxon == t])
shape1 <- sapply(taxa_group1, function(t) master$Symbol[master$Taxon == t])
inlet1 <- sapply(taxa_group1, function(t) master$Inlet[master$Taxon == t])

taxa_group2 <- c("M. rufus", "M. rufus", "M. berthae", "M. mittermeieri",
                 "M. margotmarshae", "M. mamiratra", "M. margotmarshae", "M. arnholdi",
                 "M. boraha", "M. jollyae", "M. jollyae", "M. marohita",
                 "M. macarthurii", "M. manitatra", "M. manitatra", "M. manitatra",
                 "M. ganzhorni", "M. ganzhorni", "M. murinus (central)", "M. ravelobensis",
                 "M. ravelobensis", "M. bongolavensis")
color2 <- sapply(taxa_group2, function(t) master$Code[master$Taxon == t])
shape2 <- sapply(taxa_group2, function(t) master$Symbol[master$Taxon == t])
inlet2 <- sapply(taxa_group2, function(t) master$Inlet[master$Taxon == t])

rects <- data.frame(
  xmin = c(0, 9.5, 13.5, 15.5, 19.5),
  xmax = c(3.5, 10.5, 14.5, 18.5, Inf)
)

## Plot
pC <- ggplot(ibd_all) + 
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "#f6f6f6") +
    geom_errorbar(aes(x=Pair, y=log(Mean1), ymin=log(quantile5_lower1), ymax=log(quantile95_upper1)), color = "black", width=.1, size=.1, position = position_nudge(x = -0.225)) +
    theme_classic() +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size=10), axis.title=element_text(size=10), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks = element_line(linewidth = 0.8), legend.position = "none", aspect.ratio=3, plot.margin=grid::unit(c(0,0,0,0)+0.75, "cm")) +
    geom_vline(xintercept = c(1.5, 2.5, 4.5, 5.5, 6.5, 7.5, 8.5, 11.5, 12.5, 16.5, 17.5, 20.5, 21.5), color="lightgrey") +
    geom_vline(xintercept = c(3.5, 9.5, 10.5, 13.5, 14.5, 15.5, 18.5, 19.5, Inf), color="black") +
    geom_hline(yintercept=log(0.497), linetype ="dashed", size =0.8, color="#332288") +
    geom_hline(yintercept=log(0.303), linetype ="dashed", size =0.8, color="#882255") +
    geom_hline(yintercept=log(0.439), linetype ="dashed", size =0.8, color="#ED9AC4") +
    geom_hline(yintercept=c(-2,-1,0,1,2),color="lightgrey") +
    geom_hline(yintercept=Inf,color="black") +
    geom_errorbar(aes(x=Pair, y=log(Mean1), ymin=log(quantile5_lower1), ymax=log(quantile95_upper1)), color = "black", width=.1, size=0.5, position = position_nudge(x = -0.225)) +
    geom_star(aes(x=Pair, y=log(Mean1)), starshape=rev(shape1), size=2.5, fill=rev(color1), position = position_nudge(x = -0.225)) +
    geom_star(aes(x=Pair, y=log(Mean1)), starshape=rev(inlet1), size=1.25, fill='white', position = position_nudge(x = -0.225)) +
    geom_errorbar(aes(x=Pair, y=log(Mean2), ymin=log(quantile5_lower2), ymax=log(quantile95_upper2)), color = "black", width=.1, size=.5, position = position_nudge(x = 0.225)) +
    geom_star(aes(x=Pair, y=log(Mean2)), starshape=rev(shape2), size=2.5, fill=rev(color2), position = position_nudge(x = 0.225)) +
    geom_star(aes(x=Pair, y=log(Mean2)), starshape=rev(inlet2), size=1.25, fill='white', position = position_nudge(x = 0.225)) +
    labs(y = "log NRMSE") +
    scale_y_continuous(breaks = seq(-2.0, 2.0, by = 1.0), expand = c(0, 0), labels = scales::number_format(accuracy = 0.1)) 
pC <- pC  + coord_flip(ylim = c(-2.05,2.25)) 


#################################################
## D) Plot genealogical divergence index (gdi) ##
#################################################

## Read gdi data
gdi_all <- read_excel("gdi.xlsx", sheet="gdi_conv")

## Create custom ordering
gdi_all$Pair <- factor(gdi_all$Pair, levels=gdi_all$Pair[order(gdi_all$Order, decreasing = TRUE)])

## Plot
pD <- ggplot(gdi_all) +  
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "#f6f6f6") +
    geom_errorbar(aes(x=Pair, y=Mean1, ymin=HPD95_lower1, ymax=HPD95_upper1), color = "black", width=0.1, size=0.5, position = position_nudge(x = -0.225)) +
    theme_classic() + 
    theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size=10), axis.title=element_text(size=10), axis.title.y = element_blank(), axis.ticks = element_line(linewidth = 0.8), 
        legend.position = "none", aspect.ratio=3, axis.text.y=element_blank(), plot.margin=grid::unit(c(0,0,0,0)+0.75, "cm")) +
    geom_vline(xintercept = c(1.5,2.5,4.5,5.5,6.5,7.5,8.5,11.5,12.5,16.5,17.5,20.5,21.5), color="lightgrey") +
    geom_vline(xintercept = c(3.5, 9.5,10.5,13.5,14.5,15.5,18.5,19.5,Inf), color="black") +
    geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8,1),color="lightgrey") +
    geom_hline(yintercept=c(0.25), linetype="dashed", size=0.8) +
    geom_hline(yintercept=Inf,color="black") +
    geom_errorbar(aes(x=Pair, y=Mean1, ymin=HPD95_lower1, ymax=HPD95_upper1), color = "black", width=.1, size=.5, position = position_nudge(x = -0.225)) +
    geom_star(aes(x=Pair, y=Mean1), starshape=rev(shape1), size=2.5, fill=rev(color1), position = position_nudge(x = -0.225)) +
    geom_star(aes(x=Pair, y=Mean1), starshape=rev(inlet1), size=1.25, fill='white', position = position_nudge(x = -0.225)) +
    geom_errorbar(aes(x=Pair, y=Mean2, ymin=HPD95_lower2, ymax=HPD95_upper2), color = "black", width=.1, size=.5, position = position_nudge(x = 0.225)) +
    geom_star(aes(x=Pair, y=Mean2), starshape=rev(shape2), size=2.5, fill=rev(color2), position = position_nudge(x = 0.225)) +
    geom_star(aes(x=Pair, y=Mean2), starshape=rev(inlet2), size=1.25, fill='white', position = position_nudge(x = 0.225)) +
    labs(y = expression(italic("gdi"))) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0))
pD <- pD + coord_flip(ylim = c(-0.05,1.05)) 


##########################################
## E) Plot morphometric differentiation ##
##########################################

## Load morphometric data
morpho_all <- read_excel("morphoENM.xlsx", sheet="morpho")

## Create custom ordering
morpho_all$Pair <- factor(morpho_all$Pair, levels=morpho_all$Pair[order(morpho_all$Order, decreasing = TRUE)])
morpho_all$Mean <- as.double(morpho_all$Mean)

## Assign colors
taxa <- c("ber-ruf", "myo-ruf", "myo-ber", "leh-mit", "mam-marg", "sam-mam",
                "sam-marg", "sp1-arn", "sim-bor", "mar-jol", "ger-jol", "ger-mar",
                "jon-mac", "gan-man", "mur-man", "mur-man", "mur-gan", "mur-gan",
                "mur-gan", "bon-rav", "dan-rav", "dan-bon")
colors <- sapply(taxa, function(t) master$Code[master$Taxon == t])

## Plot
pE <- ggplot(morpho_all) +  
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "#f6f6f6") +
    geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8,1),color="lightgrey") +
    geom_hline(yintercept=Inf,color="black") +
    geom_errorbar(aes(x=Pair, y=1-Mean, ymin=1-upper95, ymax=1-lower95), color = "black", width=.1, size=.1) +
    theme_classic() + 
    theme(panel.grid.major.x = element_blank(), axis.text = element_text(size=10), axis.title=element_text(size=10), axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.8), legend.position = "none", aspect.ratio=3, axis.text.y=element_blank(), plot.margin=grid::unit(c(0,0,0,0)+0.75, "cm")) +
    geom_vline(xintercept = c(1.5,2.5,4.5,5.5,6.5,7.5,8.5,11.5,12.5,16.5,17.5,20.5,21.5), color="lightgrey") +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.483, ymax = 0.566, alpha = 0.025, fill = "#332288") +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.633, ymax = 0.804, alpha = 0.025, fill = "#882255") +
    geom_vline(xintercept = c(3.5, 9.5,10.5,13.5,14.5,15.5,18.5,19.5,Inf), color="black") +
    geom_errorbar(aes(x=Pair, y=1-Mean, ymin=1-upper95, ymax=1-lower95), color = "black", width=.1, size=.5) +
    geom_point(aes(x=Pair, y=1-Mean), size = 2.5, color="black") +
    geom_point(aes(x=Pair, y=1-Mean), size = 2, color=rev(colors)) +
    labs(y = "Differentiation") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0))
pE <- pE + coord_flip(ylim = c(-0.05,1.05))



##############################################
## F) Plot ecological niche differentiation ##
##############################################
 
# Load ENM data
ENM_all <- read_excel("morphoENM.xlsx", sheet="enm")

## Create custom ordering
ENM_all$Pair <- factor(ENM_all$Pair, levels=ENM_all$Pair[order(ENM_all$Order, decreasing = TRUE)])
ENM_all$schoenersD <- as.double(ENM_all$schoenersD)

## Plot
pF <- ggplot(ENM_all) + 
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "#f6f6f6") +
    geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8,1),color="lightgrey") +
    geom_hline(yintercept=Inf,color="black") +
    geom_errorbar(aes(x=Pair, y=1-schoenersD, ymin=1-upper95, ymax=1-lower95), color = "black", width=.1, size=.1) +
    theme_classic() + theme(panel.grid.major.x = element_blank(), axis.text = element_text(size=10), axis.title=element_text(size=10), axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.8), legend.position = "none", aspect.ratio=3, axis.text.y=element_blank(), plot.margin=grid::unit(c(0,0,0,0)+0.75, "cm")) +
    geom_vline(xintercept = c(1.5,2.5,4.5,5.5,6.5,7.5,8.5,11.5,12.5,16.5,17.5,20.5,21.5), color="lightgrey") +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.574, ymax = 0.698, alpha = 0.025, fill = "#332288") +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.336, ymax = 0.457, alpha = 0.025, fill = "#882255") +
    geom_vline(xintercept = c(3.5, 9.5,10.5,13.5,14.5,15.5,18.5,19.5,Inf), color="black") +
    geom_errorbar(aes(x=Pair, y=1-schoenersD, ymin=1-upper95, ymax=1-lower95), color = "black", width=.1, size=.5) +
    geom_point(aes(x=Pair, y=1-schoenersD), size = 2.5, color="black") +
    geom_point(aes(x=Pair, y=1-schoenersD), size = 2, color=rev(colors)) +
    labs(y = "Differentiation") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0, 0))
pF <- pF + coord_flip(ylim = c(-0.05,1.05))

svg("Fig_S3_CDEF.svg",height=10,width=4*8/3)
plot_grid(pC, pD, pE, pF, align = "h", ncol = 4, rel_heights = c(1/4, 1/4, 1/4, 1/4))
dev.off()
