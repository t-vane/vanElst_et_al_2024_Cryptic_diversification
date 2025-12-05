### Script to generate Supplementary Figures S1-S5 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ape")
library("dplyr")
library("ggnewscale")
library("ggplot2")
library("ggtree")
library("phytools")

## Set root and outgroup samples
root <- c("Cheirogaleuscro_106189", "Cheirogaleusmaj_106245", "Cheirogaleusmed_106354")
out <- c("Cheirogaleuscro_106189", "Cheirogaleusmaj_106245", "Cheirogaleusmed_106354","Mirzazaz_DLC2316m","Mirzazaz_DLC319m","Mirzazaz_DLC323f")

## Read coloration file
colors_prelim <- as.data.frame(read.csv("colors.csv", header = TRUE))

## Loop over five different thresholds of missing data
thresholds <- c(0.05, 0.25, 0.5, 0.75, 0.95)

for th in thresholds {
    # Read trees
    tree1_unrooted <- read.tree(paste0("allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss",th,".noinv.treefile"))
    tree2_unrooted <- read.tree(paste0("final.maxmiss",th,".thinned.svdq.tre"))

    # Root and ladderize
    tree1 <- ladderize(ape::root(tree1_unrooted, root))
    tree2 <- ladderize(ape::root(tree2_unrooted, root))

    # Drop outgroups to improve visibility
    tree1 <- drop.tip(tree1, tip = out)
    tree2 <- drop.tip(tree2, tip = out)

    # Set coloring
    colors1 <- colors_prelim %>% arrange(factor(individual, levels = tree1$tip.label))
    colors2 <- colors_prelim %>% arrange(factor(individual, levels = tree2$tip.label))
    colors_lines <- rbind(colors1,colors2)
    coloration <- colors_lines$color

    # Generate plot objects and save associated data
    p1 <- ggtree(tree1) + 
        geom_tree(size = 1.2) + 
        geom_treescale() + 
        geom_tiplab(align=TRUE, linetype='dashed', hjust=-0.05, size = 7, color=colors1$color) + 
        xlim_tree(0.25) + 
        geom_text2(aes(subset=(label!="100/100" & !isTip), label = label, size = 5, hjust = -0.08))
    d1 <- p1$data
    p2 <- ggtree(tree2) + 
        geom_treescale() + 
        geom_nodelab() + 
        geom_tiplab(align=TRUE, linetype='dashed', hjust=-0.05, size = 5, color=colors2$color) + 
        xlim_tree(0.25)
    d2 <- p2$data

    # Reverse x-axis and set offset for right tree
    d2$x <- d2$x/200
    d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.1

    # Bind data frames
    dd <- bind_rows(d1, d2) %>%  filter(!is.na(label))
    # Filter out non-tips
    dd <- dd[dd$isTip == TRUE,]
    # Set dd$x to desired position for left and right tree
    dd$x[1:208] <- 0.09
    dd$x[209:416] <- dd$x[209:416] - 0.025

    # Final plot
    pdf(paste0("Figs_S1-S5_",th, ".pdf"), height = 60, width = 40)
    p1 + 
        geom_tree(data=d2, size = 1.2) + 
        ggnewscale::new_scale_fill() + 
        geom_tiplab(data=d2, hjust=1.1, size = 7, color=colors2$color) +
        geom_text2(data=d2, aes(subset=(label!=100 & !isTip), label = label, size = 5, hjust = 1.65)) + 
        geom_line(aes(x,y, group=label, color=coloration), data=dd) + 
        scale_color_manual(values=c("#00FF00"="#00FF00", "#FF0000"="#FF0000", "#DF73FF"="#DF73FF", "#38A800"="#38A800",
                                     "#FFAA00"="#FFAA00", "#0070FF"="#0070FF", "#BED2FF"="#BED2FF", "#73DFFF"="#73DFFF",
                                     "#FFFF73"="#FFFF73", "#FFBEE8"="#FFBEE8", "#D1FF73"="#D1FF73", "#A80000"="#A80000",
                                     "#FF7F7F"="#FF7F7F", "#E6E600"="#E6E600", "#00A9E6"="#00A9E6", "#001DE6"="#001DE6",
                                     "#AAFF00"="#AAFF00", "#E64C00"="#E64C00")) + 
        theme(legend.position = "none")
    dev.off()
}


