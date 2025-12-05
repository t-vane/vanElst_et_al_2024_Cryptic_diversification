### Script to generate Supplementary Figure S6 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ape")
library("ggtree")
library("patchwork")
library("phytools")

## Set root and outgroup samples
root <- c("Cheirogaleuscro_106189", "Cheirogaleusmaj_106245", "Cheirogaleusmed_106354")
out <- c("Cheirogaleuscro_106189", "Cheirogaleusmaj_106245", "Cheirogaleusmed_106354","Mirzazaz_DLC2316m","Mirzazaz_DLC319m","Mirzazaz_DLC323f")

## Set renaming map
list_gsub <- read.csv(sep = ";", text="
old; new
Marnholdi; italic('M. arnholdi')                 
Msp1; paste(italic('M.'), ' sp. 1')      
Mmamiratra; italic('M. mamiratra')                  
Mmargotmarshae; italic('M. margotmarshae')                 
Msambiranensis; italic('M. sambiranensis')
Mtavaratra; italic('M. tavaratra')               
Mberthae; italic('M. berthae')
Mrufus; italic('M. rufus')
Mmyoxinus; italic('M. myoxinus')
Mleh_Ambatovy; paste(italic('M. lehilahytsara'),' (Ambatovy)')
Mleh_Ankafobe; paste(italic('M. lehilahytsara'),' (Ankafobe)')
Mleh_Ambavala; paste(italic('M. lehilahytsara'), ' (Ambavala)')
Mleh_Riamalandy; paste(italic('M. lehilahytsara'), ' (Riamalandy)')
Mmittermeieri; italic('M. mittermeieri')
Mtanosi; italic('M. tanosi')
Mboraha; italic('M. boraha')
Msimmonsi_N; paste(italic('M. simmonsi'), ' (north)')
Msimmonsi_C; paste(italic('M. simmonsi'), ' (south)')
Mgerpi; italic('M. gerpi')
Mjollyae; italic('M. jollyae')
Mmarohita; italic('M. marohita')
Mjonahi; italic('M. jonahi')
Mmacarthurii; italic('M. macarthurii')
Mbongolavensis; italic('M. bongolavensis')
Mravelobensis; italic('M. ravelobensis')
Mdanfossi; italic('M. danfossi')
Mganzhorni; italic('M. ganzhorni')
Mmanitatra; italic('M. manitatra')
Mmurinus_S; paste(italic('M. murinus'), ' (south)')
Mmurinus_C; paste(italic('M. murinus'), ' (central)')
Mmurinus_N; paste(italic('M. murinus'), ' (north)')
Mgriseorufus; italic('M.griseorufus')
")

## Loop over five different thresholds of missing data
thresholds <- c(0.05, 0.25, 0.5, 0.75, 0.95)
plots <- list()

for th in thresholds {
    # Read tree
    tree_unrooted <- read.tree(paste0("final.maxmiss", th, ".thinned.speciesAssignments.svdq.tre"))

    # Root and ladderize
    tree <- ladderize(ape::root(tree_unrooted, root))

    # Drop outgroups to improve visibility
    tree <- drop.tip(tree1, tip = out)

    # Change tip labels
    for(x in 1:nrow(list_gsub)) {
        tree$tip.label <- gsub(list_gsub[x,"old"],list_gsub[x,"new"], tree$tip.label)
    }

    # Plot
    p <- ggtree(tree) + 
        xlim_tree(25) + 
        geom_tree(size = 0.5) + 
        theme(legend.position = "none") +
        geom_tiplab(align=TRUE, parse = TRUE, linetype='dashed', size = 5) + 
        geom_text2(aes(subset=(label!=100 & !isTip), label = label, size = 5, hjust = -0.2))

    plots[[s]] <- p
}

pdf("Fig_S6.pdf", width = 18, height = 18)
(plots$0.05 | plots$0.25 | plots$0.5) /
(plots$0.75 | plots$0.95)
dev.off()