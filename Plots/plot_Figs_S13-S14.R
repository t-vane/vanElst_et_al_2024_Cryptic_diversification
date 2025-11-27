### Script to generate Supplementary Figures S13-S14 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library("ape")
library("deeptime")
library("ggplot2")
library("ggtree")
library("ggstar")
library("phytools")
library("treeio")

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

## Loop over two different versions of posterior age distribution (with and without uncertainty in estimates)
version <- c("nouncertainty", "uncertainty")

for v in version {
    # Read and clean tree
    tree <- read.beast(paste0("BPP_mean_FigTree_", v, ".tre"))
    tree <- drop.tip(tree, tip = c("MlehilahytsaraS"))
    
    # Rename tips
    for(x in 1:nrow(list_gsub)) {
        tree@phylo$tip.label <- gsub(list_gsub[x,"old"],list_gsub[x,"new"], tree@phylo$tip.label)
    }

    # Set node ages and letters
    q <- ggtree(tree)
    node_ages <- round(max(q$data$x) - q$data$x, digits=0)
    node_letters <- c(rep(NA, 30),letters[1:25])

    # Plot
    pdf(paste0("Figs_S13-S14_", v, ".pdf"), height = 10, width = 15.5)
    p <- ggtree(tree) + 
        geom_tree() + 
        geom_tiplab(align=TRUE, parse = TRUE, linetype='dashed', size = 5) + 
        xlim_tree(200) +
        geom_vline(xintercept = c(-2250,-2000,-1750,-1500,-1250,-1000,-750,-500,-250), color="dimgray",linetype = "dashed") + 
        geom_range("height_95_HPD", color='red', size=3, alpha=.4) +  
        geom_text2(aes(subset=(!isTip), label=node_ages), color = "black", hjust = -0.22, size=5) +
        geom_text2(aes(subset=(node > 30), label=node_letters), color = "dimgray", hjust = 2, vjust = -0.5, size = 5) + 
        theme_tree2() + xlab(label = "Thousand years ago (ka)") + 
        theme(axis.text.x = element_text(size = 15),axis.title.x = element_text(size = 15)) +
        scale_x_continuous(breaks=seq(-2250,0,250), labels=abs(seq(-2250,0,250))) 
        # Revert scale
        revts(p)
    dev.off()
}
