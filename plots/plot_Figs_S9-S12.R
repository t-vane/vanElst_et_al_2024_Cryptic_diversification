### Script to generate Supplementary Figures S9-S12 of van Elst et al. (2025), Nature Ecology & Evolution (https://doi.org/10.1038/s41559-024-02547-w)

library(ggplot2)
library(gridExtra)
library(scales)

## Read posterior distributions
x1 <- read.table("mcmc.final.run1.txt",header=TRUE,sep="\t")
x2 <- read.table("mcmc.final.run2.txt",header=TRUE,sep="\t")
x3 <- read.table("mcmc.final.run3.txt",header=TRUE,sep="\t")
x4 <- read.table("mcmc.final.run4.txt",header=TRUE,sep="\t")

## Generate parameter labels
param_labs <- c("Sample", 
    expression(italic("\u03B8")[italic("M. arnholdi")]),expression(italic("\u03B8")[italic("M. berthae")]),expression(italic("\u03B8")[italic("M. bongolavensis")]),expression(italic("\u03B8")[italic("M. boraha")]),
    expression(italic("\u03B8")[italic("M. danfossi")]),expression(italic("\u03B8")[italic("M. ganzhorni")]),expression(italic("\u03B8")[italic("M. gerpi")]),expression(italic("\u03B8")[italic("M. griseorufus")]),
    expression(italic("\u03B8")[italic("M. jollyae")]),expression(italic("\u03B8")[italic("M. jonahi")]),expression(italic("\u03B8")[paste(italic("M. lehilahytsara"), " (north)")]),expression(italic("\u03B8")[paste(italic("M. lehilahytsara"), " (south)")]),
    expression(italic("\u03B8")[italic("M. macarthurii")]),expression(italic("\u03B8")[italic("M. mamiratra")]),expression(italic("\u03B8")[italic("M. manitatra")]),expression(italic("\u03B8")[italic("M. margotmarshae")]),
    expression(italic("\u03B8")[italic("M. marohita")]),expression(italic("\u03B8")[italic("M. mittermeieri")]),expression(italic("\u03B8")[paste(italic("M. murinus"), " (north)")]),expression(italic("\u03B8")[paste(italic("M. murinus"), " (central)")]),
    expression(italic("\u03B8")[italic("M. myoxinus")]),expression(italic("\u03B8")[italic("M. ravelobensis")]),expression(italic("\u03B8")[italic("M. rufus")]),expression(italic("\u03B8")[italic("M. sambiranensis")]),
    expression(italic("\u03B8")[italic("M. simmonsi")]),expression(italic("\u03B8")[paste(italic("M."), " sp. 1")]),expression(italic("\u03B8")[italic("M. tanosi")]),expression(italic("\u03B8")[italic("M. tavaratra")]),
    expression(italic("\u03B8")[italic("Mirza zaza")]),expression(italic("\u03B8")[Root~all]),expression(italic("\u03B8")[paste("Root ",italic("Microcebus"))]),
    
    expression(italic("\u03B8")[Node~a]),expression(italic("\u03B8")[Node~b]),expression(italic("\u03B8")[Node~c]),expression(italic("\u03B8")[Node~d]),expression(italic("\u03B8")[Node~e]),expression(italic("\u03B8")[Node~f]),
    expression(italic("\u03B8")[Node~g]),expression(italic("\u03B8")[Node~h]),expression(italic("\u03B8")[Node~i]),expression(italic("\u03B8")[Node~j]),expression(italic("\u03B8")[Node~k]),expression(italic("\u03B8")[Node~l]),
    expression(italic("\u03B8")[Node~m]),expression(italic("\u03B8")[Node~n]),expression(italic("\u03B8")[Node~o]),expression(italic("\u03B8")[Node~p]),expression(italic("\u03B8")[Node~q]),expression(italic("\u03B8")[Node~r]),
    expression(italic("\u03B8")[Node~s]),expression(italic("\u03B8")[Node~t]),expression(italic("\u03B8")[Node~u]),expression(italic("\u03B8")[Node~v]),expression(italic("\u03B8")[Node~w]),expression(italic("\u03B8")[Node~x]),
    expression(italic("\u03B8")[Node~y]),expression(italic("\u03B8")[Node~z]),expression(italic("\u03C4")[Root~all]),expression(italic("\u03C4")[paste("Root ",italic("Microcebus"))]),expression(italic("\u03C4")[Node~a]),expression(italic("\u03C4")[Node~b]),
    expression(italic("\u03C4")[Node~c]),expression(italic("\u03C4")[Node~d]),expression(italic("\u03C4")[Node~e]),expression(italic("\u03C4")[Node~f]),expression(italic("\u03C4")[Node~g]),expression(italic("\u03C4")[Node~h]),
    expression(italic("\u03C4")[Node~i]),expression(italic("\u03C4")[Node~j]),expression(italic("\u03C4")[Node~k]),expression(italic("\u03C4")[Node~l]),expression(italic("\u03C4")[Node~m]),expression(italic("\u03C4")[Node~n]),
    expression(italic("\u03C4")[Node~o]),expression(italic("\u03C4")[Node~p]),expression(italic("\u03C4")[Node~q]),expression(italic("\u03C4")[Node~r]),expression(italic("\u03C4")[Node~s]),expression(italic("\u03C4")[Node~t]),
    expression(italic("\u03C4")[Node~u]),expression(italic("\u03C4")[Node~v]),expression(italic("\u03C4")[Node~w]),expression(italic("\u03C4")[Node~x]),expression(italic("\u03C4")[Node~y]),expression(italic("\u03C4")[Node~z]), 

    "Log likelihood")
                  
## Plot likelihood distribution
# Generate and bind data frames
d1 <- data.frame(value = x1[,"lnL"], variable = "chain 1")
d2 <- data.frame(value = x2[,"lnL"], variable = "chain 2")
d3 <- data.frame(value = x3[,"lnL"], variable = "chain 3")
d4 <- data.frame(value = x4[,"lnL"], variable = "chain 4")
dat <- rbind(d1,d2,d3,d4)
dat$variable <- factor(dat$variable, c("chain 1", "chain 2", "chain 3", "chain 4"))

# Plot
p <- ggplot(dat, aes(x = variable, y = value)) + 
    geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
    scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac")) + 
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + 
    labs(title=param_labs[which(colnames(x1) == "lnL")]) +
    scale_y_continuous(label=comma)

png("Fig_S9.png",res=300,height=5*300,width=8*300)
grid.arrange(p87,ncol=1,nrow=1)
dev.off()

## Plot tau posterior distribution
# Initialize variables
tau_cols <- grep("^tau_", colnames(x1), value = TRUE)
mnh12 <- c()
mnh34 <- c()
mnhCounter <- 1
plots <- list()

# Loop over tau columns
for (i in tau_cols){
    # Generate and bind data frames
    d1 <- data.frame(value = x1[,i], variable = "chain 1")
    d2 <- data.frame(value = x2[,i], variable = "chain 2")
    d3 <- data.frame(value = x3[,i], variable = "chain 3")
    d4 <- data.frame(value = x4[,i], variable = "chain 4")
    dat <- rbind(d1,d2,d3,d4)
    dat$variable <- factor(dat$variable, c("chain 1", "chain 2", "chain 3", "chain 4"))

    # Plot violins
    p <- ggplot(dat, aes(x = variable, y = value)) + 
        eom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
        scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac")) + 
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + 
        labs(title=param_labs[which(colnames(x1) == i)])
    plots[[length(plots)+1]] <- p
  
    # Generate median node height lists
    mnh12[mnhCounter] <- median(c(x1[,i],x2[,i]))
    mnh34[mnhCounter] <- median(c(x3[,i],x4[,i]))
    mnhCounter <- mnhCounter + 1
}
# Plot node heights
mnhData <- as.data.frame(cbind(mnh12,mnh34))
p <- ggplot(data=mnhData,aes(x=mnh12,y=mnh34,alpha=0.7)) + 
    geom_point(size=3) + 
    geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + 
    theme(legend.position="none") + 
    labs(title="Median Node Heights",x="chains 1 + 2",y="chains 3 + 4") + 
    xlim(0.0,0.005) + 
    ylim(0.0,0.005)
plots[[length(plots)+1]] <- p

# Plot all
png("Fig_S10.png",res=300,height=10*300,width=8*300)
do.call(grid.arrange, c(plots, ncol=4))
dev.off()


## Plot theta posterior distribution
# Initialize variables
theta_cols <- grep("^theta_", colnames(x1), value = TRUE)
mnh12 <- c()
mnh34 <- c()
mnhCounter <- 1
plots <- list()

# Loop over theta columns
for (i in tau_cols){
    # Generate and bind data frames
    d1 <- data.frame(value = x1[,i], variable = "chain 1")
    d2 <- data.frame(value = x2[,i], variable = "chain 2")
    d3 <- data.frame(value = x3[,i], variable = "chain 3")
    d4 <- data.frame(value = x4[,i], variable = "chain 4")
    dat <- rbind(d1,d2,d3,d4)
    dat$variable <- factor(dat$variable, c("chain 1", "chain 2", "chain 3", "chain 4"))

    # Plot violins
    p <- ggplot(dat, aes(x = variable, y = value)) + 
        eom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
        scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac")) + 
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + 
        labs(title=param_labs[which(colnames(x1) == i)])
    plots[[length(plots)+1]] <- p
  
    # Generate median node height lists
    mnh12[mnhCounter] <- median(c(x1[,i],x2[,i]))
    mnh34[mnhCounter] <- median(c(x3[,i],x4[,i]))
    mnhCounter <- mnhCounter + 1
}
# Plot node heights
mnhData <- as.data.frame(cbind(mnh12,mnh34))
p <- ggplot(data=mnhData,aes(x=mnh12,y=mnh34,alpha=0.7)) + 
    geom_point(size=3) + 
    geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + 
    theme(legend.position="none") + 
    labs(title="Median Node Heights",x="chains 1 + 2",y="chains 3 + 4") + 
    xlim(0.0,0.005) + 
    ylim(0.0,0.005)
plots[[length(plots)+1]] <- p

# Plot all
png("Fig_S11.png",res=300,height=10*300,width=8*300)
do.call(grid.arrange, c(plots[1:29], ncol=4, nrow=8))
dev.off()

png("Fig_S12.png",res=300,height=10*300,width=8*300)
do.call(grid.arrange, c(plots[30:length(plots)], ncol=4, nrow=8))
dev.off()
