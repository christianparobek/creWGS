## An R script to plot the output of Gubbins
## Started 23 November, 2015
## Modified for CRE 25 April, 2016

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

library(ape)
library(ade4)
library(stringr)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)


################################################
######### DEFINE SOME USEFUL FUNCTIONS #########
################################################

sample.detector <- function(df, index) {
  index <- str_pad(index, 2, pad = "0") # pad a leading zero if necessary
  sample <- (paste("CRE", index, sep = "")) # slap an A on the front to match my Ab names
  row_indices <- str_detect(df[,4], sample) # detect each row with that sample
  return(df[row_indices,]) # return a dataframe with data for that row
} # given a df & index (samplename), return rows matching that name 


################################################
########## READ & PROCESS TREE DATA ############
################################################

ec_tree <- read.tree("../ec/gubbins/ec.final_tree.tre")
kp_tree <- read.tree("../kp/gubbins/kp.final_tree.tre")
  # read in the tree

ec_tree_log <- ec_tree
kp_tree_log <- kp_tree
  # duplicate tree structure

ec_tree_log$edge.length <- log1p(ec_tree$edge.length)
kp_tree_log$edge.length <- log1p(kp_tree$edge.length)
  # log transform the edge lengths

plot(ec_tree_log, type = "phylogram", use.edge.length = TRUE, x.lim = c(-5, 50))
plot(kp_tree_log, type = "phylogram", use.edge.length = TRUE, x.lim = c(-5, 50))
  # one way of plotting

ec_newick <- write.tree(ec_tree_log) 
kp_newick <- write.tree(kp_tree_log) 
  # convert to nwk

ec_phylog <- newick2phylog(ec_newick)
kp_phylog <- newick2phylog(kp_newick) 
  # convert to phylog (an ade4 object)


################################################
########## READ & PROCESS EMBL DATA ############
################################################

## Start with recombination_predictions.embl
## In gedit, edit the embl file until left with
## something in this form:
## 3584801  3599125	colour="4"	taxa="CRE149"
## where it's feature start, feature end, color, taxa name

ec_data <- read.table("../ec/gubbins/tmp2", header = FALSE)
kp_data <- read.table("../kp/gubbins/tmp2", header = FALSE)

blanktheme <- theme(line = element_blank(),
                    text = element_blank(), 
                    panel.border = element_blank(), 
                    legend.position="none",
                    plot.margin=unit(c(0.2,0.1,0,0.05), "cm"),
                    panel.background = element_rect(fill = "white"))

ec_l = list() # make an empty list
for (i in 001:156) {
  
  ec_subset <- sample.detector(ec_data, i)
  
  if (as.logical(dim(ec_subset))[1] == TRUE) {
    p <- ggplot(ec_subset, aes(xmin=V1, xmax=V2, ymin=1, ymax=2, color = factor(V3))) + 
      geom_rect(size = 0.2) + blanktheme + scale_x_continuous(limits = c(0, 5000000)) + 
      scale_colour_manual(drop=TRUE, limits = levels(ec_subset$V3), values = c("#d7191c", "#2b83ba"))
  } else {
    p <- ggplot() + geom_blank() + blanktheme + scale_x_continuous(limits = c(0, 5000000))
  }
  ec_l[[i]] <- p
} # for each index, get the data and make a graph


# sample.detector() isn't perfect because when I search for 015 it will also catch 155
# so the quick fix is just change 015 and 014 to something that it didn't recognize
# here I'm using 001


###################################
######## COMBINED PLOTTING ########
###################################

svg("gubbins.svg", width = 8, height = 8)

vp1 <- viewport(x=0.04,y=0.35,width=0.6, height=0.65, just = c("left", "bottom"))
vp2 <- viewport(x=0.04,y=0,width=0.6, height=0.35, just = c("left", "bottom"))

vp3 <- viewport(x=0.45,y=0.41,width=0.4, height=0.57, just = c("left", "bottom"))
vp4 <- viewport(x=0.45,y=0.045,width=0.4, height=0.29, just = c("left", "bottom"))

pushViewport(vp1)
par(new=TRUE, fig=gridFIG())
plot(ec_phylog)
upViewport()

pushViewport(vp2)
par(new=TRUE, fig=gridFIG())
plot(kp_phylog)
axis(1, at=log1p(seq(0,1000, by=100)),line = 3.5, labels = c(0, 100, rep("", 9)), cex.axis = 0.8, mgp=c(3, .3, 0), cex.axis = 0.5)
upViewport()

grid.arrange(ec_l[[001]], ec_l[[001]], ec_l[[149]], ec_l[[013]], ec_l[[114]], ec_l[[012]],
             ec_l[[155]], ec_l[[119]], ec_l[[009]], ec_l[[005]], ec_l[[010]], ec_l[[002]],
             ec_l[[001]], ec_l[[004]], ec_l[[003]],
             ncol = 1, vp = vp3, newpage = FALSE)

grid.arrange(l[[006]], l[[150]], l[[011]], l[[007]], l[[008]], l[[156]], 
             l[[152]], ncol = 1, vp = vp4, newpage = FALSE)
  ## Run gubbins_for_kp.r, make sure it points to the right kp files in the kp/gubbins/ directory
  ## Ideally should be integrated into this script, above, but couldn't make it work quickly

grid.text("A", x = 0.03, y = 0.98, gp=gpar(fontsize=17, col="black"))
grid.text("B", x = 0.03, y = 0.33, gp=gpar(fontsize=17, col="black"))

dev.off()
