## An R script to plot the output of Gubbins
## Started 23 November, 2015
## Modified for CRE 25 April, 2016

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

library(ape)
library(ade4)
#library(genoPlotR)
library(stringr)
library(ggplot2)
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

tree <- read.tree("ec.final_tree.tre")
  # read in the tree

tree_log <- tree
  # duplicate tree structure

tree_log$edge.length <- log1p(tree$edge.length)
  # log transform the edge lengths

plot(tree_log, type = "phylogram", use.edge.length = TRUE, x.lim = c(-5, 50))
  # one way of plotting, but labsls not even

my_newick <- write.tree(tree_log) 
  # convert to nwk

my_phylog <- newick2phylog(my_newick) 
  # convert to phylog (an ade4 object)

svg("Ec_tree.svg", width = 8, height = 7)
plot(my_phylog)
  # plot the tree
axis(1, at=log1p(seq(0,10000, by=10)),line = 3.5, labels = c(0, 10, rep("", 999)), cex.axis = 0.8, mgp=c(3, .3, 0))
  # add a log-scaled x axis, labeling the first two ticks only
dev.off()

################################################
########## READ & PROCESS EMBL DATA ############
################################################

## Start with recombination_predictions.embl
## In gedit, edit the embl file until left with
## something in this form:
## 3584801  3599125	colour="4"	taxa="CRE149"
## where it's feature start, feature end, color, taxa name

data <- read.table("tmp2", header = FALSE)

blanktheme <- theme(line = element_blank(),
                    text = element_blank(), 
                    panel.border = element_blank(), 
                    legend.position="none",
                    plot.margin=unit(c(0.2,0.1,0,0.05), "cm"),
                    panel.background = element_rect(fill = "white"))

l = list() # make an empty list
for (i in 001:155) {
  
  data_subset <- sample.detector(data, i)
  
  if (as.logical(dim(data_subset))[1] == TRUE) {
    p <- ggplot(data_subset, aes(xmin=V1, xmax=V2, ymin=1, ymax=2, color = factor(V3))) + 
      geom_rect(size = 0.2) + blanktheme + scale_x_continuous(limits = c(0, 5000000)) + 
      scale_colour_manual(drop=TRUE, limits = levels(data$V3), values = c("#d7191c", "#2b83ba"))
  } else {
    p <- ggplot() + geom_blank() + blanktheme + scale_x_continuous(limits = c(0, 5000000))
  }

  l[[i]] <- p
} # for each index, get the data and make a graph

# sample.detector() isn't perfect because when I search for 015 it will also catch 155
# so the quick fix is just change 015 and 014 to something that it didn't recognize
# here I'm using 001

svg("Ec_tracks.svg", width = 3, height = 6)

grid.arrange(l[[001]], l[[001]], l[[149]], l[[013]], l[[114]], l[[012]],
             l[[155]], l[[119]], l[[009]], l[[005]], l[[010]], l[[002]],
             l[[001]], l[[004]], l[[003]],
             ncol = 1)

dev.off()




##################################################
## Failed attempts at combining tree and tracks ##
##################################################

library(grid)
library(gridBase)
library(gridExtra)
layout(matrix(2, ncol = 2))
plot(my_phylog)

# First base plot
plot(1:10)


par(mfrow=c(1,2))

plot.new()              ## suggested by @Josh
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <-plotViewport(c(1.8,1,0,1))







library(grid)
library(gridBase)
library(gridExtra)


layout(matrix(c(1,3, 2,3, 4,3), nrow = 3, ncol = 2, byrow = TRUE))

# First base plot
plot(1:10)

# second base plot 
frame()
# Grid regions of current base plot (ie from frame)
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
# Table grob
grob <-  tableGrob(iris[1:2,1:2])  
grid.draw(grob)

popViewport(3)

# third base plot
plot(1:10)

# fourth
frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)  
grid.draw(grob)
popViewport(3)
