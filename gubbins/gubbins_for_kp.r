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
library(gridExtra)

################################################
######### DEFINE SOME USEFUL FUNCTIONS #########
################################################

sample.detector <- function(df, index) {
  index <- str_pad(index, 3, pad = "0") # pad a leading zero if necessary
  sample <- (paste("CRE", index, sep = "")) # slap a CRE on the front to match my CRE names
  row_indices <- str_detect(df[,4], sample) # detect each row with that sample
  return(df[row_indices,]) # return a dataframe with data for that row
} # given a df & index (samplename), return rows matching that name 


################################################
########## READ & PROCESS TREE DATA ############
################################################

tree <- read.tree("kp.final_tree.tre")
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

svg("Kp_tree.svg", width = 8, height = 2.8)
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
for (i in 001:156) {
  
  data_subset <- sample.detector(data, i)
  
  if (as.logical(dim(data_subset))[1] == TRUE) {
    p <- ggplot(data_subset, aes(xmin=V1, xmax=V2, ymin=1, ymax=2, color = factor(V3))) + 
      geom_rect(size = 0.1) + blanktheme + scale_x_continuous(limits = c(0, 5250000)) + 
      scale_colour_manual(drop=TRUE, limits = levels(data$V3), values = c("#d7191c", "#2b83ba"))
  } else {
    p <- ggplot() + geom_blank() + blanktheme + scale_x_continuous(limits = c(0, 5250000))
  }

  l[[i]] <- p
} # for each index, get the data and make a graph

# sample.detector() isn't perfect because when I search for 015 it will also catch 155
# so the quick fix is just change 015 and 014 to something that it didn't recognize
# here I'm using 001

svg("Kp_tracks.svg", width = 3, height = 2.4)

grid.arrange(l[[006]], l[[150]], l[[011]], l[[007]], l[[008]], l[[156]], l[[152]], ncol = 1)

dev.off()
