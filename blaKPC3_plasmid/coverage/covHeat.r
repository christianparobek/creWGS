## make heatmap of coverage data across blaKPC3 plasmid
## started 9 May 2016


##############################
###### USEFUL FUNCTIONS ######
##############################

library(ggplot2)
library(reshape2)
library(stringr)


slide.func <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}
  # sliding window transform


##############################
######## READ IN DATA ########
##############################

files <- list.files(path="./genomecov/")

setwd("genomecov")

sites <- read.table(files[1], header=FALSE, sep="\t")[,2] 
  # read in site / position names

df <- do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,3]))
  # read in just the data column of each file

##############################
######## APPLY WINDOW ########
##############################

data <- df
window <- 100
step = 100

new_df <- apply(df, 2, slide.func, window, step)


##############################
########## TRANSFORM #########
##############################

mat <- replace(new_df, new_df < 25, 0)
mat2 <- as.data.frame(replace(mat, mat >= 25, 1))
  # define cutoff for coloring

mat3 <- cbind(1:nrow(mat2), mat2)
  # add indices

names(mat3) <- append("index", str_extract(files, "......"))
  # add names to the mat3 object
  # remove the file extensions before appending

#########
### Pretty order
#########
keeps <- c("index", "CRE156", "CRE152", "CRE150", "CRE008", "CRE006", "CRE011", "CRE007", "CRE009", "CRE005", "CRE010",
           "CRE012", "CRE014", "CRE013", "CRE114", "CRE001", "CRE002", "CRE003", "CRE015", "CRE004", "CRE149", "CRE119", "CRE155")
mat4 <- mat3[keeps]
names(mat4) <- c("index", "Kp156", "Kp152", "Kp150", "Kp008", "Kp006", "Kp011", "Kp007", "Ec009", "Ec005", "Ec010",
                 "Ec012", "Ec014", "Ec013", "Ec114", "Ec001", "Ec002", "Ec003", "Ec015", "Ec004", "Ec149", "Ec119", "Ec155")
  # rename according to species

#########
### Paper-consistent order
#########
keeps <- c("index", "CRE156", "CRE152", "CRE150", "CRE149", "CRE114", "CRE015", "CRE014", "CRE013", "CRE012", "CRE011", "CRE010", "CRE009",
           "CRE008", "CRE007", "CRE006", "CRE005", "CRE004", "CRE003", "CRE002", "CRE001")
mat4 <- mat3[keeps]
names(mat4) <- c("index", "Kp156", "Kp152", "Kp150", "Ec149", "Ec114", "Ec015", "Ec014", "Ec013", "Ec012", "Kp011", "Ec010", "Ec009",
                 "Kp008", "Kp007", "Kp006", "Ec005", "Ec004", "Ec003", "Ec002", "Ec001")
  # rename according to species

mat4_m <- melt(mat4, id.var = "index")
  # melt it for ggplot plotting

##############################
############ PLOT ############
##############################

tiff("../heatmap.tiff", compress = "lzw", units = "in", height = 9, width = 6, res = 300)

ggplot(mat4_m, aes(index, variable, group = variable)) +
  geom_blank() +
  geom_rect(aes(xmin=345, xmax=445, ymin = 0.5, ymax = 20.5), fill = "gray70", inherit.aes = FALSE) +
  geom_rect(aes(xmin=364, xmax=373, ymin = 0.5, ymax = 20.5), fill = "gray40", inherit.aes = FALSE) + 
  geom_tile(aes(fill = as.factor(value)), height = 0.8) + 
  labs(x = "Position along blaKPC3 Reference (kb)", y = "Clinical Isolate") + 
  theme(axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.margin = unit(c(2, 2, 1, 1), "line")) + 
  scale_fill_manual(breaks = c(0, 1), values = c("#a6cee3", "#1f78b4")) + 
  scale_x_continuous(breaks = c(1, 100, 200, 300, 400, 473), 
                     labels = c(1, 10, 20, 30, 40, 47.3), 
                     limits = c(1, 473), expand = c(0, 0)) +
  geom_segment(aes(x = 473, y = 1, xend = 500, yend = 8))

dev.off()


##############################
###### PERCENT COVERAGE ######
##############################

apply(mat4, 2, sum)/nrow(mat4)

  