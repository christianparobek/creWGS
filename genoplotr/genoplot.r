## Gene map of Tn2-Tn4401
## CRE project



###########################
####### IMPORT LIBS #######
###########################

library(genoPlotR)
library(stringr)
library(grid)



###########################
####### DEFINE FUNCS ######
###########################

getter <- function(samp_name){
  
  ## subset just to cre01
  s <- data[data$sample_name == samp_name,]
  
  sample_df <- data.frame(name = c("blaTEM1b", "tnpR_tn2", "tnpA_tn2_left", "tnpR_tn4401", "tnpA_tn4401", "istA_tn4401", "istB_tn4401", "blaKPC3", "tnpA_tn4401_2", "tnpA_tn2_right"),
                          start = c(s$blaTEM1b_begin,
                                    s$tnpR_Tn2_begin,
                                    s$tnpA_Tn2_left_begin,
                                    s$tnpR_Tn4401b_begin, 
                                    s$tnpA_Tn4401b_begin, 
                                    s$istA_Tn4401b_begin, 
                                    s$istB_Tn4401b_begin, 
                                    s$blaKPC3_begin, 
                                    s$tnpA_2_tn4401b_begin,
                                    s$tnpA_Tn2_right_begin),
                          end = c(s$blaTEM1b_end,
                                  s$tnpR_Tn2_end,
                                  s$tnpA_Tn2_left_end,
                                  s$tnpR_Tn4401b_end, 
                                  s$tnpA_Tn4401b_end, 
                                  s$istA_Tn4401b_end, 
                                  s$istB_Tn4401b_end, 
                                  s$blaKPC3_end, 
                                  s$tnpA_2_tn4401b_end,
                                  s$tnpA_Tn2_right_end),
                          strand = c(-1,-1,1,-1,1,1,1,1,-1,1), lwd = 0.5,
                          col = c("#ca0020", "#92c5de", "#92c5de","#92c5de", "#92c5de", "#92c5de", "#92c5de", "#ca0020", "#92c5de", "#92c5de"))
  
  sample_seg <- dna_seg(sample_df)
  
  return(sample_seg)
}
  # This gets start and end coordinates for arrows

comp_getter <- function(df, upper_samp, lower_samp){
  
  u_s <- with(df[names == upper_samp,], c(utlb, utle+1, utrb))
  # get the upper starts
  u_e <- with(df[names == upper_samp,], c(utle, utrb-1, utre))
  # get the upper ends
  l_s <- with(df[names == lower_samp,], c(utlb, utle+1, utrb))
  # get the lower starts
  l_e <- with(df[names == lower_samp,], c(utle, utrb-1, utre))
  # get the lower ends
  
  out_df <- data.frame(start1 = u_s, end1 = u_e, start2 = l_s, end2 = l_e, col = c("gray", "#f4a582", "gray"))
  out_comp <- comparison(out_df)
  return(out_comp)
}
  # This gets coordinates for track comparisons


###########################
######## READ DATA ########
###########################

data <- read.table("final_positions_7_11_fix.txt", header = TRUE, sep = ",")


##########################################
######## EXTRACT THE SEGMENT DATA ########
##########################################

# define the tn2 segment
tn2 <- data.frame(name = c("blaTEMb", "tnpR", "tnpA"), 
                  start = c(150, 1191, 1912), 
                  end = c(1006, 1748, 4917), 
                  strand = c(-1, -1, 1), lwd = 0.5,
                  col = c("#ca0020", "#92c5de", "#92c5de"))
tn2_seg <- dna_seg(tn2)

## Get the easy ones
cre001_seg <- getter("CRE001"); cre002_seg <- getter("CRE002")
cre003_seg <- getter("CRE003"); cre004_seg <- getter("CRE004")
cre005_seg <- getter("CRE005"); cre007_seg <- getter("CRE007")
cre008_seg <- getter("CRE008"); cre009_seg <- getter("CRE009")
cre010_seg <- getter("CRE010"); cre011_seg <- getter("CRE011")
cre012_seg <- getter("CRE012"); cre013_seg <- getter("CRE013")
cre014_seg <- getter("CRE014"); cre015_seg <- getter("CRE015")
cre114_seg <- getter("CRE114")

## Get the difficult ones
## for CRE006 and CRE149
cre006_df <- data.frame(name = c("tnpR_tn4401", "tnpA_tn4401", "istA_tn4401", "istB_tn4401", "blaKPC3", "tnpA_tn4401_2"),
                        start = c(3111,4936,8072,9094,10260,11514), end = c(4826,7965,9097,9873,11141,12833),
                        strand = c(-1,1,1,1,1,-1), lwd = 0.5, col = c("#92c5de", "#92c5de", "#92c5de", "#92c5de", "#ca0020", "#92c5de"))
cre006_seg <- dna_seg(cre006_df)

## for CRE149
cre149_df <- data.frame(name = c("tnpR_tn2", "tnpA_tn2_left", "tnpR_tn4401", "tnpA_tn4401", "istA_tn4401", "istB_tn4401", "blaKPC3", "tnpA_tn4401_2", "tnpA_tn2_right"),
                        start = c(1251,1789,2740,4565,7701,8723,9889,11143,12824), end = c(1621,2685,4451,7590,8722,9498,10766,12458,14924),
                        strand = c(-1,1,-1,1,1,1,1,-1,1), lwd = 0.5, col = c("#92c5de", "#92c5de", "#92c5de", "#92c5de", "#92c5de", "#92c5de", "#ca0020", "#92c5de", "#92c5de"))
cre149_seg <- dna_seg(cre149_df)
# indels were causing funkiness when I tried doing by hand, so just used CRE07 coords, but adjusted tnpR_tn2 down to the correct size for CRE149

## for CRE152, and CRE156
cre152_df <- data.frame(name = c("tnpA_tn2_left", "tnpR_tn4401", "tnpA_tn4401", "istA_tn4401", "istB_tn4401", "blaKPC3", "tnpA_tn4401_2"),
                        start = c(1789,2740,4565,7701,8723,9889,11143), end = c(2685,4451,7590,8722,9498,10766,12458),
                        strand = c(1,-1,1,1,1,1,-1), lwd = 0.5, col = c("#92c5de", "#92c5de", "#92c5de", "#92c5de", "#92c5de", "#ca0020", "#92c5de"))
cre152_seg <- dna_seg(cre152_df)

## for CRE150
cre150_seg <- cre152_seg[-1,]
cre150_seg[cre150_seg$name == "blaKPC3",]$col <- "#fee090"

###########################################
######## GET COMPARISON TRACK DATA ########
###########################################

## comparison boundary extractor
names <- data$sample_name
utlb <- data$blaTEM1_begin # upper_tn2_left_begin
utle <- data$tnpA_Tn2_left_end # upper_tn2_left_end
utrb <- data$tnpA_Tn2_right_begin # upper_tn2_right_begin
utre <- data$tnpA_Tn2_right_end # upper_tn2_right_end
comps <- data.frame(names, utlb, utle, utrb, utre)


tn2_cre001 <- data.frame(start1 = c(150, 2811), end1 = c(2810, 4948), start2 = c(30235, 42904), end2 = c(32895, 45011), col = "gray")
tn2_cre001_comp <- comparison(tn2_cre001)

cre001_cre002 <- comp_getter(comps, "CRE001", "CRE002")
cre002_cre003 <- comp_getter(comps, "CRE002", "CRE003")
cre003_cre004 <- comp_getter(comps, "CRE003", "CRE004")
cre004_cre005 <- comp_getter(comps, "CRE004", "CRE005")
## do cre005_cre006 by hand (below)
## do cre006_cre007 by hand (below)
cre005_cre008 <- comp_getter(comps, "CRE005", "CRE008")
cre008_cre009 <- comp_getter(comps, "CRE008", "CRE009")
cre009_cre010 <- comp_getter(comps, "CRE009", "CRE010")
cre010_cre010 <- comp_getter(comps, "CRE010", "CRE010") # using 10 instead of 11 because of 11's artifactual dup
cre010_cre012 <- comp_getter(comps, "CRE010", "CRE012") # using 10 instead of 11 because of 11's artifactual dup
cre012_cre013 <- comp_getter(comps, "CRE012", "CRE013")
cre013_cre014 <- comp_getter(comps, "CRE013", "CRE014")
cre014_cre015 <- comp_getter(comps, "CRE014", "CRE015")
cre015_cre114 <- comp_getter(comps, "CRE015", "CRE114")
cre114_cre149 <- comp_getter(comps, "CRE114", "CRE149")


## cre005_cre006 comparison track
cre005_cre006_df <- data.frame(start1 = cre005_seg[cre005_seg$name == "tnpA_tn2_left",]$end, 
                               end1 = cre005_seg[cre005_seg$name == "tnpA_tn2_right",]$start,
                               start2 = cre006_seg$start[1], end2 = rev(cre006_seg$end)[1], col = "#f4a582")
cre005_cre006 <- comparison(cre005_cre006_df)

## cre006_cre007 comparison track
## but actually cre006_cre005 because of artifactual dup in 7
cre006_cre005_df <- data.frame(start1 = cre006_seg$start[1], end1 = rev(cre006_seg$end)[1], 
                               start2 = cre005_seg[cre005_seg$name == "tnpA_tn2_left",]$end, 
                               end2 = cre005_seg[cre005_seg$name == "tnpA_tn2_right",]$start, col = "#f4a582")
cre006_cre005 <- comparison(cre006_cre005_df)

## Fix the 114_149 comparison
cre114_cre149 <- comp_getter(comps, "CRE114", "CRE149")
cre114_cre149$end2[1] <- cre149_seg[cre149_seg$name == "tnpR_tn4401",]$start - 1
cre114_cre149$start2[2] <- cre149_seg[cre149_seg$name == "tnpR_tn4401",]$start
cre114_cre149$end2[2] <- cre149_seg[cre149_seg$name == "tnpA_tn2_right",]$start - 1
cre114_cre149$start2[3] <- cre149_seg[cre149_seg$name == "tnpA_tn2_right",]$start
cre114_cre149$end2[3] <- cre149_seg[cre149_seg$name == "tnpA_tn2_right",]$end

## Make 149_150 comparison
cre149_cre150_df <- data.frame(start1 = cre149_seg[cre149_seg$name == "tnpA_tn2_left",]$end, 
                               end1 = cre149_seg[cre149_seg$name == "tnpA_tn2_right",]$start,
                               start2 = cre152_seg$end[1], end2 = rev(cre152_seg$end)[1], col = "#f4a582")
cre149_cre150 <- comparison(cre149_cre150_df)

## Make 150_152 comparison
cre150_cre152_df <- data.frame(start1 = cre152_seg$end[1], 
                               end1 = rev(cre152_seg$end)[1],
                               start2 = cre152_seg$end[1], end2 = rev(cre152_seg$end)[1], col = "#f4a582")
cre150_cre152 <- comparison(cre150_cre152_df)

###########################
####### PLOTTING IT #######
###########################


dna_segs <- list(tn2_seg, cre001_seg, cre002_seg, cre003_seg, cre004_seg, cre005_seg, cre006_seg, cre005_seg, cre008_seg, cre009_seg, cre010_seg, cre010_seg, cre012_seg, cre013_seg, cre014_seg, cre015_seg, cre114_seg, cre149_seg, cre150_seg, cre152_seg, cre152_seg)
comparisons <- list(tn2_cre001_comp, cre001_cre002, cre002_cre003, cre003_cre004, cre004_cre005, cre005_cre006, cre006_cre005, cre005_cre008, cre008_cre009, cre009_cre010, cre010_cre010, cre010_cre012, cre012_cre013, cre013_cre014, cre014_cre015, cre015_cre114, cre114_cre149, cre149_cre150, cre150_cre152, cre150_cre152)
names(dna_segs) <- c("Tn2", "Ec01", "Ec02", "Ec03", "Ec04", "Ec05", "Kp06", "Kp07", "Kp08", "Ec09", "Ec10", "Kp11", "Ec12", "Ec13", "Ec14", "Ec15", "Ec114", "Ec149", "Kp150", "Kp152", "Kp156") 


svg("genoplot.svg", height = 8, width = 5)

plot_gene_map(dna_segs = dna_segs, comparisons = comparisons, offsets = c(4500,0,0,0,0,0,2800,0,0,0,0,0,0,0,0,0,0,1500,2800,1800,1800), plot_new = FALSE)
grid.text("blaTEM-1", x = 0.18, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpR", x = 0.24, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpR", x = 0.33, y = 0.924, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("tnpA", x = 0.45, y = 0.924, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("istA", x = 0.57, y = 0.924, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("istB", x = 0.622, y = 0.924, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("blaKPC-3", x = 0.69, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpA", x = 0.78, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpA", x = 0.87, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))

dev.off()


## Notes
# cre06 and cre150 are same length wrt tn2 / 4401
# cre152 and 156 are pretty much identical

# cre08 may or may not have some indels in the 5' end of tn2
# looking at IGV pileup of cre08 vs cre01, seems like the "deletion" in the 5' end of tn2
# has pretty similar coverage to the tn4401 element, consistent with the same tn2-tn4401 element
# clearly doesn't have the same entire plasmid
# cre008 also has an intact blaTEM, which `cre008_seg` misrepresents, so just use cre005 for this
# hold on though, cre008 clearly has snps relative to the others, so 

# so just use `cre007_seg` as a stand-in
# cre007 and cre011 are identical to the others from beginning of blaTEM through end of tn2 element...
# with the exception of that repeated sequence in between blaKPC and tnpA_tn4401_2...
# which I determined to be a de novo assembly artifact (see lab notebook).
# So use cre005 and cre010 as standins for these.



##############################
####### MAKING PLAMIDS #######
##############################
library(ggplot2)

data <- read.table("pblaKPC-3_coords", header = TRUE, sep = ",")

data$width <- data$end - data$start # add a width column

data$index <- nrow(data):1

p <- ggplot(data, aes(xmin = start, xmax = end, ymin = index - 0.01, ymax = index + 0.01)) + labs(x = "Plasmid Length") +
  geom_rect(fill = "black") + theme_bw() + scale_y_continuous(breaks = length(data$sample):1, labels = data$sample) +
  geom_rect(aes(xmin=3066, xmax=13058, ymin = index - 0.3, ymax = index + 0.3), fill = "#f4a582") +
  geom_rect(aes(xmin=3110, xmax=4818, ymin = index - 0.2, ymax = index + 0.2), fill = "#92c5de", color = "black", size = 0.3) +
  geom_rect(aes(xmin=4936, xmax=7957, ymin = index - 0.2, ymax = index + 0.2), fill = "#92c5de", color = "black", size = 0.3) +
  geom_rect(aes(xmin=8071, xmax=9088, ymin = index - 0.2, ymax = index + 0.2), fill = "#92c5de", color = "black", size = 0.3) +
  geom_rect(aes(xmin=9094, xmax=9861, ymin = index - 0.2, ymax = index + 0.2), fill = "#92c5de", color = "black", size = 0.3) +
  geom_rect(aes(xmin=10256, xmax=11129, ymin = index - 0.2, ymax = index + 0.2), fill = "#ca0020", color = "black", size = 0.3) +
  geom_rect(aes(xmin=11390, xmax=12698, ymin = index - 0.2, ymax = index + 0.2), fill = "#92c5de", color = "black", size = 0.3) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title.x = element_text(size = rel(0.8)))
## Coordinates are based of Ec04.
## Harcoded coordinates here and in input data file are approximate.

##############################
####### COMBINED PLOT PLAMIDS #######
##############################

svg("plasmid_map.svg", width = 7, height = 7)

grid.newpage()

print(p, vp = viewport(x=0, width=0.4, height=0.95,
                       just="left", name="left"))

pushViewport(viewport(x=1, y = 0.52, width=0.6, height=0.921,
                      just= "right", name="right"))

plot_gene_map(dna_segs = dna_segs, comparisons = comparisons, offsets = c(4500,0,0,0,0,0,2800,0,0,0,0,0,0,0,0,0,0,1500,2800,1800,1800), plot_new = FALSE)
grid.text("blaTEM-1", x = 0.20, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpR", x = 0.25, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpR", x = 0.34, y = 0.9175, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("tnpA", x = 0.46, y = 0.9175, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("istA", x = 0.575, y = 0.9175, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("istB", x = 0.625, y = 0.9175, rot = 0, just = c("left", "bottom"), gp=gpar(fontsize=7, col="black"))
grid.text("blaKPC-3", x = 0.71, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpA", x = 0.78, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))
grid.text("tnpA", x = 0.87, y = 0.94, rot = 45, just = c("left", "bottom"), gp=gpar(fontsize=8, col="black"))

upViewport(1)

grid.text("A", x = 0.02, y = 0.97, gp=gpar(fontsize=20, col="black"))
grid.text("B", x = 0.39, y = 0.97, gp=gpar(fontsize=20, col="black"))


dev.off()

