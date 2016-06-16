## Read in data
ec <- read.table("ec/coverage/coverage_ec.txt", header = FALSE)
kp <- read.table("kp/coverage/coverage_kp.txt", header = FALSE)

## Order data acendingly
ec <- ec[with(ec, order(V2)), ]
## Add a third column with numbers in it
ec$V6 <- 1:length(ec$V1)

## Order data acendingly
kp <- kp[with(kp, order(V2)), ]
## Add a third column with numbers in it
kp$V6 <- 1:length(kp$V1)



##Plot the coverage graph
tiff(file="cov.tiff", width = 4, height = 5, units = "in", res = 600, compression = "lzw")

## set up plot parameters
par(mfrow = c(2,1), mar = c(2,5,3,1))

plot(ec$V2 ~ ec$V6, axes=FALSE, xlab="", ylab="", ylim=c(0,1), col="black", pch=20, type = "n")
points(ec$V2 ~ ec$V6, col="#80cdc1", pch=20)
points(ec$V3 ~ ec$V6, col="#018571", pch=20)
points(ec$V4 ~ ec$V6, col="#dfc27d", pch=20)
points(ec$V5 ~ ec$V6, col="#a6611a", pch=20)
axis(1, at=1:length(ec$V1), labels=ec$V1, cex.axis=.7, las=3, cex = 0.5)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), las = 2, cex.axis = 0.7)
legend(9,0.5, legend=c("5x Coverage", "10x Coverage", "25x Coverage", "50x Coverage"), col=c("#80cdc1", "#018571", "#dfc27d", "#a6611a"), pch=20, cex = 0.6)

par(mar = c(3,5,2,1))
plot(kp$V2 ~ kp$V6, axes=FALSE, xlab="", ylab="", ylim=c(0,1), col="black", pch=20, type = "n")
points(kp$V2 ~ kp$V6, col="#80cdc1", pch=20)
points(kp$V3 ~ kp$V6, col="#018571", pch=20)
points(kp$V4 ~ kp$V6, col="#dfc27d", pch=20)
points(kp$V5 ~ kp$V6, col="#a6611a", pch=20)
axis(1, at=1:length(kp$V1), labels=kp$V1, cex.axis=.7, las=3, cex = 0.5)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), las = 2, cex.axis = 0.7)

mtext("Fraction of Genome at Fold Coverage", side = 2, outer = TRUE, line = -2, cex = 0.8)

dev.off()
