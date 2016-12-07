
baccount =  read.csv("../data/microscopy/final_cell_count_post.csv", header = TRUE)
str(baccount)
names(baccount)
summary(baccount)
attach(baccount)
shapiro.test(Residence.Time..hrs.)# p < .05, R.T. times are normally distributed
logRT= log(Residence.Time..hrs.)#Graph region cannot appropriately capture the wide gradient of R.T. values. Transformation makes this clearer.
shapiro.test(logRT)# p < .05, data is parametric
#Use Pearson's Product Moment for correlation.

#R.T. Vs. Total Cell Count. Logarithmic Transformation 
logcount=log(DAPI.Count)
logcell=log(Cell_ml)
shapiro.test(logcell)# p < .05 Data is parametric.
shapiro.test(logcount)#p > .05. Data is nonparametric. 
#The data was transformed to linearize the graph.
#As opposed to transformation, Use Kendall's to test correlation.



png("../figures/Density_RT.png", height = 1200, width = 1200, res = 192)
plot(logRT, logcell, 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n",
     cex = 1.2, pch = 16, xlim = log10(c(0.1, 1000)), ylim = c(4,10))
axis(side = 1, 
     labels = c(.1, 1, 10, 100, 1000),
     at = log10(c(.1, 1, 10, 100, 1000)), cex.lab = 1.5, lwd.ticks = 2)
axis(side = 2,
     las = 1, cex.lab = 1.5, lwd.ticks = 2)
axis(side = 3, at = log10(c(.1, 1, 10, 100, 1000)), labels = F, lwd.ticks = 2)
axis(side = 4, labels = F, lwd.ticks = 2)
box(lwd = 2)
abline(lm(logcell ~ logRT), lwd = 2)
mtext("Residence Time (Hours)", line = 3, side = 1, cex = 1.5)
mtext(expression("Bacterial Density" ~ (log[10] ~ cells ~ ml^{-1}) ), 
      line = 2.25, side = 2, cex = 1.2)
dev.off()

png("../figures/Activity_RT.png", height = 1200, width = 1200, res = 192)
plot(logRT, X..Active, 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n",
     cex = 1.2, pch = 16, ylim = c(0,100), xlim = log10(c(.1, 1000)))
axis(side = 1, 
     labels = c(.1, 1, 10, 100, 1000),
     at = log10(c(.1, 1, 10, 100, 1000)), cex.lab = 1.5, lwd.ticks = 2)
axis(side = 2,
     las = 1, cex.lab = 1.5, lwd.ticks = 2)
axis(side = 3, at = log10(c(.1, 1, 10, 100, 1000)), labels = F, lwd.ticks = 2)
axis(side = 4, labels = F, lwd.ticks = 2)
box(lwd = 2)
abline(lm(X..Active ~ logRT), lwd = 2)
mtext("Residence Time (Hours)", line = 3, side = 1, cex = 1.5)
mtext(expression("Percent Activity"), 
      line = 2.25, side = 2, cex = 1.5)
dev.off()