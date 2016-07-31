#Residence Time Chemostat Experiment
#Single Cell Heterogeneity Correlation Analysis
#Jaylen C. Beatty, July 18th, 2016, Indiana University, Lennon Lab REU
baccount =  read.csv(file.choose(), header = TRUE)
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

plot(logRT, logcell, xlab= "Ln Residence Time (hours)", ylab= "Ln Total Number of Cells (cells/mL)", main = "Population Sizes decrease with Residence Time",  cex= .8, col= "blue")
abline(lm(logcell ~ logRT))
reglogtotalnum= lm(formula= logcount ~ logRT)
summary(reglogtotalnum)#y = - .3629x + 3.904
summary(reglogtotalnum)$r.squared
text(x= ('4.00'), y= 9, paste("r-squared = ", round(summary(reglogtotalnum)$r.squared, 3)))
text(x= ('4.00'), y= 8, labels ="y= -.3629x + 3.904")
cor.test(logcell, logRT, method= "pearson")# p < .05, r=  -.8345473





plot(logRT, logcount, xlab= "Ln Residence Time (hours)", ylab= "Ln Total Number of Cells Per Transect", main = "Population Sizes decrease with Residence Time",  cex= .8, col= "blue")
abline(lm(logcount ~ logRT))
reglogtotalnum= lm(formula= logcount ~ logRT)
summary(reglogtotalnum)#y = - .3629x + 3.904
summary(reglogtotalnum)$r.squared #y= -.47242x + 4.59627
text(x= ('4.00'), y= 5, paste("r-squared = ", round(summary(reglogtotalnum)$r.squared, 3)))
text(x= ('4.00'), y= 4, labels= "p= 2.2e-16")
text(x= ('4.00'), y= 4.5, labels ="y= -.3629x + 3.904")
#The paste function denotes the text to be some variable from some operation.
#In this case, it is the r^2 value.
cor.test(logcount, logRT, method= "kendall")# p < .05, tau=  -.58322
#Suggests that there is a statistically significant correlation between residence time
#and total cell quantity, and that this relationship is negatively correlated
residlogtotalnum= resid(reglogtotalnum)
summary(residlogtotalnum)
plot(logRT, residlogtotalnum, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))
#The regression from Trial 1 showed a strong, linear relationship. The Trial 2 data seems
#to have complicated this relationship. If the patterns observed in Trial 1 held in Trial 2,
#then chemostat 9's community would have possessed the highest cell count, as opposed to chemostat 6's.
#Re-sampling of chemostat 9 may be necessary. It may also be necessary to standardize cells per mL
#because each stained sample was filtered down in different volumes. 

#R.T Vs. % Activity
shapiro.test(X..Active)# p < .05, data is parametric, use pearson's product moment to test correlation
plot(logRT, X..Active, xlab= "Ln Residence Time (hours)", ylab= "% of Active Cells Per Transect", main = "% Activity Declines with Residence Time",  cex= .8, col= "red")
abline(lm(X..Active ~ logRT))
regactive = lm(formula = (X..Active ~ logRT))
summary(regactive)#y = -3.15x + 54.3021
text(x= ('4.00'), y= 80, paste("r-squared = ", round(summary(regactive)$r.squared, 3)))
text(x= ('4.00'), y=10, labels = "y = -3.15x + 54.3021 ")
cor.test(X..Active, logRT, method= "pearson")
#p < .05, suggesting that there is a correlation between % active cells and Residence Time.
# r = - .4231026, suggesting a fairly weak but nonetheless present negative correlation.
#This supports the initial hypothesis in regards to the presence of correlation and the directionality
#However, the weakness of the correlation is somewhat surprising.
residactive= resid(regactive)
summary(residactive)
plot(logRT, residactive, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))
#The relationship has remained consistent between trial 1 and 2. Based off the scatter of the residuals
#a linear correlation could be considered appropriate. It appears that respiration may be slightly higher
#in lower R.T. systems, but there's a wide range of variation in a given community.

#R.T. vs. Mortality
shapiro.test(X..Dead)# p > .05, data is non-parametric, use kendall's.
plot(logRT, X..Dead, xlab= "Ln Residence Time(hours)", ylab= "% of Dead Cells Per Transect", main= "Mortality Is Not Correlated With Residence Time", cex = .8, col= "green")
warnings()
abline(lm(X..Dead ~ logRT))
regdead= lm(formula= (X..Dead ~ logRT))
summary(regdead)#y= -1.2339x + 52.7133
text(x= ('3.90'), y= 10, paste("r-squared = ", round(summary(regdead)$r.squared, 3)))
text(x= ('0.67'), y=10, labels = "y= -1.2339x + 52.7133")
cor.test(X..Dead, logRT, method= "kendall")
# p > .05, suggesting that there is no correlation between mortality and residence time. 
regdead= lm(formula= (X..Dead ~ logRT))
summary(regdead)
residdead= resid(regdead)
summary(residdead)
plot(logRT, residdead, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))
#Inputting the data from Trial 2 has yielded a lack of statistical significance between mortality and R.T.
#This implies that mortality across the R.T. gradient is fairly consistent when all other factors are controlled
#This could simply be indicative of the ability for the wastewater community to survive in OECD media;
#the high or low resource conditions observed in the experiment don't alter their survival significantly.

#R.T. vs. Dormancy
shapiro.test(X..Dormant)# p < .05, data is parametric.
plot(logRT, X..Dormant, xlab= "Ln Residence Time(hours)", ylab= "% of Inactive Cells per Transect", main= "Dormancy increases with Residence Time ", cex = .8)
regdorm= lm(formula= (X..Dormant ~ logRT))
summary(regdorm)#y= 3.15x + 45.6979
abline(lm(X..Dormant ~ logRT))
text(x= ('4.00'), y= 90, paste("r-squared = ", round(summary(regdorm)$r.squared, 3)))
text(x= ('4.00'), y=14, labels = "y= 3.15x + 45.679")
cor.test(X..Dormant, logRT, method= "pearson")
#p < .05, suggesting that a correlation exists between Dormancy Levels and Residence Time.
#r = .4231026, suggesting a weak positive correlation between the two variables
#This supports the initial hypothesis in all regards except the overall strength of the correlation.
#It's unsurprising that this relationship exists considering that it is the inverse of the activity level relationship.
regdorm= lm(formula= (X..Dormant ~ logRT))
residdorm= resid(regdorm)
summary(residdorm)
plot(logRT, residdorm, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))

#Overall it appears that the predictions drawn surrounding the single cell heterogeneity variables 
#were supported. The higher R.T. chemostats had fewer cells, fewer deaths, fewer active members and more dormant members
#and vice versa for the lower R.T. chemostats. 
#The fact that the strength of the relationship between R.T., activity levels, mortality and dormancy is quite interesting

#It must be stated that the levels of activity, dormancy and cell death based upon SYTOX and CTC staining pattersn
#Could be a direct of the sampling. The samples were removed from aeration, resource input, and washout when the stained solution is
#prepared and stored. This return to an infinite residence time could induce dormancy mechanisms in member of the lower R.T. chemostats
#thereby causing the analysis to underestimate the number of active cells and overestimate the number of dormant cells in regards to the original
#chemostat community. Additionally, the dyes sit in solution while the samples are stored both on the filter and in the liquid proper.
#There is likely residual CTC and SYTOX that could be staining cells that have reactivated as a result of the sample preparation,
#and cells that have died inbetween the original SYTOX incubation period and the microscopy analysis. This latter factor is especially notable
#because many of the cells, especially those form the higher R.T. samples are co-stained with SYTOX and CTC, suggesting they were active at one point and died later.

#The mortality factor's variation could be indicative of how the high stress on both ends of the gradient
#can induce cell death. Alternatively, perhaps the percentages were higher in the faster chemostats because of the sheer
#number of cells, or perhaps because r-strategists who grow and die quickly are more available in the faster chemostats

#The activity level's and dormancy level's seemingly low variation may be more of a result of the loss of the high flow rate
#which may have prompted some members of the community to go dormant. It's also possible that the results observed
#in the lower R.T. chemostats are influenced by the prevalence of biofilm production that was unique to the lower R.T. chemostats
#I'd expect most members of the films themselves to be predominantly active, but perhaps that's not the case for the
#planktonic cells who were present in the stained samples. It's also entirely plausible that the original predictions
#overestimaed how much of a liability inactivity in fast systems would have actually been. 

#Based on the way that mortality, activity and dormancy correlated with the trial 1 chemostats
#It's possible that when acquring the "middle value" chemostats from trial 2 the correlations will break down
#ultimately suggesting that R.T. isn't that influential in altering these particular variables.

logestimate= lm(DAPI.Count ~ logRT, data= baccount)
coef(logestimate)
plot(Residence.Time..hrs., DAPI.Count, xlab= "Log Residence Time (hours)", ylab= "Total Number of Cells Per Transect", main = "Population Sizes decrease with Residence Time",  cex= .8, col= "blue")
curve((-14.55014*log(x)+76.21212), add=TRUE, col = 2)

#With regards to the inclusion of the trial 2 data, the relationship observed between % activity and % dormancy
#remains consistent. There is a weak correlation between residence time and these factors, but 
#there is a wide range of variance; a notable percentage of the populations in the high R.T. chemostats remained active
#while a notable percentage in the low R.T. chemostats remained dormant. 

#It must be noted that in regards to CTC counts, many cells, more so in the higher R.T. chemostats
#appeared co-stained with SYTOX because the two stains fluoresce at similar wavelengths. Using the "Green Fire Blue" Look-up Table in FIJI helped to parse the two stains apart
#in trial 2, but many were legitimately stained with both based on the image modifications under that lookup table
#This suggests that many of the cells who remained active or reactivated themselves died off, and these deaths appeared
#more prevalent in the higher R.T. systems, though accurate counts were not collected.

#The relationship between total # of cells and residence time remains the same, but the model 
#has become less predicative of the relationship between this variable and residence time, likely due to the strange
#results derived from trial 2. The comparisons between very high and very low R.T. environments still hold
#firm in the context of cell number comparisons. The middle values may be more complex. 
#It is also possible that the variation in the patterns observed between trial 1 and 2 are due
#to variation within the source communities, which implies that residence time will alter community sizes
#differently depending on the fine-scale composition of the original community and the environmental conditions

#Mortality has been shown to not be correlated with residence time based on the inclusion of Trial 2
#data. The prior prediction of high stress on both ends of gradient could also serve to explain this lack
#of significant correlation. If both conditions force cells to undergo different sorts of adaptations
#to survive, then mortality may remain consistent across the gradient. It's also possibly that the temporal
#and concentration-based dynamics associated with residence time as an environmental variable is 
#less influential than the biotic variables associated with the cells in the community.