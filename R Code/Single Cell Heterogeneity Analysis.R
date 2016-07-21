#Residence Time Chemostat Experiment
#Single Cell Heterogeneity Correlation Analysis
#Jaylen C. Beatty, July 18th, 2016, Indiana University, Lennon Lab REU
baccount =  read.csv(file.choose(), header = TRUE)
str(baccount)
names(baccount)
summary(baccount)
attach(baccount)
shapiro.test(Residence.Time..hrs.)# p < .05, R.T. times are normally distributed
shapiro.test(logRT)# p < .05 
shapiro.test(DAPI.Count)# p < .05, Total cell counts normally distributed
#Use Pearson's Product Moment for correlation.

#R.T. Vs. Total Cell Count
logRT= log(Residence.Time..hrs.)#Graph region cannot appropriately capture the wide gradient of R.T. values. Transformation makes this clearer.
logcount=log(DAPI.Count)
plot(logRT, DAPI.Count, xlab= "Log Residence Time (hours)", ylab= "Total Number of Cells Per Transect", main = "Population Sizes decrease with Residence Time",  cex= .8, col= "blue")
text(x= ('3.00'), y=200, labels='r= -.8408713')
text(x= ('3.00'), y=155, labels ="y= -.3643x + 118.6097")
abline(lm(DAPI.Count ~ logRT))
regtotalnum= lm(formula= DAPI.Count ~ logRT)
#Linear Equation: -.3643x + 118.6097

#Check Residuals to determine appropriateness of linear correlation
residtotalnum=resid(regtotalnum)
summary(residtotalnum)
plot(logRT, residtotalnum, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))
cor.test(DAPI.Count, logRT, method= "pearson")
#Because p < .05, the null hypothesis is refuted, suggesting there is a correlation between 
# Communal cell quantity and residence time. The r value is equal to -.8408713, suggesting a 
# strong negative correlation. This fits the prediction of the original hypothesis.
summary(lm(DAPI.Count ~ Residence.Time..hrs.))
#The linear equation for this correlation is: -.3637x + 118.7109

#R.T Vs. % Activity
shapiro.test(X..Active)# p < .05, data is parametric, use pearson's product moment to test correlation
plot(logRT, X..Active, xlab= "Log Residence Time (hours)", ylab= "% of Active Cells Per Transect", main = "% Active Cells across Residence Time Gradient",  cex= .8, col= "red")
abline(lm(X..Active ~ logRT))
text(x= ('3.00'), y=80, labels='r= -.4231026')
text(x= ('3.00'), y=70, labels ="y= -3.284x + 54.475")
cor.test(X..Active, logRT., method= "pearson")
#p < .05, suggesting that there is a correlation between % active cells and Residence Time.
# r = - .4231026, suggesting a fairly weak but nonetheless present negative correlation.
#This supports the initial hypothesis in regards to the presence of correlation and the directionality
#However, the weakness of the correlation is somewhat surprising.
regactive = lm(formula = (X..Active ~ logRT))
summary(regactive)
residactive= resid(regactive)
summary(residactive)
plot(logRT, residactive, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))

#R.T. vs. Mortality
shapiro.test(X..Dead)# p < .05, data is parametric, use pearson's product moment.
plot(logRT, X..Dead, xlab= "Log Residence Time(hours)", ylab= "% of Dead Cells Per Transect", main= "Mortality Across Residence Time Gradient", cex = .8, col= "green")
warnings()
abline(lm(X..Dead ~ logRT))
text(x= ('3.00'), y=40, labels='r= -.3165323')
text(x= ('3.00'), y=30, labels ="y= -2.1611x + 57.7660")
cor.test(X..Dead, logRT, method= "pearson")
# p < .05, suggesting that there is a correlation between mortality and residence time. 
#r = -.3165323, suggesting a weak negative correlation similar to that of the relationship
#observed between R.T. and activity levels. 
regdead= lm(formula= (X..Dead ~ logRT))
summary(regdead)
residdead= resid(regdead)
summary(residdead)
plot(logRT, residdead, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))


#R.T. vs. Dormancy
shapiro.test(X..Dormant)# p < .05, data is parametric.
plot(logRT, X..Dormant, xlab= "Log Residence Time(hours)", ylab= "% of Inactive Cells per Transect", main= "Dormancy Levels Across Residence Time Gradient", cex = .8)
abline(lm(X..Dormant ~ logRT))
text(x= ('3.00'), y=80, labels='r= .4231026')
text(x= ('3.00'), y=70, labels ="y= 3.284x + 45.525")
cor.test(X..Dormant, logRT, method= "pearson")
#p < .05, suggesting that a correlation exists between Dormancy Levels and Residence Time.
#r = .4231026, suggesting a weak positive correlation between the two variables
#This supports the initial hypothesis in all regards except the overall strength of the correlation.
#It's unsurprising that this relationship exists considering that it is essentially the inverse of the activity level relationship.
regdorm= lm(formula= (X..Dormant ~ logRT))
summary(regdorm)
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
curve((-21.57748*log(x)+126.95155), add=TRUE, col = 2)

#Logarithmic Transformation of Total Cell # Graph
plot(logRT, logcount, xlab= "Log Residence Time (hours)", ylab= "Log Total Number of Cells Per Transect", main = "Population Sizes decrease with Residence Time",  cex= .8, col= "blue")
abline(lm(logcount ~ logRT))
reglogtotalnum= lm(formula= logcount ~ logRT)
summary(reglogtotalnum)#y= -.47242x + 4.59627
text(x= ('4.00'), y= 5, labels='r= -.922685')
text(x= ('4.00'), y=4, labels ="y= -.47242x + 4.59627")
cor.test(logcount, logRT, method= "pearson")# p < .05, r -.922685
plot(logRT, residlogtotalnum, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot.", abline(0,0))
residlogtotalnum= resid(reglogtotalnum)
summary(residlogtotalnum)
plot(logRT, residlogtotalnum, xlab= "Residence Time", ylab= "Residual", main = "Resid. Plot for R.T.", abline(0,0))
