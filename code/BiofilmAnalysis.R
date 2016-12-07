#Residence Time Chemostat Experiment
#Biofilm Production Correlation Analysis
#Jaylen C. Beatty, July 19th, 2016, Indiana University, Lennon Lab REU
film=read.csv(file.choose(), header = TRUE)
str(film)
names(film)
summary(film)
attach(film)
shapiro.test(RawAbsorb)#p > .05, data nonparametric. Use Spearman's or transform
#R.T. values are normally distributed.
logrt= log(Restime)#Graph region cannot appropriately capture the wide gradient of R.T. values. Transformation makes this clearer.
logabsorb=log(RawAbsorb)
shapiro.test(logabsorb)
plot(logrt, ModAbsorb, xlab= "Ln Residence Time (hours)", ylab= "Absorbance Readings", main = "Biofilm Production Is Not Correlated with Res. Time",  cex= .8)
text(x= ('3.00'), y=.30, labels='p > .05')
(text(x=('3.00'), y= .25, labels= 'r= -.032133'))
abline(lm(ModAbsorb ~ logrt))
regfilm= lm(formula= logabsorb ~ logrt)
cor.test(ModAbsorb, logrt, method="pearson")# p > .05, suggesting there is no correlation.
summary(regfilm)

#test data with controls subtracted
shapiro.test(ModAbsorb)# p > .05, data nonparametric
plot(logrt, ModAbsorb)
abline(lm(ModAbsorb ~ logrt))
cor.test(ModAbsorb, logrt, method="pearson")
regfilm= lm(formula= ModAbsorb ~ logrt)
summary(regfilm)

recabsorb= (1/ModAbsorb)
plot(logrt, recabsorb)
shapiro.test(recabsorb)
abline(lm(recabsorb ~ logrt))
sqrtabsorb= sqrt(ModAbsorb)
shapiro.test(sqrtabsorb)
sqrtrestime= sqrt(Restime)
shapiro.test(sqrtrestime)

allelofilm= aov(ModAbsorb ~ Restime, data= film)
summary(allelofilm)
TukeyHSD(allelofilm)
boxplot(ModAbsorb ~ Chemostat, data= film, main= "Biofilm Production Comparison" , xlab= "Community Label", ylab= " Absorbance Readings", cex.lab= .7)
