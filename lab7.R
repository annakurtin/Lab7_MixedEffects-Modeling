## ----setup, include=FALSE-----------------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE)
r <- getOption("repos")
r["CRAN"] <- "https://ftp.osuosl.org/pub/cran/"
options(repos = r)


## ----load packages, include=FALSE---------------------------------------------------
#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("ggplot2", "sandwich", "lme4", "lmtest", "merTools", "ResourceSelection", "GGally", "Hmisc", "plotrix", "pander", "lattice", "jstats", "sjstats", "jtools")
#run function to install packages
ipak(packages)


## -----------------------------------------------------------------------------------
wolfkde2 <- read.csv("Data/wolfkde.csv", header=TRUE, sep = ",", na.strings="NA", dec=".")
wolfkde3 <-na.omit(wolfkde2)
wolfkde3$usedFactor <-as.factor(wolfkde3$usedFactor)

head(wolfkde3)
table(wolfkde2$pack, wolfkde2$used)


## -----------------------------------------------------------------------------------
top.env <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3)
pander(summary(top.env))
ggcoef(top.env, exclude_intercept = TRUE)


## ----message=FALSE------------------------------------------------------------------
top.env.bv <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3, subset=pack== "Bow Valley")
pander(summary(top.env.bv))

## but subset by packs
top.env.rd <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3, subset=pack== "Red Deer")
pander(summary(top.env.rd))


## -----------------------------------------------------------------------------------
top.env.mixed <- glmer(used~Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn+(1|pack), data=wolfkde3,family=binomial(link="logit"))
summary(top.env.mixed)


## -----------------------------------------------------------------------------------
head(wolfkde3)
wolfkde3$Elevation2_sc <-scale(wolfkde3$Elevation2)
hist(wolfkde3$Elevation2)
hist(wolfkde3$Elevation2_sc)
plot(wolfkde3$Elevation2, wolfkde3$Elevation2_sc)
summary(wolfkde3$Elevation2)
summary(wolfkde3$Elevation2_sc)

wolfkde3$DistFromHighHumanAccess2_sc <- scale(wolfkde3$DistFromHighHumanAccess2)
plot(wolfkde3$DistFromHighHumanAccess2, wolfkde3$DistFromHighHumanAccess2_sc, )


## -----------------------------------------------------------------------------------
top.env2 <- glm(used ~ Elevation2_sc + DistFromHighHumanAccess2_sc + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3)
pander(summary(top.env2))


## -----------------------------------------------------------------------------------
sd(wolfkde3$Elevation2)


## -----------------------------------------------------------------------------------
exp(top.env2$coefficients[2])


## -----------------------------------------------------------------------------------
top.env.mixed2 <- glmer(used~Elevation2_sc + DistFromHighHumanAccess2_sc + openConif+modConif+closedConif+mixed+herb+shrub+water+burn+(1|pack), data=wolfkde3,family=binomial(link="logit"))
summary(top.env.mixed2)


## -----------------------------------------------------------------------------------
AIC(top.env2, top.env.mixed2, top.env.rd, top.env.bv)


## -----------------------------------------------------------------------------------
# Bring in data
elk <- read.table("Data/lab7_elk_migrant.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
head(elk)
elk$elkuidF <- as.factor(elk$elkuid)

# get to know our data
table(elk$elkuid, elk$year)
table(elk$elkuid, elk$used)


## -----------------------------------------------------------------------------------
ggplot(elk, aes(x=utmx, y = utmy, color = elkuidF)) + geom_point()


## -----------------------------------------------------------------------------------
elk.used <- subset(elk, elk$used == 1)
ggplot(elk.used, aes(x=utmx, y = utmy, color = elkuidF)) + geom_point()
write.csv(elk.used, "Data/lab8_elk_migrant_used.csv")  ## we might want to use this later, like in Lab 8


## -----------------------------------------------------------------------------------
hist(elk$ctotrisk)
hist(log(elk$ctotrisk))

summary(elk[,30:31]) ## So 1579 NA's predation risk values, and 11 NA's for total herbaceous vegetation
length(elk$ctotrisk)


## -----------------------------------------------------------------------------------
elk2 <- elk[complete.cases(elk[30:31]), ]
summary(elk2)
length(elk2$ctotrisk)


## -----------------------------------------------------------------------------------
elk2$ctotrisk[elk2$ctotrisk>1]=1
table(elk2$elkuid, elk2$year)
table(elk2$elkuid, elk2$used)

# Compute sample size for each elk
n = tapply(elk$idn, elk$elkuid,length)
n


## -----------------------------------------------------------------------------------
wolf = tapply(na.omit(elk2$ctotrisk), elk2$elkuid[which((elk2$ctotrisk!="NA")==TRUE)],mean)
wolf
hist(wolf)


## -----------------------------------------------------------------------------------
forage = tapply(na.omit(elk2$totalherb), elk$elkuid[which((elk2$totalherb!="NA")==TRUE)],mean)
forage
hist(forage)


## -----------------------------------------------------------------------------------
elk2$totalherb_sc <- scale(elk2$totalherb)
elk2$ctotrisk_sc <- scale(elk2$ctotrisk)
elk2$ctotrisk2_sc <- scale(elk2$ctotrisk2)
elk2$riskforage_sc <- scale(elk2$riskforage)
elk2$for2_sc <- scale(elk2$for2)
elk2$risk2_sc <- scale(elk2$risk2)

##### AGain, just to double check what scale is doing
plot(elk2$ctotrisk_sc, elk2$ctotrisk)


## -----------------------------------------------------------------------------------
# Fitting best model(s)
forage = glm(used~totalherb, data=elk2,family=binomial(link="logit"))
risk = glm(used~ctotrisk, data=elk2,family=binomial(link="logit"))
forANDrisk = glm(used~totalherb+ctotrisk, data=elk2,family=binomial(link="logit"))
forrisk = glm(used~totalherb+ctotrisk+ctotrisk*totalherb, data=elk2,family=binomial(link="logit"))

AIC(forage, risk, forANDrisk, forrisk)


## -----------------------------------------------------------------------------------
summary(forrisk)
ggcoef(forrisk, exclude_intercept = TRUE)


## -----------------------------------------------------------------------------------
forrisk_sc = glm(used~totalherb_sc+ctotrisk_sc+ctotrisk_sc*totalherb_sc, data=elk2,family=binomial(link="logit"))
summary(forrisk_sc)
ggcoef(forrisk_sc, exclude_intercept = TRUE)


## -----------------------------------------------------------------------------------
hist(elk2$ctotrisk)


## -----------------------------------------------------------------------------------
mean(elk2$ctotrisk)


## -----------------------------------------------------------------------------------
# Calculate some summary statistics for forage
hist(elk$totalherb)
hist(log(elk$totalherb))
quantile(elk$totalherb,na.rm=TRUE)
mean(elk$totalherb,na.rm=TRUE)

herb.lo = 5
herb.med = 15
herb.hi = 50


## -----------------------------------------------------------------------------------
predrisk = seq(0,1,0.01)
pred.lo = predict(forrisk,type="response", newdata = data.frame(totalherb=herb.lo, ctotrisk=predrisk))
pred.med = predict(forrisk,type="response", newdata = data.frame(totalherb=herb.med, ctotrisk=predrisk))
pred.hi = predict(forrisk,type="response", newdata = data.frame(totalherb=herb.hi, ctotrisk=predrisk))

# Make plot
plot(elk2$ctotrisk,elk2$used, xlab="Risk", ylab="Pr(Use)")
lines(predrisk,pred.lo, lty=1)
lines(predrisk,pred.med, lty=2)
lines(predrisk,pred.hi, lty=3)
legend(x=0.7,y=0.95, legend=c("Observed","Low Forage","Medium Forage","High Forage"), pch=c(1,rep(-1,3)),lty=c(-1,1:3),bty="n")


## -----------------------------------------------------------------------------------
elk2$usedF <- as.factor(elk2$used)
ggplot(elk2, aes(x=ctotrisk, y = used)) + stat_smooth(method="glm", method.args = list(family="binomial"))
ggplot(elk2, aes(x=totalherb, y = used)) + stat_smooth(method="glm", method.args = list(family="binomial"))
ggplot(elk2, aes(x=riskforage, y = used)) + geom_rug() + stat_smooth(method="glm", method.args = list(family="binomial"))


## -----------------------------------------------------------------------------------
elk2$forage.cat  <- as.factor(as.numeric(cut2(elk2$totalherb, g=3)))
elk2$forage.cat  <- cut2(elk2$totalherb, g=3)
ggplot(elk2, aes(x=ctotrisk, y = used, fill = forage.cat)) + stat_smooth(method="glm", method.args = list(family="binomial"))


## -----------------------------------------------------------------------------------
elk2$risk.cat  <- cut2(elk2$ctotrisk, g=3)
ggplot(elk2, aes(x=totalherb, y = used, fill = risk.cat)) + stat_smooth(method="glm", method.args = list(family="binomial"))


## -----------------------------------------------------------------------------------
mep(forrisk_sc)


## -----------------------------------------------------------------------------------
forrisk2 <-vcovHC(forrisk, elk$elkuid, type = "const", sandwhich = TRUE)
forrisk2

coeftest(forrisk, forrisk2)


## -----------------------------------------------------------------------------------
forriskFI = glm(used~totalherb+ctotrisk+ctotrisk*totalherb+elkuidF, data=elk,family=binomial(link="logit"))
summary(forriskFI)
ggcoef(forriskFI, exclude_intercept = TRUE)


## -----------------------------------------------------------------------------------
table(elk2$elkuid, elk2$used)


## -----------------------------------------------------------------------------------
AIC(forriskFI, forage, risk, forANDrisk, forrisk)


## -----------------------------------------------------------------------------------
elkids = sort(unique(elk2$elkuid))
modelnames = paste("elk",elkids)
models = list()

for (i in 1:length(elkids)){
	models[[i]]=glm(used~totalherb_sc+ctotrisk_sc+ctotrisk_sc*totalherb_sc, data=elk2,subset = elkuid==elkids[i], family=binomial(link="logit"))

}

names(models)=modelnames
#lapply(models,summary) #### Note I supressed this because it just spits out 17 models, 1 for each elk

# This creates a dataframe with the beta coefficients for each model/elk
coefs = as.data.frame(lapply(models, coef))
coefs

#Calculate means for each of the coefficients
mean(as.numeric(coefs[1,]))
mean(as.numeric(coefs[2,]))
mean(as.numeric(coefs[3,]))
mean(as.numeric(coefs[4,]))


## -----------------------------------------------------------------------------------
par(mfrow=c(2,2))
hist(as.numeric(coefs[1,]), main="intercept",breaks=10)
hist(as.numeric(coefs[2,]), main="Forage",breaks=10)
hist(as.numeric(coefs[3,]), main ="Risk",breaks=10)
hist(as.numeric(coefs[4,]), main="Forage*Risk",breaks=10)


## -----------------------------------------------------------------------------------
fr.ri = glmer(used~totalherb_sc+ctotrisk_sc+ctotrisk_sc*totalherb_sc+(1|elkuid), data=elk2,family=binomial(link="logit"), verbose=FALSE)
summary(fr.ri)


## -----------------------------------------------------------------------------------
confint(fr.ri, method = "Wald")


## -----------------------------------------------------------------------------------
ggcoef(fr.ri)


## -----------------------------------------------------------------------------------
fr.ri2 = glmer(used~totalherb_sc+ctotrisk_sc+ctotrisk_sc*totalherb_sc+(1|elkuid), data=elk2,family=binomial(link="logit"), verbose=FALSE, nAGQ = 10)
summary(fr.ri2)


## -----------------------------------------------------------------------------------
# Useful functions
fixef(fr.ri) 
ranef(fr.ri)


## -----------------------------------------------------------------------------------
coef(fr.ri) 


## -----------------------------------------------------------------------------------
B0 <- cbind(as.numeric(coefs[1,]),ranef(fr.ri)$elkuid[,1]+fixef(fr.ri)[1])
rownames(B0)=rownames(ranef(fr.ri)$elkuid)
colnames(B0) = c("Two-Step","Random Effects")
str(B0)
B0
par(mfrow=c(1,1))
plot(B0[,1], B0[, 2])
abline(lm(B0[,1] ~ B0[, 2]))


## -----------------------------------------------------------------------------------
# Make a histogram using the package plotrix
multhist(list(B0[,1],B0[,2]), xlab = "Intercept Coefficient",ylab="Frequency (# Elk)", col = c("gray","tan"))
legend("topright", legend = colnames(B0), fill = c("gray","tan"), bty = "n")


## -----------------------------------------------------------------------------------
fr.rc = glmer(used~totalherb_sc+ctotrisk_sc+totalherb_sc*ctotrisk_sc+(ctotrisk_sc|elkuid), data=elk2,family=binomial(link="logit"), verbose=FALSE)

fixef(fr.rc) # This is the fixed effects coefficients
ranef(fr.rc) # These are the random effects, which in this model is just (1|elkuid), so, one coefficient for each individual elk


## -----------------------------------------------------------------------------------
coef(fr.rc)
##### Note here that the coefficient for predation risk is allowed to vary
summary(fr.rc)


## -----------------------------------------------------------------------------------
lattice::dotplot(ranef(fr.rc, condVar = TRUE))


## -----------------------------------------------------------------------------------
B0.rc = cbind(as.numeric(coefs[1,]),coef(fr.rc)$elkuid[,1])
rownames(B0.rc)=rownames(ranef(fr.ri)$elkuid)
colnames(B0.rc) = c("Two-Step","Random Effects")
B.risk = cbind(as.numeric(coefs[3,]),coef(fr.rc)$elkuid[,3])
rownames(B.risk)=rownames(ranef(fr.ri)$elkuid)
colnames(B.risk) = c("Two-Step","Random Effects")

## lets look at the Intercepts
B0.rc


## -----------------------------------------------------------------------------------
B.risk


## -----------------------------------------------------------------------------------
plot(B.risk[,1], B.risk[,2])
abline(lm(B.risk[,1] ~ B.risk[, 2]))


## -----------------------------------------------------------------------------------
# Make histogram of betas
par(mfrow = c(1,2))
multhist(list(B0.rc[,1],B0.rc[,2]), xlab = "Intercept Coefficient",ylab="Frequency (# Elk)", col = c("gray","tan"), ylim=c(0,10))
legend(2.4,10, legend = colnames(B0), fill = c("gray","tan"), bty = "n")
multhist(list(B.risk[,1],B.risk[,2]), xlab = "Risk Coefficient",ylab="Frequency (# Elk)", col = c("gray","tan"),ylim = c(0,10))
legend(2.4,10, legend = colnames(B0), fill = c("gray","tan"), bty = "n")


## -----------------------------------------------------------------------------------
AIC(forriskFI, forage, risk, forANDrisk, forrisk, fr.ri, fr.rc)


## -----------------------------------------------------------------------------------
anova(fr.rc, forrisk, test = "Chisq") ## put GLMM first


## -----------------------------------------------------------------------------------
anova(fr.rc, fr.ri, test = "Chisq") ## put GLMM first


## -----------------------------------------------------------------------------------
jtools::summ(fr.rc, confint = TRUE, digits =2)


## ---- warning = FALSE---------------------------------------------------------------
performance::icc(fr.rc)


## ---- warning = FALSE---------------------------------------------------------------
performance::icc(fr.ri)


## -----------------------------------------------------------------------------------
fr.rc2 = glmer(used~totalherb+ctotrisk+totalherb*ctotrisk+(ctotrisk|elkuid), data=elk2,family=binomial(link="logit"), verbose=FALSE)
summary(fr.rc2)


## -----------------------------------------------------------------------------------
elk.c = elk2[which((elk2$ctotrisk!="NA"&elk2$totalherb!="NA")==TRUE),]


## -----------------------------------------------------------------------------------
par(mfrow = c(1,1))
plot(elk.c$ctotrisk, elk.c$used,xlab="Risk",ylab="Pr(Use)")
elkids = sort(unique(elk.c$elkuid))
ltypes = rep(1:6,each=3)
lwide = rep(seq(1,3,1),each = 6)
colors = rep(c("red","black","seagreen", "violet", "tan", "orange"),3)

# Begin loop
for (i in elkids){
  # To plot predictions you need to create new data for just that elk
  dat = as.data.frame(model.matrix(terms(fr.rc2),elk.c)[elk.c$elkuid==i,])
  dat$totalherb = mean(dat$totalherb) # Use the mean forage for an elk
  dat[,colnames(dat)=="totalherb:ctotrisk"] = dat$totalherb*dat$ctotrisk # Recalculate the interaction term with the mean forage
  dat$pred.prob = plogis(as.matrix(dat)%*%t(as.matrix(coef(fr.rc2)$elkuid[which(elkids==i),]))) # Use matrix algebra to get prediction based on coefficients for each individual elk
  ord = order(dat$ctotrisk) # store the order we want to plot in
  # Draw a line for you prediction
  lines(dat$ctotrisk[ord], dat$pred.prob[ord], lty=ltypes[which(elkids==i)],lwd = lwide[which(elkids==i)],col=colors[which(elkids==i)])
}

legend("right", legend = c("Observed", paste("Elk ",elkids)), pch=c(1,rep(-1,length(elkids))),lty = c(-1,ltypes[1:length(elkids)]),lwd = c(-1,lwide[1:length(elkids)]),col = c(1,colors[1:length(elkids)]), bty = "n")



## -----------------------------------------------------------------------------------
ggplot(elk2, aes(x=ctotrisk, y = used, colour = elkuidF)) + stat_smooth(method="glm", method.args = list(family="binomial"))


## -----------------------------------------------------------------------------------
elk2$naive.pred <-predict(forrisk, type = "response")
elk2$fr.rc.pred <- predict(fr.rc, type = "response")
hist(elk2$fr.rc.pred)


## -----------------------------------------------------------------------------------
elk2$fr.rc.pred2 <- predict(fr.rc, re.form = NA, type = "response")
summary(elk2$fr.rc.pred2)


## -----------------------------------------------------------------------------------
elk2$fr.rc.pred3 <- predict(fr.rc, re.form = ~(1|elkuid) , type = "response")
summary(elk2$fr.rc.pred3)
hist(elk2$fr.rc.pred3)


## -----------------------------------------------------------------------------------
plot(elk2$fr.rc.pred2, elk2$fr.rc.pred)


## -----------------------------------------------------------------------------------
ggpairs(elk2[46:49])


## -----------------------------------------------------------------------------------
ggplot(elk2, aes(utmx, utmy, col = naive.pred)) + geom_point(size=2) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')


## -----------------------------------------------------------------------------------
ggplot(elk2, aes(utmx, utmy, col = fr.rc.pred)) + geom_point(size=2) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')


## -----------------------------------------------------------------------------------
ggplot(elk2, aes(utmx, utmy, col = fr.rc.pred2)) + geom_point(size=2) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')


## -----------------------------------------------------------------------------------
ggplot(elk2, aes(utmx, utmy, col = fr.rc.pred3)) + geom_point(size=2) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')


## -----------------------------------------------------------------------------------
boot.fr.rc <- bootMer(fr.rc, FUN = function(x) as.numeric(logLik(x), nsim = 100))
boot.fr.rc$mle


## -----------------------------------------------------------------------------------
preds <- predictInterval(fr.rc, newdata = elk2, which = "fixed", n.sims = 99, stat = "median", type = "probability")
hist(preds)
head(preds)


## -----------------------------------------------------------------------------------
preds <- predictInterval(fr.rc, newdata = elk2, which = "random", n.sims = 99, stat = "median", type = "probability")
head(preds)


## ----eval=FALSE, include=FALSE------------------------------------------------------
## knitr::purl(input = "README.Rmd", output = "lab7.R", documentation = 1)

