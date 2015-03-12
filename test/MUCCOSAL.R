####################################################################
## Analysis of Mississippi EGG and GET (Western Blot) data
## 11/05/14
####################################################################

load("egg_get.rda")

save.image("egg_get.rda")

#setwd("C:\\Users\\Christy\\Dropbox\\Consulting\\EGG Abell")
#setwd("C:\\Users\\cmpink02\\Dropbox\\Consulting\\EGG Abell")

#data <- read.csv("Mississippi_GET.csv", header = TRUE)

data <- read.csv("all_data.csv", header = TRUE)

####################################################################
## Create the FAR values
####################################################################
data$cut.far <- data$cutfreq/data$cutamp
data$cut5far <- data$cut5freq/data$cut5amp
data$mucfar <- data$mucfreq/data$mucamp
data$muc5far <- data$muc5freq/data$muc5amp


####################################################################
## Check to see what happens with permutations of the beta values.
####################################################################
beta1 <- function (x, y, n, sub){
  beta.obs.boot <- numeric(n)
  for(i in 1:n) {
    beta.obs.boot[i] <- coef(lm(sample(y, length(y)) ~ x))[2]
  }
  summary(beta.obs.boot)
  plot(density(beta.obs.boot), lwd=3, col="steelblue",
       main=c("Re-sampling betas. True beta is: ", round(coef(lm(y ~ x))[2], 2)), sub=sub)
}


beta2 <- function (x, y, n, sub) {
  ibm.lm <- lm(y~x)
  ibm.fit <- fitted(ibm.lm)
  ibm.resid <- resid(ibm.lm)
  beta.resid.boot <- numeric(n)
  for(i in 1:n) {
    beta.resid.boot[i] <- coef(lm(
      ibm.fit + ibm.resid[sample(y, length(y))] ~ x))[2]
  }
  plot(density(beta.resid.boot), lwd=3, col="steelblue",
       main=c("Re-sampling beta residuals for amplitude. True beta is: ", round(coef(lm(y ~ x))[2], 2)), sub=sub)
} 


set.seed(1321654)
n<-1000

###### BASELINE and 5-day
# Frequency 
beta1(data$cutfreq, with(data, cut5freq - cutfreq), n, "Cutaneous EGG Frequency")
beta2(data$cutfreq, with(data, cut5freq - cutfreq), n, "Cutaneous EGG Frequency")
beta1(data$mucfreq, with(data, muc5freq - mucfreq), n, "Mucosal EGG Frequency")
beta2(data$mucfreq, with(data, muc5freq - mucfreq), n, "Mucosal EGG Frequency")
# Amplitude
beta1(data$cutamp, with(data, cut5amp - cutamp), n, "Cutaneous EGG Amplitude")
beta2(data$cutamp, with(data, cut5amp - cutamp), n, "Cutaneous EGG Amplitude")
beta1(data$mucamp, with(data, muc5amp - mucamp), n, "Mucosal EGG Amplitude")
beta2(data$mucamp, with(data, muc5amp - mucamp), n, "Mucosal EGG Amplitude")
# FAR
beta1(data$cut.far, with(data, cut5far - cut.far), n, "Cutaneous EGG FAR")
beta2(data$cut.far, with(data, cut5far - cut.far), n, "Cutaneous EGG FAR")
beta1(data$mucfar, with(data, muc5far - mucfar), n, "Mucosal EGG FAR")
beta2(data$mucfar, with(data, muc5far - mucfar), n, "Mucosal EGG FAR")


####################################################################
## Check to see if there is a correlation between the cutaneous data
## and the mucosal data.
####################################################################

###### BASELINE
cor.test(data$cutfreq, data$mucfreq, use="complete")
cor.test(data$cutamp, data$mucamp, use="complete")
cor.test(data$cut.far, data$mucfar, use="complete")

###### 5-Day follow-up
cor.test(data$cut5freq, data$muc5freq, use="complete")
cor.test(data$cut5amp, data$muc5amp, use="complete")
cor.test(data$cut5far, data$muc5far, use="complete")


####################################################################
## Check to see if there is a correlation between the cutaneous data
## and the mucosal data baseline and change from baseline values.
####################################################################

###### Cutaneous
cor.test(data$cutfreq, with(data, cut5freq - cutfreq), use="complete")
cor.test(data$cutamp, with(data, cut5amp - cutamp), use="complete")
cor.test(data$cut.far, with(data, cut5far - cut.far), use="complete")

###### Mucosal
cor.test(data$mucfreq, with(data, muc5freq - mucfreq), use="complete")
cor.test(data$mucamp, with(data, muc5amp - mucamp), use="complete")
cor.test(data$mucfar, with(data, muc5far - mucfar), use="complete")


####################################################################
## Create plots of the baseline versus change in values for mucosal
## and cutaneous for each of the three values: freq, amp, FAR.
####################################################################
library(RColorBrewer)
basechgplt <- function (x, y, xlab){
  cols<-brewer.pal(10, "BrBG")
  lm.fit <- lm(y ~ x)
  lo.fit <- loess(y ~ x)
  xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
  ypred <- predict(lo.fit, data.frame(x = xseq))
  plot(x, y, xlab = xlab, ylab = "Difference from Baseline", col = cols[9])
  abline(lm.fit, col=cols[10], lwd = 2)
  lines(xseq, ypred, col = cols[8])
  legend("topright", c("Linear", "Loess"), lty = 1, col = cols[c(10,8)], bty = "n")
}


###### Cutaneous
basechgplt(data$cutfreq, with(data, cut5freq - cutfreq), "Baseline EGG Frequency")
basechgplt(data$cutamp, with(data, cut5amp - cutamp), "Baseline EGG Amplitude")
basechgplt(data$cut.far, with(data, cut5far - cut.far), "Baseline EGG FAR")

###### Mucosal
basechgplt(data$mucfreq, with(data, muc5freq - mucfreq), "Baseline Mucosal Frequency")
basechgplt(data$mucamp, with(data, muc5amp - mucamp), "Baseline Mucosal Amplitude")
basechgplt(data$mucfar, with(data, muc5far - mucfar), "Baseline Mucosal FAR")





####################################################################
## First Analysis
## Look at EGG difference pre/post in freq, amp, and ratio (far)
####################################################################
## NOTE - Need to combine plots into 1x3 figure ...

## Frequency
x <- data$mucfreq
y <- with(data, muc5freq - mucfreq)

cor(x, y, use = "complete")
##  -0.6661728
lm.fit <- lm(y ~ x)
summary(lm.fit)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  4.25143    0.32854   12.94   <2e-16 ***
## x           -0.78081    0.05789  -13.49   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 2.022 on 228 degrees of freedom
##   (343 observations deleted due to missingness)
## Multiple R-squared:  0.4438,  Adjusted R-squared:  0.4413
## F-statistic: 181.9 on 1 and 228 DF,  p-value: < 2.2e-16

lo.fit <- loess(y ~ x)
xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
ypred <- predict(lo.fit, data.frame(x = xseq))

## Plot of frequency changes
## pdf("Freq.pdf")
## png("Freq.png", width = 7, height = 7, units = "in", res = 600)
#####par(mfrow = c(1, 3))
#####par(mar = c(5, 4, 2, 2))

plot(x, y, xlab = "Baseline Mucosal Frequency", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")

## dev.off()


## Comparing baseline and follow-up densities
## pdf("EGG Frequency Density.pdf")
## idx.rm <- which(is.na(data$cutfreq))
## plot(density(data$cutfreq[-idx.rm]), xlab = "EGG Frequency",
##      main = "Density of EGG Frequency")
## idx.rm <- which(is.na(data$Temp.Cut.Freq))
## lines(density(data$Temp.Cut.Freq[-idx.rm]), col = 2)
## legend("topright", c("Baseline", "Follow-up"), col = 1:2, lty = 1)
## dev.off()
## t = 1.5681, df = 119, p-value = 0.1195


## y <- with(data, Temp.Cut.Freq - cutfreq)

####################################################################
## Amplitude
####################################################################

x <- data$mucamp
y <- data$muc5Amp - data$mucamp

n<-1000
beta.obs.boot <- numeric(n)
for(i in 1:n) {
  beta.obs.boot[i] <- coef(lm(sample(y, length(y)) ~ x))[2]
}
summary(beta.obs.boot)
plot(density(beta.obs.boot), lwd=3, col="steelblue",
     main=c("Re-sampling betas for amplitude. True beta is: ", round(coef(lm(y ~ x))[2], 2)))
#abline(v=coef(lm(y ~ x))[2], lwd=3, col='gold')

ibm.lm <- lm(y~x)
ibm.fit <- fitted(ibm.lm)
ibm.resid <- resid(ibm.lm)
beta.resid.boot <- numeric(n)
for(i in 1:n) {
  beta.resid.boot[i] <- coef(lm(
    ibm.fit + ibm.resid[sample(y, length(y))] ~ x))[2]
}
plot(density(beta.resid.boot), lwd=3, col="steelblue",
     main=c("Re-sampling beta residuals for amplitude. True beta is: ", round(coef(lm(y ~ x))[2], 2)))

#x1<-matrix(rep(x, 1000), 1000, length(x), byrow=TRUE)


## plot(x, y)
cor(x, y, use = "complete")
## -0.7465946
lm.fit <- lm(y ~ x)
summary(lm.fit)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  0.08596    0.00863   9.961   <2e-16 ***
## x           -0.71642    0.04294 -16.683   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.1063 on 221 degrees of freedom
##   (350 observations deleted due to missingness)
## Multiple R-squared:  0.5574,	Adjusted R-squared:  0.5554
## F-statistic: 278.3 on 1 and 221 DF,  p-value: < 2.2e-16

## abline(lm.fit)
lo.fit <- loess(y ~ x)
xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
ypred <- predict(lo.fit, data.frame(x = xseq))
## lines(xseq, ypred, col = 2)

## Plot of Amplitude changes
## pdf("Amplitude.pdf")
## png("Amplitude.png", width = 7, height = 7, units = "in", res = 600)
plot(x, y, xlab = "Baseline Mucosal Amplitude", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")
## dev.off()




####################################################################
## FAR (freq / amp ratio)
####################################################################

data$mucfar <- data$mucfreq/data$mucamp
x <- data$mucfar
## stange neg number
#idx <- which(data$cut.far < 0)
##     bTSS  cutdate cutfreq cutamp   cut.far tempdate temp.yes  mucdate
## 308   NA 10/27/09   16.67  -0.09 -185.2222 10/29/09        1 10/29/09
## Probably a mistake - should be zero?

data$muc5FAR <- data$muc5Freq/data$muc5Amp
y <- data$muc5FAR - data$mucfar
## plot(x, y)
cor(x, y, use = "complete")
## -0.58748
lm.fit <- lm(y ~ x)
summary(lm.fit)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 60.33900    7.07076   8.534 2.28e-15 ***
## x           -0.68815    0.06376 -10.792  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 68.35 on 221 degrees of freedom
##   (350 observations deleted due to missingness)
## Multiple R-squared:  0.3451,	Adjusted R-squared:  0.3422
## F-statistic: 116.5 on 1 and 221 DF,  p-value: < 2.2e-16


## abline(lm.fit)
lo.fit <- loess(y ~ x)
xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
ypred <- predict(lo.fit, data.frame(x = xseq))
## lines(xseq, ypred, col = 2)

## Plot of Amplitude changes
## pdf("FAR.pdf")
## png("FAR.png", width = 7, height = 7, units = "in", res = 600)
plot(x, y, xlab = "Baseline Mucosal FAR", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")

dev.off()