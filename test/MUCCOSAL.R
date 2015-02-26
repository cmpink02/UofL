####################################################################
## Analysis of Mississippi EGG and GET (Western Blot) data
## 11/05/14
####################################################################

load("egg_get.rda")

save.image("egg_get.rda")

#setwd("C:\\Users\\Christy\\Dropbox\\Consulting\\EGG Abell")
setwd("C:\\Users\\cmpink02\\Dropbox\\Consulting\\EGG Abell")

#data <- read.csv("Mississippi_GET.csv", header = TRUE)
data <- read.csv("muc_all.csv", header = TRUE)

### Use this to re-create the original output
# data<-data[which(data$Study == "Original"),]

### Use this to exclude the weird value
# data$cutamp[data$cutamp < 0] <- NA
####################################################################
## First Analysis
## Look at EGG difference pre/post in freq, amp, and ratio (far)
####################################################################


## NOTE - Need to combine plots into 1x3 figure ...

## Frequency
x <- data$mucfreq
y <- with(data, muc5Freq - mucfreq)
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
par(mfrow = c(1, 3))
par(mar = c(5, 4, 2, 2))

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