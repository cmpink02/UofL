####################################################################
## Analysis of Mississippi EGG and GET (Western Blot) data
## 11/05/14
####################################################################

load("egg_get.rda")

save.image("egg_get.rda")

setwd("C:\\Users\\Christy\\Dropbox\\Consulting\\EGG Abell")

data <- read.csv("Mississippi_GET.csv", header = TRUE)

####################################################################
## First Analysis
## Look at EGG difference pre/post in freq, amp, and ratio (far)
####################################################################


## NOTE - Need to combine plots into 1x3 figure ...

## Frequency
x <- data$cutfreq
y <- with(data, Temp.Cut.Freq - cutfreq)
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
## Multiple R-squared:  0.4438,	Adjusted R-squared:  0.4413
## F-statistic: 181.9 on 1 and 228 DF,  p-value: < 2.2e-16

lo.fit <- loess(y ~ x)
xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
ypred <- predict(lo.fit, data.frame(x = xseq))

## Plot of frequency changes
## pdf("Freq.pdf")
## png("Freq.png", width = 7, height = 7, units = "in", res = 600)
par(mfrow = c(1, 3))
par(mar = c(5, 4, 2, 2))

plot(x, y, xlab = "Baseline EGG Frequency", ylab = "Difference from Baseline")
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

x <- data$cutamp
y <- data$Temp.Cut.Amp - data$cutamp
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
plot(x, y, xlab = "Baseline EGG Amplitude", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")
## dev.off()

####################################################################
## FAR (freq / amp ratio)
####################################################################

data$cut.far <- data$cutfreq/data$cutamp
x <- data$cut.far
## stange neg number
idx <- which(data$cut.far < 0)
##     bTSS  cutdate cutfreq cutamp   cut.far tempdate temp.yes  mucdate
## 308   NA 10/27/09   16.67  -0.09 -185.2222 10/29/09        1 10/29/09
## Probably a mistake - should be zero?

data$Temp.Cut.FAR <- data$Temp.Cut.Freq/data$Temp.Cut.Amp
y <- data$Temp.Cut.FAR - data$cut.far
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
## pdf("FAR.pdf")x <- data$cutfreq
y <- with(data, Temp.Cut.Freq - cutfreq)
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
## Multiple R-squared:  0.4438,	Adjusted R-squared:  0.4413
## F-statistic: 181.9 on 1 and 228 DF,  p-value: < 2.2e-16

lo.fit <- loess(y ~ x)
xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = 0.1)
ypred <- predict(lo.fit, data.frame(x = xseq))

## Plot of frequency changes
## pdf("Freq.pdf")
## png("Freq.png", width = 7, height = 7, units = "in", res = 600)
par(mfrow = c(1, 3))
par(mar = c(5, 4, 2, 2))

plot(x, y, xlab = "Baseline EGG Frequency", ylab = "Difference from Baseline")
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

x <- data$cutamp
y <- data$Temp.Cut.Amp - data$cutamp
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
plot(x, y, xlab = "Baseline EGG Amplitude", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")
## dev.off()

####################################################################
## FAR (freq / amp ratio)
####################################################################

data$cut.far <- data$cutfreq/data$cutamp
x <- data$cut.far
## stange neg number
idx <- which(data$cut.far < 0)
##     bTSS  cutdate cutfreq cutamp   cut.far tempdate temp.yes  mucdate
## 308   NA 10/27/09   16.67  -0.09 -185.2222 10/29/09        1 10/29/09
## Probably a mistake - should be zero?

data$Temp.Cut.FAR <- data$Temp.Cut.Freq/data$Temp.Cut.Amp
y <- data$Temp.Cut.FAR - data$cut.far
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
plot(x, y, xlab = "Baseline EGG FAR", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")

## png("FAR.png", width = 7, height = 7, units = "in", res = 600)
plot(x, y, xlab = "Baseline EGG FAR", ylab = "Difference from Baseline")
abline(lm.fit)
lines(xseq, ypred, col = 2)
legend("topright", c("Linear", "Loess"), lty = 1, col = 1:2, bty = "n")

dev.off()


## Need to tidy up figures but otherwise am done

## Overall characteristics of data
## DM
table(data$dm)
sum(is.na(data$dm))  ## 22
sum(is.na(data$male)) ## 22

table(data$dm)



tab <- xtabs(~xgrp + dm, data = data)
chisq.test(tab)
## X-squared = 4.6906, df = 3, p-value = 0.1959

## Gender
tab <- xtabs(~xgrp + male, data = data)
chisq.test(tab)
## X-squared = 5.5459, df = 3, p-value = 0.1359

## Age
summary(lm(data$age ~ xgrp))

dim(data)
sum(is.na(data$cutamp))
sum(is.na(x))
sum(is.na(y))


####################################################################
## 11/18/2014
## Updates for EGG changes  vs. baseline
####################################################################

## 1. Need to divide baseline scores into quartiles and calculate
##      differences per quartile

## Frequency
x <- data$cutfreq
xgrp <- cut(x, summary(x)[c(1, 2, 3, 5, 6)], include.lowest = TRUE)
table(xgrp)
## [0.2,3.6] (3.6,4.6] (4.6,6.2]  (6.2,18]
##       104       115        92       102
sum(is.na(x))  ## 160, ok
y <- with(data, Temp.Cut.Freq - cutfreq)
tapply(y, xgrp, function(x) median(x, na.rm = TRUE))
## [0.2,3.6] (3.6,4.6] (4.6,6.2]  (6.2,18]
##      1.27      0.46     -0.63     -1.50
tapply(y, xgrp, function(x) mean(x, na.rm = TRUE))
 ## [0.2,3.6]  (3.6,4.6]  (4.6,6.2]   (6.2,18]
 ## 1.8860000  0.7423529  0.2814894 -2.0175000
res <- tapply(y, xgrp, function(x) quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
res <- do.call(rbind, res)
res <- round(res, 2)
## Use this ...
## NOTE - Need to tack on MEANS as well since t-test is based on these ...
means <- tapply(y, xgrp, function(x) mean(x, na.rm = TRUE))
means <- round(means, 2)

## Look at paired t-tests for each group
p.t <- tapply(y, xgrp, function(x) t.test(x, mu = 0)$p.value)
##    [0.2,3.6]    (3.6,4.6]    (4.6,6.2]     (6.2,18]
## 1.640512e-07 1.334924e-05 4.402724e-01 1.918245e-06

## Anova
anova <- summary(lm(y ~ xgrp))
names(anova)
## F-statistic:  29.5 on 3 and 226 DF,  p-value: 3.936e-16

## Associations of baseline EGG groups with patient characteristics
## DM
tab <- xtabs(~xgrp + dm, data = data)
chisq.test(tab)
## X-squared = 4.6906, df = 3, p-value = 0.1959

## Gender
tab <- xtabs(~xgrp + male, data = data)
chisq.test(tab)
## X-squared = 5.5459, df = 3, p-value = 0.1359

## Age
summary(lm(data$age ~ xgrp))
## F-statistic: 1.536 on 3 and 409 DF,  p-value: 0.2045

## Delayed emptying
data$delayed.empty <- ifelse(data$blgetsol4h > 10, 1, 0)  ## emptying later
xtabs(~delayed.empty + b.delayed, data = data)  ## nearly uniform agreement ...
tab <- xtabs(~xgrp + delayed.empty, data = data)
chisq.test(tab)
## X-squared = 9.4743, df = 3, p-value = 0.02361
round(prop.table(tab, margin = 1), 2)
##            delayed.empty
## xgrp           0    1
##   [0.2,3.6] 0.27 0.73
##   (3.6,4.6] 0.36 0.64
##   (4.6,6.2] 0.48 0.52
##   (6.2,18]  0.43 0.57

tab <- table(data$delayed.empty)
prop.table(tab)


prop.trend.test(tab[,2], rowSums(tab))
## X-squared = 6.7972, df = 1, p-value = 0.00913

cor(x, data$delayed.empty, use="complete")
cor.test(x, data$blgetsol4h, use="complete")

egg.lab <- c("Frequency", rep("", 3), "Amplitude", rep("", 3), "FAR", rep("", 3))
baseline <- c(levels(xgrp), rep(NA, 8))
diff <- c(paste(means, ", ", res[,2], " (", res[,1], ", ", res[,3], ")", sep = ""), rep("", 8))
p.t <- c(p.t, rep(NA, 8))
N <- c(table(xgrp), rep(NA, 8))

res.df <- data.frame(EGG = egg.lab, Baseline = baseline, N = N, Difference = diff, Pvalue = p.t,
                     stringsAsFactors = FALSE)

####################################################################
## Amplitude
####################################################################

idx <- which(data$cutamp < 0)
data$cutamp[idx]

## remove observation < 0
x <- data$cutamp[-idx]
y <- data$Temp.Cut.Amp[-idx] - data$cutamp[-idx]

xgrp <- cut(x, summary(x)[c(1, 2, 3, 5, 6)], include.lowest = TRUE)
table(xgrp)
sum(is.na(x))  ## 160, ok
res <- tapply(y, xgrp, function(x) quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
res <- do.call(rbind, res)
res <- round(res, 2)
## Use this ...

means <- tapply(y, xgrp, function(x) mean(x, na.rm = TRUE))
means <- round(means, 2)


## Look at paired t-tests for each group
p.t <- tapply(y, xgrp, function(x) t.test(x, mu = 0)$p.value)
round(p.t, 3)
## [0.01,0.04] (0.04,0.07] (0.07,0.13]  (0.13,2.1]
##       0.000       0.000       0.193       0.004
p.w <- tapply(y, xgrp, function(x) wilcox.test(x, mu = 0)$p.value)
## [0.01,0.04] (0.04,0.07] (0.07,0.13]  (0.13,2.1]
##       0.000       0.000       0.534       0.001

## Same result ...


## Anova
anova <- summary(lm(y ~ xgrp))
anova
## F-statistic:  12.7 on 3 and 219 DF,  p-value: 1.098e-07

## Look at paired t-tests for each group
p.t <- tapply(y, xgrp, function(x) t.test(x, mu = 0)$p.value)
round(p.t, 3)
## [0.01,0.04] (0.04,0.07] (0.07,0.13]  (0.13,2.1]
##       0.000       0.000       0.193       0.004

## Anova
anova <- summary(lm(y ~ xgrp))
anova
## F-statistic: 25.07 on 3 and 218 DF,  p-value: 5.676e-14

## Associations of baseline EGG groups with patient characteristics
## DM
tab <- xtabs(~xgrp + dm, data = data[-idx, ])
chisq.test(tab)
## X-squared = 10.3111, df = 3, p-value = 0.0161
prop.table(tab, margin = 1)
## try trend test - probably not significant ...
prop.trend.test(tab[,2], rowSums(tab))
## X-squared = 1.7326, df = 1, p-value = 0.1881

## Gender
tab <- xtabs(~xgrp + male, data = data[-idx, ])
chisq.test(tab)
## X-squared = 5.4836, df = 3, p-value = 0.1396

## Age
summary(lm(data$age[-idx] ~ xgrp))
## F-statistic: 0.3266 on 3 and 407 DF,  p-value: 0.8061

## Delayed emptying
tab <- xtabs(~xgrp + delayed.empty, data = data[-idx, ])
chisq.test(tab)
## X-squared = 1.5608, df = 3, p-value = 0.6683
round(prop.table(tab, margin = 1), 2)
##              delayed.empty
## xgrp             0    1
##   [0.01,0.04] 0.43 0.57
##   (0.04,0.07] 0.36 0.64
##   (0.07,0.13] 0.38 0.62
##   (0.13,2.1]  0.34 0.66

prop.trend.test(tab[,2], rowSums(tab))
## X-squared = 1.171, df = 1, p-value = 0.2792
## Neither significant ...


## res.df <- data.frame(EGG = egg.lab, Baseline = baseline, N = N, Difference = diff, Pvalue = p.t)

res.df$Baseline[5:8] <- levels(xgrp)
res.df$Difference[5:8] <- paste(means, ", ", res[,2], " (", res[,1], ", ", res[,3], ")", sep = "")
res.df$Pvalue[5:8] <- p.t
res.df$N[5:8] <- table(xgrp)


####################################################################
## FAR (freq / amp ratio)
####################################################################

data$cut.far <- data$cutfreq/data$cutamp
idx <- which(data$cut.far < 0)
x <- data$cut.far[-idx]
## stange neg number
data$Temp.Cut.FAR <- data$Temp.Cut.Freq/data$Temp.Cut.Amp
y <- data$Temp.Cut.FAR[-idx] - data$cut.far[-idx]


xgrp <- cut(x, summary(x)[c(1, 2, 3, 5, 6)], include.lowest = TRUE)
table(xgrp)
sum(is.na(x))  ## 160, ok
res <- tapply(y, xgrp, function(x) quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
res <- do.call(rbind, res)
res <- round(res, 2)
## Use this ...

means <- tapply(y, xgrp, function(x) mean(x, na.rm = TRUE))
means <- round(means, 2)


## Look at paired t-tests for each group
p.t <- tapply(y, xgrp, function(x) t.test(x, mu = 0)$p.value)
round(p.t, 3)
## [0.01,0.04] (0.04,0.07] (0.07,0.13]  (0.13,2.1]
##       0.000       0.000       0.193       0.004

## Anova
anova <- summary(lm(y ~ xgrp))
anova
## F-statistic: 25.07 on 3 and 218 DF,  p-value: 5.676e-14

## Associations of baseline EGG groups with patient characteristics
## DM
tab <- xtabs(~xgrp + dm, data = data[-idx, ])
chisq.test(tab)
## X-squared = 5.9869, df = 3, p-value = 0.1123

## Gender
tab <- xtabs(~xgrp + male, data = data[-idx, ])
chisq.test(tab)
## X-squared = 10.5449, df = 3, p-value = 0.01446
prop.table(tab, margin = 1)
## try trend test - probably not significant ...
prop.trend.test(tab[,2], rowSums(tab))
## X-squared = 0.7219, df = 1, p-value = 0.3955

## Age
summary(lm(data$age[-idx] ~ xgrp))
## F-statistic: 0.3018 on 3 and 407 DF,  p-value: 0.8241

## Delayed emptying
tab <- xtabs(~xgrp + delayed.empty, data = data[-idx, ])
chisq.test(tab)
## X-squared = 3.6282, df = 3, p-value = 0.3045

## res.df <- data.frame(EGG = egg.lab, Baseline = baseline, N = N, Difference = diff, Pvalue = p.t)
res.df$Baseline[9:12] <- levels(xgrp)
res.df$Difference[9:12] <- paste(means, ", ", res[,2], " (", res[,1], ", ", res[,3], ")", sep = "")
res.df$Pvalue[9:12] <- p.t
res.df$N[9:12] <- table(xgrp)

res.df
write.csv(res.df, file = "EGG_results.csv")

####################################################################
## New Items 11/18/14
## 1. Can you describe: Age, Sex, Diagnosis (it will be mainly DM or non DM=ID)?
## 2. Can you look at correlation of BASELine EGG WITH SX, and GET values (delayed >10% of 4 ours)
## 3. Lastly can you do graphs of change in sx and change in GET, just like the EGG, to see if similar?
####################################################################

## bget4h > 10, tget4h??
vars <- c("age", "male", "idiopathic", "dm")
summary(data[,vars])
summary(data$age)
##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
## 12.00   38.00   47.00   47.39   57.00   86.00      23

idx <- which(is.na(data$age))

table(data$male[-idx])
##   0   1
## 441 110
sum(is.na(data$male))  ## 22
prop.table(table(data$male[-idx]))
##        0        1
## 0.800363 0.199637

table(data$dm[-idx])
##   0   1
## 382 169
sum(is.na(data$dm))  ## 22
prop.table(table(data$dm[-idx]))
##         0         1
## 0.6932849 0.3067151

sum(is.na(data$age))  ## 23
dim(data)
nrow(data) - 22
xtabs(~idiopathic + dm, data = data)
## no overlap except when both zero
table(data$postsur)  ## 38 1's

xtabs(~delayed.empty + b.delayed, data = data)
table(data$delayed.empty)
prop.table(table(data$delayed.empty))
sum(is.na(data$delayed.empty))  ## 120

####################################################################
## Graphs for change in baseline emptying?
####################################################################

summary(data$blgetsol4h)
summary(data$tget4h)
summary(data$tget4h - data$blgetsol4h)
x <- data$blgetsol4h
y <- data$tget4h - data$blgetsol4h
## Plot vs. baseline?
pdf("Change in gastric emptying.pdf")
plot(x, y, xlab = "Baseline percent gastric emptying at 4 hours",
     ylab = "Difference post-pre gastric emptying at 4 hours")
dev.off()

## just an artifact of constraints?
xx <- runif(100, 0, 100)
yy <- runif(100, 0, 100)
pdf("Random differences.pdf")
plot(xx, yy-xx, xlab = "Random baseline",
     ylab = "Difference between random baseline and follow-up")
dev.off()




## Well, pretty strong negative correlation here ...
cor.test(x, y)
## t = -10.0432, df = 378, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.5348567 -0.3756943
## sample estimates:
##        cor
## -0.4589499

table(data$b.delayed)
xtabs(~delayed.empty + b.delayed, data = data)


####################################################################
## Biomarkers from Western Blot (GET)
####################################################################

## GET variables
## First number is protein and following numbers are bands
## scl70	p100
## ssa	p60	p52
## rpn	p68	p34	p22
## pribo	p38
## AntiCen	p1718
## Ku	p86	p80	p66
## ssb	p48	p43
## Sm	p29	p28	p16
## Jo1	p56
## UnDef	 UnDef
## Composite Variables
## classic without undefined	GIBS (GI Blotting Score)	Total groups of antibodies

data$scl70 <- ifelse(data$p100 > 0, 1, 0)
table(data$scl70)
##   0   1
## 120  11
data$ssa <- ifelse(data$p60 > 0 | data$p52 > 0, 1, 0)

data$rnp <- ifelse(data$p68 > 0 | data$p34 > 0 | data$p22 > 0, 1, 0)
data$rnp.tot <- data$p68 + data$p34 + data$p22
table(data$rnp.tot)
## GIBS 2, 5, and 15 based on 1, 2, and 3 bands present
data$rnp.gibs <- ifelse(data$rnp.tot == 0, 0,
                        ifelse(data$rnp.tot > 0 & data$rnp.tot <=4, 2,
                               ifelse(data$rnp.tot > 4 & data$rnp.tot <=7, 5, 15)))
table(data$rnp.gibs)


data$pribo <- ifelse(data$p38 > 0, 1, 0)
data$AntiCen <- ifelse(data$p1718 > 0, 1, 0)
data$Ku <- ifelse(data$p86 > 0 | data$p80 > 0 | data$p66 > 0, 1, 0)
data$ssb <- ifelse(data$p48 > 0 | data$p43 > 0, 1, 0)
data$Sm <- ifelse(data$p29 > 0 | data$p28 > 0 | data$p16 > 0, 1, 0)
data$Jo1 <- ifelse(data$p56 > 0, 1, 0)
data$undef01 <- ifelse(data$UnDef > 0 | data$UnDef.1 > 0, 1, 0)
## Everyone has one or the other so not useful ...
table(data$UnDef)  ## This one is actually empty
table(data$UnDef.1)  ## Total number of undefined bands it seems ...

## Check out distribution of these
vars <- c("ssa", "rnp", "pribo", "AntiCen", "Ku", "ssb", "Sm", "Jo1")
apply(data[,vars], 2, table)
##   ssa rnp pribo AntiCen  Ku ssb  Sm Jo1
## 0  37  43    99     128  54  65  90 107
## 1 186 183    67       1 155 125 102   9



data$rnp.tot <- data$p68 + data$p34 + data$p22
table(data$rnp.tot)
##  0  3  4  6  7  9 10 11 12
## 43 54  4 18  5  2  1  3  1

## rnp association with symptom scores
vars2 <- c("bV", "bN", "bANORES", "bBL", "bAP", "bTSS")
cor(data$rnp.tot, data[,vars2], use = "complete")
##             bV        bN   bANORES       bBL          bAP     bTSS
## [1,] 0.1189888 0.1414821 0.3110871 0.2441019 -0.007742139 0.275678


cor.test(data$rnp.tot, data$bTSS)
## data:  data$rnp.tot and data$bTSS
## t = 2.6025, df = 83, p-value = 0.01096
cor.test(data$rnp.tot, data$bANORES)
## t = 2.9641, df = 82, p-value = 0.003973
cor.test(data$rnp.tot, data$bN)
## t = 1.2999, df = 83, p-value = 0.1972

vars2 <- c("bV", "bN", "bANORES", "bBL", "bAP", "bTSS")
cor(data$rnp.gibs, data[,vars2], use = "complete")
##              bV         bN   bANORES       bBL         bAP   bTSS
## [1,] 0.03226235 0.05599655 0.2683114 0.1624263 -0.02885454 0.1678

tab <- xtabs(~rnp.gibs + b.delayed, data = data)
round(prop.table(tab, margin = 1), 2)
chisq.test(tab)

table(data$b.delayed)
prop.table(table(data$b.delayed))


## Remove AntiCen before calculating t-tests below since only 1
vars <- c("ssa", "rnp", "pribo", "Ku", "ssb", "Sm", "Jo1")

## Symptom Scores
## bV	baseline vomiting score	0-4
## bN	baseline nausea score	0-4
## bANORES	baseline anorexia score	0-4
## bBL	baseline bloating score	0-4
## bAP	baseline abdominal pain score	0-4
## bTSS	baseline total symptom score	0-20

## Correlation between symptom scores and GET variables ...
## Initially try just running a bunch of t-tests?
## Big problem with multiple comparisons though ...
vars2 <- c("bV", "bN", "bANORES", "bBL", "bAP", "bTSS")
pvals <- matrix(NA, nrow = length(vars), ncol = length(vars2))
rownames(pvals) <- vars
colnames(pvals) <- vars2
for (i in 1:length(vars)) {  ## bands
    for (j in 1:length(vars2)) {  ## symptom scores
        res <- try(t.test(data[,vars2[j]] ~ data[,vars[i]]))
        if (class(res) == "try-error")
            next
        pvals[i,j] <- res$p.value
    }
}
round(pvals, 3)


apply(data[,vars2], 2, function(x) tapply(x, data$pribo, function(x) mean(x, na.rm = TRUE)))


library(multtest)
padj <- mt.rawp2adjp(pvals, proc=c("BH", "BY"))
padj.bh <- padj$adjp[order(padj$index),"BH"]
##             rawp        BH        BY
## [1,] 0.004977675 0.1782341 0.7711733
## [2,] 0.008487340 0.1782341 0.7711733


## Nothing passes multiple comparisons adjustment ...
padj.by <- padj$adjp[order(padj$index),"BY"]
res <- data.frame(nams, numbers, hrs.text, pvals, padj.bh, padj.by)


round(pvals, 3)

## NOTE - Didn't use BASELINE emptying below ...


## Gastric Emptying
## tget4h	post tempGES gastric emptying forth hour	%
## NOTE - This is DELAYED if >10% (or should be <10%?)

## NOTE - Here looked at POST instead of BASELINE by mistake ...
summary(data$tget1h)
summary(data$tget4h)
data$delayed.empty <- ifelse(data$tget4h > 10, 1, 0)  ## emptying later
table(data$delayed.empty)
##   0   1
## 150 263
## Most have delayed emptying

pvals2 <- numeric(length(vars))
for (i in 1:length(pvals2)) {
    res <- chisq.test(data[,vars[i]], data$delayed.empty)
    pvals2[i] <- res$p.value
}
names(pvals2) <- vars
round(pvals2, 2)
##  0.5396787 0.5189102 0.8214060 0.9991661 0.9721638 0.9279907 0.4685504
## Nothing associated ...




