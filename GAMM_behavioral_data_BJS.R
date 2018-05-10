# Clear variables from the memory
closeAllConnections()
rm(list=ls())

############################################################
# (1) Set environment and load dataset and libraries
############################################################
## (1.1) Set environment
readPath <- setwd( "C:/Users/sarah/OneDrive - Temple University/PhD Temple 2015-2020/PhD Research/Video analysis/Behavioral data")

## (1.2) Load .csv dataset
Data <- read.csv("Behavioral_data_by_arousal_Survivors_BJS.csv")

## (1.3) Load libraries
library(lme4)
library(lmerTest)
library(fitdistrplus)
library(car)
library(ggplot2)
library(mgcv)
library(itsadug)
library(AER)

############################################################
# (2) Format conversions
############################################################
## (2.1) Date and time formating
Data$First_movement <- as.POSIXct(as.character(Data$First_movement), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Arousal_start <- as.POSIXct(as.character(Data$Arousal_start), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Active_start <- as.POSIXct(as.character(Data$Active_start), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Arousal_end <- as.POSIXct(as.character(Data$Arousal_end), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Last_movement <- as.POSIXct(as.character(Data$Last_movement), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Date_died <- as.POSIXct(as.character(Data$Date_died), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")

## (2.2) Factor formating
Data$Bat_ID <- factor(Data$Bat_ID)
Data$Group <- factor(Data$Group)

############################################################
# (3) Fit generalized additive mixed models  
############################################################
## (3.1) GAMM with  response variable Grooming_TotalDur_after
# (3.1.1) Add 1 to grooming duration values in dataset (to remove zeros)
Data_Grooming_TotalDur_after <- Data
Data_Grooming_TotalDur_after$Grooming_TotalDur_after <- (Data_Grooming_TotalDur_after$Grooming_TotalDur_after + 1)

# (3.1.2) Plot Grooming_TotalDur_after as a function of Day
plot(Grooming_TotalDur_after ~ Day, pch = 20, data = Data_Grooming_TotalDur_after)

# (3.1.3) Plot Grooming_TotalDur_after as a function of Day and Group
ggplot(data = Data_Grooming_TotalDur_after, aes(x = Day, y = Grooming_TotalDur_after, colour = Group, group = Group)) +  
       geom_point(aes(col = Data_Grooming_TotalDur_after$Group)) +
       geom_smooth(se = T)

# (3.1.4) Run GAMM models with normal distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm1 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), 
              data = Data_Grooming_TotalDur_after)
summary(gamm1$gam) # gam style summary of fitted model
anova(gamm1$gam) 
# Inspect model residuals:
check_resid(gamm1$gam)
gam.check(gamm1$gam)
plot(gamm1$gam, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day smoothing term:
gamm2 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day), 
              random = list(Bat_ID =~ 1),
              data = Data_Grooming_TotalDur_after)
summary(gamm2$gam)
# Inspect model residuals:
check_resid(gamm2$gam)
gam.check(gamm2$gam)
plot(gamm2$gam, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm3 <- gamm(Grooming_TotalDur_after ~ OFGroup 
              + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after)
summary(gamm3$gam)
# Inspect model residuals:
check_resid(gamm3$gam)
gam.check(gamm3$gam)
plot(gamm3$gam, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm1$lme, gamm2$lme, gamm3$lme)

# (3.1.5) Run GAMM models with gamma distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm11 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), 
              data = Data_Grooming_TotalDur_after,
              family = Gamma(),
              niterPQL = 100)
summary(gamm11$gam)
# Inspect model residuals:
check_resid(gamm11$gam)
gam.check(gamm11$gam)
plot(gamm11$gam, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day smoothing term:
gamm12 <- gamm(Grooming_TotalDur_after ~ Group 
             + s(Day), 
             random = list(Bat_ID =~ 1),
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gamm12$gam)
# Inspect model residuals:
check_resid(gamm12$gam)
gam.check(gamm12$gam)
plot(gamm12$gam, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm13 <- gamm(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
             random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gamm13$gam)
# Inspect model residuals:
check_resid(gamm13$gam)
gam.check(gamm13$gam)
plot(gamm13$gam, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm11$lme, gamm12$lme, gamm13$lme)

# (3.1.6) Run GAMM models with negative binomial distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm21 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group),
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = negbin(c(1,10)),
               niterPQL = 100)
summary(gamm21$gam)
anova(gamm21$gam)
# Inspect distribution of model residuals:
check_resid(gamm21$gam)
gam.check(gamm21$gam)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gamm21$gam, type = "pearson")^2)                     # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-(length(coef(gamm21$gam))+1)     # estimated resid df (N-p)
resid.ssq/resid.df                                                              # ratio should be approx 1

plot(gamm21, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Day smoothing term:
gamm22 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day), # interaction term for factor (Group)
               random = list(Bat_ID =~ 1),
               data = Data_Grooming_TotalDur_after,
               family = negbin(c(1,10)))
summary(gamm22$gam)
# Inspect model residuals:
check_resid(gamm22$gam)
gam.check(gamm22$gam)
plot(gamm22$gam, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm23 <- gamm(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
             random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
             data = Data_Grooming_TotalDur_after,
             family = negbin(c(1,10)))
summary(gamm23$gam)
# Inspect model residuals:
check_resid(gamm23$gam)
gam.check(gamm23$gam)
plot(gamm23$gam, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm21, gamm22, gamm23)