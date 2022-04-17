knitr::opts_chunk$set(echo = TRUE)
rm(list = ls())
library(gdata)




sur <- read.table("http://web1.sph.emory.edu/dkleinb/allDatasets/surv2datasets/addicts.dat", skip=19)

colnames(sur) <- c("ID", "Clinic", "Status", "Survival", "Prison", "Methodone")

head(sur)
sur <- as.data.frame(sur)
sur <- sur[,-1]


dev_off()
hist(x=sur$Methodone,pch = 19,col ="lightgray",xlab="Methodone dose",boxwex = 0.8,medlwd = 0.5,border="blue",outcex=0.7)
# Clinic is a factor variable, everything else is numeric.

sur$Clinic <- as.factor(sur$Clinic)
sur$Status <- as.numeric(sur$Status)
sur$Survival <- as.numeric(sur$Survival)
sur$Prison <- as.numeric(sur$Prison)
sur$Methodone <- as.numeric(sur$Methodone)

# Creation of a new variable in order to modelize the Kaplan-Meier Curve with Methadone data
sur$Dosage[sur$Methodone < 55] <- "Low"
sur$Dosage[(sur$Methodone <= 70) & (sur$Methodone >=55)] <- "Medium"
sur$Dosage[sur$Methodone > 70] <- "High"
sur$Dosage <- ordered(sur$Dosage, levels = c("Low", "Medium", "High")) # Reorder since it is an ordinal factor
head(sur)

# Check sur types:
str(sur)

library(survival)
library(ggplot2)
library(ggfortify)
library(ggsurvplot)
library(survminer)

head(sur)

#Définition des différents modèles
model_fit_Clinic <- survfit(Surv(Survival, Status) ~ Clinic, data = sur)

model_fit_Prison <- survfit(Surv(Survival, Status) ~ Prison, data = sur)

model_fit_Dosage <- survfit(Surv(Survival, Status) ~ Dosage, data = sur)

#Log-Rank test for comparing different group
library("coin") 
 

survdiff(Surv(Survival, Status) ~ Prison, data = sur)
logrank_test(Surv(Survival, Status) ~ Prison, data = sur, distribution = "exact")

survdiff(Surv(Survival, Status) ~ Dosage, data = sur)

logrank_test(Surv(Survival, Status) ~ Dosage, data = sur, distribution = "exact")

#Basic KM curves without results table
ggsurvplot(model_fit_Clinic, data = sur)
ggsurvplot(model_fit_Prison, data = sur, pval = TRUE,conf.int = T)
ggsurvplot(model_fit_Dosage, data = sur, pval = TRUE,conf.int = T)

#KM with grid
autoplot(model_fit_ur, pval = TRUE,conf.int = T)Clinic) + labs(x = "\n Survival Time (Days)", y = "Survival Probabilities \n", title = "Survival Time of \n Methadone Patients \n") 
autoplot(model_fit_Prison) + labs(x = "\n Survival Time (Days)", y = "Survival Probabilities \n", title = "Survival Time of \n Methadone Patients \n") 
autoplot(model_fit_Dosage) + labs(x = "\n Survival Time (Days)", y = "Survival Probabilities \n", title = "Survival Time of \n Methadone Patients \n") 


library("RTCGA.clinical")
library(survminer)

#Full KM curves with numbers at risk
#Clinic
ggsurvplot(
  model_fit_Clinic,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

#Prison
ggsurvplot(
  model_fit_Prison,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

#Dosage
ggsurvplot(
  model_fit_Dosage,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1000),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

ggsurv <- ggsurvplot(
  model_fit_Dosage,                     # survfit object with calculated statistics.
  data = sur,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF","#008000"),
  xlim = c(0,500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Low", "Medium","High")    # change legend labels.
)
ggsurv




#We model the hazard function using Cox regression with all variables

all_model <- coxph(Surv(Survival, Status) ~ Clinic + Prison + Dosage, data = sur)
summary(all_model)

ggforest(all_model)
sur

ci <- confint(GBSG2.coxph)
exp(cbind(coef(GBSG2.coxph), ci))["horThyes",] 

GBSG2.zph <- cox.zph(all_model)
GBSG2.zph 
plot(GBSG2.zph, var = "Methodone")


res <- residuals(all_model,type = "martingale")
plot(res ~ Methodone, data = sur, ylim = c(-2.5, 2.5), pch=".", ylab = "Martingale residuals")
abline(h=0, lty=2) 


fit2=coxph(Surv(time,delta)~wtime+factor(dtype)+factor(gtype)+score,data=hodg,method='breslow')
resid(fit2,type='martingale')
plot(hodg$wtime, resid(fit2),
     xlab="Waiting Time to Transplant (months)", ylab="Martingale Residuals",
     main='Figure 11.4 on page 361')
lines(lowess(hodg$wtime, resid(fit2)),col='red')

#Model investigation and diagnostic
names(all_model)


ggcoxdiagnostics(all_model, type = "deviance",ggtheme = ggplot2::theme_bw())

