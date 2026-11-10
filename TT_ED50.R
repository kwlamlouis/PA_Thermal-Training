#clear workspace
rm(list=ls())

#set language to English
Sys.setenv(LANG = "en")

####<FVFM>####
####TT####
####1)Data####
library(tidyverse) #for "group_by", "summarize"
library(ggplot2)

fvfm <- read_csv("TT_ED50.csv")

fvfm$Mother <-as.factor(fvfm$Mother)
fvfm$colony_id <-as.factor (fvfm$colony_id)
fvfm$Age <-as.factor (fvfm$Age)
fvfm$Tx <-as.factor (fvfm$Tx)
fvfm$Combo <-as.factor (fvfm$Combo)

summary(fvfm)
head(fvfm)

####3) Plots####
#Summary (need Treatments' avg, SD, SE)
fvfm_summary<-  
  fvfm%>%
  group_by(Mother, Age, Tx, Combo, temp)%>%
  summarise(mean = mean(fvfm_1), sd=sd(fvfm_1))

#Use dose curve approach; code used is based on approach described by Evensen et al. 2022
#Evensen NR, VoA.26stra CR, Fine M, Perna G, Buitrago-Lopez C, Cardenas A, Banc-Prandi G, Rowe K, Barshis DJ.
#(2022) Empirically derived thermal threshA.26ds of four coral species along the Red Sea using a portable and 
#standardized experimental approach. Coral Reefs 41(2): 239-52. [https://doi.org/10.1007/s00338-022-02233-y]

#Treatment subsets

#Adult_26
A.26_fvfm <-
  fvfm_summary %>%
  filter(Combo == "A_26")
A.26_fvfm

#Adult_29
A.29_fvfm <-
  fvfm_summary %>%
  filter(Combo == "A_29")
A.29_fvfm

#Recruit_26
R.26_fvfm <-
  fvfm_summary %>%
  filter(Combo == "R_26")
R.26_fvfm

#Recruit_29
R.29_fvfm <-
  fvfm_summary %>%
  filter(Combo == "R_29")
R.29_fvfm

####*Adult 26 plots####
library(drc) #for "drm"
library(Rmisc) #for "CI"

#Run population-level fit
A.26_pop <- drm(mean ~ temp, data = A.26_fvfm, fct = LL.3())
summary(A.26_pop)
plot(A.26_pop)

#Run individual colony fits
A.26_indiv <- drm(mean ~ temp, data = A.26_fvfm, curveid=Mother, fct = LL.3())
summary(A.26_indiv)
plot(A.26_indiv)

#extract coeffs by colony, then compute 95% CIs
A.26_genocoeffs_50<-data.frame(ED(A.26_indiv, c(50)))
A.26_coeff_mean<-mean(A.26_genocoeffs_50$Estimate)
A.26_coeff_mean
#36.58099

A.26_summary<-data.frame(CI(A.26_genocoeffs_50$Estimate, ci=0.95))
A.26_coeff_lower<-A.26_summary[3,]
A.26_coeff_upper<-A.26_summary[1,]

#______________________________________________
A.26_genocoeffs_5<-data.frame(ED(A.26_indiv, c(5)))
A.26_coeff_mean_5<-mean(A.26_genocoeffs_5$Estimate)
A.26_coeff_mean_5
#31.38261

A.26_genocoeffs_95<-data.frame(ED(A.26_indiv, c(95)))
A.26_coeff_mean_95<-mean(A.26_genocoeffs_95$Estimate)
A.26_coeff_mean_95
#42.65199
#________________________________________________


####*Adult 29 plots####
#Run population-level fit
A.29_pop <- drm(mean ~ temp, data = A.29_fvfm, fct = LL.3())
summary(A.29_pop)
plot(A.29_pop)

#Run individual colony fits
A.29_indiv <- drm(mean ~ temp, data = A.29_fvfm, curveid=Mother, fct = LL.3())
summary(A.29_indiv)
plot(A.29_indiv)

#extract coeffs by cA.26ony, then compute 95% CIs
A.29_genocoeffs_50<-data.frame(ED(A.29_indiv, c(50)))
A.29_coeff_mean<-mean(A.29_genocoeffs_50$Estimate)
A.29_coeff_mean
#37.28684

A.29_summary<-data.frame(CI(A.29_genocoeffs_50$Estimate, ci=0.95))
A.29_coeff_lower<-A.29_summary[3,]
A.29_coeff_upper<-A.29_summary[1,]


####*Recruit 26 plots####
#Run population-level fit
R.26_pop <- drm(mean ~ temp, data = R.26_fvfm, fct = LL.3())
summary(R.26_pop)
plot(R.26_pop)

#Run individual cA.26ony fits
R.26_indiv <- drm(mean ~ temp, data = R.26_fvfm, curveid=Mother, fct = LL.3())
summary(R.26_indiv)
plot(R.26_indiv)

#extract coeffs by cA.26ony, then compute 95% CIs
R.26_genocoeffs_50<-data.frame(ED(R.26_indiv, c(50)))
R.26_coeff_mean<-mean(R.26_genocoeffs_50$Estimate)
R.26_coeff_mean
#36.65604

R.26_summary<-data.frame(CI(R.26_genocoeffs_50$Estimate, ci=0.95))
R.26_coeff_lower<-R.26_summary[3,]
R.26_coeff_upper<-R.26_summary[1,]



####*Recruit 29 plots####
#Run population-level fit
R.29_pop <- drm(mean ~ temp, data = R.29_fvfm, fct = LL.3())
summary(R.29_pop)
plot(R.29_pop)

#Run individual cA.26ony fits
R.29_indiv <- drm(mean ~ temp, data = R.29_fvfm, curveid=Mother, fct = LL.3())
summary(R.29_indiv)
plot(R.29_indiv)

#extract coeffs by cA.26ony, then compute 95% CIs
R.29_genocoeffs_50<-data.frame(ED(R.29_indiv, c(50)))
R.29_coeff_mean<-mean(R.29_genocoeffs_50$Estimate)
R.29_coeff_mean
#37.56425

R.29_summary<-data.frame(CI(R.29_genocoeffs_50$Estimate, ci=0.95))
R.29_coeff_lower<-R.29_summary[3,]
R.29_coeff_upper<-R.29_summary[1,]


####*Calculate ED50s####
library(reshape2) #for "melt"

#combine cA.26ony-ED50s into data frame for statistical analysis

ED50s <- data.frame(cbind(A.26_genocoeffs_50[,1],A.29_genocoeffs_50[,1],R.26_genocoeffs_50[,1], R.29_genocoeffs_50[,1]))


ED50s <- ED50s %>% 
  dplyr::rename(A.26= X1,
                A.29=X2,
                R.26=X3,
                R.29=X4
  )


ED50s$colony<-as.factor(1:nrow(ED50s))
str(ED50s)

ED50s_long<-melt(ED50s, id="colony")

ED50s_long<-ED50s_long %>% 
  dplyr::rename(Treatment= variable,
                ED50=value)


####*ED50 plot####

PA_coeff_means<-data.frame(A.26_coeff_mean, A.29_coeff_mean, R.26_coeff_mean, R.29_coeff_mean)
PA_coeff_lowers<-data.frame(A.26_coeff_lower,  A.29_coeff_lower, R.26_coeff_lower, R.29_coeff_lower)
PA_coeff_uppers<-data.frame(A.26_coeff_upper,  A.29_coeff_upper, R.26_coeff_lower, R.29_coeff_upper)

A.26_preddata = data.frame(temp = seq(26,39, length.out = 400))
A.26_pred = as.data.frame(predict(A.26_pop, newdata = A.26_preddata, interval = 'confidence'))
A.26_preddata = data.frame(A.26_preddata, mean = A.26_pred$Prediction, Lower = A.26_pred$Lower, Upper = A.26_pred$Upper)

A.29_preddata = data.frame(temp = seq(26,39, length.out = 400))
A.29_pred = as.data.frame(predict(A.29_pop, newdata = A.29_preddata, interval = 'confidence'))
A.29_preddata = data.frame(A.29_preddata, mean =A.29_pred$Prediction, Lower = A.29_pred$Lower, Upper = A.29_pred$Upper)

R.26_preddata = data.frame(temp = seq(26,39, length.out = 400))
R.26_pred = as.data.frame(predict(R.26_pop, newdata = R.26_preddata, interval = 'confidence'))
R.26_preddata = data.frame(R.26_preddata, mean = R.26_pred$Prediction, Lower = R.26_pred$Lower, Upper = R.26_pred$Upper)

R.29_preddata = data.frame(temp = seq(26,39, length.out = 400))
R.29_pred = as.data.frame(predict(R.29_pop, newdata = R.29_preddata, interval = 'confidence'))
R.29_preddata = data.frame(R.29_preddata, mean =R.29_pred$Prediction, Lower = R.29_pred$Lower, Upper = R.29_pred$Upper)


combo_plot <- ggplot() +
  geom_jitter(data = fvfm_summary, aes(x = temp, y = mean, color = Combo), size = 1.5, width = 0.25) +
  scale_x_continuous(limits=c(26,40), breaks=c(27, 30, 33, 36, 39)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = A.26_preddata, aes(x = temp, y = mean), color = 'steelblue', show.legend = FALSE) +
  geom_ribbon(data = A.26_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'steelblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = A.26_coeff_mean), color = 'steelblue', show.legend = FALSE, size =0.75) +
  
  geom_line(data = A.29_preddata, aes(x = temp, y = mean), color = 'red2', show.legend = FALSE) +
  geom_ribbon(data = A.29_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red2', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = A.29_coeff_mean), color = 'red2', show.legend = FALSE, size =0.75) +
  
  geom_line(data = R.26_preddata, aes(x = temp, y = mean), color = 'lightblue', show.legend = FALSE) +
  geom_ribbon(data = R.26_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'lightblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = R.26_coeff_mean), color = 'lightblue', show.legend = FALSE, size =0.75) +
  
  geom_line(data = R.29_preddata, aes(x = temp, y = mean), color = 'violet', show.legend = FALSE) +
  geom_ribbon(data = R.29_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'violet', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = R.29_coeff_mean), color = 'violet', show.legend = FALSE, size =0.75) +
  
  ggtitle("") +
  scale_color_manual(values=c("steelblue", "red2", "lightblue", "violet")) +
  ylab("Maximum quantum yield (Fv/Fm)") +
  xlab("Temperature (°C)") +
  theme_classic()+
  
  theme(text=element_text(size=14,  family="sans"))

combo_plot

#Ratio: 500 x 700

####4) Analysis####
library(lmerTest)
library(emmeans)

#compare ED50s (ANOVA approach)

ED50_lm <- lm(ED50 ~ Treatment, data = ED50s_long)
Anova(ED50_lm, type = "III")

#check assumptions
plot(fitted(ED50_lm),residuals(ED50_lm))
hist(residuals(ED50_lm))
qqnorm(residuals(ED50_lm))

#posthoc 
ED50.posthoc1 <- emmeans(ED50_lm, pairwise ~ Treatment, type="response")
summary(ED50.posthoc1)
# contrast    estimate    SE df t.ratio p.value
# A.26 - R.29  -0.9833 0.253 16  -3.884  0.0065 **
# R.26 - R.29  -0.9082 0.253 16  -3.587  0.0119 *
#P value adjustment: tukey method for comparing a family of 4 estimates 

####_________####
####Ren####
####1) Libraries####

library(tidyverse)
library(ggplot2)
library(drc)
library(lmerTest)
library(emmeans)
library(Rmisc)
library(reshape2)
library(car)
library(patchwork)


####2) Data####
fvfm <- read_csv("Ren ED50/fvfm_data20250811.csv")
summary(fvfm)


fvfm$treatment <-as.factor (fvfm$treatment)
fvfm$temp <-as.numeric(fvfm$temp)
fvfm$colony_id <-as.factor (fvfm$colony_id)
fvfm$nubbin_id <-as.factor (fvfm$nubbin_id)
fvfm$Date <- as.Date(fvfm$Date, format ='%m/%d/%Y')

summary(fvfm)
head(fvfm)


####3) Plots####
#Combine "fvfm1" and "fvfm2" into one column
fvfm_long <- 
  fvfm %>%
  pivot_longer(
    cols = c(fvfm1, fvfm2),
    values_to = "fvfm"
  )

#Date       treatment colony_id nubbin_id  temp fvfm_avg     SD     SE name   fvfm
#2025-07-03 OL        PA1       PA1-1      28.5    0.734 0.0134 0.0095 fvfm1 0.724
#2025-07-03 OL        PA1       PA1-1      28.5    0.734 0.0134 0.0095 fvfm2 0.743

#Summary (need Treatments' avg, SD, SE)
fvfm_long2<-  
  fvfm_long%>%
  group_by(Date,treatment, temp)%>%
  summarise(mean = mean(fvfm), sd=sd(fvfm), n= n(), se=sd/sqrt(n))

#set "temp" as factor here 
fvfm_long2$temp<-as.factor(fvfm_long2$temp)
head(fvfm_long2)

#Draw plot
ggplot(fvfm_long2, aes(x=Date, y=mean, colour=temp))+ 
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  facet_wrap(~treatment)+
  scale_color_manual(values = c("lightblue", "blue", "darkblue", "lightpink", "violet", "red2", "red4"))+
  labs(y="fvfm")+
  theme_bw()+
  theme(text=element_text(size=13,  family="serif"), axis.text.x=element_text(angle=90, vjust=0.5),legend.position ="right")


#Use dose curve approach; code used is based on approach described by Evensen et al. 2022
#Evensen NR, Voolstra CR, Fine M, Perna G, Buitrago-Lopez C, Cardenas A, Banc-Prandi G, Rowe K, Barshis DJ.
#(2022) Empirically derived thermal thresholds of four coral species along the Red Sea using a portable and 
#standardized experimental approach. Coral Reefs 41(2): 239-52. [https://doi.org/10.1007/s00338-022-02233-y]

#Treatment subsets

#OL
OL_fvfm <-
  fvfm %>%
  filter(treatment == "OL")
OL_fvfm

#FTS
FTS_fvfm <-
  fvfm %>%
  filter(treatment == "FTS")
FTS_fvfm

#RAS
fvfm_RAS <-
  fvfm %>%
  filter(treatment == "RAS")
fvfm_RAS


####*OL plots####

#Run population-level fit
OL_pop <- drm(fvfm_avg ~ temp, data = OL_fvfm, fct = LL.3())
summary(OL_pop)
plot(OL_pop)

#Run individual colony fits
OL_indiv <- drm(fvfm_avg ~ temp, data = OL_fvfm, curveid=colony_id, fct = LL.3())
summary(OL_indiv)
plot(OL_indiv)

#extract coeffs by colony, then compute 95% CIs
OL_genocoeffs_50<-data.frame(ED(OL_indiv, c(50)))
OL_coeff_mean<-mean(OL_genocoeffs_50$Estimate)
OL_coeff_mean
#37.28

OL_summary<-data.frame(CI(OL_genocoeffs_50$Estimate, ci=0.95))
OL_coeff_lower<-OL_summary[3,]
OL_coeff_upper<-OL_summary[1,]


####*FTS plots####

#Run population-level fit
FTS_pop <- drm(fvfm_avg ~ temp, data = FTS_fvfm, fct = LL.3())
summary(FTS_pop)
plot(FTS_pop)

#Run individual colony fits
FTS_indiv <- drm(fvfm_avg ~ temp, data = FTS_fvfm, curveid=colony_id, fct = LL.3())
summary(FTS_indiv)
plot(FTS_indiv)

#extract coeffs by colony, then compute 95% CIs
FTS_genocoeffs_50<-data.frame(ED(FTS_indiv, c(50)))
FTS_coeff_mean<-mean(FTS_genocoeffs_50$Estimate)
FTS_coeff_mean
#34.38

FTS_summary<-data.frame(CI(FTS_genocoeffs_50$Estimate, ci=0.95))
FTS_coeff_lower<-FTS_summary[3,]
FTS_coeff_upper<-FTS_summary[1,]



####*RAS plots####

#Run population-level fit
RAS_pop <- drm(fvfm_avg ~ temp, data = fvfm_RAS, fct = LL.3())
summary(RAS_pop)
plot(RAS_pop)

#Run individual colony fits
RAS_indiv <- drm(fvfm_avg ~ temp, data = fvfm_RAS, curveid=colony_id, fct = LL.3())
summary(RAS_indiv)
plot(RAS_indiv)

#extract coeffs by colony, then compute 95% CIs
RAS_genocoeffs_50<-data.frame(ED(RAS_indiv, c(50)))
RAS_coeff_mean<-mean(RAS_genocoeffs_50$Estimate)
RAS_coeff_mean
#34.50

RAS_summary<-data.frame(CI(RAS_genocoeffs_50$Estimate, ci=0.95))
RAS_coeff_lower<-RAS_summary[3,]
RAS_coeff_upper<-RAS_summary[1,]



####*Calculate ED50s####
#combine colony-ED50s into data frame for statistical analysis

ED50s <- data.frame(cbind(OL_genocoeffs_50[,1],FTS_genocoeffs_50[,1],RAS_genocoeffs_50[,1]))
                         

ED50s <- ED50s %>% 
  dplyr::rename(OL= X1,
                FTS=X2,
                RAS=X3,
                )


ED50s$colony_id<-as.factor(1:nrow(ED50s))
str(ED50s)

library(reshape2)
ED50s_long<-melt(ED50s, id="colony_id")

ED50s_long<-ED50s_long %>% 
  dplyr::rename(Treatment= variable,
                ED50=value)


####*ED50 plot####

PA_coeff_means<-data.frame(OL_coeff_mean, FTS_coeff_mean, RAS_coeff_mean)
PA_coeff_lowers<-data.frame(OL_coeff_lower,  FTS_coeff_lower, RAS_coeff_lower)
PA_coeff_uppers<-data.frame(OL_coeff_upper,  FTS_coeff_upper, RAS_coeff_upper)

OL_preddata = data.frame(temp = seq(27,39, length.out = 400))
OL_pred = as.data.frame(predict(OL_pop, newdata = OL_preddata, interval = 'confidence'))
OL_preddata = data.frame(OL_preddata, fvfm_avg = OL_pred$Prediction, Lower = OL_pred$Lower, Upper = OL_pred$Upper)

FTS_preddata = data.frame(temp = seq(27,39, length.out = 400))
FTS_pred = as.data.frame(predict(FTS_pop, newdata = FTS_preddata, interval = 'confidence'))
FTS_preddata = data.frame(FTS_preddata, fvfm_avg =FTS_pred$Prediction, Lower = FTS_pred$Lower, Upper = FTS_pred$Upper)

RAS_preddata = data.frame(temp = seq(27,39, length.out = 400))
RAS_pred = as.data.frame(predict(RAS_pop, newdata = RAS_preddata, interval = 'confidence'))
RAS_preddata = data.frame(RAS_preddata, fvfm_avg = RAS_pred$Prediction, Lower = RAS_pred$Lower, Upper = RAS_pred$Upper)


combo_plot <- ggplot() +
  geom_jitter(data = fvfm, aes(x = temp, y = fvfm_avg, color = treatment), size = 1.5, width = 0.25) +
  scale_x_continuous(limits=c(25,40), breaks=c(27, 30, 33, 36, 39)) +
  scale_y_continuous(limits=c(-0.05, 0.85), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = OL_preddata, aes(x = temp, y = fvfm_avg), color = 'darkred', show.legend = FALSE) +
  geom_ribbon(data = OL_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkred', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = OL_coeff_mean), color = 'darkred', show.legend = FALSE, size =0.75) +
  
  geom_line(data = FTS_preddata, aes(x = temp, y = fvfm_avg), color = 'steelblue', show.legend = FALSE) +
  geom_ribbon(data = FTS_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'steelblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = FTS_coeff_mean), color = 'steelblue', show.legend = FALSE, size =0.75) +
  
  geom_line(data = RAS_preddata, aes(x = temp, y = fvfm_avg), color = 'coral', show.legend = FALSE) +
  geom_ribbon(data = RAS_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'coral', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = RAS_coeff_mean), color = 'coral', show.legend = FALSE, size =0.75) +

  ggtitle("") +
  scale_color_manual(values=c("darkred", "steelblue", "coral")) +
  ylab("Maximum quantum yield (Fv/Fm)") +
  xlab("Temperature (°C)") +
  theme_classic()+

  theme(text=element_text(size=14,  family="sans"))
combo_plot

####4) Analysis####

#compare ED50s (ANOVA approach)

ED50_lm <- lm(ED50 ~ Treatment, data = ED50s_long)

#check assumptions
plot(fitted(ED50_lm),residuals(ED50_lm))
hist(residuals(ED50_lm))
qqnorm(residuals(ED50_lm))

#posthoc 
ED50.posthoc1 <- emmeans(ED50_lm, pairwise ~ Treatment, type="response")
summary(ED50.posthoc1)

# contrast  estimate    SE df t.ratio p.value
# OL - FTS     2.904 0.445  6   6.531  0.0015 **
# OL - RAS     2.785 0.445  6   6.263  0.0019 **
# FTS - RAS   -0.119 0.445  6  -0.268  0.9614


####_________####

####Hu####
####1) Libraries####
library(readr)
library(tidyverse)
library(drc)
library(lmerTest)
library(emmeans)
library(Rmisc)
library(reshape2)
library(car)
library(patchwork)



####2) Data####
#Part II (measuring each species separately)
fvfm_Hu <- read_csv("C:/Users/lamkw/Downloads/3. R test/1. TT/Hu ED50/fvfm_data20250817_Hu.csv")
summary(fvfm_Hu)

fvfm_Hu$Period <- as.factor(fvfm_Hu$Period)
fvfm_Hu$species <-as.factor (fvfm_Hu$species)
fvfm_Hu$temp <-as.numeric(fvfm_Hu$temp)
fvfm_Hu$colony_id <-as.factor (fvfm_Hu$colony_id)
fvfm_Hu$nubbin_id <-as.factor (fvfm_Hu$nubbin_id)

summary(fvfm_Hu)
head(fvfm_Hu)

#Filter out "CBASS"
fvfm_Hu <-
  fvfm_Hu%>%
  filter(Period =="CBASS")

####3) Plots####

#quick look at data (post-experiment; March 31st)
fvfm_Hu_boxplot_post <-
  ggplot(fvfm_Hu, aes(x=species, y=post_fvfm_Hu, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  labs(title="B. CBASS Fv/Fm: post-experiment")+
  theme_classic()

POST <- fvfm_Hu_boxplot_post +  theme(text=element_text(size=16,  family="sans"))

fvfm_Hu_boxplot_post + facet_wrap(~site)

fvfm_Hu_COMBO <- PRE / POST 
fvfm_Hu_COMBO 

#Use dose curve approach; code used is based on approach described by Evensen et al. 2022
#Evensen NR, VoTurbstra CR, Fine M, Perna G, Buitrago-Lopez C, Cardenas A, Banc-Prandi G, Rowe K, Barshis DJ.
#(2022) Empirically derived thermal threshTurbds of four coral species along the Red Sea using a portable and 
#standardized experimental approach. Coral Reefs 41(2): 239-52. [https://doi.org/10.1007/s00338-022-02233-y]


#species subsets

#Turb
Turb_fvfm_Hu <-
  fvfm_Hu %>%
  filter(species == "Turb")

Turb_fvfm_Hu


#Ech
fvfm_Hu_Ech <-
  fvfm_Hu %>%
  filter(species == "Ech")

fvfm_Hu_Ech


#Por
fvfm_Hu_Por <-
  fvfm_Hu %>%
  filter(species == "Por")
fvfm_Hu_Por



####*Turb plots####

#Run population-level fit
Turb_pop <- drm(fvfm_avg ~ temp, data = Turb_fvfm_Hu, fct = LL.3())
summary(Turb_pop)
plot(Turb_pop)

#Run individual colony fits
#using "shared-slope logistic model" to stabilize ED₅₀ estimation.
#
Turb_indiv <- drm(
  fvfm_avg ~ temp,
  data = Turb_fvfm_Hu,
  curveid = colony_id,
  fct = LL.3(names = c("b","d","e")),
  pmodels = list(
    ~ 1,               # b : common (intercept-only)
    ~ colony_id,       # d : varies by colony
    ~ colony_id        # e : varies by colony (ED50 per colony)
  ))
summary(Turb_indiv)
plot(Turb_indiv)

#____________________________________________________________________
#Problem with this model:
Turb_indiv0 <- drm(fvfm_avg ~ temp, data = Turb_fvfm_Hu, curveid=colony_id, fct = LL.3())
summary(Turb_indiv0)
#警告訊息：於 sqrt(diag(varMat))： 產生了 NaNs
#This is due to sharp decline from 36C to 39C
#For all colonies (TR1–TR5), Fv/Fm stays flat (~0.6) from 30–36 °C, then drops 
#very sharply between 36–39 °C -> meaning the model sees almost no “lower tail”.
#So to fit a 3-parameter log-logistic (LL.3), drc has to make the transition extremely 
#steep to connect the upper plateau (~0.6) and the final point (~0.1)
#-> that in turn causes the variance–covariance matrix to become unstable, giving the error
plot(Turb_indiv0)
Turb_genocoeffs_50<-data.frame(ED(Turb_indiv0, c(50)))
#_____________________________________________________________________


#extract coeffs by colony, then compute 95% CIs
Turb_genocoeffs_50<-data.frame(ED(Turb_indiv, c(50)))
Turb_coeff_mean<-mean(Turb_genocoeffs_50$Estimate)
Turb_coeff_mean
#37.78

Turb_summary<-data.frame(CI(Turb_genocoeffs_50$Estimate, ci=0.95))
Turb_coeff_lower<-Turb_summary[3,]
Turb_coeff_upper<-Turb_summary[1,]



####*Ech plots####

#Run population-level fit
TxEch_pop <- drm(fvfm_avg ~ temp, data = fvfm_Hu_Ech, fct = LL.3())
summary(TxEch_pop)
plot(TxEch_pop)

#Run individual cTurbony fits
TxEch_indiv <- drm(fvfm_avg ~ temp, data = fvfm_Hu_Ech, curveid=colony_id, fct = LL.3())
summary(TxEch_indiv)
plot(TxEch_indiv)

#extract coeffs by cTurbony, then compute 95% CIs
TxEch_genocoeffs_50<-data.frame(ED(TxEch_indiv, c(50)))
TxEch_coeff_mean<-mean(TxEch_genocoeffs_50$Estimate)
TxEch_coeff_mean
#37.12

TxEch_summary<-data.frame(CI(TxEch_genocoeffs_50$Estimate, ci=0.95))
TxEch_coeff_lower<-TxEch_summary[3,]
TxEch_coeff_upper<-TxEch_summary[1,]


####*Por plots####

#Run population-level fit
TxPor_pop <- drm(fvfm_avg ~ temp, data = fvfm_Hu_Por, fct = LL.3())
summary(TxPor_pop)
plot(TxPor_pop)

#Run individual cTurbony fits
TxPor_indiv <- drm(fvfm_avg ~ temp, data = fvfm_Hu_Por, curveid=colony_id, fct = LL.3())
summary(TxPor_indiv)
plot(TxPor_indiv)

#extract coeffs by cTurbony, then compute 95% CIs
TxPor_genocoeffs_50<-data.frame(ED(TxPor_indiv, c(50)))
TxPor_coeff_mean<-mean(TxPor_genocoeffs_50$Estimate)
TxPor_coeff_mean
#36.46

TxPor_summary<-data.frame(CI(TxPor_genocoeffs_50$Estimate, ci=0.95))
TxPor_coeff_lower<-TxPor_summary[3,]
TxPor_coeff_upper<-TxPor_summary[1,]



####*Calculate ED50s####
#combine cTurbony-ED50s into data frame for statistical analysis

ED50s <- data.frame(cbind(Turb_genocoeffs_50[,1],TxEch_genocoeffs_50[,1],TxPor_genocoeffs_50[,1]))


ED50s <- ED50s %>% 
  dplyr::rename(Turb= X1,
                TxEch=X2,
                TxPor=X3,
  )


ED50s$colony_id<-as.factor(1:nrow(ED50s))
str(ED50s)

library(reshape2)
ED50s_long<-melt(ED50s, id="colony_id")

ED50s_long<-ED50s_long %>% 
  dplyr::rename(species= variable,
                ED50=value)


####*ED50 plot####

PA_coeff_means<-data.frame(Turb_coeff_mean, TxEch_coeff_mean, TxPor_coeff_mean)
PA_coeff_lowers<-data.frame(Turb_coeff_lower,  TxEch_coeff_lower, TxPor_coeff_lower)
PA_coeff_uppers<-data.frame(Turb_coeff_upper,  TxEch_coeff_upper, TxPor_coeff_upper)

Turb_preddata = data.frame(temp = seq(27,39, length.out = 400))
Turb_pred = as.data.frame(predict(Turb_pop, newdata = Turb_preddata, interval = 'confidence'))
Turb_preddata = data.frame(Turb_preddata, fvfm_avg = Turb_pred$Prediction, Lower = Turb_pred$Lower, Upper = Turb_pred$Upper)

TxEch_preddata = data.frame(temp = seq(27,39, length.out = 400))
TxEch_pred = as.data.frame(predict(TxEch_pop, newdata = TxEch_preddata, interval = 'confidence'))
TxEch_preddata = data.frame(TxEch_preddata, fvfm_avg =TxEch_pred$Prediction, Lower = TxEch_pred$Lower, Upper = TxEch_pred$Upper)

TxPor_preddata = data.frame(temp = seq(27,39, length.out = 400))
TxPor_pred = as.data.frame(predict(TxPor_pop, newdata = TxPor_preddata, interval = 'confidence'))
TxPor_preddata = data.frame(TxPor_preddata, fvfm_avg = TxPor_pred$Prediction, Lower = TxPor_pred$Lower, Upper = TxPor_pred$Upper)


combo_plot <- ggplot() +
  geom_jitter(data = fvfm_Hu, aes(x = temp, y = fvfm_avg, color = species), size = 1.5, width = 0.25) +
  scale_x_continuous(limits=c(25,40), breaks=c(27, 30, 33, 36, 39)) +
  scale_y_continuous(limits=c(-0.01, 0.85), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = Turb_preddata, aes(x = temp, y = fvfm_avg), color = 'darkred', show.legend = FALSE) +
  geom_ribbon(data = Turb_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkred', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = Turb_coeff_mean), color = 'darkred', show.legend = FALSE, size =0.75) +
  
  geom_line(data = TxEch_preddata, aes(x = temp, y = fvfm_avg), color = 'steelblue', show.legend = FALSE) +
  geom_ribbon(data = TxEch_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'steelblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = TxEch_coeff_mean), color = 'steelblue', show.legend = FALSE, size =0.75) +
  
  geom_line(data = TxPor_preddata, aes(x = temp, y = fvfm_avg), color = 'coral', show.legend = FALSE) +
  geom_ribbon(data = TxPor_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'coral', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = TxPor_coeff_mean), color = 'coral', show.legend = FALSE, size =0.75) +
  
  ggtitle("") +
  scale_color_manual(values=c("steelblue", "coral", "darkred")) +
  ylab("Maximum quantum yield (Fv/Fm)") +
  xlab("Temperature (°C)") +
  theme_classic()+
  
  theme(text=element_text(size=14,  family="sans")) +
  
  geom_text(data = PA_coeff_means, aes(label=round(Turb_coeff_mean, digits = 2)), x = 31, y = 0.2, show.legend = FALSE, color = 'darkred') +
  geom_text(data = PA_coeff_means, aes(label=round(TxEch_coeff_mean, digits = 2)), x = 31, y = 0.3, show.legend = FALSE, color = 'steelblue') +
  geom_text(data = PA_coeff_means, aes(label=round(TxPor_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'coral') 

combo_plot

####4) Analysis####

#compare ED50s (ANOVA approach)

ED50_lm <- lm(ED50 ~ species, data = ED50s_long)


#check assumptions
plot(fitted(ED50_lm),residuals(ED50_lm))
hist(residuals(ED50_lm))
qqnorm(residuals(ED50_lm))

#posthoc 
ED50.posthoc1 <- emmeans(ED50_lm, pairwise ~ species, type="response")
summary(ED50.posthoc1)
# contrast      estimate    SE df t.ratio p.value
# Turb - TxPor     1.318 0.441 12   2.987  0.0284*

#_________________________________________________________
##Using original code for analysis
#spoiler: same result.
Turb_genocoeffs_50_0<-data.frame(ED(Turb_indiv0, c(50)))

#ED50
ED50s_0 <- data.frame(cbind(Turb_genocoeffs_50_0[,1],TxEch_genocoeffs_50[,1],TxPor_genocoeffs_50[,1]))
ED50s_0 <- ED50s_0 %>% 
  dplyr::rename(Turb= X1,
                TxEch=X2,
                TxPor=X3,
  )


ED50s_0$colony_id<-as.factor(1:nrow(ED50s_0))

ED50s_long_0<-melt(ED50s_0, id="colony_id")

ED50s_long_0<-ED50s_long_0 %>% 
  dplyr::rename(species= variable,
                ED50=value)

#Stats
ED50_lm0 <- lm(ED50 ~ species, data = ED50s_long_0)

#check assumptions
plot(fitted(ED50_lm0),residuals(ED50_lm0))
hist(residuals(ED50_lm0))
qqnorm(residuals(ED50_lm0))
summary(ED50_lm0)
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   38.0930     0.3280 116.152  < 2e-16 ***
#speciesTxEch  -0.9693     0.4638  -2.090  0.05859 .  
#speciesTxPor  -1.6353     0.4638  -3.526  0.00418 ** 


#double check with ANOVA, type = III
library(car)
Anova(ED50_lm0, type = "III")
#            Sum Sq Df    F value  Pr(>F)    
#(Intercept) 7255.4  1 13491.2373 < 2e-16 ***
#species        6.8  2     6.2868 0.01356 *  

#posthoc 
ED50.posthoc0 <- emmeans(ED50_lm0, pairwise ~ species, type="response")
summary(ED50.posthoc0)

