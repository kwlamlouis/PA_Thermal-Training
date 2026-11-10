Sys.setlocale("LC_TIME", "C")
#Set working directory
setwd("C:/Users/lamkw/Downloads/1. 珊瑚實驗室/5. Thermal Training")

#------------------------------------------------------------------------------#
####[Library]####
library(ggplot2)
library(ggh4x)
library(ggpattern)
library(ggridges)
library(lubridate)
library(ggpubr) #Line plot
library(extrafont)
windowsFonts(Arial = windowsFont("Arial"))

library(nlme) #lme
library(glmmTMB) #nb
library(DHARMa)  
library(lmerTest) #lmer
library(emmeans)
library(car)  # For Levene's Test, "Anova"

library(survival)
library(survminer)

#[CBASS x ED50]
library(reshape2) #"melt"
library(drc) #for "drm"
library(Rmisc) #for "CI"

library(tidyverse)
library(dplyr)

#------------------------------------------------------------------------------#
#_____________________________________________________________________
####[Reproductive Output]#### 
RO <- read_csv("1. RO.csv")

RO$Date <- as.Date(RO$Date , format='%m/%d/%Y')
RO$Time <- factor(RO$Time, levels=c("Pre", "Post"))
RO$Lunar_Day <-as.numeric(RO$Lunar_Day)
RO$Source <- as.factor(RO$Source)
RO$Cult_Tx <-as.factor(RO$Cult_Tx)
RO$Tank <-as.factor(RO$Tank)
RO$larvae  <-as.integer(RO$larvae)
RO_no.NA<- na.omit(RO)


#Organize data: 
#Sum of RO of the same nubbins in the same tank
#Ex. Tank 1: 2.1_1, 2.1_2, 2.1_3 -> 2.1 sum of RO
RO_sum <- 
  RO_no.NA %>%
  group_by(Time, Cult_Tx, Tank, Source, Nubbin.ID) %>%
  summarise(larvae = sum(larvae),  .groups = "drop")


####1. Stats####
RO_GLMM <- glmmTMB(larvae ~ Cult_Tx*Time +(1 | Nubbin.ID) + (1|Tank),
                   contrasts = list(Cult_Tx = "contr.sum", Time = "contr.sum"),
                   data = RO_sum, 
                   family = nbinom2)

RO_GLMM.simulationOutput <- simulateResiduals(fittedModel = RO_GLMM, plot = TRUE)
testDispersion(RO_GLMM)

summary(RO_GLMM)
# Groups    Name        Variance  Std.Dev. 
# Nubbin.ID (Intercept) 1.459e+00 1.208e+00
# Tank      (Intercept) 1.543e-12 1.242e-06
# Number of obs: 201, groups:  Nubbin.ID, 108; Tank, 6

Anova(RO_GLMM, type = "III") 
#             Chisq Df Pr(>Chisq)
#(Intercept)      0  1     0.9969
#Cult_Tx          0  1     0.9975
#Time             0  1     0.9973
#Cult_Tx:Time     0  1     0.9981


##Post hoc
RO.posthoc <- emmeans(RO_GLMM, pairwise ~ Cult_Tx|Time)
#Time = Pre:
# contrast              estimate      SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29     2.57   0.411 Inf   6.272 <0.0001 ***

#Time = Post:
# contrast              estimate      SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29    12.07 220.000 Inf   0.055  0.9562

RO.posthoc2 <- emmeans(RO_GLMM, pairwise ~ Time|Cult_Tx)
#Cult_Tx = 26:
#Cult_Tx = 26:
# contrast   estimate     SE  df z.ratio p.value
# Pre - Post     3.71   0.538 Inf   6.889 <0.0001 ***

#Cult_Tx = 29:
# contrast   estimate     SE  df z.ratio p.value
# Pre - Post    13.20 220.000 Inf   0.060  0.9520


####1.1 PRE only####
RO_PRE.sum <-
  RO_sum %>%
  filter(Time != "Post") %>%
  droplevels()


RO_PRE_GLMM <- glmmTMB(larvae ~ Cult_Tx + (1|Tank),
                       data = RO_PRE.sum, 
                       family = nbinom2)

RO_PRE_GLMM.simulationOutput <- simulateResiduals(fittedModel = RO_PRE_GLMM, plot = TRUE)
testDispersion(RO_PRE_GLMM.simulationOutput)
RO_PRE_GLMM$fit$convergence  # Normally 0
RO_PRE_GLMM$sdr$pdHess      # Should be TRUE

summary(RO_PRE_GLMM)
#Random effects:
#Conditional model:
# Groups Name        Variance  Std.Dev.
# Tank   (Intercept) 2.079e-09 4.56e-05
#Number of obs: 108, groups:  Tank, 6
#Dispersion parameter for nbinom2 family (): 0.644 

#Conditional model:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   1.2095     0.1852   6.531 6.55e-11 ***
#Cult_Tx29    -2.4904     0.3602  -6.914 4.71e-12 ***


####1.2 26C Tx only####
RO_26C.sum <-
  RO_sum %>%
  filter(Cult_Tx == "26") %>%
  droplevels()


RO_26C_GLMM <- glmmTMB(larvae ~ Time + (1|Tank) + (1|Nubbin.ID),
                       data = RO_26C.sum, 
                       family = nbinom2)

RO_26C_GLMM.simulationOutput <- simulateResiduals(fittedModel = RO_26C_GLMM, plot = TRUE)
testDispersion(RO_26C_GLMM.simulationOutput)
RO_26C_GLMM$fit$convergence  # Normally 0
RO_26C_GLMM$sdr$pdHess      # Should be TRUE

summary(RO_26C_GLMM)
#Conditional model:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   0.8380     0.2782   3.013  0.00259 ** 
#TimePost     -3.6230     0.5452  -6.646 3.01e-11 ***


####2. Plot####
####2.1 RO box plot####
RO_sum2 <- 
  RO_sum %>%
  mutate(Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

RO_boxplot <- 
  ggplot(data=RO_sum2, aes(x = Cult_Tx, y=larvae, fill=Cult_Tx)) + 
  geom_boxplot(width=.35,alpha=.25,outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values=c("26°C" ="steelblue", "29°C" ="red"))+
  labs(y = "Larvae count", x="Thermal treatment")+
  ylim(0,50)+
  facet_grid(~Time)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
RO_boxplot
#Total reproductive output of nine nubbins originating from the same parental colony 
#across the three cultur tanks and acclimated to either 26°C (blue) or 29°C (orange) 
#for two years before (PRE) and two months after the heat stress experiment (POST). 
#Each point represents the total number of larvae produced by the replicated nubbins.
#Boxes indicate the median and interquartile range.
#_____________________________________________________#

##[Supplementary Figures]##
####2.2 Individual nubbins' RO & RT####
RO_no.NA2 <- 
  RO_no.NA %>%
  mutate(Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

# Plot with independent nested facets
Unified_Plot <- 
  ggplot(data = RO_no.NA2, aes(x = Date, y = larvae)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(y="Larvae count")+
  
  # Uses ggh4x to drop empty panels within rows/columns dynamically
  facet_nested(
    Cult_Tx + Tank ~ Time + Source, 
    scales = "free_x", 
    nest_line = TRUE # Adds a visual indicator grouping the sub-headers
  ) +
  
  theme_light() +
  theme(
    text = element_text(size = 11, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(color = "black", size = 9),
    ggh4x.facet.nestline = element_line(color = "grey50") # Style for nesting line
  )
Unified_Plot



####2.3 Proportional RO - colony contributions####
RO_colony_Proportion <-
  RO_sum %>%
  group_by(Time, Cult_Tx) %>%
  mutate(Release_Rate = sum(larvae > 0) / n(), Proportion= larvae/sum(larvae)* Release_Rate,
         Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

RO.Proportion_Plot <-
  ggplot(RO_colony_Proportion, aes(Cult_Tx, Proportion, fill=Source))+
  geom_col()+
  facet_wrap(~Time)+
  scale_y_continuous(labels=scales::percent, limits = c(0, 1))+
  theme_classic()+
  labs(y = "Contribution to total reproductive output", x = "", fill = "Nubbin ID") +
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
RO.Proportion_Plot


####[10% Tissue loss>]####
####1. Stats####
# Read your data
Tissue_loss <- read_csv("3. TT_ Tissue_loss.csv")
head(TT_Tissue_loss)

#For each fragment, get its event time and status
#Organize data to "until the event (tissue loss) was observed"
tt_summary <- Tissue_loss %>%
  group_by(Colony.ID, Nubbin.ID, Cult_Tx) %>%
  summarize(
    time   = if (any(Tissue_loss == 1)) {                #If the coral ever had tissue loss, 
      min(postTTday[Tissue_loss == 1])                  #we take the earliest day (min(postTTday)) it happened
    } else {
      max(postTTday)                                   # or else, we censored (i.e., not observing the event) at the last day of the observation period (i.e., we record the last day the sample was observed)
    },
    status = as.integer(any(Tissue_loss == 1))         
  ) %>%
  ungroup()

tt_summary %>% count(Cult_Tx)

# Check your condensed data has exactly one row per ID
stopifnot(nrow(tt_summary) == length(unique(tt_summary$Nubbin.ID)))

# 3. Create the survival object and fit
#"Kaplan-Meier survival analysis was performed to compare tissue loss between culture temperatures. The log-rank test indicated no significant difference (Chi2 = 1.27, p = 0.26)."
Tissue_obj <- Surv(time = tt_summary$time, event = tt_summary$status)
Tissue_km_fit   <- survfit(Tissue_obj ~ Cult_Tx, data = tt_summary)
surv_diff <- survdiff(Surv(time, status) ~ Cult_Tx, data = tt_summary)
# Chisq= 1.3  on 1 degrees of freedom, p= 0.3


####2. Plot####
my_theme <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"))

TissueLoss.Plot <- ggsurvplot(
  Tissue_km_fit,
  data = tt_summary,
  risk.table = FALSE,  pval = TRUE,  conf.int = TRUE,
  ggtheme = my_theme,
  legend = "right",
  xlab = "Days of heat stress exposure",
  ylab = "Cumulative probability of maintaining < 10% tissue loss",
  legend.title = "Thermal treatment",
  legend.labs = c("26℃", "29℃"),
  palette = c("#0072B2", "#FF0000")
)

#The survival curve will show the proportion of coral fragments without tissue loss over time.
#If the curve for 29C remains higher for longer, it supports your hypothesis that long-term exposure to 29?XC increased thermal tolerance.
#The p-value (log-rank test) tells you whether the difference between groups is statistically significant.



####[Survival]####
####1. Stats####
#Each nubbin (n=108) only has one record in the datasheet (i.e., it will always be 108 columns)
#Because we're recording when ("Date") the "event" occurred (either dead or alive ).
TT_SR <- read_csv("4. Updated_TT_SR.csv")
head(TT_SR)

TT_SR$Cult_Tx <-as.factor(TT_SR$Cult_Tx)
TT_SR$Colony.ID <-as.factor(TT_SR$Colony.ID)
TT_SR$Nubbin.ID <-as.factor(TT_SR$Nubbin.ID)
TT_SR$Tank <-as.factor(TT_SR$Tank)
TT_SR$Mortality <- as.integer(TT_SR$Mortality )

# 3. Create the survival object and fit
surv_obj <- Surv(time = TT_SR$Post_day , event = TT_SR$Mortality)
surv_km_fit   <- survfit(surv_obj ~ Cult_Tx, data = TT_SR)

#log-rank test: another stats to double check
survdiff(surv_obj ~ Cult_Tx, data = TT_SR)
#            N Observed Expected        (O-E)^2/E   (O-E)^2/V
#Cult_Tx=26 54        9        7.97          0.132     0.267
#Cult_Tx=29 54        7        8.03          0.132     0.267
#Chisq= 0.3  on 1 degrees of freedom, p= 0.6 


####2. Plot####
my_theme <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"))

SR.Plot <- ggsurvplot(
  surv_km_fit,
  data = TT_SR,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = TRUE,
  ggtheme = my_theme,
  legend = "right",
  xlab = "Days of post-heat stress recovery",
  ylab = "Cumulative survival probability",
  legend.title = "Thermal Treatment",
  legend.labs = c("26℃", "29℃"),
  palette = c("#0072B2", "#FF0000")
)



####[FvFm]####
TT_FvFm <- read_csv("2. TT_FvFm.csv")

TT_FvFm$Session <- as.factor(TT_FvFm$Session)
TT_FvFm$Date <-as.Date(TT_FvFm$Date, format ='%m/%d/%Y')
TT_FvFm$Timeline <- factor(TT_FvFm$Timeline, levels=c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13"))
TT_FvFm$Cult_Tx <-as.factor(TT_FvFm$Cult_Tx)
TT_FvFm$Colony.ID <-as.factor(TT_FvFm$Colony.ID)
TT_FvFm$Nubbin.ID <-as.factor(TT_FvFm$Nubbin.ID)
TT_FvFm$FvFm.reps <- as.factor(TT_FvFm$FvFm.reps)
TT_FvFm$Tank <-as.factor(TT_FvFm$Tank)
TT_FvFm$FvFm <-as.numeric(TT_FvFm$FvFm)

# 1. Remove missing Fv/Fm values
FvFm_clean <- TT_FvFm %>%
  filter(!is.na(FvFm)) %>%
  droplevels()

# Check if there are any "na" values (case-insensitive) in the FvFm column
any(is.na(FvFm_clean$FvFm))

# 2. Calculate mean Fv/Fm from technical replicates
FvFm_summary <- FvFm_clean %>%
  group_by(
    Session, Timeline, Cult_Tx,
    Colony.ID, Tank, Nubbin.ID
  ) %>%
  summarise(
    Mean_FvFm = mean(FvFm),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  droplevels()

# Separate the phases
FvFm_PRE <- FvFm_summary %>%
  filter(Session == "Pre_stress") %>%
  droplevels()

FvFm_PHS <- FvFm_summary %>%
  filter(Session == "HS",
         Timeline != "T10") %>%
  droplevels()

FvFm_POST <- FvFm_summary %>%
  filter(Session == "Recovery") %>%
  droplevels()

#### 1. Stats####
#### 1.1 PRE ####
# Model 1: Equal residual variance
PRE_equal <- glmmTMB(
  Mean_FvFm ~ Cult_Tx +
    (1 | Colony.ID) +
    (1 | Tank),
  data = FvFm_PRE,
  family = stats::gaussian()
)

# Model 2: Different residual variance by treatment
PRE_unequal <- glmmTMB(
  Mean_FvFm ~ Cult_Tx +
    (1 | Colony.ID) +
    (1 | Tank),
  dispformula = ~ Cult_Tx,
  data = FvFm_PRE,
  family = stats::gaussian()
)

AIC(PRE_equal, PRE_unequal)
#            df       AIC
#PRE_equal    5 -561.2122
#PRE_unequal  6 -582.2899 <- better

simulationOutput.PRE_unequal <- simulateResiduals(fittedModel = PRE_unequal, plot = T)
plotResiduals(PRE_unequal,  FvFm_PRE$Cult_Tx)

summary(PRE_unequal)
#Random effects:
#Conditional model:
# Groups    Name        Variance  Std.Dev.
# Colony.ID (Intercept) 1.389e-04 0.011787  --> meaning including "Colony.ID" is important
# Tank      (Intercept) 1.509e-05 0.003885
# Residual                     NA       NA
#Number of obs: 108, groups:  Colony.ID, 12; Tank, 6

#Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.719963   0.005458  131.90  < 2e-16 ***
#Cult_Tx29   -0.027290   0.008064   -3.38 0.000714 ***

#Dispersion model:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -4.6759     0.1040  -44.95  < 2e-16 ***
#Cult_Tx29     0.7391     0.1489    4.96 6.87e-07 ***
# => This tests whether the model's remaining variation differs between treatments.
#  The positive coefficient indicates greater residual variability in the 29 °C group.


#### 1.2 PHS ####
## Use "Tank-level summaries" ##
PHS_tank <- FvFm_PHS %>%
  group_by(Cult_Tx, Tank, Timeline) %>%
  summarise(
    Mean_FvFm = mean(Mean_FvFm),
    n_nubbins = n(),
    .groups = "drop"
  )

# Select T1 to T6
PHS_tank_early <- droplevels(
  subset(PHS_tank, Timeline %in% paste0("T", 1:6))
)

PHS_m12_early <- lmer(
  Mean_FvFm ~ Cult_Tx * Timeline +
    (1 | Tank),
  data = PHS_tank_early,
  REML = TRUE
)

simulationOutput.PHS_m12.early <- simulateResiduals(fittedModel = PHS_m12_early, plot = T)
plotResiduals(simulationOutput.PHS_m12.early,  form = PHS_tank_early$Cult_Tx)
plotResiduals(simulationOutput.PHS_m12.early,  form = PHS_tank_early$Timeline)
plotResiduals(simulationOutput.PHS_m12.early,  form = PHS_tank_early$Tank)
#no good

# Check whether the tank-specific residual patterns remain after accounting for the differences already estimated by (1 | Tank).
# Simulate residuals conditional on the fitted tank effects
simulationOutput.PHS_m12_early <- simulateResiduals(PHS_m12_early,
  n = 1000,
  re.form = NULL,
  plot = TRUE
)
#Good

VarCorr(PHS_m12_early)
# Groups   Name        Std.Dev.
# Tank     (Intercept) 0.025951
# Residual             0.011869

summary(PHS_m12_early)
#Random effects:
# Groups   Name        Variance  Std.Dev.
# Tank     (Intercept) 0.0006735 0.02595 
# Residual             0.0001409 0.01187 

anova(PHS_m12_early,
      type = 3,
      ddf = "Satterthwaite"
)
#Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#Cult_Tx          0.00005344 5.3444e-05     1     4  0.3794 0.57127  
#Timeline         0.00155287 3.1058e-04     5    20  2.2045 0.09427 .
#Cult_Tx:Timeline 0.00047164 9.4328e-05     5    20  0.6696 0.65105 


#### 1.3 POST####
# Remove nubbin.ID because it has nearly zero estimated differences among individual nubbins
POST_glm2 <- glmmTMB(
  Mean_FvFm ~ Cult_Tx * Timeline +
    (1 | Colony.ID) +
    (1 | Tank),
  family = stats::gaussian(),
  dispformula = ~ Timeline,
  contrasts = list(Cult_Tx = "contr.sum", Timeline = "contr.sum"),
  data = FvFm_POST,
)

simulationOutput.POST.gml2 <- simulateResiduals(fittedModel = POST_glm2, plot = T)
plotResiduals(simulationOutput.POST.gml2,  form = FvFm_POST$Cult_Tx)
plotResiduals(simulationOutput.POST.gml2,  form = FvFm_POST$Timeline)
plotResiduals(simulationOutput.POST.gml2,  form = FvFm_POST$Tank)

Anova(POST_glm2, type = "III")
#                     Chisq Df Pr(>Chisq)    
#Cult_Tx             0.9974  1  0.3179368    
#Timeline           74.9614  2  < 2.2e-16 ***
#Cult_Tx:Timeline   14.9949  2  0.0005545 ***

POST_glm2.posthoc <- emmeans(POST_glm2, pairwise ~ Cult_Tx|Timeline)
#Timeline = T11:
# contrast              estimate     SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29  -0.0249 0.0195 Inf  -1.280  0.2006

#Timeline = T12:
# contrast              estimate     SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   0.0336 0.0176 Inf   1.906  0.0566

#Timeline = T13:
# contrast              estimate     SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   0.0424 0.0253 Inf   1.677  0.0936

POST_glm2.posthoc <- emmeans(POST_glm2, pairwise ~ Timeline|Cult_Tx)
#Cult_Tx = 26:
# contrast  estimate     SE  df z.ratio p.value
# T11 - T12  -0.0333 0.0115 Inf  -2.903  0.0103 *
# T11 - T13   0.0590 0.0173 Inf   3.406  0.0019 **
# T12 - T13   0.0923 0.0163 Inf   5.669 <0.0001 ***

#Cult_Tx = 29:
# contrast  estimate     SE  df z.ratio p.value
# T11 - T12   0.0252 0.0114 Inf   2.212  0.0692
# T11 - T13   0.1263 0.0170 Inf   7.428 <0.0001 ***
# T12 - T13   0.1011 0.0160 Inf   6.319 <0.0001 ***

####2. Plot####
FvFm_clean <-
FvFm_clean%>%
  mutate(Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

TT_FvFm_Plot<- 
  ggline(FvFm_clean, x = "Timeline", y = "FvFm", color = "Cult_Tx",
         add = c("mean_se")) +
  scale_color_manual(values=c("26°C" ="steelblue", "29°C" ="red"))+
  labs(y = "Fv/Fm", x = "Timepoint", color = "Thermal treatment") +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Arial")) +
  theme(legend.title = element_text(size=10))+
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) 
TT_FvFm_Plot

####[Red index]####
RGB <- read_csv("7. RGB.csv")

RGB$Batch<-as.factor(RGB$Batch)
RGB$Timeline <- as.factor(RGB$Timeline)
RGB$Check_date <- as.Date(RGB$Check_date, format='%m/%d/%Y')
RGB$Post_day <-as.factor(RGB$Post_day)
RGB$Tank<- as.factor(RGB$Tank)
RGB$Cult_Tx<- as.factor(RGB$Cult_Tx)
RGB$Colony.ID<- as.factor(RGB$Colony.ID)
RGB$Nubbin.ID<- as.factor(RGB$Nubbin.ID)
head(RGB)


####1. Stats####
##Method 1: Lmer: using "Post_day" as factor
RGB$Post_day <-as.factor(RGB$Post_day)

Red_stats_allTHREE.lmer <- lmer(Coral_D1_D5_D1 ~ Cult_Tx*Post_day + (1|Colony.ID)+ (1|Nubbin.ID) + (1|Tank), 
                                data = RGB)

Red_stats_allTHREE.lmer.simulationOutput <- simulateResiduals(fittedModel = Red_stats_allTHREE.lmer, plot = TRUE)
testDispersion(Red_stats_allTHREE.lmer)

anova(Red_stats_allTHREE.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                  Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#Cult_Tx          0.00185 0.00185     1   4.317  0.0798  0.790553    
#Post_day         0.92929 0.46465     2 184.779 20.0174 1.353e-08 ***
#Cult_Tx:Post_day 0.31999 0.15999     2 184.779  6.8927  0.001297 ** 


##Posthoc
Red_stats_allTHREE.lmer_posthoc1 <- emmeans(Red_stats_allTHREE.lmer, pairwise~Cult_Tx|Post_day)
#Post_day = 0:
# contrast              estimate    SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29   0.1302 0.177 4.47   0.734  0.4996

#Post_day = 7:
# contrast              estimate    SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29   0.0499 0.178 4.50   0.281  0.7915

#Post_day = 14:
# contrast              estimate    SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29  -0.0310 0.178 4.51  -0.174  0.8691

Red_stats_allTHREE.lmer_posthoc2 <- emmeans(Red_stats_allTHREE.lmer, pairwise~Post_day|Cult_Tx)
#Cult_Tx = 26:
#  contrast               estimate     SE  df t.ratio p.value
# Post_day0 - Post_day7    0.1046 0.0307 198   3.401  0.0023 **
# Post_day7 - Post_day14   0.1135 0.0315 192   3.605  0.0012 **

#Cult_Tx = 29:
# contrast               estimate     SE  df t.ratio p.value
# Post_day0 - Post_day7    0.0242 0.0301 194   0.803  0.7016
# Post_day7 - Post_day14   0.0326 0.0310 193   1.051  0.5455


####2. Plot####
#Data: mean Standard_Red/Tx
RGB_summary <-
  RGB%>%
  group_by(Timeline, Post_day, Cult_Tx)%>%
  summarise(mean.std_red=mean(Coral_D1_D5_D1), sd.std_red=sd(Coral_D1_D5_D1), n=n(), SE = sd.std_red/sqrt(n))

RGB_summary2<-
  RGB_summary%>%
  mutate(Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

Stand_Red.plot <-
  ggplot(RGB_summary2, aes(x=as.factor(RGB_summary2$Post_day), y=mean.std_red, colour=Cult_Tx))+
  geom_line(aes(group=Cult_Tx))+
  geom_point()+
  geom_errorbar(aes(ymin=mean.std_red-SE, ymax=mean.std_red+SE), width=.2)+
  scale_color_manual(values=c("steelblue", "red"))+
  theme_classic()+
  labs(x="Days since first valid RGB measurement", y="Standardized red index value", colour="Thermal treatment")+
  ylim(0, 1.5)+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Stand_Red.plot



####[Calcification rate]####
TT_DryWeight <- read_csv("5. TT_DryWeight.csv")

TT_DryWeight$Section <- as.factor(TT_DryWeight$Section)
TT_DryWeight$Timeline <- as.factor(TT_DryWeight$Timeline)
TT_DryWeight$Cult_Tx <- as.factor(TT_DryWeight$Cult_Tx)
TT_DryWeight$Tank <- as.factor(TT_DryWeight$Tank)
TT_DryWeight$Nubbin.ID <- as.factor(TT_DryWeight$Nubbin.ID)
TT_DryWeight$Date <-as.Date(TT_DryWeight$Date, format ='%m/%d/%Y')
TT_DryWeight$Growth <-as.numeric(TT_DryWeight$Growth)
TT_DryWeight$Dry_weight <- as.numeric(TT_DryWeight$Dry_weight)
head(TT_DryWeight)

nrow(TT_DryWeight)

####1. Stats####
# Combine the t-distribution with the dispersion correction
# "t_family" is specifically designed to handle those "heavy tails" (outliers/extremes) that are triggering the KS test error, while keeping negatives intact.
# "dispformula = ~Timeline": as the experiment progressed, the variance in growth increased, which is a biological finding in itself
Cal_analysis <- TT_DryWeight %>%
  filter(Timeline != "Tm4")

GR_glm_allTHREE <- glmmTMB(Growth ~ Cult_Tx * Timeline + (1|Parental.ID)+ (1|Nubbin.ID) + (1|Tank), 
                           family = t_family(),           # Handles the heavy tails/outliers
                           dispformula = ~ Timeline,      # To account for heteroscedasticity, the dispersion parameter was modeled as a function of the experimental timepoint
                           contrasts = list(Cult_Tx = "contr.sum", Timeline = "contr.sum"),
                           data = Cal_analysis)

GR_glm_allTHREE.lmer.simulationOutput <- simulateResiduals(fittedModel = GR_glm_allTHREE, plot = TRUE)
testDispersion(GR_glm_allTHREE.lmer.simulationOutput)

GR_glm_allTHREE_anova <- Anova(GR_glm_allTHREE, type = "III", test.statistic = "Chisq")
#                   Chisq Df Pr(>Chisq) 
#Cult_Tx           1.0251  1     0.3113    
#Timeline         95.1793  4  < 2.2e-16 ***
#Cult_Tx:Timeline 33.6462  4  8.806e-07 ***

GR_glm_allTHREE_posthoc <- emmeans(GR_glm_allTHREE, pairwise~Cult_Tx|Timeline)
#Timeline = Tm1:
# contrast              estimate    SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   -0.406 0.251 Inf  -1.620  0.1053

#Timeline = Tm2:
# contrast              estimate    SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29    0.272 0.258 Inf   1.056  0.2909

#Timeline = Tm3:
# contrast              estimate    SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   -0.304 0.325 Inf  -0.937  0.3488

#Timeline = Tm5:
# contrast              estimate    SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   -0.329 0.260 Inf  -1.268  0.2049

#Timeline = Tm6:
# contrast              estimate    SE  df z.ratio p.value
# Cult_Tx26 - Cult_Tx29   -0.477 0.281 Inf  -1.695  0.0900

GR_glm_allTHREE_posthoc2 <- emmeans(GR_glm_allTHREE, pairwise~Timeline|Cult_Tx)
#Cult_Tx = 26:
#Tm1 - Tm2   0.0211 0.0955 Inf   0.222  0.9995
#Tm2 - Tm3   0.0487 0.2140 Inf   0.227  0.9994
#Tm5 - Tm6  -0.4422 0.1270 Inf  -3.476  0.0046 **

#Cult_Tx = 29:
#Tm1 - Tm2   0.6995 0.0840 Inf   8.325 <0.0001 ***
#Tm2 - Tm3  -0.5279 0.1260 Inf  -4.199  0.0003 ***
#Tm5 - Tm6  -0.5899 0.1280 Inf  -4.625 <0.0001 ***

####1. Plot####
#Note: Remove "T0" because "growth rate" = T1 - T0 ..etc; hence, "GR inputs" start from "T1"
TT_DryWeight2<-
  TT_DryWeight%>%
  filter(!is.na(Growth) & Timeline!="T0")%>%
  mutate(Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")))

GR.plot<-
  ggline(TT_DryWeight2, x = "Date", y = "Growth", color = "Cult_Tx",
         add = c("mean_se"))+
  labs(y = "Calcification rate (mg/g/day)", x = "", color = "Thermal treatment") +
  scale_color_manual(values=c("26°C" ="steelblue", "29°C" ="red"))+
  theme_classic()+
  theme(text=element_text(size=12,  family="Arial"), legend.position ="right")+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
GR.plot



####[Symbiont density]####
Symb <- read_csv("8. TT symb count.csv")

Symb$Section<-as.factor(Symb$Section)
Symb$Tank<- as.factor(Symb$Tank)
Symb$Cult_Tx<- as.factor(Symb$Cult_Tx)
Symb$Colony.ID<- as.factor(Symb$Colony.ID)
Symb$Nubbin.ID<- as.factor(Symb$Nubbin.ID)
Symb$`Ten^5_cell_cm2` <- as.numeric(Symb$`Ten^5_cell_cm2`)
head(Symb)

# 1. Remove missing Fv/Fm values
Symb_clean <- Symb %>%
  filter(!is.na(`Ten^5_cell_cm2`)) %>%
  droplevels()

# Check if there are any "na" values (case-insensitive) 
any(is.na(Symb_clean$`Ten^5_cell_cm2`))

####1. Stats####
# Calculate the average for each fragment (ID)
Symb_no.NA.avg <- 
  Symb_clean %>%
  group_by(Section, Cult_Tx, Tank, Colony.ID, Nubbin.ID)%>% 
  summarise(mean_density = mean(`Ten^5_cell_cm2`, na.rm = TRUE), .groups = 'drop')

Sym.lmer_Tank.ID <- lmer(
  mean_density ~ Section* Cult_Tx + (1|Tank), 
  data = Symb_no.NA.avg)

Sym.lmer_Tank.ID.simulationOutput <- simulateResiduals(fittedModel = Sym.lmer_Tank.ID, plot = T)
testDispersion(Sym.lmer_Tank.ID.simulationOutput)

anova(Sym.lmer_Tank.ID)
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
#Section         395.50  395.50     1 188.084  3.2081 0.07488 
#Cult_Tx          75.69   75.69     1   4.018  0.6140 0.47690  
#Section:Cult_Tx 571.43  571.43     1 188.084  4.6352 0.03260 *

#Post hoc
Sym.lmer_Tank.ID_posthoc1 <- emmeans(Sym.lmer_Tank.ID, pairwise~ Section | Cult_Tx)
#Cult_Tx = 26:
# contrast                 estimate   SE  df t.ratio p.value
# End of PHS - End of Post    0.576 2.25 188   0.257  0.7978

#Cult_Tx = 29:
# contrast                 estimate   SE  df t.ratio p.value
# End of PHS - End of Post   -6.283 2.26 188  -2.781  0.0060 **

Sym.lmer_Tank.ID_posthoc2 <- emmeans(Sym.lmer_Tank.ID, pairwise~ Cult_Tx | Section)
#Section = End of PHS:
# contrast              estimate   SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29     8.62 6.78 4.41   1.271  0.2668

#Section = End of Post:
# contrast              estimate   SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29     1.76 6.84 4.55   0.257  0.8082


####2. Plot####
Symb_no.NA.ave_PlotUse <- 
  Symb_no.NA.avg %>%
  mutate(
    Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")),
    Section = factor(Section, levels = c("End of PHS", "End of Post"))
  )

Symb_boxplot <- 
  ggplot(data=Symb_no.NA.ave_PlotUse, aes(x = Cult_Tx, y=mean_density, fill=Cult_Tx)) + 
  geom_boxplot(width=.35,alpha=.25,outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values=c("26°C" ="steelblue", "29°C" ="red"))+
  scale_y_continuous(limits=c(0,100), breaks = c(0,20, 40, 60, 80, 100))+
  labs(x = "", y = expression(paste("Algal symbiont count (", 10^5, " cells/cm"^2, ")")))+
  facet_grid(~Section)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Symb_boxplot



####[Chl a]####
Chla <- read_csv("9. TT_Chla.csv")

Chla$Section<-as.factor(Chla$Section)
Chla$Tank<- as.factor(Chla$Tank)
Chla$Cult_Tx<- as.factor(Chla$Cult_Tx)
Chla$Colony.ID<- as.factor(Chla$Colony.ID)
Chla$Nubbin.ID<- as.factor(Chla$Nubbin.ID)
Chla$Chla <- as.numeric(Chla$Chla)

# 1. Remove missing Fv/Fm values
Chla_clean <- Chla %>%
  filter(!is.na(Chla)) %>%
  droplevels()

# Check if there are any "na" values (case-insensitive) in the FvFm column
any(is.na(Chla_clean$Chla))


####1. Stats####
Chla_noN.lmer <- lmer(
  Chla ~ Section * Cult_Tx + (1 | Nubbin.ID) + (1|Tank), 
  data = Chla_clean
)

#assumption check
testDispersion(Chla_noN.lmer)
simulateResiduals(fittedModel = Chla_noN.lmer, plot = T)
VarCorr(Chla_noN.lmer)

anova(Chla_noN.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Section         1.1607  1.1607     1 86.639 17.4708 6.947e-05 ***
#Cult_Tx         0.0038  0.0038     1  3.940  0.0575    0.8224    
#Section:Cult_Tx 3.4358  3.4358     1 86.639 51.7174 2.145e-10 ***

#Posthoc
Chla_posthoc1 <- emmeans(Chla_noN.lmer, pairwise~Section|Cult_Tx)
Chla_posthoc1
#Cult_Tx = 26:
# contrast                 estimate     SE  df t.ratio p.value
# End of PHS - End of Post    0.424 0.0533 101   7.955 <0.0001 ***

#Cult_Tx = 29:
# contrast                 estimate     SE  df t.ratio p.value
# End of PHS - End of Post   -0.112 0.0524  97  -2.147  0.0343 *

Chla_posthoc2 <- emmeans(Chla_noN.lmer, pairwise~Cult_Tx|Section)
Chla_posthoc2
#Section = End of PHS:
# contrast              estimate     SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29    0.284 0.0734 6.51   3.867  0.0071 **

#Section = End of Post:
# contrast              estimate     SE   df t.ratio p.value
# Cult_Tx26 - Cult_Tx29   -0.253 0.0766 7.63  -3.302  0.0116 *


####2. Plot####
Chla_no.NA_PlotUse <- 
  Chla_no.NA %>%
  mutate(
    Cult_Tx = factor(paste0(Cult_Tx, "°C"), levels = c("26°C", "29°C")),
    Section = factor(Section, levels = c("End of PHS", "End of Post"))
  )

Chla_no.NA_boxplot <- 
  ggplot(data=Chla_no.NA_PlotUse, aes(x = Cult_Tx, y=Chla, fill=Cult_Tx)) + 
  geom_boxplot(width=.35,alpha=.25,outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values=c("26°C" ="steelblue", "29°C" ="red"))+
  labs(x = "", y = expression(paste("Chlorophyll a (", mu, "g/cm"^2, ")")))+
  facet_grid(~Section)+
  scale_y_continuous(limits=c(0,3.0), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3.0))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text=element_text(size=12,  family="Arial"))+
  theme(axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin = margin(t=10)))
Chla_no.NA_boxplot



####[CBASS X ED50]####
fvfm <- read_csv("TT_ED50.csv")

fvfm$Mother <-as.factor(fvfm$Mother)
fvfm$colony_id <-as.factor (fvfm$colony_id)
fvfm$Age <-as.factor (fvfm$Age)
fvfm$Tx <-as.factor (fvfm$Tx)
fvfm$Combo <-as.factor (fvfm$Combo)
fvfm$Tank <- as.factor(fvfm$Tank)

summary(fvfm)
head(fvfm)

####1. Plots####
#Combo subsets
#Adult_26
A.26_fvfm <-
  fvfm %>%
  filter(Age =="Adult", Tx =="26")
A.26_fvfm

#Adult_29
A.29_fvfm <-
  fvfm %>%
  filter(Age =="Adult", Tx =="29")
A.29_fvfm

#Recruit_26
R.26_fvfm <-
  fvfm %>%
  filter(Age =="Recruit", Tx =="26")
R.26_fvfm

#Recruit_29
R.29_fvfm <-
  fvfm %>%
  filter(Age =="Recruit", Tx =="29")
R.29_fvfm

####*Adult 26 plots####
#Run population-level fit
A.26_pop <- drm(fvfm_mean ~ temp, data = A.26_fvfm, fct = LL.3())
summary(A.26_pop)
plot(A.26_pop)

#Run individual colony fits
A.26_indiv <- drm(fvfm_mean ~ temp, data = A.26_fvfm, curveid=Mother, fct = LL.3())
summary(A.26_indiv)
plot(A.26_indiv)

#extract coeffs by colony, then compute 95% CIs
A.26_genocoeffs_50<-data.frame(ED(A.26_indiv, c(50)))
A.26_coeff_fvfm_mean<-mean(A.26_genocoeffs_50$Estimate)
A.26_coeff_fvfm_mean
#36.58099

A.26_summary<-data.frame(CI(A.26_genocoeffs_50$Estimate, ci=0.95))
A.26_coeff_lower<-A.26_summary[3,]
A.26_coeff_upper<-A.26_summary[1,]


####*Adult 29 plots####
#Run population-level fit
A.29_pop <- drm(fvfm_mean ~ temp, data = A.29_fvfm, fct = LL.3())
summary(A.29_pop)
plot(A.29_pop)

#Run individual colony fits
A.29_indiv <- drm(fvfm_mean ~ temp, data = A.29_fvfm, curveid=Mother, fct = LL.3())
summary(A.29_indiv)
plot(A.29_indiv)

#extract coeffs by colony, then compute 95% CIs
A.29_genocoeffs_50<-data.frame(ED(A.29_indiv, c(50)))
A.29_coeff_fvfm_mean<-mean(A.29_genocoeffs_50$Estimate)
A.29_coeff_fvfm_mean
#37.28684

A.29_summary<-data.frame(CI(A.29_genocoeffs_50$Estimate, ci=0.95))
A.29_coeff_lower<-A.29_summary[3,]
A.29_coeff_upper<-A.29_summary[1,]


####*Recruit 26 plots####
#Run population-level fit
R.26_pop <- drm(fvfm_mean ~ temp, data = R.26_fvfm, fct = LL.3())
summary(R.26_pop)
plot(R.26_pop)

#Run individual colony fits
R.26_indiv <- drm(fvfm_mean ~ temp, data = R.26_fvfm, curveid=Mother, fct = LL.3())
summary(R.26_indiv)
plot(R.26_indiv)

#extract coeffs by colony, then compute 95% CIs
R.26_genocoeffs_50<-data.frame(ED(R.26_indiv, c(50)))
R.26_coeff_fvfm_mean<-mean(R.26_genocoeffs_50$Estimate)
R.26_coeff_fvfm_mean
#36.65604

R.26_summary<-data.frame(CI(R.26_genocoeffs_50$Estimate, ci=0.95))
R.26_coeff_lower<-R.26_summary[3,]
R.26_coeff_upper<-R.26_summary[1,]


####*Recruit 29 plots####
#Run population-level fit
R.29_pop <- drm(fvfm_mean ~ temp, data = R.29_fvfm, fct = LL.3())
summary(R.29_pop)
plot(R.29_pop)

#Run individual colony fits
R.29_indiv <- drm(fvfm_mean ~ temp, data = R.29_fvfm, curveid=Mother, fct = LL.3())
summary(R.29_indiv)
plot(R.29_indiv)

#extract coeffs by colony, then compute 95% CIs
R.29_genocoeffs_50<-data.frame(ED(R.29_indiv, c(50)))
R.29_coeff_fvfm_mean<-mean(R.29_genocoeffs_50$Estimate)
R.29_coeff_fvfm_mean
#37.56425

R.29_summary<-data.frame(CI(R.29_genocoeffs_50$Estimate, ci=0.95))
R.29_coeff_lower<-R.29_summary[3,]
R.29_coeff_upper<-R.29_summary[1,]


####*Calculate ED50s####
#combine colony-ED50s into data frame for statistical analysis

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
  dplyr::rename(Combo= variable,
                ED50=value)


# Extract separate Age and Tx columns from the Combo column strings (e.g., "A.26")
ED50s_long <- ED50s_long %>%
  mutate(
    Age = case_when(
      grepl("^A", Combo) ~ "Adult",
      grepl("^R", Combo) ~ "Recruit"
    ),
    Tx = case_when(
      grepl("26$", Combo) ~ "26",
      grepl("29$", Combo) ~ "29"
    )
  )

# Ensure the new variables are recognized as factors
ED50s_long$Age <- as.factor(ED50s_long$Age)
ED50s_long$Tx <- as.factor(ED50s_long$Tx)


####*ED50 plot####

PA_coeff_fvfm_means<-data.frame(A.26_coeff_fvfm_mean, A.29_coeff_fvfm_mean, R.26_coeff_fvfm_mean, R.29_coeff_fvfm_mean)
PA_coeff_lowers<-data.frame(A.26_coeff_lower,  A.29_coeff_lower, R.26_coeff_lower, R.29_coeff_lower)
PA_coeff_uppers<-data.frame(A.26_coeff_upper,  A.29_coeff_upper, R.26_coeff_lower, R.29_coeff_upper)

A.26_preddata = data.frame(temp = seq(30,39, length.out = 400))
A.26_pred = as.data.frame(predict(A.26_pop, newdata = A.26_preddata, interval = 'confidence'))
A.26_preddata = data.frame(A.26_preddata, fvfm_mean = A.26_pred$Prediction, Lower = A.26_pred$Lower, Upper = A.26_pred$Upper)

A.29_preddata = data.frame(temp = seq(30,39, length.out = 400))
A.29_pred = as.data.frame(predict(A.29_pop, newdata = A.29_preddata, interval = 'confidence'))
A.29_preddata = data.frame(A.29_preddata, fvfm_mean =A.29_pred$Prediction, Lower = A.29_pred$Lower, Upper = A.29_pred$Upper)

R.26_preddata = data.frame(temp = seq(30,39, length.out = 400))
R.26_pred = as.data.frame(predict(R.26_pop, newdata = R.26_preddata, interval = 'confidence'))
R.26_preddata = data.frame(R.26_preddata, fvfm_mean = R.26_pred$Prediction, Lower = R.26_pred$Lower, Upper = R.26_pred$Upper)

R.29_preddata = data.frame(temp = seq(30,39, length.out = 400))
R.29_pred = as.data.frame(predict(R.29_pop, newdata = R.29_preddata, interval = 'confidence'))
R.29_preddata = data.frame(R.29_preddata, fvfm_mean =R.29_pred$Prediction, Lower = R.29_pred$Lower, Upper = R.29_pred$Upper)

combo_plot <- ggplot() +
  geom_jitter(data = fvfm, aes(x = temp, y = fvfm_mean, color = Combo), size = 1.5, width = 0.25) +
  scale_x_continuous(limits=c(30,40), breaks=c(30, 33, 36, 39)) +
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = A.26_preddata, aes(x = temp, y = fvfm_mean), color = 'steelblue', show.legend = FALSE) +
  geom_ribbon(data = A.26_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'steelblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_fvfm_means, aes(xintercept = A.26_coeff_fvfm_mean), color = 'steelblue', show.legend = FALSE, linewidth =0.75) +
  
  geom_line(data = A.29_preddata, aes(x = temp, y = fvfm_mean), color = 'red2', show.legend = FALSE) +
  geom_ribbon(data = A.29_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red2', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_fvfm_means, aes(xintercept = A.29_coeff_fvfm_mean), color = 'red2', show.legend = FALSE, linewidth =0.75) +
  
  geom_line(data = R.26_preddata, aes(x = temp, y = fvfm_mean), color = 'lightblue', show.legend = FALSE) +
  geom_ribbon(data = R.26_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'lightblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_fvfm_means, aes(xintercept = R.26_coeff_fvfm_mean), color = 'lightblue', show.legend = FALSE, linewidth =0.75) +
  
  geom_line(data = R.29_preddata, aes(x = temp, y = fvfm_mean), color = 'violet', show.legend = FALSE) +
  geom_ribbon(data = R.29_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'violet', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_fvfm_means, aes(xintercept = R.29_coeff_fvfm_mean), color = 'violet', show.legend = FALSE, linewidth =0.75) +
  
  ggtitle("") +
  scale_color_manual(values=c("steelblue", "red2", "lightblue", "violet")) +
  labs(y = "Fv/Fm", x = "Temperature (°C)", color = "Life-stage x Thermal treatment") +
  theme_classic()+
  theme(text = element_text(size = 12, family = "Arial"))
combo_plot


####2. Analysis####
#compare ED50s (ANOVA approach)

ED50_lm <- lm(ED50 ~ Age*Tx, 
              contrasts = list(Age = "contr.sum", Tx = "contr.sum"),
              data = ED50s_long)
Anova(ED50_lm, type = "III")
#            Sum Sq Df    F value  Pr(>F)    
#Age             0.2  1 9.6900e-01 0.3395777    
#Tx              3.3  1 2.0323e+01 0.0003575 ***
#Age:Tx          0.1  1 3.1950e-01 0.5797536 
#Residuals      2.6 16 

#check assumptions
simulateResiduals(fittedModel = ED50_lm, plot = T)
testDispersion(ED50_lm)

#posthoc 
ED50.posthoc1 <- emmeans(ED50_lm, pairwise ~ Age*Tx, type="response")
summary(ED50.posthoc1) 
# contrast                    estimate    SE df t.ratio p.value
# Adult Tx26 - Recruit Tx26    -0.0751 0.253 16  -0.296  0.9906
# Adult Tx26 - Adult Tx29      -0.7059 0.253 16  -2.788  0.0575
# Adult Tx26 - Recruit Tx29    -0.9833 0.253 16  -3.884  0.0065 **
# Recruit Tx26 - Adult Tx29    -0.6308 0.253 16  -2.492  0.0995
# Recruit Tx26 - Recruit Tx29  -0.9082 0.253 16  -3.587  0.0119 *
# Adult Tx29 - Recruit Tx29    -0.2774 0.253 16  -1.096  0.6971


## [Analysis 2] 30C vs 33C vs 36C vs 39C
fvfm_summary$temp <- as.factor(fvfm_summary$temp)
Dose.curve_lmer <- lmer(mean ~ Age*Tx*temp + (1|Mother) + (1|Tank), 
                        contrasts = list(Age = "contr.sum", Tx = "contr.sum"),
                        data = fvfm_summary)

testDispersion(Dose.curve_lmer)
Dose.curve_simulateResiduals <- simulateResiduals(fittedModel = Dose.curve_lmer, plot = T)
plotResiduals(Dose.curve_simulateResiduals,  form = fvfm_summary$Mother)
plotResiduals(Dose.curve_simulateResiduals,  form = fvfm_summary$Tank)

Anova(Dose.curve_lmer, type = "III")
#                Chisq Df Pr(>Chisq)    
#Age          172.540  1  < 2.2e-16 ***
#Tx            15.561  1  7.987e-05 ***
#temp        1615.154  3  < 2.2e-16 ***
#Age:Tx        26.107  1  3.231e-07 ***
#Age:temp      92.138  3  < 2.2e-16 ***
#Tx:temp       45.033  3  9.103e-10 ***
#Age:Tx:temp   20.582  3  0.0001286 ***


#posthoc 
Dose.curve_lmer.posthoc1 <- emmeans(Dose.curve_lmer, ~ Age * Tx * temp)
pairs(Dose.curve_lmer.posthoc1, by = c("Age", "Tx"), adjust = "tukey")
#Age = Adult, Tx = 26
# contrast          estimate     SE df t.ratio p.value
# temp30 - temp33   0.0206 0.0197 41.4   1.048  0.7226
# temp33 - temp36   0.2769 0.0197 41.4  14.082 <0.0001 ***
# temp36 - temp39   0.2200 0.0197 41.4  11.190 <0.0001 ***

#Age = Recruit, Tx = 26:
# contrast          estimate     SE df t.ratio p.value
# temp30 - temp33   0.0004 0.0197 41.4   0.020  1.0000
# temp33 - temp36   0.2443 0.0197 41.4  12.424 <0.0001 ***
# temp36 - temp39   0.1879 0.0197 41.4   9.559 <0.0001 ***

#Age = Adult, Tx = 29:
# contrast        estimate     SE df t.ratio p.value
# temp30 - temp33   0.0245 0.0197 41.4   1.248  0.6006
# temp33 - temp36   0.1923 0.0197 41.4   9.779 <0.0001 ***
# temp36 - temp39   0.3080 0.0197 41.4  15.665 <0.0001 ***

#Age = Recruit, Tx = 29:
# contrast        estimate     SE df t.ratio p.value
# temp30 - temp33  -0.0008 0.0197 41.4  -0.041  1.0000
# temp33 - temp36   0.1375 0.0197 41.4   6.995 <0.0001 ***
# temp36 - temp39   0.1670 0.0197 41.4   8.494 <0.0001 ***

pairs(Dose.curve_lmer.posthoc1, by = c("temp"), adjust = "tukey")
#temp = 30:
# contrast                    estimate     SE df t.ratio p.value
# Adult Tx26 - Recruit Tx26     0.1098 0.0193 83   5.675 <0.0001 ***
# Adult Tx26 - Adult Tx29      -0.0159 0.0193 83  -0.824  0.8432
# Adult Tx26 - Recruit Tx29     0.2337 0.0193 83  12.078 <0.0001 ***
# Recruit Tx26 - Adult Tx29    -0.1257 0.0193 83  -6.499 <0.0001 ***
# Recruit Tx26 - Recruit Tx29   0.1239 0.0193 83   6.402 <0.0001 ***
# Adult Tx29 - Recruit Tx29     0.2496 0.0193 83  12.901 <0.0001 ***

#temp = 33:
# contrast                    estimate     SE df t.ratio p.value
# Adult Tx26 - Recruit Tx26     0.0896 0.0193 83   4.631 <0.0001 ***
# Adult Tx26 - Adult Tx29      -0.0120 0.0193 83  -0.620  0.9253
# Adult Tx26 - Recruit Tx29     0.2123 0.0193 83  10.971 <0.0001 ***
# Recruit Tx26 - Adult Tx29    -0.1016 0.0193 83  -5.251 <0.0001 ***
# Recruit Tx26 - Recruit Tx29   0.1227 0.0193 83   6.340 <0.0001 ***
# Adult Tx29 - Recruit Tx29     0.2243 0.0193 83  11.592 <0.0001 ***

#temp = 36:
# contrast                    estimate     SE df t.ratio p.value
# Adult Tx26 - Recruit Tx26     0.0570 0.0193 83   2.946  0.0213 *
# Adult Tx26 - Adult Tx29      -0.0966 0.0193 83  -4.993 <0.0001 ***
# Adult Tx26 - Recruit Tx29     0.0729 0.0193 83   3.770  0.0017 **
# Recruit Tx26 - Adult Tx29    -0.1536 0.0193 83  -7.939 <0.0001 ***
# Recruit Tx26 - Recruit Tx29   0.0159 0.0193 83   0.824  0.8432
# Adult Tx29 - Recruit Tx29     0.1695 0.0193 83   8.763 <0.0001 ***

#temp = 39:
# contrast                    estimate     SE df t.ratio p.value
# Adult Tx26 - Recruit Tx26     0.0249 0.0193 83   1.289  0.5727
# Adult Tx26 - Adult Tx29      -0.0086 0.0193 83  -0.445  0.9705
# Adult Tx26 - Recruit Tx29     0.0199 0.0193 83   1.030  0.7322
# Recruit Tx26 - Adult Tx29    -0.0335 0.0193 83  -1.733  0.3132
# Recruit Tx26 - Recruit Tx29  -0.0050 0.0193 83  -0.258  0.9939
# Adult Tx29 - Recruit Tx29     0.0285 0.0193 83   1.475  0.4572


####____________________####
####[Save as HD]####
#Check current location
getwd()

png("GR.png", units="in", width=9, height=7, res=600)
GR.plot
dev.off()

#Save data
write.csv(as.data.frame(tt_summary), "Tissue.loss_Results.csv", row.names = TRUE)

  