

############################
#load packages
############################
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggstatsplot")
install.packages("dplyr")
install.packages("lme4")
install.packages("car")
install.packages("emmeans")
install.packages("nlme")
install.packages("marginaleffects")
install.packages("piecewiseSEM")
install.packages("rstantools")
install.packages("multcomp")
install.packages("treemapify")
install.packages("relaimpo")
install.packages("r2glmm")
install.packages("patchwork")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("gridExtra")
install.packages("MuMIn")
install.packages("boot")
install.packages("corrr") ## for PCA
install.packages("ggcorrplot") ## for PCA
install.packages("FactoMineR") ## for PCA
install.packages("factoextra")
install.packages("ggfortify")
install.packages("lavaan")
install.packages("patchwork")
install.packages("cowplot")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(nlme)
library(marginaleffects)
library(piecewiseSEM)
library(rstantools)
library(multcomp)
library(treemapify)
library(relaimpo)
library(r2glmm)
library(patchwork)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(MuMIn)
library(boot)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(lavaan)
library(patchwork)
library(cowplot)

################## load data ##########################################################################################

stat <- read.csv("./data/data_C3.csv")
stat
names(stat)
nrow(stat) # 5293
table(stat$year)

hist(sqrt(stat$ratio_CN_topsoil))
hist(log(stat$ratio_CN_topsoil))

hist(sqrt(stat$phosphorus))
hist(log(stat$phosphorus))

hist(sqrt(stat$aridity))
hist(log(stat$aridity))

stat$sqrtMI <- sqrt(stat$aridity) # to meet normality

stat$AI <- 1-stat$aridity ### aridity index
hist(stat$AI)

stat$N_fix <- ifelse(stat$N_fix == "Non_fixers", "OP", 
                             ifelse(stat$N_fix == "Fixers", "N2FP", stat$N_fix))

################## Subset the variables we need  ######################################
###################################################################################################
stat_subset <-  subset(stat, select = c(lat, lon, species,Genus,N_fix, myco_type, Economy, big_D13, AI,
                                        nitrogen_topsoil,soc_topsoil, ratio_CN_topsoil, ratio_NC_topsoil, 
                                        phosphorus, aridity, sqrtMI,ln_N_topsoil,
                                        lnP, ln_CN_topsoil, ln_NC_topsoil,chi,beta,lnbeta, pre_mean, tmp))

colSums(is.na(stat_subset)) 
stat_subset <-  na.omit(stat_subset)## remove NA
nrow(stat_subset) # 4953

###################lnbeta outliers following MAD method ##############
beta_mean <- mean(stat_subset$beta, na.rm = T)
beta_mean #424.05

lnbeta_median <- median(stat_subset$lnbeta, na.rm = T)
lnbeta_median #5.29
lnbeta_mean <- mean(stat_subset$lnbeta, na.rm = T)
lnbeta_mean #5.39
lnbeta_mad <- mad(stat_subset$lnbeta, na.rm = T) # median absolute deviation
lnbeta_mad # 0.94

stat_clean <- subset(stat_subset, lnbeta < lnbeta_median + 3 * lnbeta_mad &
                       lnbeta > lnbeta_median - 3 * lnbeta_mad)

nrow(stat_clean) # 4867
table(stat_clean$Economy)
length(table(stat_clean$species))


######################## Exclude NM #####################################################
stat_no_NM = subset (stat_clean, Economy!="NM")
table(stat_no_NM$Economy) ## 1140 mining, 3605 scavenging 
length(unique(stat_no_NM$Genus))
length(unique(stat_no_NM$species)) #2599
nrow(stat_no_NM) # 4745
names(stat_no_NM) # 4745


################### dividing by Tropical, Temperate, Boreal ####################################
bins <- c(-Inf,-30, 30,Inf)
labels <- c("Temperate", "Tropical", "Temperate")
# Combine the original data with the new latitude category variable
stat_no_NM$biome <- cut(stat_no_NM$lat, breaks = bins, labels = labels, include.lowest = TRUE)


#################### Linear mixed effects model for beta #######################################
########################## H1: beta decreases with P and N increase, and with aridity decrease
########################## H2: These relationships are affected by mycorrhizal association,  N fixation affects only N/C-beta relationship

lmer_beta <- lmer(lnbeta~ lnP + ln_NC_topsoil + Economy + sqrtMI + 
                          sqrtMI*Economy + ln_NC_topsoil*Economy + lnP*Economy 
                          + N_fix + N_fix*ln_NC_topsoil+
                          (1|Genus), data=stat_no_NM)
Anova(lmer_beta) 
vif(lmer_beta) ## High VIF caused by the interaction term between N_fix and ln_NC_topsoil with no interaction found 
summary(lmer_beta)

################# removing the interaction term N-fix * ln NC_topsoil 
lmer_beta_dropp <- lmer(lnbeta~ lnP + ln_NC_topsoil + Economy + sqrtMI + 
                    sqrtMI*Economy + ln_NC_topsoil*Economy + lnP*Economy + N_fix + 
                    (1|Genus), data=stat_no_NM)
Anova(lmer_beta_dropp) 
vif(lmer_beta_dropp) ## all VIF under 10
r.squaredGLMM(lmer_beta_dropp) ## R2 0.42
############# export result for table S1 #####################
Beta_model_S1 <- data.frame(Var = c('lnP','lnN/C', 'Myco_NAS', 'GAI', 'N fixation',
                                 'GAI X MAS', 'ln N/C X MAS', 'ln P X MAS'
                                 ))
Beta_model_S1$df <- as.matrix(Anova(lmer_beta_dropp))[1:8, 2]

Beta_model_S1$Slope <- c(
  summary(emtrends(lmer_beta_dropp, ~lnP, var = "lnP"))[1, 2],
  summary(emtrends(lmer_beta_dropp, ~ln_NC_topsoil, var = "ln_NC_topsoil"))[1, 2],
  NA,
  summary(emtrends(lmer_beta_dropp, ~sqrtMI, var = "sqrtMI"))[1, 2],
  NA, NA, NA,NA)

Beta_model_S1$SE <- c(
  summary(emtrends(lmer_beta_dropp, ~lnP, var = "lnP"))[1, 3],
  summary(emtrends(lmer_beta_dropp, ~ln_NC_topsoil, var = "ln_NC_topsoil"))[1, 3],
  NA,
  summary(emtrends(lmer_beta_dropp, ~sqrtMI, var = "sqrtMI"))[1, 3],
  NA, NA, NA,NA)

Beta_model_S1$p <- as.matrix(Anova(lmer_beta_dropp))[1:8, 3]

Beta_model_S1$VIF <- as.matrix(vif(lmer_beta_dropp))[1:8, 1]

Beta_model_S1
write.csv(Beta_model_S1, "./output/Beta_model_S1.csv")

################### Keep the original model ############################
lmer_beta <- lmer(lnbeta~ lnP + ln_NC_topsoil + Economy + sqrtMI + 
                    sqrtMI*Economy + ln_NC_topsoil*Economy + lnP*Economy + N_fix + N_fix*ln_NC_topsoil +
                    (1|Genus), data=stat_no_NM)
Anova(lmer_beta) 
vif(lmer_beta)
summary(lmer_beta)
AIC(lmer_beta) # 11881
r.squaredGLMM(lmer_beta) ## R2 0.42
shapiro.test(residuals(lmer_beta))

qqnorm(residuals(lmer_beta))
qqline(residuals(lmer_beta))

residuals <- resid(lmer_beta)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## good
plot(fitted(lmer_beta), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity :OK

densityPlot(residuals(lmer_beta))

options(lmerTest.pbkrtest.limit = 6000) ### when we have large dataset
emm_options(pbkrtest.limit = 6000) 

comp_slopes_P <- cld(emtrends(lmer_beta, ~ Economy, var = 'lnP')) ## comp slopes for lnP
comp_slopes_P ## no sig # between slopes
test(emtrends(lmer_beta, pairwise~Economy, "lnP")) 

comp_slopes_N <- cld(emtrends(lmer_beta, ~ Economy, var = 'ln_NC_topsoil')) ## comp slopes for N/C
comp_slopes_N ## mining different from scavenging 

test(emtrends(lmer_beta, pairwise~Economy, "ln_NC_topsoil")) # p-values slopes comparison
test(emtrends(lmer_beta, pairwise~N_fix, "ln_NC_topsoil")) # p-values slopes comparison

comp_slopes_aridity <- cld(emtrends(lmer_beta, ~ Economy, var = 'sqrtMI')) ## comp slopes for aridity
comp_slopes_aridity ## no sig # between slopes
test(emtrends(lmer_beta, pairwise~Economy, "sqrtMI")) # p-values slopes comparison


P_test <- test(emtrends(lmer_beta, ~ Economy, var = "lnP")) # test if slopes are different from 0
P_test

N_test <- test(emtrends(lmer_beta, ~ Economy, var = "ln_NC_topsoil")) # test if slopes are different from 0
N_test

aridity_test <- test(emtrends(lmer_beta, ~ Economy, var = "sqrtMI")) # test if slopes are different from 0
aridity_test

N_fix_test <- test(emtrends(lmer_beta, ~ N_fix, var = "ln_NC_topsoil")) # test if slopes are different from 0
N_fix_test


comp_slopes_Nfix <- cld(emtrends(lmer_beta, ~ N_fix, var = 'ln_NC_topsoil')) ## comp slopes for aridity
comp_slopes_Nfix ## no sig # between slopes

cld(emmeans(lmer_beta, ~Economy))

###################### Export tables from the lmer model #####################################
###########################################################################################
Beta_model <- data.frame(Var = c('lnP','lnN/C', 'Myco_NAS', 'GAI', 'N fixation',
                                 'GAI X MAS', 'ln N/C X MAS', 'ln P X MAS', 
                                 'ln N/C X N fixation'))
Beta_model$df <- as.matrix(Anova(lmer_beta))[1:9, 2]

Beta_model$Slope <- c(
  summary(emtrends(lmer_beta, ~lnP, var = "lnP"))[1, 2],
  summary(emtrends(lmer_beta, ~ln_NC_topsoil, var = "ln_NC_topsoil"))[1, 2],
  NA,
  summary(emtrends(lmer_beta, ~sqrtMI, var = "sqrtMI"))[1, 2],
  NA, NA, NA,NA,NA)

Beta_model$SE <- c(
  summary(emtrends(lmer_beta, ~lnP, var = "lnP"))[1, 3],
  summary(emtrends(lmer_beta, ~ln_NC_topsoil, var = "ln_NC_topsoil"))[1, 3],
  NA,
  summary(emtrends(lmer_beta, ~sqrtMI, var = "sqrtMI"))[1, 3],
  NA, NA, NA,NA,NA)

Beta_model$p <- as.matrix(Anova(lmer_beta))[1:9, 3]

Beta_model$VIF <- as.matrix(vif(lmer_beta))[1:9, 1]

Beta_model

write.csv(Beta_model, "./output/Beta_model.csv")

Beta_model_comp <- data.frame(Var = c('lnP mining', 'lnP scavenging', 
                                 'lnN/C mining', 'lnN/C scavenging',
                                 'GAI mining', 'GAI scavenging',
                                  'lnN/C Fixers', 'lnN/C Non Fixers'))
Beta_model_comp$Slope <- c(
                       summary(emtrends(lmer_beta, ~Economy, var = "lnP"))[1, 2],
                       summary(emtrends(lmer_beta, ~Economy, var = "lnP"))[2, 2],
                       summary(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[1, 2],
                       summary(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[2, 2],
                       summary(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[1, 2],
                       summary(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[2, 2],
                       summary(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[1, 2],
                       summary(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[2, 2]
                       )
Beta_model_comp$SE <- c(
                    summary(emtrends(lmer_beta, ~Economy, var = "lnP"))[1, 3],
                    summary(emtrends(lmer_beta, ~Economy, var = "lnP"))[2, 3],
                    summary(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[1, 3],
                    summary(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[2, 3],
                    summary(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[1, 3],
                    summary(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[2, 3],
                    summary(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[1, 3],
                    summary(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[2, 3]
                      )
Beta_model_comp$p <- c(
                test(emtrends(lmer_beta, ~Economy, var = "lnP"))[1, 6],
                test(emtrends(lmer_beta, ~Economy, var = "lnP"))[2, 6],
                test(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[1, 6],
                test(emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"))[2, 6],
                test(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[1, 6],
                test(emtrends(lmer_beta, ~Economy, var = "sqrtMI"))[2, 6],
                test(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[1, 6],
                test(emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"))[2, 6]
                )
Beta_model_comp$group <- c(
  cld(object = emtrends(lmer_beta, ~Economy, var = "lnP"),adjust = "Tukey",Letters = letters, alpha = 0.05)[2, 7],
  cld(object = emtrends(lmer_beta, ~Economy, var = "lnP"),adjust = "Tukey",Letters = letters, alpha = 0.05)[1, 7],
  cld(object = emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"),adjust = "Tukey",Letters = letters, alpha = 0.05)[2, 7],
  cld(object = emtrends(lmer_beta, ~Economy, var = "ln_NC_topsoil"),adjust = "Tukey",Letters = letters, alpha = 0.05)[1, 7],
  cld(object = emtrends(lmer_beta, ~Economy, var = "sqrtMI"),adjust = "Tukey",Letters = letters, alpha = 0.05)[1, 7],
  cld(object = emtrends(lmer_beta, ~Economy, var = "sqrtMI"),adjust = "Tukey",Letters = letters, alpha = 0.05)[2, 7],
  cld(object = emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"),adjust = "Tukey",Letters = letters, alpha = 0.05)[1, 7],
  cld(object = emtrends(lmer_beta, ~N_fix, var = "ln_NC_topsoil"),adjust = "Tukey",Letters = letters, alpha = 0.05)[2, 7]
)
Beta_model_comp

write.csv(Beta_model_comp, "./output/Beta_model_comp.csv")
############ Extract slopes for N, P, and aridity for Economy + N fix ########################################### 
Slope_N_mining <- N_test [1, 2] 
Slope_N_scavenging <- N_test [2, 2] 

Slope_P_mining <- P_test [1, 2] 
Slope_P_scavenging <- P_test [2, 2] 

Slope_AI_mining <- aridity_test [1, 2] 
Slope_AI_scavenging <- aridity_test [2, 2] 

Slope_N_fix <- N_fix_test [1, 2] 
Slope_N_Nonfix<- N_fix_test [2, 2] 

#################################################################################################

############ Extract Intercepts for N, P, and aridity ########################################### 
intercepts_N <- summary(emmeans(lmer_beta, ~Economy, var ='ln_NC_topsoil', at = list(ln_NC_topsoil = 0)))
intercepts_P <- summary(emmeans(lmer_beta, ~Economy, var ='lnP', at = list(lnP = 0)))
intercepts_AI <- summary(emmeans(lmer_beta, ~Economy, var ='sqrtMI', at = list(aridity = 0)))

int_N_mining <- intercepts_N [1, 2] 
int_N_scavenging <- intercepts_N [2, 2] 

int_P_mining <- intercepts_P [1, 2] 
int_P_scavenging <- intercepts_P [2, 2] 

int_AI_mining <- intercepts_AI [1, 2] 
int_AI_scavenging <- intercepts_AI [2, 2] 


intercepts_Nfix <- summary(emmeans(lmer_beta, ~N_fix, var ='ln_NC_topsoil', at = list(ln_NC_topsoil = 0)))
int_N_fix <- intercepts_Nfix [1, 2] 
int_N_Nonfix <- intercepts_Nfix [2, 2] 

############################### Plots without intervals of confidence ###################################

############## regression lines for NC ratio for each Economy strategy ###################################
mining = subset(stat_no_NM, Economy =="Mining")
N_seq_mining <- seq(min(mining$ln_NC_topsoil, na.rm = T), max(mining$ln_NC_topsoil, na.rm = T), 0.001)
N_seq_mining
N_trend_mining <- int_N_mining + N_seq_mining * Slope_N_mining
N_trend_mining <- as.data.frame(cbind(N_seq_mining, N_trend_mining))
N_trend_mining$Economy= "Mining"


Scavenging = subset(stat_no_NM, Economy =="Scavenging")
N_seq_scav<- seq(min(Scavenging$ln_NC_topsoil, na.rm = T), max(Scavenging$ln_NC_topsoil, na.rm = T), 0.001)
N_seq_scav
N_trend_scavenging <- int_N_scavenging + N_seq_scav * Slope_N_scavenging
N_trend_scavenging <- as.data.frame(cbind(N_seq_scav, N_trend_scavenging))
N_trend_scavenging$Economy= "Scavenging"

N_test
NC_plot_bis <- (NC_plot_bis <- ggplot(data = stat_no_NM, aes(x = ln_NC_topsoil, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen","Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging","Mining")) +
      
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 3) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 3) +
      
    geom_line(data = N_trend_scavenging, aes(x = N_seq_scav , y = N_trend_scavenging ), col = 'lightseagreen', lwd = 2, alpha = 0.8) +
    geom_line(data = N_trend_mining, aes(x = N_seq_mining, y = N_trend_mining), col = 'tan3', lwd = 2, alpha = 0.8) +
  
      theme(legend.position ="right",
            legend.justification = "top",
          axis.title.y = element_text(size = 50, colour = 'black'),
          axis.title.x = element_text(size = 50, colour = 'black'),
          axis.text.x = element_text(size = 50, colour = 'black'),
          axis.text.y = element_text(size = 50, colour = 'black'),
          axis.line = element_line(size = 2, colour = "black"), 
          axis.ticks = element_line(size = 1, colour = "black"), 
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          #panel.grid.major = element_line(colour = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size=40),
          legend.background = element_rect(fill = "white", color = "white", size = 0.8), 
          legend.key = element_rect(color = "white", fill = "white"))+
         # legend.position =c(0.99, 0.99),
          #legend.justification = c(1, 1)) +
    ylab(expression('ln ' * italic('β'))) +
    xlab(expression('ln ' * 'N/C ratio')))+
  guides(fill = guide_legend(
    title= NULL,
    title.position = "top",
    override.aes = list(size = 5),
    keywidth = unit(1, "cm"),
    keyheight = unit(1.5, "cm"),
    title.hjust = 0.5,
    nrow = 2,
    byrow = TRUE
  ))

NC_plot_bis

################################################################################################################

############## regression lines for P  for each Economy strategy ###################################
P_seq_mining <- seq(min(mining$lnP, na.rm = T), max(mining$lnP, na.rm = T), 0.001)
P_seq_mining
P_trend_mining <- int_P_mining + P_seq_mining * Slope_P_mining
P_trend_mining <- as.data.frame(cbind(P_seq_mining, P_trend_mining))
P_trend_mining$Economy= "Mining"


P_seq_scav<- seq(min(Scavenging$lnP, na.rm = T), max(Scavenging$lnP, na.rm = T), 0.001)
P_seq_scav
P_trend_scavenging <- int_P_scavenging + P_seq_scav * Slope_P_scavenging
P_trend_scavenging <- as.data.frame(cbind(P_seq_scav, P_trend_scavenging))
P_trend_scavenging$Economy= "Scavenging"


P_test
P_plot_bis <- (P_plot_bis <- ggplot(data = stat_no_NM, aes(x = lnP, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging", "Mining")) +
      #scale_shape_manual(values = c("Temperate" = 21, "Tropical" = 24)) + 
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 3) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 3) +
    
    geom_line(data = P_trend_scavenging, aes(x = P_seq_scav , y = P_trend_scavenging ), col = 'lightseagreen', lwd = 2, alpha = 1.5) +
    geom_line(data = P_trend_mining, aes(x = P_seq_mining, y = P_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +
    
    theme(legend.position ="right",
          legend.justification = "top",
          axis.title.y = element_text(size = 50, colour = 'black'),
          axis.title.x = element_text(size = 50, colour = 'black'),
          axis.text.x = element_text(size = 50, colour = 'black'),
          axis.text.y = element_text(size = 50, colour = 'black'),
          axis.line = element_line(size = 2, colour = "black"), 
          axis.ticks = element_line(size = 1, colour = "black"), 
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          #panel.grid.major = element_line(colour = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size=40),
          legend.background = element_rect(fill = "white", color = "white", size = 0.8), 
          legend.key = element_rect(color = "white", fill = "white")) +
          #legend.position =c(1.1, 1.1),
          #legend.justification = c(1, 1)) +
    ylab(expression('ln ' * italic('β'))) +
      xlab(expression('ln ' * P[i]))) +
  guides(fill = guide_legend(
    title= NULL,
    title.position = "top",
    override.aes = list(size = 5),
    keywidth = unit(1, "cm"),
    keyheight = unit(1.5, "cm"),
    title.hjust = 0.5,
    nrow = 2,
    byrow = TRUE
  ))

P_plot_bis

######################################################################################################
############## regression lines for aritidy  for each Economy strategy ###################################
mining = subset(stat_no_NM, Economy =="Mining")
AI_seq_mining <- seq(min(mining$sqrtMI, na.rm = T), max(mining$sqrtMI, na.rm = T), 0.001)
AI_seq_mining
AI_trend_mining <- int_AI_mining + AI_seq_mining * Slope_AI_mining
AI_trend_mining <- as.data.frame(cbind(AI_seq_mining, AI_trend_mining))
AI_trend_mining$Economy= "Mining"


Scavenging = subset(stat_no_NM, Economy =="Scavenging")
AI_seq_scav<- seq(min(Scavenging$sqrtMI, na.rm = T), max(Scavenging$sqrtMI, na.rm = T), 0.001)
AI_seq_scav
AI_trend_scavenging <- int_AI_scavenging + AI_seq_scav * Slope_AI_scavenging
AI_trend_scavenging <- as.data.frame(cbind(AI_seq_scav, AI_trend_scavenging))
AI_trend_scavenging$Economy= "Scavenging"

aridity_test
MI_plot_bis <- (MI_plot_bis <- ggplot(data = stat_no_NM, aes(x = sqrtMI, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging", "Mining")) +
      
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 3) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 3) +
    
    geom_line(data = AI_trend_scavenging, aes(x = AI_seq_scav , y = AI_trend_scavenging ), col = 'lightseagreen',  lwd = 2, alpha = 1.5) +
    geom_line(data = AI_trend_mining, aes(x = AI_seq_mining, y = AI_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +

    scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
    theme(legend.position ="none",
          #legend.justification = "top",
          axis.title.y = element_text(size = 50, colour = 'black'),
          axis.title.x = element_text(size = 50, colour = 'black'),
          axis.text.x = element_text(size = 50, colour = 'black'),
          axis.text.y = element_text(size = 50, colour = 'black'),
          axis.line = element_line(size = 2, colour = "black"), 
          axis.ticks = element_line(size = 1, colour = "black"), 
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          #panel.grid.major = element_line(colour = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size=40),
          legend.background = element_rect(fill = "white", color = "white", size = 0.8), 
          legend.key = element_rect(color = "white", fill = "white"))+
          #legend.position =c(0.98, 0.98),
          #legend.justification = c(1, 1)) +
    ylab(expression('ln ' * italic('β'))) +
      xlab(expression(bold("\u221A") * MI))+
  guides(fill = guide_legend(
    title= NULL,
    #title.position = "top",
    override.aes = list(size = 5),
    keywidth = unit(1, "cm"),
    keyheight = unit(1.5, "cm"),
    title.hjust = 0.5,
    nrow = 2,
    byrow = TRUE
  )))

MI_plot_bis

#####################################################################################################
############## regression lines for NC ratio for each N fixation strategies ###################################
stat_no_NM$N_fix
N2FP = subset(stat_no_NM, N_fix =="N2FP")
OP = subset(stat_no_NM, N_fix =="OP")

fixers_seq <- seq(min(N2FP$ln_NC_topsoil, na.rm = T), max(N2FP$ln_NC_topsoil, na.rm = T), 0.001)
fixers_seq
N_trend_fixers <- int_N_fix + fixers_seq * Slope_N_fix
N_trend_fixers <- as.data.frame(cbind(fixers_seq, N_trend_fixers))
N_trend_fixers$N_fix= "N2FP"

Nonfixers_seq <- seq(min(OP$ln_NC_topsoil, na.rm = T), max(OP$ln_NC_topsoil, na.rm = T), 0.001)
Nonfixers_seq
N_trend_Nonfixers <- int_N_Nonfix + Nonfixers_seq * Slope_N_Nonfix
N_trend_Nonfixers <- as.data.frame(cbind(Nonfixers_seq, N_trend_Nonfixers))
N_trend_Nonfixers$N_fix= "OP"


N_fix_test
Nfix_Plot_bis <- (Nfix_Plot_bis <- ggplot(data = stat_no_NM, aes(x = ln_NC_topsoil, y = lnbeta, fill=N_fix)) +
    
    scale_fill_manual(values = c("N2FP" = "red3", "OP" = "black")) +
      
      geom_jitter(data = subset(stat_no_NM, N_fix == "N2FP"), 
                  aes(fill = N_fix), shape = 24, color = "red3", 
                  alpha = 0.5, size = 3) +
      geom_jitter(data = subset(stat_no_NM, N_fix == "OP"), 
                  aes(fill = N_fix), shape = 21, color = "black", 
                  alpha = 0.5, size = 3) +
      
    geom_line(data = N_trend_fixers, aes(x = fixers_seq , y = N_trend_fixers ), col = 'red3', lwd = 2, alpha = 0.8) +
    geom_line(data = N_trend_Nonfixers, aes(x = Nonfixers_seq, y = N_trend_Nonfixers), col = 'black', lwd = 2, alpha = 0.8) +
   
      
    #scale_x_continuous(breaks = seq(-3.5, 0, by = 0.5)) +
      
    theme(legend.position ="right",
          legend.justification ="top",
          axis.title.y = element_text(size = 50, colour = 'black'),
          axis.title.x = element_text(size = 50, colour = 'black'),
          axis.text.x = element_text(size = 50, colour = 'black'),
          axis.text.y = element_text(size = 50, colour = 'black'),
          axis.line = element_line(size = 2, colour = "black"), 
          axis.ticks = element_line(size = 1, colour = "black"), 
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          #panel.grid.major = element_line(colour = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size=40),
          legend.background = element_rect(fill = "white", color = "white", size = 0.8), 
          legend.key = element_rect(color = "white", fill = "white"))+ 
          #legend.position =c(1, 1),
          #legend.justification = c(1, 1)) +
    ylab(expression('ln ' * italic('β'))) +
    xlab(expression('ln ' * 'N/C ratio')))+
guides(fill = guide_legend(
  title= NULL,
  title.position = "top",
  override.aes = list(size = 4),
  keywidth = unit(1, "cm"),
  keyheight = unit(1.5, "cm"),
  nrow = 2,
  byrow = TRUE))
 
Nfix_Plot_bis

###############SEM######################################################################
########################################################################################
table(stat_no_NM$N_fix)
########### Set categorical as levels   #########################################################
stat_no_NM$Nfix_level = 0 ## creation of a column with Ntrt= 0 
stat_no_NM$Nfix_level[stat_no_NM$N_fix == 'N2FP'] = 1 
stat_no_NM$Nfix_level[stat_no_NM$N_fix == 'OP'] = 2

stat_no_NM$Economy_level = 0 ## creation of a column with Ntrt= 0 
stat_no_NM$Economy_level[stat_no_NM$Economy == 'Mining'] = 1
stat_no_NM$Economy_level[stat_no_NM$Economy == 'Scavenging'] = 2

stat_no_NM$ln_N_topsoil = log(stat_no_NM$nitrogen_topsoil)
stat_no_NM$ln_C_topsoil = log(stat_no_NM$soc_topsoil)

############# Subset #########################
Temperate = subset(stat_no_NM,biome == "Temperate")
Tropical = subset(stat_no_NM,biome == "Tropical")

nrow(Temperate) # 2481
table(Temperate$Economy) # 872 mining, 1609 scavenging 

nrow(Tropical) # 2264
table(Tropical$Economy) # 286 mining, 1996 scavenging 

Mining = subset(stat_no_NM,Economy == "Mining")
Scavenging = subset(stat_no_NM,Economy == "Scavenging")

arid =  subset (stat_no_NM, aridity < 0.65)
humid =  subset (stat_no_NM, aridity >= 0.65)

########### SEM beta model all points #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = stat_no_NM,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = stat_no_NM,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI + Nfix_level + Economy_level ,random = ~1|Genus, data = stat_no_NM)
)
summary(beta_sem)
plot(beta_sem)

########### SEM beta model Temperate #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = Temperate,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = Temperate,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI + Nfix_level + Economy_level ,random = ~1|Genus, data = Temperate)
)
summary(beta_sem)
plot(beta_sem)

########### SEM beta model Tropical #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = Tropical,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = Tropical,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI + Nfix_level + Economy_level ,random = ~1|Genus, data = Tropical)
)
summary(beta_sem)
plot(beta_sem)

########### SEM beta model scavenging #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = Scavenging,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = Scavenging,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI + Nfix_level ,random = ~1|Genus, data = Scavenging)
)
summary(beta_sem)
plot(beta_sem)

########### SEM beta model Mining #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = Mining,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = Mining,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI + Nfix_level ,random = ~1|Genus, data = Mining)
)
summary(beta_sem)
plot(beta_sem)

##################### Merge plots for figures 3, 4, and 5 ########################################

Fig3 <- plot_grid(NC_plot_bis, Nfix_Plot_bis, ncol = 2, 
                         align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                         label_x = c(0.12, 0.12))
Fig3


ggsave("fig/Fig3.tiff", Fig3, 
       width = 60, height = 30, units = "cm", dpi = 1000, type = "cairo")

NC_plot_bis.g <- ggplotGrob(NC_plot_bis)
Nfix_Plot_bis.g <- ggplotGrob(Nfix_Plot_bis)

Fig3_bis <- cbind(NC_plot_bis.g,Nfix_Plot_bis.g,size = 'max')

jpeg(filename = "Fig3_bis.jpeg", 
     width = 20, height = 10, units = 'in', res = 800)
grid.newpage()
grid.draw(Fig3_bis)
grid.text("(a)", x = 0.01, y = 0.95, just = "left", gp = gpar(fontsize = 30, fontface = "bold"))  # Top-left
grid.text("(b)", x = 0.55, y = 0.95, just = "left", gp = gpar(fontsize = 30, fontface = "bold"))  # Top-right
dev.off()


Fig4 <- plot_grid(MI_plot_bis, P_plot_bis, ncol = 2, 
                  align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                  label_x = c(0.12, 0.12))
Fig4

ggsave("fig/Fig4.tiff", Fig4, 
       width = 60, height = 30, units = "cm", dpi = 1000, type = "cairo")


MI_plot_bis.g <- ggplotGrob(MI_plot_bis)
P_plot_bis.g <- ggplotGrob(P_plot_bis)

Fig4_bis <- cbind(MI_plot_bis.g,P_plot_bis.g ,size = 'max')

jpeg(filename = "Fig4_bis.jpeg", 
     width = 20, height = 10, units = 'in', res = 800)
grid.newpage()
grid.draw(Fig4_bis)
grid.text("(a)", x = 0.1, y = 0.95, just = "left", gp = gpar(fontsize = 30, fontface = "bold"))  # Top-left
grid.text("(b)", x = 0.53, y = 0.95, just = "left", gp = gpar(fontsize = 30, fontface = "bold"))  # Top-right
dev.off()





################################ Figures Supplementary Information ################################################
###################################################################################################


################################ Fig S1 and S2 ################################################
###################################################################################################

install.packages("rnaturalearth")
install.packages("rnaturalearthdata")

library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(plotbiomes)
library(gridExtra)
library(grid)

whittaker_plot <- whittaker_base_plot() +
  guides(color = guide_legend(title = NULL)) 
whittaker_plot

whittaker_plot <- whittaker_base_plot() +
  theme(legend.position = c(0.15, 0.84),
        legend.text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 30, colour = 'black'),
        axis.text.y = element_text(size = 30, colour = 'black'))
whittaker_plot 

stat_no_NM$pre_mean_cm = stat_no_NM$pre_mean/10
stat_no_NM$tmp
names(stat_no_NM)
beta_diagram <- subset(stat_no_NM,pre_mean_cm>0 ) ## these points where pre_mean =0 are in fact NA
nrow(beta_diagram) # 4428
names(beta_diagram)
beta_diagram$Economy
table(beta_diagram$N_fix)
length(table(beta_diagram$species)) # 2406

beta_diagram$N_fix <- ifelse(beta_diagram$N_fix == "Non Fixers", "OP", 
                             ifelse(beta_diagram$N_fix == "Fixers", "N2FP", beta_diagram$N_fix))


myco_color_palette <- c("Scavenging" = "magenta", "Mining" = "blue")
N_shape_palette <- c("N2FP" = 8,  # Filled circle
                     "OP" = 16)  # Filled triangle


whittaker_beta <- whittaker_plot +
  geom_point(data = beta_diagram, aes(x = tmp, y = pre_mean_cm, color=Economy, shape = N_fix),size = 2.5)+
  scale_color_manual(values = myco_color_palette)+
  scale_shape_manual(values = N_shape_palette)+
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 1, override.aes = list(size = 5))) +
  theme(legend.position = c(0.15, 0.68),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 30, colour = 'black'),
        axis.text.y = element_text(size = 30, colour = 'black'))+
  guides(fill = guide_legend(
    title= NULL,
    title.position = "top",
    override.aes = list(size = 8),
    keywidth = unit(1, "cm"),
    keyheight = unit(1, "cm"),
    #nrow = 20,
    byrow = TRUE))
whittaker_beta

ggsave("Fig/whittaker_beta.tiff", whittaker_beta, 
       width = 45, height = 25, units = "cm", dpi = 800, type = "cairo")

names(beta_diagram)

#################### World map ##############################################
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create world map with points
Worldmap_beta <- ggplot(data = world) +
  geom_sf(fill = "white", color = "black") + # Set the map background to white
  geom_point(data = beta_diagram, aes(x = lon, y = lat, color = Economy, shape = N_fix), size=3) +
  scale_color_manual(values = myco_color_palette) +
  scale_shape_manual(values = N_shape_palette) +
  theme(
    legend.position = c(0.11, 0.30),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"), # Set overall background to white
    panel.grid.major = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.x = element_text(size = 30, colour = 'black'),
    axis.text.y = element_text(size = 30, colour = 'black')
  ) +
  guides(
    color = guide_legend(title = NULL), # Remove legend title for color (Economy)
    shape = guide_legend(title = NULL), # Remove legend title for shape (N_fix)
    fill = guide_legend(
      title = NULL,
      title.position = "top",
      override.aes = list(size = 8),
      keywidth = unit(1, "cm"),
      keyheight = unit(1, "cm"),
      byrow = TRUE
    )
  )
Worldmap_beta

ggsave("Fig/Worldmap_beta.tiff", Worldmap_beta, 
       width = 45, height = 25, units = "cm", dpi = 800, type = "cairo")


################### proportions ########################### 
# Calculate the overall proportion of scavenging and mining in Economy, and OP vs N2FP
proportions_economy <- beta_diagram %>%
  group_by(Economy) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

proportions_economy #### Mining: 24%, scavenging 76%

proportions_Nfix<- beta_diagram %>%
  group_by(N_fix) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))
proportions_Nfix #### N2FP: 8.6%, scavenging 91.4%


# Calculate the proportion of OP and N2FP within each Economy group
proportions_nfix_by_economy <- beta_diagram %>%
  group_by(Economy, N_fix) %>%
  summarise(count = n()) %>%
  group_by(Economy) %>%
  mutate(proportion = count / sum(count))

proportions_nfix_by_economy
### Mining-N2FP : 5.8%
### Mining-OP: 94.2%
### Scavenging-N2FP:9.5%
### Mining-N2FP:90.5%


################################ Fig S3, S4 and S5 ################################################
###################################################################################################


############################## MI vs big Delta 13C #################################
MI_big_D13<- (MI_big_D13 <- ggplot(data = stat_no_NM, aes(x = sqrtMI, y = big_D13, fill=Economy)) + 
                
                scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                                  breaks = c("Scavenging", "Mining"),
                                  labels = c("Scavenging", "Mining")) +
                
                geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                            aes(fill = Economy), shape = 24, color = "tan3", 
                            alpha = 0.5, size = 3) +
                geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                            aes(fill = Economy), shape = 21, color = "lightseagreen", 
                            alpha = 0.5, size = 3) +
                
                #geom_line(data = AI_trend_scavenging, aes(x = AI_seq_scav , y = AI_trend_scavenging ), col = 'lightseagreen',  lwd = 2, alpha = 1.5) +
                #geom_line(data = AI_trend_mining, aes(x = AI_seq_mining, y = AI_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +
                
                scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
                theme(legend.position ="none",
                      #legend.justification = "top",
                      axis.title.y = element_text(size = 50, colour = 'black'),
                      axis.title.x = element_text(size = 50, colour = 'black'),
                      axis.text.x = element_text(size = 50, colour = 'black'),
                      axis.text.y = element_text(size = 50, colour = 'black'),
                      axis.line = element_line(size = 2, colour = "black"), 
                      axis.ticks = element_line(size = 1, colour = "black"), 
                      #panel.background = element_rect(fill = 'white', colour = 'black'),
                      #panel.grid.major = element_line(colour = "white"),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank(),
                      legend.text = element_text(size=30),
                      legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
                      legend.key = element_rect(color = "white", fill = "white"))+
                #legend.position =c(0.98, 0.98),
                #legend.justification = c(1, 1)) +
                ylab(expression('∆'* italic('13')['C'] * ' (‰)')) +
                xlab(expression(bold("\u221A") * MI))+
                guides(fill = guide_legend(
                  title= NULL,
                  #title.position = "top",
                  override.aes = list(size = 5),
                  keywidth = unit(1, "cm"),
                  keyheight = unit(1.5, "cm"),
                  title.hjust = 0.5,
                  nrow = 2,
                  byrow = TRUE
                )))
MI_big_D13
##################################################################################################################

############################## ln_NC_topsoil vs big Delta 13C #################################
NC_big_D13<- (NC_big_D13 <- ggplot(data = stat_no_NM, aes(x = ln_NC_topsoil, y = big_D13, fill=Economy)) + 
                
                scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                                  breaks = c("Scavenging", "Mining"),
                                  labels = c("Scavenging", "Mining")) +
                
                geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                            aes(fill = Economy), shape = 24, color = "tan3", 
                            alpha = 0.5, size = 3) +
                geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                            aes(fill = Economy), shape = 21, color = "lightseagreen", 
                            alpha = 0.5, size = 3) +
                
                #geom_line(data = AI_trend_scavenging, aes(x = AI_seq_scav , y = AI_trend_scavenging ), col = 'lightseagreen',  lwd = 2, alpha = 1.5) +
                #geom_line(data = AI_trend_mining, aes(x = AI_seq_mining, y = AI_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +
                
                scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
                theme(legend.position ="none",
                      #legend.justification = "top",
                      axis.title.y = element_text(size = 50, colour = 'black'),
                      axis.title.x = element_text(size = 50, colour = 'black'),
                      axis.text.x = element_text(size = 50, colour = 'black'),
                      axis.text.y = element_text(size = 50, colour = 'black'),
                      axis.line = element_line(size = 2, colour = "black"), 
                      axis.ticks = element_line(size = 1, colour = "black"), 
                      #panel.background = element_rect(fill = 'white', colour = 'black'),
                      #panel.grid.major = element_line(colour = "white"),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank(),
                      legend.text = element_text(size=30),
                      legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
                      legend.key = element_rect(color = "white", fill = "white"))+
                #legend.position =c(0.98, 0.98),
                #legend.justification = c(1, 1)) +
                ylab(expression('∆'* italic('13')['C'] * ' (‰)')) +
                xlab(expression('ln ' * 'N/C ratio'))+
                guides(fill = guide_legend(
                  title= NULL,
                  #title.position = "top",
                  override.aes = list(size = 5),
                  keywidth = unit(1, "cm"),
                  keyheight = unit(1.5, "cm"),
                  title.hjust = 0.5,
                  nrow = 2,
                  byrow = TRUE
                )))
NC_big_D13
##################################################################################################################

############################## lnP vs big Delta 13C #################################
P_big_D13<- (P_big_D13 <- ggplot(data = stat_no_NM, aes(x = lnP, y = big_D13, fill=Economy)) + 
               
               scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                                 breaks = c("Scavenging", "Mining"),
                                 labels = c("Scavenging", "Mining")) +
               
               geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                           aes(fill = Economy), shape = 24, color = "tan3", 
                           alpha = 0.5, size = 3) +
               geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                           aes(fill = Economy), shape = 21, color = "lightseagreen", 
                           alpha = 0.5, size = 3) +
               
               #geom_line(data = AI_trend_scavenging, aes(x = AI_seq_scav , y = AI_trend_scavenging ), col = 'lightseagreen',  lwd = 2, alpha = 1.5) +
               #geom_line(data = AI_trend_mining, aes(x = AI_seq_mining, y = AI_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +
               
               scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
               theme(legend.position ="none",
                     #legend.justification = "top",
                     axis.title.y = element_text(size = 50, colour = 'black'),
                     axis.title.x = element_text(size = 50, colour = 'black'),
                     axis.text.x = element_text(size = 50, colour = 'black'),
                     axis.text.y = element_text(size = 50, colour = 'black'),
                     axis.line = element_line(size = 2, colour = "black"), 
                     axis.ticks = element_line(size = 1, colour = "black"), 
                     #panel.background = element_rect(fill = 'white', colour = 'black'),
                     #panel.grid.major = element_line(colour = "white"),
                     panel.border = element_blank(),
                     panel.grid = element_blank(),
                     panel.background = element_blank(),
                     legend.text = element_text(size=30),
                     legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
                     legend.key = element_rect(color = "white", fill = "white"))+
               #legend.position =c(0.98, 0.98),
               #legend.justification = c(1, 1)) +
               ylab(expression('∆'* italic('13')['C'] * ' (‰)')) +
               xlab(expression('ln ' * P[i])) +
               guides(fill = guide_legend(
                 title= NULL,
                 #title.position = "top",
                 override.aes = list(size = 5),
                 keywidth = unit(1, "cm"),
                 keyheight = unit(1.5, "cm"),
                 title.hjust = 0.5,
                 nrow = 2,
                 byrow = TRUE
               )))
P_big_D13
##################################################################################################################
stat_no_NM
############################# lnP vs big Delta 13C #################################
beta_big_D13<- (beta_big_D13 <- ggplot(data = stat_no_NM, aes(x = big_D13, y = lnbeta, fill=Economy)) + 
                  
                  scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                                    breaks = c("Scavenging", "Mining"),
                                    labels = c("Scavenging", "Mining")) +
                  
                  geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                              aes(fill = Economy), shape = 24, color = "tan3", 
                              alpha = 0.5, size = 4) +
                  geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                              aes(fill = Economy), shape = 21, color = "lightseagreen", 
                              alpha = 0.5, size = 4) +
                  
                  #geom_line(data = AI_trend_scavenging, aes(x = AI_seq_scav , y = AI_trend_scavenging ), col = 'lightseagreen',  lwd = 2, alpha = 1.5) +
                  #geom_line(data = AI_trend_mining, aes(x = AI_seq_mining, y = AI_trend_mining), col = 'tan3', linetype="dashed", lwd = 1.5, alpha = 0.8) +
                  
                  scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
                  theme(legend.position ="none",
                        #legend.justification = "top",
                        axis.title.y = element_text(size = 50, colour = 'black'),
                        axis.title.x = element_text(size = 50, colour = 'black'),
                        axis.text.x = element_text(size = 50, colour = 'black'),
                        axis.text.y = element_text(size = 50, colour = 'black'),
                        axis.line = element_line(size = 2, colour = "black"), 
                        axis.ticks = element_line(size = 1, colour = "black"), 
                        #panel.background = element_rect(fill = 'white', colour = 'black'),
                        #panel.grid.major = element_line(colour = "white"),
                        panel.border = element_blank(),
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        legend.text = element_text(size=30),
                        legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
                        legend.key = element_rect(color = "white", fill = "white"))+
                  #legend.position =c(0.98, 0.98),
                  #legend.justification = c(1, 1)) +
                  ylab(expression('ln ' * italic('β'))) +
                  xlab(expression('∆'* italic('13')['C'] * ' (‰)')) +
                  
                  guides(fill = guide_legend(
                    title= NULL,
                    #title.position = "top",
                    override.aes = list(size = 5),
                    keywidth = unit(1, "cm"),
                    keyheight = unit(1.5, "cm"),
                    title.hjust = 0.5,
                    nrow = 2,
                    byrow = TRUE
                  )))
beta_big_D13
##################################################################################################################

############################# lnP vs big Delta 13C #################################
Chi_beta <- (Chi_beta <- ggplot(data = stat_no_NM, aes(x = chi, y = lnbeta, fill=Economy)) + 
                  
                  scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                                    breaks = c("Scavenging", "Mining"),
                                    labels = c("Scavenging", "Mining")) +
               
                  geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                              aes(fill = Economy), shape = 24, color = "tan3", 
                              alpha = 0.5, size = 4) +
                  geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                              aes(fill = Economy), shape = 21, color = "lightseagreen", 
                              alpha = 0.5, size = 4) +
               
                  scale_x_continuous(labels = function(x) ifelse(x %% 1 == 0, as.character(round(x)), as.character(x))) +
                  theme(legend.position ="right",
                        legend.justification = "top",
                        axis.title.y = element_text(size = 50, colour = 'black'),
                        axis.title.x = element_text(size = 50, colour = 'black'),
                        axis.text.x = element_text(size = 50, colour = 'black'),
                        axis.text.y = element_text(size = 50, colour = 'black'),
                        axis.line = element_line(size = 2, colour = "black"), 
                        axis.ticks = element_line(size = 1, colour = "black"), 
                        #panel.background = element_rect(fill = 'white', colour = 'black'),
                        #panel.grid.major = element_line(colour = "white"),
                        panel.border = element_blank(),
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        legend.text = element_text(size=30),
                        legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
                        legend.key = element_rect(color = "white", fill = "white"))+
                  #legend.position =c(0.98, 0.98),
                  #legend.justification = c(1, 1)) +
                  xlab(expression(chi['isotopes'] * " (Pa Pa"^{-1} * ")"))+
               ylab(expression('ln ' * italic('β'))) +
                  
                  guides(fill = guide_legend(
                    title= NULL,
                    #title.position = "top",
                    override.aes = list(size = 5),
                    keywidth = unit(1, "cm"),
                    keyheight = unit(1.5, "cm"),
                    title.hjust = 0.5,
                    nrow = 2,
                    byrow = TRUE
                  )))
Chi_beta
##################################################################################################################

##################### Merge plots ########################################
merged_FigS4 <- plot_grid(beta_big_D13,Chi_beta, ncol = 2, 
                          align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                          label_x = c(0.12, 0.12))
merged_FigS4


merged_FigS5 <- plot_grid(NC_big_D13,MI_big_D13, P_big_D13,ncol = 3, 
                            align = "vh", labels = c("(a)", "(b)", "(c)"), label_size = 30, 
                            label_x = c(0.12, 0.12, 0.12))
merged_FigS5

ggsave("fig/merged_FigS4.tiff", merged_FigS4, 
       width = 80, height = 30, units = "cm", dpi = 800, type = "cairo")

ggsave("fig/merged_FigS5.tiff", merged_FigS5, 
       width = 80, height = 30, units = "cm", dpi = 1000, type = "cairo")


########################## Figure S6 ##################################################################

####################### lnP depending on biomes ####################################
P_lmer <- lmer(lnP~ biome +
                 (1|Genus),
               data = stat_no_NM)
Anova(P_lmer)
summary(P_lmer)
emmeans_results_P <- emmeans(P_lmer, pairwise ~ biome)

test_P <- cld(emmeans(P_lmer, ~biome))
test_P [2, 7]

summary(emmeans_results_P)
P_letters <- data.frame(x = c(1, 2),
                        biome = c("Temperate", "Tropical"),
                        lnP = c(100,200),
                        y = c(6, 3), 
                        group = c(test_P[2,7],test_P[1,7]))

P_letters$letter[P_letters$group == "  2"] <- "a"
P_letters$letter[P_letters$group == " 1 "] <- "b"

############# Box plot lnP-Biomes #######################################################
biome_P <- ggplot(data = stat_no_NM, 
                  aes(x = biome, y = lnP, fill=biome)) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) + 
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y = element_text(size = 60, colour = 'black'),
        axis.title.x = element_text(size = 40, colour = 'black'),
        axis.text.x = element_text(size = 40, colour = 'black'),
        axis.text.y = element_text(size = 60, colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white")) +
  geom_boxplot(outlier.color = NA) +
  geom_text(data = P_letters, aes(x = x, y = y, label = letter), size = 16) +
  labs(fill = "Biomes") +
  ylab(expression('ln ' * 'P')) +
  xlab('Biomes')
biome_P

####################### lnNC depending on biomes ####################################
NC_lmer <- lmer(ln_NC_topsoil ~ biome +
                  (1|Genus),
                data = stat_no_NM)
Anova(NC_lmer)
summary(NC_lmer)
emmeans_results_NC <- emmeans(NC_lmer, pairwise ~ biome)

test_N <- cld(emmeans(NC_lmer, ~biome))
test_N [2, 7]

summary(emmeans_results_NC)
N_letters <- data.frame(x = c(1, 2),
                        biome = c("Temperate", "Tropical"),
                        ln_NC_topsoil = c(70,70),
                        y = c(-1.2, -1.8), 
                        group = c(test_N[2,7],test_N[1,7]))

N_letters$letter[N_letters$group == "  2"] <- "a"
N_letters$letter[N_letters$group == " 1 "] <- "b"

############# Box plot lnNC-Biomes #######################################################
biome_NC <- ggplot(data = stat_no_NM, 
                   aes(x = biome, y = ln_NC_topsoil, fill=biome)) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) + 
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y = element_text(size = 60, colour = 'black'),
        axis.title.x = element_text(size = 40, colour = 'black'),
        axis.text.x = element_text(size = 40, colour = 'black'),
        axis.text.y = element_text(size = 60, colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white")) +
  geom_boxplot(outlier.color = NA) +
  geom_text(data = N_letters, aes(x = x, y = y, label = letter), size = 16) +
  labs(fill = "Biomes") +
  ylab(expression('ln ' * 'N/C ratio'))+
  xlab('Biomes')

biome_NC

library(gridExtra)
library(grid)

merged_FigS6 <- plot_grid(biome_NC, biome_P, ncol = 2, 
                       align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                       label_x = c(0.90, 0.90))
merged_FigS6


ggsave("fig/merged_FigS6 .tiff", merged_FigS6 , 
       width = 60, height = 30, units = "cm", dpi = 1000, type = "cairo")
