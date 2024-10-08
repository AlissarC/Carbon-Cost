
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

stat <- read.csv("../data/data_C3.csv")
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


################## Subset the variables we need  ######################################
###################################################################################################
stat_subset <-  subset(stat, select = c(species,Genus,N_fix, myco_type, Economy, ratio_NC_topsoil, 
                                        phosphorus, aridity, sqrtMI,
                                        lnP, ln_NC_topsoil,lnbeta, pre_mean, tmp))

colSums(is.na(stat_subset)) 
stat_subset <-  na.omit(stat_subset)## remove NA
nrow(stat_subset) # 4953

###################lnbeta outliers following MAD method ##############

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

# change the name of non fixers 
stat_no_NM$N_fix[stat_no_NM$N_fix == "Non_fixers"] <- "Non Fixers"

length(unique(stat_no_NM$species)) #2599

#################### Linear mixed effects model for beta #######################################
########################## H1: beta decreases with P and N increase, and with aridity decrease
########################## H2: These relationships are affected by mycorrhizal association,  N fixation affects only N/C-beta relationship

lmer_beta <- lmer(lnbeta~ lnP + ln_NC_topsoil + Economy + sqrtMI + 
                             sqrtMI*Economy + ln_NC_topsoil*Economy + lnP*Economy + N_fix + N_fix*ln_NC_topsoil +
                          (1|Genus), data=stat_no_NM)
Anova(lmer_beta) 
vif(lmer_beta) ## High VIF caused by the interaction term between N_fix and ln_NC_topsoil with no interaction found 


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
NC_plot_bis <- (NC_beta_plot_bis <- ggplot(data = stat_no_NM, aes(x = ln_NC_topsoil, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen","Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging","Mining")) +
      
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 2) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 2) +
      
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
          legend.text = element_text(size=30),
          legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
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
P_plot_bis <- (P_beta_plot_bis <- ggplot(data = stat_no_NM, aes(x = lnP, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging", "Mining")) +
      
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 2) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 2) +
    
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
          legend.text = element_text(size=30),
          legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
          legend.key = element_rect(color = "white", fill = "white")) +
          #legend.position =c(1.1, 1.1),
          #legend.justification = c(1, 1)) +
    ylab(expression('ln ' * italic('β'))) +
    xlab(expression('ln ' * 'P'))) +
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
MI_plot_bis <- (MI_beta_plot_bis <- ggplot(data = stat_no_NM, aes(x = sqrtMI, y = lnbeta, fill=Economy)) + 
    
    scale_fill_manual(values = c("Scavenging" = "lightseagreen", "Mining" = "tan3"), 
                      breaks = c("Scavenging", "Mining"),
                      labels = c("Scavenging", "Mining")) +
      
      geom_jitter(data = subset(stat_no_NM, Economy == "Mining"), 
                  aes(fill = Economy), shape = 24, color = "tan3", 
                  alpha = 0.5, size = 2) +
      geom_jitter(data = subset(stat_no_NM, Economy == "Scavenging"), 
                  aes(fill = Economy), shape = 21, color = "lightseagreen", 
                  alpha = 0.5, size = 2) +
    
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
          legend.text = element_text(size=30),
          legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
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
Fixers = subset(stat_no_NM, N_fix =="Fixers")
Non_Fixers = subset(stat_no_NM, N_fix =="Non Fixers")

fixers_seq <- seq(min(Fixers$ln_NC_topsoil, na.rm = T), max(Fixers$ln_NC_topsoil, na.rm = T), 0.001)
fixers_seq
N_trend_fixers <- int_N_fix + fixers_seq * Slope_N_fix
N_trend_fixers <- as.data.frame(cbind(fixers_seq, N_trend_fixers))
N_trend_fixers$N_fix= "Fixers"

Nonfixers_seq <- seq(min(Non_Fixers$ln_NC_topsoil, na.rm = T), max(Non_Fixers$ln_NC_topsoil, na.rm = T), 0.001)
Nonfixers_seq
N_trend_Nonfixers <- int_N_Nonfix + Nonfixers_seq * Slope_N_Nonfix
N_trend_Nonfixers <- as.data.frame(cbind(Nonfixers_seq, N_trend_Nonfixers))
N_trend_Nonfixers$N_fix= "Non Fixers"


N_fix_test
Nfix_Plot_bis <- (Nfix_plot_bis <- ggplot(data = stat_no_NM, aes(x = ln_NC_topsoil, y = lnbeta, fill=N_fix)) +
    
    scale_fill_manual(values = c("Fixers" = "red3", "Non Fixers" = "black")) +
      
      geom_jitter(data = subset(stat_no_NM, N_fix == "Fixers"), 
                  aes(fill = N_fix), shape = 24, color = "red3", 
                  alpha = 0.5, size = 2) +
      geom_jitter(data = subset(stat_no_NM, N_fix == "Non Fixers"), 
                  aes(fill = N_fix), shape = 21, color = "black", 
                  alpha = 0.5, size = 2) +
      
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
          legend.text = element_text(size=30),
          legend.background = element_rect(fill = "white", color = "white", size = 0.5), 
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

##################### Merge plots ########################################

merged_plot_N <- plot_grid(NC_plot_bis, Nfix_Plot_bis, ncol = 2, 
                              align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                              label_x = c(0.12, 0.12))
merged_plot_N


merged_MI_P <- plot_grid(MI_plot_bis, P_plot_bis, ncol = 2, 
                           align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                           label_x = c(0.12, 0.12))
merged_MI_P


ggsave("data/merged_plot_N.tiff", merged_plot_N, 
       width = 60, height = 30, units = "cm", dpi = 800, type = "cairo")

ggsave("data/merged_MI_P.tiff", merged_MI_P, 
       width = 60, height = 30, units = "cm", dpi = 800, type = "cairo")

###############SEM######################################################################
########################################################################################

########### Set categorical as levels   #########################################################
stat_no_NM$Nfix_level = 0 ## creation of a column with Ntrt= 0 
stat_no_NM$Nfix_level[stat_no_NM$N_fix == 'Fixers'] = 1 
stat_no_NM$Nfix_level[stat_no_NM$N_fix == 'Non Fixers'] = 2

stat_no_NM$Economy_level = 0 ## creation of a column with Ntrt= 0 
stat_no_NM$Economy_level[stat_no_NM$Economy == 'Mining'] = 1
stat_no_NM$Economy_level[stat_no_NM$Economy == 'Scavenging'] = 2

table(stat_no_NM$Economy_level)


########### SEM beta model all points #########################################################
beta_sem <- psem(
  lme1 <- lme(ln_NC_topsoil~sqrtMI , data = stat_no_NM,random = ~1|Genus), 
  lme2 <- lme(lnP~sqrtMI, data = stat_no_NM,random = ~1|Genus),
  lme3 <- lme(lnbeta~lnP+ ln_NC_topsoil + sqrtMI +
             Nfix_level + Economy_level, 
                random = ~1|Genus, data = stat_no_NM, method = "ML"),
  data = stat_no_NM
)
summary(beta_sem)
plot(beta_sem)

r_squared_lme1 <- r.squaredGLMM(lme1) # 0.33
r_squared_lme2 <- r.squaredGLMM(lme2) # 0.29
r_squared_lme3 <- r.squaredGLMM(lme3) # 0.41

line.thick <- data.frame(
  summary(beta_sem)$coefficients,
  line.thickness = abs(summary(
    beta_sem)$coefficients$Std.Estimate) * 16.67) %>%
  mutate(line.thickness = round(line.thickness, digits = 2)) %>%
  dplyr::select(-Var.9)
line.thick


