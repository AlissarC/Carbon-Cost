
library(tidyverse)

final_data <- read.csv("./data/data_clean_beta.csv")
final_data
names(final_data)
nrow(final_data) #5876

####### calculus of nitrogen, soc, and C/N means for all layers ############################## 

nitrogen <- final_data [grep("nitrogen_", colnames(final_data))] #selection of columns containing nitrogen
nitrogen
final_data$nitrogen_sum = rowSums(nitrogen, na.rm = T)
final_data$nitrogen_mean = rowMeans(nitrogen, na.rm = T)
names(final_data)
final_data$nitrogen_mean

soc <- final_data [grep("soc_", colnames(final_data))] #selection of columns containing soc
soc
final_data$soc_sum = rowSums(soc, na.rm = T)
final_data$soc_mean = rowMeans(soc, na.rm = T)

final_data$ratio_CN_mean <- final_data$soc_mean/(final_data$nitrogen_mean)
#######################################################################################
names(final_data)
############################# ratios C/N for each layer and transformation of the aridity index  #############################

final_data <- final_data %>% 
  mutate(ratio_CN_0_5cm = soc_0.5cm/(nitrogen_0.5cm)) %>%  
  mutate(ratio_CN_15_30cm = soc_15.30cm/(nitrogen_15.30cm)) %>% 
  mutate(ratio_CN_30_60cm = soc_30.60cm/(nitrogen_30.60cm)) %>% 
  mutate(aridity = awi_pm_sr_yr/10000)
final_data

############ calculation of the mean of some soil variables for top soil : 0-5cm, 5-15cm ###########
nitrogen_topsoil <- rowMeans(final_data[, c("nitrogen_0.5cm", "nitrogen_5.15cm")])
soc_topsoil <- rowMeans(final_data[, c("soc_0.5cm", "soc_5.15cm")])

final_data$nitrogen_topsoil <- nitrogen_topsoil
final_data$soc_topsoil <- soc_topsoil


final_data$ratio_CN_topsoil <- final_data$soc_topsoil/final_data$nitrogen_topsoil
final_data$ratio_NC_topsoil <- final_data$nitrogen_topsoil/final_data$soc_topsoil

hist(final_data$ratio_NC_topsoil)
hist(final_data$ratio_CN_topsoil)

###################################################################################################

############ calculation of the mean of some soil variables for sub soil : 15-30cm, 30-60cm ###########
nitrogen_subsoil <- rowMeans(final_data[, c("nitrogen_15.30cm", "nitrogen_30.60cm")])
soc_subsoil <- rowMeans(final_data[, c("soc_15.30cm", "soc_30.60cm")])


final_data$nitrogen_subsoil <- nitrogen_subsoil
final_data$soc_subsoil <- soc_subsoil

final_data$ratio_CN_subsoil <- final_data$soc_subsoil/final_data$nitrogen_subsoil
final_data$ratio_NC_subsoil <- final_data$nitrogen_subsoil/final_data$soc_subsoil

hist(final_data$ratio_CN_subsoil)
###################################################################################################


##################### calculation for the first 30 cm of the soil ################
nitrogen_0_30cm <- rowMeans(final_data[, c("nitrogen_0.5cm", "nitrogen_5.15cm", "nitrogen_15.30cm" )])
final_data$nitrogen_0_30cm <- nitrogen_0_30cm

soc_0_30cm <- rowMeans(final_data[, c("soc_0.5cm", "soc_5.15cm", "soc_15.30cm" )])
final_data$soc_0_30cm <- soc_0_30cm

final_data$ratio_CN_0_30cm <- final_data$soc_0_30cm/final_data$nitrogen_0_30cm
final_data$ratio_NC_0_30cm <- final_data$nitrogen_0_30cm/final_data$soc_0_30cm

hist(final_data$ratio_CN_0_30cm)
hist(final_data$ratio_NC_0_30cm)

nrow(final_data) # 5876
names(final_data)

#### replace NA by Non-fixers
final_data$N_fix <- ifelse(is.na(final_data$N_fix), "Non_fixers", final_data$N_fix)

######### remove some rare myco_type 
final_data <- subset(final_data, myco_type != "uncertain" & myco_type != "species-specific EcM-AM or AM or NM" &myco_type != "OM" &myco_type != "species-specific: AM or rarely EcM-AM or AM")
table(final_data$myco_type)

######################### change the name of myco_groups and set them as factors ##############
final_data$myco_type <- as.factor(final_data$myco_type)
final_data$myco_type <- ifelse(final_data$myco_type == "EcM-AM", "EcMAM", as.character(final_data$myco_type))
final_data$myco_type <- ifelse(final_data$myco_type == "NM-AM", "NMAM", as.character(final_data$myco_type))
table(final_data$myco_type)
nrow(final_data) # 5821

############################# call functions to calculate beta and chi #################################################

install.packages ("R.utils") # to read source function
library (R.utils)

source("./calc_optimal_vcmax.R")

sourceDirectory("./functions", modifiedOnly = FALSE)

data_clean = final_data
################ calculate atmospheric pressure (Pa) from elevation (m)#####################

calc_patm = function(z) {
  
  kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo = 298.15   # base temperature, K (Prentice, unpublished)
  kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
  kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
  
  patm
}

data_clean$patm<- calc_patm(data_clean$z)


################ calculate  nstar (unitless relative viscosity of h2o at temperature relative to 25°C) ####################

calc_nstar = function(temp, z){ # temp in °C and z in m
  
  patm = calc_patm(z)
  
  # viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
  ns      = calc_viscosity_h2o( temp, z )  # Pa s 
  ns25    = calc_viscosity_h2o( 25, z )  # Pa s 
  nstar = ns / ns25                       # (unitless)
  
  nstar
  
}
# Apparently it call automatically calc_viscositt_h2o 
data_clean$nstar<- calc_nstar(temp = data_clean$tmp, z = data_clean$z)


##########################calculate gammastar (Pa)###############################################

calc_gammastar_pa = function(temp, z) {
  
  patm = calc_patm(z)
  rat = calc_patm(z) / calc_patm(0)
  gammastar25 = 4.332 * rat  # Pa
  Hgm=37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6 * calc_patm(0)
  O2_z = O2 * 1e-6 * calc_patm(z)
  
  temp_k = 273.15+ temp
  
  gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
  
  gStar_pa
  
}
data_clean$gammastar_pa <-calc_gammastar_pa(temp = data_clean$tmp, z = data_clean$z)

######## calcualte the Michaelis-Menton coefficient (Pa) for Rubisco from temperature##############################

calc_km_pa = function(temp, z) {
  
  patm = calc_patm(z) 
  rat = patm / calc_patm(0)
  
  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03 * rat 
  Ko25 = 28210 * rat 
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6) * patm 
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
  
}

data_clean$Km<-calc_km_pa(temp = data_clean$tmp, z = data_clean$z)


names(data_clean)


#################### calculation of chi form big Delta 13c based on Lavergne et al., 2020  #######################
data_clean$ca <- data_clean$CO2 * 1e-6 * data_clean$patm
a = 4.4
b = 28
f=12
data_clean$chi= (data_clean$big_D13 -a +(f*data_clean$gammastar_pa)/data_clean$ca)/(b-a) # Lavargne et al., 2020

names(data_clean)
data_clean["chi"]
summary(data_clean$chi)
hist(data_clean$chi)

nrow(data_clean)

data_clean$D<- data_clean$vpd*1000
data_clean$beta_num <- 1.6 * data_clean$nstar * data_clean$D * ((data_clean$chi-data_clean$gammastar_pa/data_clean$ca) ^ 2)
data_clean$beta_denom <- ((1 - data_clean$chi)^2) * (data_clean$Km + data_clean$gammastar_pa)

data_clean$beta <- data_clean$beta_num / data_clean$beta_denom

summary(data_clean$beta_num)
hist(data_clean$beta_num)
plot(data_clean$beta_num, data_clean$chi)
plot(data_clean$beta_denom, data_clean$chi)

summary(data_clean$beta_denom)
hist(data_clean$beta_denom)

summary(data_clean$beta)
hist(data_clean$beta)

nrow(data_clean) # 5821

#############################################################################


############ assign AM and NMAM to scavenging economy, EcM , ErM, and EcM-AM to mining strategy

data_clean <- data_clean %>%
  mutate(Economy = case_when(
    myco_type %in% c("AM", "NMAM") ~ "Scavenging",
    myco_type %in% c("EcM", "EcMAM", "ErM") ~ "Mining",
    myco_type %in% c("NM") ~ "NM",
    TRUE ~ "Other"  
  ))
names(data_clean)
table(data_clean$Economy) # 1211: mining, 155:NM, 4455: scavenging 

########################## Log transformation to meet normality
data_clean$lnP = log(data_clean$phosphorus)
data_clean$ln_N_topsoil = log(data_clean$nitrogen_topsoil)
data_clean$ln_N_subsoil = log(data_clean$nitrogen_subsoil)
data_clean$ln_N_mean = log(data_clean$nitrogen_mean)

data_clean$ln_CN_topsoil = log(data_clean$ratio_CN_topsoil)
data_clean$ln_CN_subsoil = log(data_clean$ratio_CN_subsoil)
data_clean$ln_CN_mean = log(data_clean$ratio_CN_mean)

data_clean$ln_NC_topsoil = log(data_clean$ratio_NC_topsoil)
data_clean$ln_NC_0_30cm = log (data_clean$ratio_NC_0_30cm )

data_clean$lnbeta= log(data_clean$beta)

##################### select only C3 species and rows when 0.1<chi<0.95 ##################
nrow(data_clean)# 5821
max(data_clean$beta)
data_C3 <- data_clean %>% 
  filter(PS_pathway == "C3") %>% 
  filter(chi>0.1&chi<0.95) 
data_C3
data_C3$biome
nrow(data_C3) # 5293 
table(data_C3$PS_pathway)

write.csv(data_C3, "./data/data_C3.csv", row.names = F)

######## data_C3: data ready for analysis, next R file: stat_analysis_Beta_Final
