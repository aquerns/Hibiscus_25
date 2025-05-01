###Code:: 9-13-23

###Author:Aleah Querns
# last modified: A Louthan 241219

rm(list = ls())
load("Data/3_vital_rate_functions_and_vital_rate_values.RData") 
load("Data/4_sensitivity values.RData")
library("reshape2")
library("MuMIn")
library("lme4")
library("dplyr")
library("ggplot2")
library(gridExtra)
library(ggplot2)
library(MuMIn)
library(popbio)
library(tidyr)
library(data.table)
library(abind)
library(gridExtra)
library(R.utils)
library("ggpubr")
library(stringr)
library(boot)

# Get LTRE contributions (Vin - Vir / SensVir)----

# first, ensure the trmts are always in the same order
if (!identical(dimnames(full_vrs_allsites)[[4]], dimnames(full_vrs_averaged)[[3]]) | 
    !identical(dimnames(full_vrs_allsites)[[4]], dimnames(sens)[[1]])){ stop("trmt order is not identical")}


LTRE_contributions <- sweep(
          # getting V-Vref; i.e., difference between level-specific rates and across-level average
                  sweep (abind(	# this abind command makes an array of level-specific rates, dimensions of 3, 4, 16
                                    full_vrs_allsites["S","sr",,], 
                                    full_vrs_allsites["M","sr",,],
                                    full_vrs_allsites["L","sr",,], 
                                    full_vrs_allsites["S","fr",,], 
                                    full_vrs_allsites["M","fr",,], 
                                    full_vrs_allsites["L","fr",,], 
                                    full_vrs_allsites["S","transtoS",,], 
                                    full_vrs_allsites["M","transtoS" ,,], 
                                    full_vrs_allsites["L","transtoS",,],
                                    full_vrs_allsites["S","transtoM",,], 
                                    full_vrs_allsites["M","transtoM", ,], 
                                    full_vrs_allsites["L","transtoM",,], 
                                    full_vrs_allsites["S","transtoL",, ], 
                                    full_vrs_allsites["M","transtoL",,], 
                                    full_vrs_allsites["L","transtoL",,],
                                    full_vrs_allsites["F","germ",, ]
                                    , along=3), 
  MARGIN=c(2,3) , 
  FUN= "-", 
  STATS= 
  cbind( # this cbind command makes a matrix of average rates, dimensions 4, 16
    full_vrs_averaged["S","sr",], 
    full_vrs_averaged["M","sr",], 
    full_vrs_averaged["L","sr",], 
    full_vrs_averaged["S","fr",], 
    full_vrs_averaged["M","fr",], 
    full_vrs_averaged["L","fr",], 
    full_vrs_averaged["S","transtoS",], 
    full_vrs_averaged["M","transtoS" ,], 
    full_vrs_averaged["L","transtoS",],
    full_vrs_averaged["S","transtoM",], 
    full_vrs_averaged["M","transtoM", ], 
    full_vrs_averaged["L","transtoM",], 
    full_vrs_averaged["S","transtoL", ], 
    full_vrs_averaged["M","transtoL",], 
    full_vrs_averaged["L","transtoL",],
    full_vrs_averaged["F","germ", ])
) # this sweep command gives an array of dimensions:  3 levels, 4 treatments,16 demographic rates, 
,
# next multiplying this difference by the sensitivity for that demographic rate and treatment
MARGIN= c(2, 3), FUN= "*", 
STATS= sens)
dimnames(LTRE_contributions)[[3]] <- c("sr1", "sr2", "sr3", 
                                       "fr1", "fr2", "fr3", 
                                       "S-S", "M-S", "L-S", 
                                       "S-M", "M-M", "L-M", 
                                       "S-L", "M-L", "L-L", "germ")
save(LTRE_contributions, file= "Data/5_LTRE_contributions.RData")

LTRE_contributions_rn<-LTRE_contributions
dimnames(LTRE_contributions_rn)[[3]] <- c("sr_1", "sr_2", "sr_3", 
                                       "fr_1", "fr_2", "fr_3", 
                                       "TranstoS_1", "TranstoS_2", "TranstoS_3", 
                                       "TranstoM_1", "TranstoM_2", "TranstoM_3", 
                                       "TranstoL_1", "TranstoL_2", "TranstoL_3", "germ_0")

allltres <- as_tibble(as.data.frame.table(LTRE_contributions_rn)) %>%
  rename(Site = Var1, HerbTrt = Var2, VR= Var3, Value = Freq) %>%
  # Separate the VR column into vital rate and size class
  separate(VR, into = c("VitalRate", "SizeClass"), sep = "_", convert = TRUE) %>%
  # Recode the size class numbers to their respective labels
  mutate(SizeClass = recode(SizeClass,
                            `0` = "Fruit",
                            `1` = "Small",
                            `2` = "Medium",
                            `3` = "Large"),
         VitalRate = recode(VitalRate, 
                            'fr' = "Fruiting",
                            "germ" = "Recruitment",
                            "sr" = "Survival",
                            "TranstoL" = "Transition to L",
                            "TranstoM" = "Transition to M",
                            "TranstoS" = "Transition to S"),
        Site = recode (Site, 
                       "C" = "Intermediate",
                       "N" = "Arid",
                       "S" = "Mesic") )

allltres$Site <-factor (allltres$Site, levels = c("Arid" , "Intermediate", "Mesic"))
allltres$SizeClass <-factor (allltres$SizeClass, levels = c("Fruit" , "Small", "Medium", "Large"))
allltres$VitalRate <-factor (allltres$VitalRate, levels = c("Survival" , "Fruiting", "Transition to S", "Transition to M", "Transition to L", "Recruitment"))


# plotting LTRE_contributions----

jpeg(filename = "all_ltre_values.jpg", width = 10, height = 10, units = "in", res = 300)

ggplot(allltres, aes(x=Site, y=Value, col=SizeClass, group=SizeClass)) +
  geom_point(size=2, alpha=0.7, aes(shape=SizeClass))+
  geom_line(size=1, alpha=0.5)+
  facet_grid(rows=vars(VitalRate), cols=vars(HerbTrt),scales="free")+
  theme_bw(base_size=16) +
  scale_color_manual(values = c("orange1", "royalblue1", "springgreen4", "firebrick2"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("LTRE Value")+ guides(col=guide_legend(title="Life Stage/Size Class"), shape=guide_legend(title="Life Stage/Size Class"))

dev.off()


