# project: Hibiscus demographic compenstaion
# This script cleans data, then provides size and time classes
# written: October 16, 2015 (Allison Louthan)
# modified: by Aleah
# packages: 
library("reshape2")
library("MuMIn")
library("lme4")
library("dplyr")
library("ggplot2")
library(gridExtra)
library(ggplot2)



#################################################################################
# Assigning size classes########################################################

hibiscusbiannual <- read.csv(file = "~/RawData/hibiscus_biannual.csv", header=TRUE)

hibiscusbiannual$level <- as.factor(hibiscusbiannual$level)
hibiscusbiannual$block <- as.factor(hibiscusbiannual$block)
hibiscusbiannual$BM <- as.numeric(hibiscusbiannual$BM)
hibiscusbiannual$trmtUHURU <- as.factor(hibiscusbiannual$trmtUHURU )
hibiscusbiannual$year <- as.factor(hibiscusbiannual$year)
hibiscusbiannual$binfruits <- as.numeric(hibiscusbiannual$binfruits)
hibiscusbiannual$fruits <- as.numeric(hibiscusbiannual$fruits)
hibiscusbiannual$binsurnext6mo <- as.numeric(hibiscusbiannual$binsurnext6mo)
hibiscusbiannual$woodystems <- as.numeric(hibiscusbiannual$woodystems)
hibiscusbiannual$height <- as.numeric(hibiscusbiannual$height)
hibiscusbiannual$plantID <- as.numeric(hibiscusbiannual$plantID)
hibiscusbiannual$BA <- as.numeric(hibiscusbiannual$BA)

# small, medium and large size class designations:
hist(hibiscusbiannual$BM, breaks= 100) # looks nice and lognormal
quantile(hibiscusbiannual$BM, c(1/3, 2/3), na.rm=TRUE) #split three quantiles

# Add S M L size classes: 0-8.705845 is small, 8.705845- 17.255058  is medium, >17.255058  is large
hibiscusbiannual$sclass1 <- as.factor(c(rep(NA,length(hibiscusbiannual$X))))
hibiscusbiannual$sclass2 <- as.factor(c(rep(NA,length(hibiscusbiannual$X))))

#sizeclass in current timestep
hibiscusbiannual <- hibiscusbiannual %>%  mutate(sclass1 = case_when(BM >= 17.255058 ~ '3',
                                                                     BM >= 8.705845 ~ '2',
                                                                     TRUE ~ '1'))
#sizeclass in next timestep
hibiscusbiannual <- hibiscusbiannual %>%  mutate(sclass2 = case_when(BMnext6mo >= 17.255058 ~ '3',
                                                                     BMnext6mo >= 8.705845 ~ '2',
                                                                     TRUE ~ '1'))

#create column for survivor values to make calculating survivorship easier #OK OR JUST USE BINSURNEXT6MO
#hibiscusbiannual1 <- hibiscusbiannual %>%  mutate(sur = case_when(surnext6mo == "sur" ~ '1',
# surnext6mo== "dd" ~ '0',
#TRUE ~ 'NA'))


options(warn= 2)
years <- c("11", "12w", "12s","13w", "13s", "14w", "14s") # 11 means 11s
trmtUHURUs <- c("CONT", "LMH", "MEGA", "MESO") # cont= open control, means all herbivores have access
# LMH means that all large mammlian herbivores are excluded
# MEGA means only megaherbivores are excluded (elephant, giraffe)
# MESO means that everything larger than a mesoherbivor is excluded 
# (includes "mega" herbivores, deer- sized herbivores)
levels <- c("S", "C", "N")# S is wet, C is intermediate, N is arid

################################################################################
#####################Create "time" column--factor of season########
levels(hibiscusbiannual$year)

#Assign timepoints to levels...
#11=0
#12W=1
#12S=2
#13W=3
#13S=4
#14W=5
#14S=6
#15W=7
hibiscusbiannual$timepoint <- c(rep(NA,length(hibiscusbiannual$plantID)))


for (i in 1:dim(hibiscusbiannual)[1] ){
  if(hibiscusbiannual$year[i]== "11"){hibiscusbiannual$timepoint[i] <- "0"}
  if(hibiscusbiannual$year[i]== "12w"){hibiscusbiannual$timepoint[i] <- "1"}
  if(hibiscusbiannual$year[i]== "12s"){hibiscusbiannual$timepoint[i] <- "2"}
  if(hibiscusbiannual$year[i]== "13w"){hibiscusbiannual$timepoint[i] <- "3"}
  if(hibiscusbiannual$year[i]== "13s"){hibiscusbiannual$timepoint[i] <- "4"}
  if(hibiscusbiannual$year[i]== "14w"){hibiscusbiannual$timepoint[i] <- "5"} 
  if(hibiscusbiannual$year[i]== "14s"){hibiscusbiannual$timepoint[i] <- "6"}
  if(hibiscusbiannual$year[i]== "15w"){hibiscusbiannual$timepoint[i] <- "7"}
}

######Checking
#are all seedlings in the small size class?
hibiscusbiannualS <- hibiscusbiannual %>% dplyr::filter(hibiscusbiannual$height<10 &hibiscusbiannual$BA <3.14)
hibiscusbiannualS$woodystems <- as.numeric(hibiscusbiannualS$woodystems)

#what are the values here?
unique(hibiscusbiannualS$woodystems)#keeping in mind not all plants apparently have measurements for woody stems

mean(hibiscusbiannualS[,"BM"], na.rm=TRUE)
mean(hibiscusbiannualS[,"woodystems"], na.rm=TRUE)

#exclude the 1 weird small plant with many woody stems

weirdplant <- hibiscusbiannual %>% dplyr::filter(hibiscusbiannual$height<10 & hibiscusbiannual$BA <3.14 & hibiscusbiannual$woodystems >3)

hibiscusbiannual  <- hibiscusbiannual %>% dplyr::filter(hibiscusbiannual$plantID!=461)

#############Remove pseudocontrols--not necessary for further analysis
#remove pseudocontrols

hibiscusbiannual<- hibiscusbiannual %>% filter(block!="N1p",
                             block!="N2p",
                             block!="S1p",
                             block!="S2p")
#rename blocks to integer numbers:easier for storing in arrays later
hibiscusbiannual <-  hibiscusbiannual %>%  mutate(block = case_when(block== "N1" ~ 1,
                                                  block== "N2" ~ 2,
                                                  block == "N3" ~3,
                                                  block == "C1" ~1,
                                                  block == "C2" ~2,
                                                  block == "C3"~3,
                                                  block == "S1"~1,
                                                  block == "S2"~2,
                                                  block == "S3"~3))

unique(hibiscusbiannual$block) #check to make sure blocks renamed

#############


#exclude rows with NAs for TrtmtUHURU--these may be data entry errors

hibiscusbiannual_noNAs <- hibiscusbiannual %>% dplyr::filter(hibiscusbiannual$trmtUHURU!= "NA")

hibiscusbiannual_noNAs <- hibiscusbiannual_noNAs[,-c(1)]
#write version with timepoints and no trt NAs--use for IPM 
write.csv(hibiscusbiannual_noNAs, file="~/Data/hibiscus_biannual_w-timepts.csv")

