###############################################################################
#########                                                        #############
#########             CORRELATIONS HEAT MAP                     ##############
###############################################################################

rm(list = ls())
library(gridExtra)
library(reshape2)
library(tidyverse)

#load in data
load("Data/6_test_for_comp_COMPRESS.RData.xz")
#############################################################################
##########        ALL HERBIVORES PRESENT                               ######
##############################################################################

###  Make data tables  ###

cors_pos<-as.data.frame.table(dfstore_pos, na.rm=F) #corrs are same across tests so technically only need to grab one df
# Rename columns for clarity
colnames(cors_pos) <- c("Treatment", "Rate1_1", "Rate1_2", "Rate2_1", "Rate2_2", "Metric", "Value")
# Create Rate1 and Rate2 by combining relevant columns
cors_pos$Rate1 <- paste(cors_pos$Rate1_1, cors_pos$Rate1_2, sep = "_")
cors_pos$Rate2 <- paste(cors_pos$Rate2_1, cors_pos$Rate2_2, sep = "_")

cors_neg<-as.data.frame.table(dfstore_neg, na.rm=F) #corrs are same across tests so technically only need to grab one df
# Rename columns for clarity
colnames(cors_neg) <- c("Treatment", "Rate1_1", "Rate1_2", "Rate2_1", "Rate2_2", "Metric", "Value")
# Create Rate1 and Rate2 by combining relevant columns
cors_neg$Rate1 <- paste(cors_neg$Rate1_1, cors_neg$Rate1_2, sep = "_")
cors_neg$Rate2 <- paste(cors_neg$Rate2_1, cors_neg$Rate2_2, sep = "_")

#Grab the relevant data depending on correlation
cors_pos<-cors_pos[,-c(2:5)]
cors_neg<-cors_neg[,-c(2:5)]

cors_wide_pos <- pivot_wider(cors_pos, names_from = Metric, values_from = Value)
cors_wide_neg <- pivot_wider(cors_neg, names_from = Metric, values_from = Value)

#create a holder to store actual P values
cors_holder <- cors_wide_pos
cors_holder$P <-NA
#grab p value from relevant correlation (pos or neg)
for (i in 1:nrow(cors_holder)){

  if (is.na(cors_holder$corr[i])){
    cors_holder$P[i] <- NA
  } else if(cors_holder$corr[i] > 0){
    cors_holder$P[i] <- cors_wide_pos$P[which(
      cors_wide_pos$Treatment== cors_holder$Treatment[i] &
        cors_wide_pos$Rate1== cors_holder$Rate1[i] &
        cors_wide_pos$Rate2== cors_holder$Rate2[i])]
  } else if (cors_holder$corr[i] < 0){
    cors_holder$P[i] <- cors_wide_neg$P[which(
      cors_wide_neg$Treatment== cors_holder$Treatment[i] &
        cors_wide_neg$Rate1== cors_holder$Rate1[i] &
        cors_wide_neg$Rate2== cors_holder$Rate2[i])]
  }
    }

df<-cors_holder    

#Make names prettier
df$Rate1<-gsub("fr", "FruitingRate", df$Rate1)
df$Rate1<-gsub("sr", "SurvivalRate", df$Rate1)
df$Rate1<-gsub("germ", "Recruitment", df$Rate1)
df$Rate1<-gsub("transtoS", "TransitiontoS", df$Rate1)
df$Rate1<-gsub("transtoM", "TransitiontoM", df$Rate1)
df$Rate1<-gsub("transtoL", "TransitiontoL", df$Rate1)

df$Rate2<-gsub("fr", "FruitingRate", df$Rate2)
df$Rate2<-gsub("sr", "SurvivalRate", df$Rate2)
df$Rate2<-gsub("germ", "Recruitment", df$Rate2)
df$Rate2<-gsub("transtoS", "TransitiontoS", df$Rate2)
df$Rate2<-gsub("transtoM", "TransitiontoM", df$Rate2)
df$Rate2<-gsub("transtoL", "TransitiontoL", df$Rate2)

df$Rate1<-gsub("_", ".", df$Rate1)
df$Rate2<-gsub("_", ".", df$Rate2)


#filter out vital rates that aren't possible but that we have placeholders for b.c of array format

gcors_filt<-df %>%
  filter(Rate1 !="S.Recruitment")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="M.Recruitment")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="L.Recruitment") 
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="S.Recruitment")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="M.Recruitment")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="L.Recruitment")

gcors_filt<-gcors_filt %>%
  filter(Rate1 !="FR.FruitingRate")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="FR.SurvivalRate")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="FR.TransitiontoS")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="FR.TransitiontoM")
gcors_filt<-gcors_filt %>%
  filter(Rate1 !="FR.TransitiontoL")

gcors_filt<-gcors_filt %>%
  filter(Rate2 !="FR.FruitingRate")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="FR.SurvivalRate")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="FR.TransitiontoS")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="FR.TransitiontoM")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="FR.TransitiontoL")

#Filter out any rates that had LTRE=0 b/c doesn't vary across N/C/S (cannot test corr)

gcors_filt<-gcors_filt %>%
  filter(Rate1 !="S.TransitiontoL")
gcors_filt<-gcors_filt %>%
  filter(Rate2 !="S.TransitiontoL")

#Make factors (give order for axes)
gcors_filt$Rate2 <- factor(gcors_filt$Rate2, levels=c(
  "S.SurvivalRate", 
  "S.TransitiontoS","S.TransitiontoM",
  "S.FruitingRate",
   "M.SurvivalRate",
  "M.TransitiontoS" ,"M.TransitiontoM","M.TransitiontoL",
  "M.FruitingRate",
  "L.SurvivalRate",
  "L.TransitiontoS", "L.TransitiontoM","L.TransitiontoL",
  "L.FruitingRate" ,
  "FR.Recruitment" ))
#Create names for axes
VEC2<- c(
  "S-Survival",
  "S-Stasis", "S-M Growth",
         "S-Fruiting Rate",
  "M-Survival",
  "M-S Regression", "M-Stasis","M-L Growth",
         "M-Fruiting Rate",
  "L-Survival",
  "L-S Regression", "L-M Regression", "L-Stasis",
        "L-Fruiting Rate", "Recruitment")
#Make factors (give order for axes)
gcors_filt$Rate1 <- factor(gcors_filt$Rate1, levels=c(
  "S.SurvivalRate", 
  "S.TransitiontoS","S.TransitiontoM",
  "S.FruitingRate",
  "M.SurvivalRate",
  "M.TransitiontoS" ,"M.TransitiontoM","M.TransitiontoL",
  "M.FruitingRate",
  "L.SurvivalRate",
  "L.TransitiontoS", "L.TransitiontoM","L.TransitiontoL",
  "L.FruitingRate" ,
  "FR.Recruitment" ))
#Create names for axes
VEC1<-  c(
  "S-Survival",
  "S-Stasis", "S-M Growth",
  "S-Fruiting Rate",
  "M-Survival",
  "M-S Regression", "M-Stasis","M-L Growth",
  "M-Fruiting Rate",
  "L-Survival",
  "L-S Regression", "L-M Regression", "L-Stasis",
  "L-Fruiting Rate", "Recruitment")

#sort replicates so we only have lower triangle
gconcat <- gcors_filt %>%
  rowwise() %>%
  mutate(pair = sort(c(Rate1, Rate2)) %>% paste(collapse = ",")) %>%
  group_by(Treatment, pair) %>%
  arrange(Rate1)%>%
  distinct(pair, .keep_all = T)


################################################################################

                 # Below: Old graph (don't use) #

################################################################################
# ##Graphing###
# A<-ggplot(gconcat %>% filter (Treatment == "CONT"), aes(x=Rate2, y=Rate1,fill=corr)) +
#   geom_tile()+
#   scale_fill_gradientn(colours = colorspace::diverge_hcl(7), na.value = "white")+
#   theme_bw(base_size=10)+
#   #add border white colour of line thickness 0.25
#   geom_tile(colour="white", size=0.25)+
#   scale_y_discrete(expand=c(0, 0), labels=VEC2)+
#   #define new breaks on x-axis
#   scale_x_discrete(expand=c(0, 0), position="bottom", labels=VEC1)+
#   #remove x and y axis labels
#   labs(x="", y="")+
#   #remove extra space
#   #theme options
#   theme(text = element_text(size=20),
#         legend.title = element_text(face = "bold", size=15),
#         #bold font for legend text
#         legend.text=element_text(size=15),
#         #set thickness of axis ticks
#         axis.ticks=element_line(size=0.4),
#         #remove plot background
#         plot.background=element_blank(),
#         panel.grid.major = element_blank(),#remove grid
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))+ 
#   geom_point(data = filter(gconcat, P<0.05,!is.na(corr), Treatment =="CONT"), color = "white", position = position_nudge(x = -0.095), shape=8)
# 
# B<-ggplot(gconcat %>% filter (Treatment == "LMH"), aes(x=Rate2, y=Rate1,fill=corr)) +
#   geom_tile()+
#   scale_fill_gradientn(colours = colorspace::diverge_hcl(7), na.value = "white")+
#   theme_bw(base_size=10)+
#   #add border white colour of line thickness 0.25
#   geom_tile(colour="white", size=0.25)+
#   scale_y_discrete(expand=c(0, 0), labels=VEC2)+
#   #define new breaks on x-axis
#   scale_x_discrete(expand=c(0, 0), position="bottom", labels=VEC1)+
#   #remove x and y axis labels
#   labs(x="", y="")+
#   #remove extra space
#   #theme options
#   theme(text = element_text(size=20),
#         legend.title = element_text(face = "bold", size=15),
#         #bold font for legend text
#         legend.text=element_text(size=15),
#         #set thickness of axis ticks
#         axis.ticks=element_line(size=0.4),
#         #remove plot background
#         plot.background=element_blank(),
#         panel.grid.major = element_blank(),#remove grid
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))+ 
#   geom_point(data = filter(gconcat, P<0.05,!is.na(corr), Treatment =="LMH"), color = "white", position = position_nudge(x = -0.095), shape=8)
# 
# C<-ggplot(gconcat %>% filter (Treatment == "MEGA"), aes(x=Rate2, y=Rate1,fill=corr)) +
#   geom_tile()+
#   scale_fill_gradientn(colours = colorspace::diverge_hcl(7), na.value = "white")+
#   theme_bw(base_size=10)+
#   #add border white colour of line thickness 0.25
#   geom_tile(colour="white", size=0.25)+
#   scale_y_discrete(expand=c(0, 0), labels=VEC2)+
#   #define new breaks on x-axis
#   scale_x_discrete(expand=c(0, 0), position="bottom", labels=VEC1)+
#   #remove x and y axis labels
#   labs(x="", y="")+
#   #remove extra space
#   #theme options
#   theme(text = element_text(size=20),
#         legend.title = element_text(face = "bold", size=15),
#         #bold font for legend text
#         legend.text=element_text(size=15),
#         #set thickness of axis ticks
#         axis.ticks=element_line(size=0.4),
#         #remove plot background
#         plot.background=element_blank(),
#         panel.grid.major = element_blank(),#remove grid
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))+ 
#   geom_point(data = filter(gconcat, P<0.05,!is.na(corr), Treatment =="MEGA"), color = "white", position = position_nudge(x = -0.095), shape=8)
# 
# D<-ggplot(gconcat %>% filter (Treatment == "MESO"), aes(x=Rate2, y=Rate1,fill=corr)) +
#   geom_tile()+
#   scale_fill_gradientn(colours = colorspace::diverge_hcl(7), na.value = "white")+
#   theme_bw(base_size=10)+
#   #add border white colour of line thickness 0.25
#   geom_tile(colour="white", size=0.25)+
#   scale_y_discrete(expand=c(0, 0), labels=VEC2)+
#   #define new breaks on x-axis
#   scale_x_discrete(expand=c(0, 0), position="bottom", labels=VEC1)+
#   #remove x and y axis labels
#   labs(x="", y="")+
#   #remove extra space
#   #theme options
#   theme(text = element_text(size=20),
#         legend.title = element_text(face = "bold", size=15),
#         #bold font for legend text
#         legend.text=element_text(size=15),
#         #set thickness of axis ticks
#         axis.ticks=element_line(size=0.4),
#         #remove plot background
#         plot.background=element_blank(),
#         panel.grid.major = element_blank(),#remove grid
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))+ 
#   geom_point(data = filter(gconcat, P<0.05,!is.na(corr), Treatment =="MESO"), color = "white", position = position_nudge(x = -0.095), shape=8)
# 
# 
# 
# 
# ###Put together plots
# 
# grid.arrange(A,B,C,D, ncol=2, nrow=2)
# #note that this is super messy because they each have their own axes--moving on....


#################################################################################
 
          # Put together plots by facet-wrapping. Better for axes #

#################################################################################
gconcat$Treatment <- factor(gconcat$Treatment, levels = c("LMH", "MESO", "MEGA", "CONT"))

jpeg(filename = "25_2_20_heatmap.jpg", width = 6, height = 5.5, units = "in", res = 300)

#1500x1500 - heatmap_1s_ipm_alltrts
#ALL HERBIVORE TRTS
ggplot(gconcat, aes(x=Rate2, y=Rate1,fill=corr)) +
  geom_tile()+
  scale_fill_gradientn(colours = colorspace::diverge_hcl(7), na.value = "white",
                       breaks=c(-1.1,-0.5, 0, 0.5, 1.1), labels= c(-1,-0.5, 0,0.5, 1))+
  theme_bw()+
  #add border white colour of line thickness 0.25
  geom_tile(colour="white", size=0.25)+
  #remove x and y axis labels
  labs(x="", y="", fill = "Correlation")+
  #remove extra space
  scale_y_discrete(expand=c(0, 0), labels= VEC2)+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0, 0), position="bottom", labels= VEC1)+
  #theme options
  theme(text = element_text(size=14, face="bold"),
        legend.title = element_text(face = "bold", size=12),
        #bold font for legend text
        legend.text=element_text(size=12),
        #set thickness of axis ticks
        axis.ticks=element_line(size=0.4),
        #remove plot background
        plot.background=element_blank(),
        panel.grid.major = element_blank(),#remove grid
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2, size=8, face="plain"),
        axis.text.y = element_text(size=8, face="plain"),
        aspect.ratio = 1) +
  geom_point(data = filter(gconcat, P<0.05,!is.na(corr)), color = "white", position = position_nudge(x = -0.095), shape=8)+
  facet_wrap(~Treatment, nrow=2, ncol=2, labeller=labeller(Treatment=c(LMH="A. LMH",
                                                             MESO="B. MESO",
                                                             MEGA= "C. MEGA",
                                                             CONT= "D. Control")
                                                      ))


dev.off()

##############################################################################

                          # Count Corrs for Processes #
                     #          Figure out what's driving   #

###############################################################################

# Define a function to classify rates
classify_rate <- function(rate) {
  case_when(
    rate %in% c("S.SurvivalRate", "M.SurvivalRate", "L.SurvivalRate") ~ "Survival",
    rate %in% c("S.TransitiontoM", "S.TransitiontoL", "M.TransitiontoL") ~ "Growth",
    rate %in% c("M.TransitiontoS", "L.TransitiontoM", "L.TransitiontoS") ~ "Regression",
    rate %in% c("S.TransitiontoS", "M.TransitiontoM", "L.TransitiontoL") ~ "Stasis",
    rate %in% c("S.FruitingRate", "M.FruitingRate", "L.FruitingRate") ~ "Fruiting",
    rate == "FR.Recruitment" ~ "Recruitment",
    TRUE ~ NA_character_  # Assign NA if it doesn't match any category
  )
}

# Apply the function to Rate1 and Rate2
gconcat <- gconcat %>%
  mutate(
    Process1 = classify_rate(Rate1),
    Process2 = classify_rate(Rate2)
  )

negsigonly<-gconcat %>% filter(P<0.05 & corr <0)


negsig_summary <- negsigonly %>%
  mutate(
    ProcessCombo = ifelse(Process1 < Process2, #Assure that process 1 is always first in combo
                          paste(Process1, Process2, sep = " - "), 
                          paste(Process2, Process1, sep = " - "))
  ) %>%
  group_by(Treatment, ProcessCombo) %>%
  summarise(N = n(), .groups = "drop")  # Sum counts for grouped combinations

# View the result
print(negsig_summary)


ggplot(negsig_summary, aes(x = reorder(ProcessCombo, N), y = N, fill = ProcessCombo)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Treatment, scales = "free_y") +  # Separate panels for each Treatment
  coord_flip() +  # Flip to make labels readable
  theme_minimal() +
  labs(
    title = "Significant Negative Correlations Between Process Combinations",
    x = "Process Combination",
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
