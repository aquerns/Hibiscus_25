
#Code: Matrix creation
#By: Aleah Querns
#Date: 9-13-23
#Purpose: Use vital rates created from IPM to 
#make matrices
# last modified: A Louthan 241219


# loading in packages and data----
rm(list = ls())
# packages
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
library(tibble)
load("Data/3_vital_rate_functions_and_vital_rate_values.RData")

###############################################################################

                        # Graph of Vital Rates #

###############################################################################

# Convert array to dataframe
df <- as.data.frame.table(full_vrs_allsites)

colnames(df) <- c("Size_Class", "Rate", "Site", "Treatment", "Value")

df$Size_Class <- gsub("F", "Fruit", df$Size_Class)
df$Size_Class <- gsub("S", "Small", df$Size_Class)
df$Size_Class <- gsub("M", "Medium", df$Size_Class)
df$Size_Class <- gsub("L", "Large", df$Size_Class)

colnames(df) <- c("Size Class", "Rate", "Site", "Treatment", "Value")


df$Rate <- gsub("sr", "Survival", df$Rate)
df$Rate <- gsub("fr", "Fruiting", df$Rate)
df$Rate <- gsub("transtoS", "Transition to Small", df$Rate)
df$Rate <- gsub("transtoM", "Transition to Medium", df$Rate)
df$Rate <- gsub("transtoL", "Transition to Large", df$Rate)
df$Rate <- gsub("germ", "Recruitment", df$Rate)

df$Site <- factor(df$Site, levels = c("N", "C", "S"))

df <- sort_by(df, df$Treatment, df$Site)

df<- df %>% filter(!is.na(Value))

df$Value <- format(df$Value, digits=3)

ggplot(data=df, aes(x=Site, y=Value))+
  geom_point(aes(col=df$`Size Class`))+
  geom_smooth(aes(group=df$`Size Class`, colour = df$'Size Class'), method="lm", se=F)+
  facet_grid(rows=vars(Rate), cols=vars(Treatment), scales="free")



###############################################################################

                         #Tables for Vital Rate Models#

###############################################################################

                        #Survival#
# Define site and treatment names
site_names <- c("C", "N", "S") #the OG naming went in a weird order
treatment_names <- c("CONT", "LMH", "MEGA", "MESO")

# Assign names to the top-level list (sites)
names(sur_mod) <- site_names

# Loop through each site and assign names to the second-level lists (treatments)
for (site in site_names) {
  names(sur_mod[[site]]) <- treatment_names
}

###############GR
# Assign names to the top-level list (sites)
names(gr_mod) <- site_names

# Loop through each site and assign names to the second-level lists (treatments)
for (site in site_names) {
  names(gr_mod[[site]]) <- treatment_names
}


###################VAR
# Assign names to the top-level list (sites)
names(var_mod) <- site_names

# Loop through each site and assign names to the second-level lists (treatments)
for (site in site_names) {
  names(var_mod[[site]]) <- treatment_names
}

####################nfr
# Assign names to the top-level list (sites)
names(nfr_mod) <- site_names

# Loop through each site and assign names to the second-level lists (treatments)
for (site in site_names) {
  names(nfr_mod[[site]]) <- treatment_names
}

####################pfr
# Assign names to the top-level list (sites)
names(pfr_mod) <- site_names

# Loop through each site and assign names to the second-level lists (treatments)
for (site in site_names) {
  names(pfr_mod[[site]]) <- treatment_names
}



#########################################################

# List of model lists
model_lists <- list(survival = sur_mod, 
                    growth = gr_mod, 
                    nfr = nfr_mod, 
                    pfr = pfr_mod, 
                    variance = var_mod)

# Initialize an empty list to store results
results_list <- list()

# Define all possible terms
terms <- c("BM", "block", "timepoint")

# Loop through each model type
for (model_name in names(model_lists)) {
  model_list <- model_lists[[model_name]]
  
  # Extract site and treatment names
  sites <- names(model_list)  
  treatments <- names(model_list[[1]])  
  
  # Loop through sites and treatments
  for (site in sites) {
    for (treatment in treatments) {
      model <- model_list[[site]][[treatment]]
      
      # Initialize fixed and random effect lists
      fixed_effects <- character(0)
      random_effects <- character(0)
      
      # Extract fixed and random effects based on model type
      if (inherits(model, "glmerMod") || inherits(model, "lmerMod")) {
        fixed_effects <- names(fixef(model))  # Fixed effects for mixed models
        random_effects <- names(ranef(model)) # Random effects for mixed models
      } else if (inherits(model, "glm") || inherits(model, "lm")) {
        fixed_effects <- names(coef(model))   # Fixed effects for glm/lm
        random_effects <- character(0)        # No random effects in glm/lm
      } else {
        next  # Skip if model type is unknown
      }
      
      # Create a named vector for presence/absence of each term
      term_presence <- setNames(rep("", length(terms)), terms)
      
      # Mark fixed effects with "F"
      term_presence[terms %in% fixed_effects] <- "F"
      
      # Mark random effects with "R"
      term_presence[terms %in% random_effects] <- "R"
      
      # Mark effects that appear in both fixed and random with "F+R"
      both_effects <- intersect(fixed_effects, random_effects)
      term_presence[both_effects] <- "F+R"
      
      # Store results as a dataframe row
      results_list[[paste(model_name, site, treatment, sep = "_")]] <- 
        c(Model = model_name, Site = site, Treatment = treatment, term_presence)
    }
  }
}

# Convert list to dataframe
results_df <- do.call(rbind, results_list)

# Convert to tibble for better formatting (optional)
results_df <- as_tibble(results_df)

# Print result
print(results_df)

library(flextable)
library(officer)
ft<- flextable(results_df)

# Export to Word document
doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_add_par("\n")  # Adds spacing after the table

print(doc, target = "models_table_output.docx")


#################################################################################

# first, get sensitiviites to vital rates for the mean matrix (averaged element-wise across N, C, and S)----
level <-c("N", "C", "S")
trmtUHURU <-c("CONT", "LMH", "MEGA", "MESO")
sizes <-c("FR", "S", "M", "L")
perturbs <-c('ambient',
            'high_sr1', 'low_sr1',
            'high_sr2', 'low_sr2',
            'high_sr3', 'low_sr3',
            'high_fr1', 'low_fr1',
            'high_fr2', 'low_fr2',
            'high_fr3', 'low_fr3',
            'high_S1', 'low_S1',
            'high_S2', 'low_S2',
            'high_S3', 'low_S3',
            'high_M1', 'low_M1',
            'high_M2', 'low_M2',
            'high_M3', 'low_M3',
            'high_L1', 'low_L1',
            'high_L2', 'low_L2',
            'high_L3', 'low_L3',
            'high_germ', 'low_germ')

#Array for storing matrices for each level/treatment/level of perturbation
matrixarray <- array(data=NA, dim=c(4,4, 4, 33),
                    dimnames = list(sizes,
                                    sizes,
                                    trmtUHURU,
                                    perturbs))

#create array for storing lambdas for each level/treatment/level of perturbation 
lamarray <- array(data=NA, dim=c(4, 33),
                  dimnames = list(
                                  trmtUHURU, #times
                                  perturbs
                  ))



for (p in 1:33){
  perturb_p = perturbs[p]
 
    for (j in 1:4){
      trmt_j = trmtUHURU[j]
      
    
  matrixpars<- full_vrs_perturb[,,trmt_j, perturb_p]
  #create filtered dataset for storage
  smallmat <- array(data=NA, dim=c(4,4), dimnames = list(sizes,
                                                         sizes)) 
smallmat[1,1] <- 0 #F-Fruit 
#technically fruits COULD contribute to F by the next survey if they mature to S, but this contribution would be so minimal we should keep it 0
smallmat[2,1] <- matrixpars[1,6] #fruit to small transition rate
  smallmat[3,1]<- 0 #fruit to medium transition rate
  smallmat[4,1]<-0 #fruit to large transition rate
  
smallmat[1,2] <- matrixpars[2,1]* matrixpars[2,2] #small to fruit transition rate: survive as a small, then fruit
  smallmat[2,2] <- matrixpars[2,1] * matrixpars[2,3] #S-S
  smallmat[3,2] <- matrixpars[2,1]* matrixpars[2,4] #S-M
  smallmat[4,2] <- matrixpars[2,1] * matrixpars[2,5] # S-L
  
   
smallmat[1,3] <- matrixpars[3,1]* matrixpars[3,2]# medium to fruit transition rate: survive as a medium, then fruit
  smallmat[2,3] <- matrixpars[3,1]* matrixpars[3,3] # M-S
  smallmat[3,3] <- matrixpars[3,1]*  matrixpars[3,4]#M-M
  smallmat[4,3] <- matrixpars[3,1]* matrixpars[3,5] #M-L
  
smallmat[1,4] <- matrixpars[4,1]* matrixpars[4,2] #large to fruit transition rate: survive as a large, then fruit
  smallmat[2,4] <- matrixpars[4,1]* matrixpars[4,3] #L-S
  smallmat[3,4] <- matrixpars[4,1]* matrixpars[4,4] #L-M
  smallmat[4,4] <- matrixpars[4,1]* matrixpars[4,5] #L-L


  #now throw the matrix into array by lvl, trt, timepoint, and perturbation level
       
  if(any(is.na(smallmat))){
    matrixarray[,,trmt_j, perturb_p] <- NA
  } else{
    matrixarray[,,trmt_j, perturb_p] <- smallmat
  }

  #NOW get lambda value for given lvl, trt, timepoint, perturbation
        
  if(any(is.na(smallmat))){
    lamarray[,j,p] <- NA
  } else{ x1<- eigen(smallmat)
        lmax <- which.max(Re(x1$values))
        lambda <- Re(x1$values[lmax])
        #Store store lambda at correct level, treatment, timepoint, perturbation level
        lamarray[j,p] <- lambda  }
 
      } # closes treatment loop
      
      
    } # closes perturbation loop
  
matrixarray_perturbed <- matrixarray
lamarray_perturbed <- lamarray
sens <-sweep(lamarray_perturbed[,-which(dimnames(lamarray_perturbed)[[2]] == "ambient")], 
             MARGIN= 1, STATS= lamarray_perturbed[, "ambient"], FUN= "-")/ # gets difference in lambdas btwn perturbed and ambient
# now, you need to divide by the change in the demographic rate, which is: 
  # perturbed rate, minus ambient rate
(
  cbind( # cbind command makes a matrix of perturbed rates
    full_vrs_perturb["S","sr" ,, 2],  # dims on the full_vrs_perturb is: sizes, rates, trmts, pertrubs2
  full_vrs_perturb["S","sr" ,, 3], 
full_vrs_perturb["M","sr" ,, 4], 
full_vrs_perturb["M","sr" ,, 5], 
full_vrs_perturb["L","sr" ,, 6], 
full_vrs_perturb["L","sr" ,, 7], 
full_vrs_perturb["S","fr" ,, 8], 
full_vrs_perturb["S","fr" ,, 9], 
full_vrs_perturb["M","fr" ,, 10], 
full_vrs_perturb["M","fr" ,, 11], 
full_vrs_perturb["L","fr" ,, 12], 
full_vrs_perturb["L","fr" ,, 13], 
full_vrs_perturb["S","transtoS" ,, 14], 
full_vrs_perturb["S","transtoS" ,, 15], 
full_vrs_perturb["M","transtoS" ,, 16], 
full_vrs_perturb["M","transtoS" ,, 17], 
full_vrs_perturb["L","transtoS" ,, 18], 
full_vrs_perturb["L","transtoS" ,, 19], 
full_vrs_perturb["S","transtoM" ,, 20], 
full_vrs_perturb["S","transtoM" ,, 21], 
full_vrs_perturb["M","transtoM" ,, 22], 
full_vrs_perturb["M","transtoM" ,, 23], 
full_vrs_perturb["L","transtoM" ,, 24], 
full_vrs_perturb["L","transtoM" ,, 25], 
full_vrs_perturb["S","transtoL" ,, 26], 
full_vrs_perturb["S","transtoL" ,, 27], 
full_vrs_perturb["M","transtoL" ,, 28], 
full_vrs_perturb["M","transtoL" ,, 29], 
full_vrs_perturb["L","transtoL" ,, 30], 
full_vrs_perturb["L","transtoL" ,, 31],
full_vrs_perturb["FR","germ" ,, 32],
full_vrs_perturb["FR","germ" ,, 33])-  
# then you minus the ambient rate
cbind( # this cbind command makes a matrix of ambient rates
  full_vrs_averaged["S","sr",], full_vrs_averaged["S","sr",], 
full_vrs_averaged["M","sr",], full_vrs_averaged["M","sr",], 
full_vrs_averaged["L","sr",], full_vrs_averaged["L","sr",], 
full_vrs_averaged["S","fr",], full_vrs_averaged["S","fr",], 
full_vrs_averaged["M","fr",], full_vrs_averaged["M","fr",], 
full_vrs_averaged["L","fr",], full_vrs_averaged["L","fr",], 
full_vrs_averaged["S","transtoS",], full_vrs_averaged["S","transtoS",], 
full_vrs_averaged["M","transtoS" ,], full_vrs_averaged["M","transtoS" ,],
full_vrs_averaged["L","transtoS",], full_vrs_averaged["L","transtoS",],
full_vrs_averaged["S","transtoM",], full_vrs_averaged["S","transtoM",], 
full_vrs_averaged["M","transtoM", ], full_vrs_averaged["M","transtoM", ], 
full_vrs_averaged["L","transtoM",], full_vrs_averaged["L","transtoM",], 
full_vrs_averaged["S","transtoL", ], full_vrs_averaged["S","transtoL", ], 
full_vrs_averaged["M","transtoL",], full_vrs_averaged["M","transtoL",], 
full_vrs_averaged["L","transtoL",],full_vrs_averaged["L","transtoL",],
full_vrs_averaged["F","germ", ], full_vrs_averaged["F","germ", ]))

# finally, average the +/ - perturbations
sens <- cbind(
  (sens[,1] + sens[,2])/2, 
  (sens[,3] + sens[,4])/2,
  (sens[,5] + sens[,6])/2,
  (sens[,7] + sens[,8])/2,
  (sens[,9] + sens[,10])/2,
  (sens[,11] + sens[,12])/2,
  (sens[,13] + sens[,14])/2,
  (sens[,15] + sens[,16])/2,
  (sens[,17] + sens[,18])/2,
  (sens[,19] + sens[,20])/2,
  (sens[,21] + sens[,22])/2,
  (sens[,23] + sens[,24])/2,
  (sens[,25] + sens[,26])/2,
  (sens[,27] + sens[,28])/2,
  (sens[,29] + sens[,30])/2,
  (sens[,31] + sens[,32])/2
) # dimensions on sensitivities should be demographic rates (in this order: small sur, med sur, lar sur, small fr, med fr, lar fr, S-S, M-S, L-S, M-S, M-M, M-L, L-S, L-M, L-L, germ), trmts

# calculate mean and CI on population growth rate using the randomly drawn coefficients-----

level<-c("N", "C", "S")
trmtUHURU<-c("CONT", "LMH", "MEGA", "MESO")
sizes <-c("FR", "S", "M", "L")
iterations <- 10000
it_vec <-1:10000

#create array for storing lambdas 
lamarray <- array(data=NA, dim=c(3, 4,  iterations),
                  dimnames = list(level,
                                  trmtUHURU,
                                  it_vec
                  ))

#holder matrix for for loop
smallmat <- array(data=NA, dim=c(4,4), dimnames = list(sizes,
                                                       sizes))
for(x in it_vec){
  iteration_x = it_vec[x]
  for (i in 1:3){
    level_i = level[i]
    for (j in 1:4){
      trmt_j = trmtUHURU[j]
      
      matrixpars<- full_vrs_random[,,level_i, trmt_j,  iteration_x]
      #create filtered dataset for storage
      
      smallmat[1,1] <- 0 #F-Fruit
      smallmat[2,1] <- matrixpars[1,6] #F-S
      smallmat[3,1]<- 0 #F-M
      smallmat[4,1]<-0 #F-L
      
      smallmat[1,2] <- matrixpars[2,1]* matrixpars[2,2] #S-Fruits
      smallmat[2,2] <- matrixpars[2,1] * matrixpars[2,3] #S-S
      smallmat[3,2] <- matrixpars[2,1]* matrixpars[2,4] #S-M
      smallmat[4,2] <- matrixpars[2,1] * matrixpars[2,5] # S-L
      
      
      smallmat[1,3] <- matrixpars[3,1]* matrixpars[3,2]# M-Fruits
      smallmat[2,3] <- matrixpars[3,1]* matrixpars[3,3] # M-S
      smallmat[3,3] <- matrixpars[3,1]*  matrixpars[3,4]#M-M
      smallmat[4,3] <- matrixpars[3,1]* matrixpars[3,5] #M-L
      
      smallmat[1,4] <- matrixpars[4,1]* matrixpars[4,2] #L-Fruit
      smallmat[2,4] <- matrixpars[4,1]* matrixpars[4,3] #L-S
      smallmat[3,4] <- matrixpars[4,1]* matrixpars[4,4] #L-M
      smallmat[4,4] <- matrixpars[4,1]* matrixpars[4,5] #L-L
      
      
      #NOW get lambda value for given lvl, trt, perturbation
      x1<- eigen(smallmat)
      lmax <- which.max(Re(x1$values))
      lambda <- Re(x1$values[lmax])
      #Store store lambda at correct level, treatment, iteration
      lamarray[i,j, x] <- lambda  
    }
    
}}

lamarray_random <- lamarray
lams_random_CIs <- apply(lamarray_random, c(1, 2), quantile, c(0.05/2, 0.5, (1-0.05/2))) # this command gives you the median and CI's on lambdas

# Convert array to a tidy dataframe
df <- as_tibble(as.data.frame.table(lams_random_CIs)) %>%
  rename(Quantile = Var1, Site = Var2, Treatment = Var3, Value = Freq) 

# Reshape to get CI columns
df_wide <- df %>%
  pivot_wider(names_from = Quantile, values_from = Value) %>%
  rename(Lower = `2.5%`, Mean = `50%`, Upper = `97.5%`)

# Print first few rows to check
head(df_wide)

#convert--North is arid, south is mesic
df_wide <- df_wide %>%
  mutate(Site = recode(Site, 
                       "N" = "Arid", 
                       "C" = "Intermediate", 
                       "S" = "Mesic"))
######## Graph of Population Growth Rates

ggplot(df_wide, aes(x = Site, y = Mean, color = Treatment, group = Treatment)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.3)) +
  labs(x = "Rainfall Level", y = expression("Population Growth Rate   (  " *lambda*")"), color = "Treatment") +
  theme_bw()+
  scale_color_manual(values=c("darkorange", "deepskyblue3", "firebrick", "forestgreen"))

###############################################################################

                 # across-population variability in lambda?----

###############################################################################

var_in_lam <- apply(
  apply(log(lamarray_random), c(2, 3),  var) , # gets across-level variance in lambda; dims are trmt, iteration
1, quantile, c(0.05/2, 0.5, (1-(0.05/2)))) # calculate CI across iternations
  
acrosssitevar <- as_tibble(as.data.frame.table(var_in_lam)) %>%
  rename(Quantile = Var1, Treatment = Var2, Value = Freq) 

# Reshape to get CI columns
acrosssitevar_wide <- acrosssitevar %>%
  pivot_wider(names_from = Quantile, values_from = Value) %>%
  rename(Lower = `2.5%`, Mean = `50%`, Upper = `97.5%`)



jpeg(filename = "varandci.jpg", width = 5, height = 4, units = "in", res = 300)
# 
#all treatments
ggplot(acrosssitevar_wide, aes(x=Treatment, y=Mean))+
  geom_point(aes(col=Treatment), size=4,  position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, col=Treatment), width=.2, linewidth=1, position=position_dodge(.9))+
  # geom_text(aes(label=Significance, y=UpperCI + 0.1), vjust=0) +
  theme_bw()+
  labs(x="Herbivore Exclusion Treatment", y=expression("Variation in Population Growth Rate (  "*lambda* ")")) +
  theme(legend.position="none")+
  scale_color_manual(values=c("darkorange", "deepskyblue3", "firebrick", "forestgreen"))

dev.off()

save(sens, file= "Data/4_sensitivity values.RData")
############ Table of sensitivity
colnames(sens) <- c("S_Survival", "M_Survival", "L_Survival",
                   "S_Fruiting", "M_Fruiting", "L_Fruiting",
                   "S_TranstoS", "M_TranstoS", "L_TranstoS", 
                   "S_TranstoM", "M_TranstoM", "L_TranstoM",
                   "S_TranstoL", "M_TranstoL", "L_TranstoL",
                   "Recruitment")
sdf <- as.data.frame.table(sens)

sdf_wide <- sdf %>%
  pivot_wider(names_from = Var1, values_from = Freq)

sdfft<-flextable (sdf_wide)
# Export to Word document
doc <- read_docx() %>%
  body_add_flextable(sdfft) %>%
  body_add_par("\n")  # Adds spacing after the table

print(doc, target = "sens_table_output.docx")


#############################################################################

 #What about confidence interval for pairwise differences? #

###########################################################################

var <- apply(lamarray_random, c(2, 3),  var) 
# Assuming your 'var' dataset has treatments as row names and each column is a replicate of variance
treatments <- rownames(var)

# Initialize list to store results
ci_results <- list()

# Iterate over pairs of treatments
for (i in 1:(length(treatments) - 1)) {
  for (j in (i + 1):length(treatments)) {
    
    # Extract the variance values for the pair of treatments
    var_t1 <- var[i,]
    var_t2 <- var[j,]
    
    # Calculate pairwise difference in variance
    diff_variance <- var_t1 - var_t2
    
    # Calculate the 95% confidence interval (2.5th and 97.5th percentiles)
    ci <- quantile(diff_variance, probs = c(0.025, 0.5, 0.975))
    
    # Store the result in the list
    ci_results[[paste(treatments[i], "vs", treatments[j])]] <- ci
  }
  }
  

# View the results
ci_results



################################################

# Create a data frame for plotting
results <- data.frame(
  comparison = character(),  # Pairwise comparisons
  median_diff = numeric(),   # Median difference in variance
  lower_ci = numeric(),      # Lower bound of the confidence interval
  upper_ci = numeric(),      # Upper bound of the confidence interval
  significance = character() # Significance indicator
)

# Iterate over pairwise comparisons and populate the results data frame
for (i in 1:(length(treatments) - 1)) {
  for (j in (i + 1):length(treatments)) {
    
    # Extract pairwise difference and confidence intervals
    var_t1 <- var[i, ]
    var_t2 <- var[j, ]
    diff_variance <- var_t1 - var_t2
    ci <- quantile(diff_variance, probs = c(0.025, 0.975))
    
    # Calculate the median difference (50th percentile)
    median_diff <- median(diff_variance)
    
    # Check if zero is within the confidence interval (to determine significance)
    significance <- ifelse(ci[1] > 0 | ci[2] < 0, "*", "ns")  # "*" for significant, "ns" for not significant
    
    # Store the results
    results <- rbind(results, data.frame(
      comparison = paste(treatments[i], "vs", treatments[j]),
      median_diff = median_diff,
      lower_ci = ci[1],
      upper_ci = ci[2],
      significance = significance
    ))
  }
}

# Plot the results

ggplot(results, aes(x = comparison, y = median_diff, color = significance)) +
  geom_point(size = 4) +  # Add points for each pairwise difference
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +  # Add error bars
  scale_color_manual(values = c("ns" = "black", "*" = "red")) +  # Color points based on significance
  labs(
    title = "Pairwise Differences in Variance Between Treatments",
    y = "Median Difference in Variance",
    x = "Pairwise Comparison"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

