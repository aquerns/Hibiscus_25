# project: Hibiscus demographic compensation
# This script cleans data, then provides size and time classes
# written: September 5, 2023 (Louthan)
# last modified: A Louthan 241220

#Code below should reflect RV ~ BM + timepoint (R) + block (R)

# modified: by Aleah
# packages & functions----
rm(list = ls())
library(reshape2)
library(MuMIn)
library(lme4)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(abind)
library(numDeriv)
library(DHARMa)
library(glmmTMB)
library(fitdistrplus)
library(logspline)
library(texreg)
library(ggpubr)
library(arm)
library(boot)
library(tidyr)
library(ggfortify)
library(kableExtra)

model.select <- "bestfit" #"global" 
# model.select should be either "bestfit" (i.e., use dredge to select the best fit coeff) or "global" (i.e., use the global model to predict vital rates)

predict_custom <-  function(model_function, binmids_function ){ 
  # this function takes a glm with block effects and averages over the block effects to get one value for the
  # size-specific vital rate values
  forpredict.glm <- as.data.frame((binmids_function))
  names(forpredict.glm) <- "BM"
  
  if (class(model_function)[1] %in% c("glmerMod", "lmerMod")) {
    if (!is.null(unique(model_function@frame$season)[1])){ # if the model contains timepoint
      forpredict.glm <- forpredict.glm %>% crossing(season=unique(model_function@frame$season)) }
    forpredict.glm$response <- predict(model_function, newdata= forpredict.glm, type= "response", re.form= NA)
    forpredict.glm <- forpredict.glm %>% group_by(BM) %>% summarise(newresponse= mean(response))
    
  } else {
    
    if (!is.null(unique(model_function$model$block)[1])){ # if the model contains block
      forpredict.glm <- forpredict.glm %>% crossing(block=unique(model_function$model$block)) } 
    if (!is.null(unique(model_function$model$timepoint)[1])){ # if the model contains timepoint
      forpredict.glm <- forpredict.glm %>% crossing(timepoint=unique(model_function$model$timepoint)) } 
    if (!is.null(unique(model_function$model$season)[1])){ # if the model contains season
      forpredict.glm <- forpredict.glm %>% crossing(season=unique(model_function$model$season)) } 
    
    forpredict.glm$response <- predict(model_function, newdata= forpredict.glm, type= "response", 
                                       re.form= NA) # because this function predicts block-specific levels let's average over the blocks
    if (!is.null(unique(model_function$model$block)[1])) {
      options(dplyr.summarise.inform = FALSE)
      forpredict.glm <- forpredict.glm %>% group_by(BM, block) %>% summarise(response= mean(response))# this works with both timepoint or season, or if timepoint & season are missing 
      options(dplyr.summarise.inform = TRUE)}
    forpredict.glm <- forpredict.glm %>% group_by(BM) %>% summarise(newresponse= mean(response))
  }
  return(forpredict.glm$newresponse)
}


# reading in data & getting formats right-----
hib_data <- read.csv("Data/hibiscus_biannual_w-timepts.csv")
hib_data$block <- factor(paste(hib_data$level, hib_data$block, sep= ""))
hib_data$level <- factor(hib_data$level, levels = c("N", "C", "S"))
hib_data$trmtUHURU<- factor(hib_data$trmtUHURU)
hib_data$timepoint<- factor(hib_data$timepoint)
hib_data$season<- factor(hib_data$season)

#get rid of unnecessary columns
hib_data <- hib_data[,!(names(hib_data) %in% c("trmtAL", "x", "y", "number", "year", "wrongspecies", "date", "measurer", "BA", "height", "woodystems",
                                               "insectherb", "herbstems", "wcover", "fcover","gcover","PA", "Arfrtsapp", "htapp", "AcinNR", "distancetoAc",
                                               "Acsp", "sclass1","sclass2", "remnotdone","BAnext6mo","heightnext6mo", "woodystemsnext6mo", "Arfrtsappnext6mo","htappnext6mo",
                                               "surnext6mo", "IDuncertain", "discardprevtransition"))]

#Set size class binmids and edges

binedges <- seq(0, max(c(hib_data$BM, hib_data$BMnext6mo), na.rm=TRUE), length.out=4) 
binmids <- binedges[1:(length(binedges)-1)]+ ((binedges[2:length(binedges)] - binedges[1:(length(binedges)-1)])/2)

#ID unique trt, levels, timepoints
trmts <- unique(hib_data$trmtUHURU)
levels <- unique(hib_data$level)
timepoints <- unique(hib_data$timepoint)


################################################################################

                            #Survival#

################################################################################

# subset data to those observations that have no NAs for objects in model
sur_subset <- hib_data %>%
  filter(!is.na(binsurnext6mo) & 
           !is.na(BM) & 
           !is.na(level) & 
           !is.na(trmtUHURU) & 
           !is.na(season) & 
           !is.na(block))

#holder for models
sur_mod <- list(list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA))


#Make models for each level and treatment
for (i in 1:length(levels)){
  for (ii in 1:length(trmts)){
    sur_subset_i_ii <- sur_subset[which(sur_subset$level == levels[i] & sur_subset$trmt == trmts[ii]),]
    
    sur_mod[[i]][[ii]] <- lme4::glmer(binsurnext6mo~ BM +
                                        (1|timepoint) +
                                        +  (1|block)
                                      , #include both season and block (which is specific to a certain level/trt) as random effects
                                      family= "binomial", 
                                      data= sur_subset_i_ii, 
                                      na.action="na.fail")# No issues with convergence
    if (model.select == "bestfit") {
      sur_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(sur_mod[[i]][[ii]], REML=F),1)[[1]]}
    
    if(isSingular(sur_mod[[i]][[ii]])) { sur_mod[[i]][[ii]] <-glm(binsurnext6mo~ BM + block +
                                                                    timepoint , 
                                                                  family= "binomial", 
                                                                  data= sur_subset_i_ii, 
                                                                  na.action="na.fail") 
    if (model.select == "bestfit") {
      sur_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(sur_mod[[i]][[ii]]),1)[[1]]}
    }
  }} # can ignore the singular warnings; they are taken care of within the loop

################################################################################

                                    #Growth rate#

################################################################################
gr_subset<- hib_data %>% 
  filter(!is.na(BMnext6mo) & 
           !is.na(BM) & 
           !is.na(level) & 
           !is.na(trmtUHURU) & 
           !is.na(season) & 
           !is.na(block))%>%
  mutate(diff_BM = BMnext6mo- BM,
         sqrtdiffbm = sqrt(diff_BM),
         logbm= log(BM),
         sqrtbm =sqrt(BM))

#holder for models
gr_mod <- list(list(NA, NA, NA, NA), 
               list(NA, NA, NA, NA), 
               list(NA, NA, NA, NA))

for (i in 1:length(levels)){
  for (ii in 1:length(trmts)){
    gr_subset_i_ii <- gr_subset[which(gr_subset$level == levels[i] & gr_subset$trmt == trmts[ii]),]
    
    gr_mod[[i]][[ii]] <- lme4::lmer( diff_BM~ BM  +
                                       (1|timepoint) +
                                       (1|block),
                                     data= gr_subset_i_ii, 
                                     REML=F,
                                     na.action="na.fail")
    if (model.select == "bestfit") {
      gr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(gr_mod[[i]][[ii]], REML=F),1)[[1]]}
    gr_mod[[i]][[ii]] <- update(gr_mod[[i]][[ii]], REML=T)
    
    if(isSingular(gr_mod[[i]][[ii]])) { gr_mod[[i]][[ii]] <-lm(diff_BM~ BM  + block +
                                                                 timepoint , 
                                                               data= gr_subset_i_ii, 
                                                               na.action="na.fail") 
    if (model.select == "bestfit") {
      gr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(gr_mod[[i]][[ii]]),1)[[1]]}
    }
  }}# can ignore the singular warnings; they are taken care of within the loop


################################################################################

                            #Variance in Growth#

################################################################################

var_subset <- 
  hib_data %>% 
  filter(
    !is.na(BM) & 
      !is.na(level) & 
      !is.na(trmtUHURU) & 
      !is.na(season) & 
      !is.na(block)) %>%
  mutate(sqrtbm = sqrt(BM))

#holder
var_mod <- list(list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA))


for (i in 1:length(levels)){
  for (ii in 1:length(trmts)){
    var_subset_i_ii <- var_subset[which(var_subset$level == levels[i] & var_subset$trmt == trmts[ii] & 
                                          !is.na(var_subset$BM) &
                                          !is.na(var_subset$BMnext6mo) &
                                          !is.na(var_subset$level) & 
                                          !is.na(var_subset$trmtUHURU) &
                                          !is.na(var_subset$season) &
                                          !is.na(var_subset$block)),]
    var_subset_i_ii$var <- NA
    var_subset_i_ii$var<-(var_subset_i_ii$BMnext6mo - (predict(gr_mod[[i]][[ii]], newdata=var_subset_i_ii)+
                                                         var_subset_i_ii$BM))^2
    
    var_mod[[i]][[ii]] <-  lme4::lmer(var ~ BM +
                                        (1|timepoint) + 
                                        (1|block), 
                                      data= var_subset_i_ii,
                                      na.action="na.fail", REML=F) 
    if (model.select == "bestfit") {
      var_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(var_mod[[i]][[ii]], REML=F),1)[[1]]}
    var_mod[[i]][[ii]] <- update(var_mod[[i]][[ii]], REML=T)
    
    if(isSingular(var_mod[[i]][[ii]])) { var_mod[[i]][[ii]] <-lm(var ~ BM + block+
                                                                   timepoint, 
                                                                 data= var_subset_i_ii,
                                                                 na.action="na.fail") 
    if (model.select == "bestfit") {
      var_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(var_mod[[i]][[ii]]),1)[[1]]}
    }
    
  }}# can ignore the singular warnings; they are taken care of within the loop

################################################################################

                              #Prob of fruiting#

################################################################################
pfr_subset <- 
  hib_data %>% 
  filter(!is.na(binfruitsnext6mo) & 
           !is.na(BM) & 
           !is.na(level) & 
           !is.na(trmtUHURU) & 
           !is.na(season) & 
           !is.na(block))

#holder
pfr_mod <- list(list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA))

for (i in 1:length(levels)){
  for (ii in 1:length(trmts)){
    pfr_subset_i_ii <- pfr_subset[which(pfr_subset$level == levels[i] & pfr_subset$trmt == trmts[ii]),]
    
    pfr_mod[[i]][[ii]] <- lme4::glmer(binfruitsnext6mo~ BM + 
                                        (1|timepoint) +(1|block), 
                                      family= "binomial", 
                                      data= pfr_subset_i_ii,
                                      na.action="na.fail")
    if (model.select == "bestfit") {
      pfr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(pfr_mod[[i]][[ii]], REML=F),1)[[1]]}
    
    if(isSingular(pfr_mod[[i]][[ii]])) { pfr_mod[[i]][[ii]] <-glm(binfruitsnext6mo~ BM + block+
                                                                    timepoint , 
                                                                  family= "binomial", 
                                                                  data= pfr_subset_i_ii,
                                                                  na.action="na.fail")
    if (model.select == "bestfit") {
      pfr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(pfr_mod[[i]][[ii]]),1)[[1]]}
    }
    
  }}# can ignore the singular warnings; they are taken care of within the loop

################################################################################

                                  #Number of fruits#

################################################################################

nfr_subset <- 
  hib_data  %>% 
  filter(!is.na(fruitsnext6mo) & 
           !is.na(BM) & 
           !is.na(level) & 
           !is.na(trmtUHURU) & 
           !is.na(season) & 
           !is.na(block) &
           binfruitsnext6mo==1&
           fruitsnext6mo!=0
  )%>%
  mutate(logfruitsperBM = log(fruitsnext6mo/BM))

#holder
nfr_mod <- list(list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA), 
                list(NA, NA, NA, NA))

for (i in 1:length(levels)){
  for (ii in 1:length(trmts)){
    nfr_subset_i_ii <- nfr_subset[which(nfr_subset$level == levels[i] & nfr_subset$trmt == trmts[ii]),]
    
    #use log fruits and log BM here
    nfr_mod[[i]][[ii]] <- lme4::lmer(logfruitsperBM~ 
                                       (1|timepoint) +(1|block), 
                                     data= nfr_subset_i_ii, na.action= "na.fail", REML=F)
    if (model.select == "bestfit") {
      nfr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(nfr_mod[[i]][[ii]], REML=F),1)[[1]]}
    nfr_mod[[i]][[ii]] <- update(  nfr_mod[[i]][[ii]], REML=T)
    
    if(isSingular(nfr_mod[[i]][[ii]])) { nfr_mod[[i]][[ii]] <-lm(logfruitsperBM~ BM + block +
                                                                   timepoint , 
                                                                 data= nfr_subset_i_ii,
                                                                 na.action="na.fail") 
    if (model.select == "bestfit") {
      nfr_mod[[i]][[ii]] <- MuMIn::get.models(MuMIn::dredge(nfr_mod[[i]][[ii]]),1)[[1]]}
    }
    
  }}# can ignore the singular warnings; they are taken care of within the loop


################################################################################
# testing assumptions of these lmer/ glmers
###############################################################################

#############################SURVIVAL############################################
# Loop through each location (1 to 3)
for (i in seq_along(sur_mod)) {
  # Loop through each treatment (1 to 4)
  for (j in seq_along(sur_mod[[i]])) {
    
    model <- sur_mod[[i]][[j]]  # Extract the model
    
    # Check if it's a valid model (glmerMod or glm)
    if (inherits(model, "glmerMod")) {
      
      cat("\nChecking residuals for GLMER Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      sim_res <- recalculateResiduals(sim_res)  # Adjust for discrete data
      
      plot(sim_res)  # Display residual plot
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      
    } else if (inherits(model, "glm") && family(model)$family == "binomial") {
      
      cat("\nChecking residuals for BINOMIAL GLM Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      sim_res <- recalculateResiduals(sim_res)  # Adjust for discrete data
      plot(sim_res)  # Display residual plot
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      
    } else {
      cat("\nSkipping model at Location:", i, "Treatment:", j, " (Not glm or glmer)\n")
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      next  # Skip to the next iteration if it's not a glm or glmer
    }
  }
}


#############################GROWTH############################################

# Loop through each location (1 to 3)
for (i in seq_along(gr_mod)) {
  # Loop through each treatment (1 to 4)
  for (j in seq_along(gr_mod[[i]])) {
    
    model <- gr_mod[[i]][[j]]  # Extract the model
    
    if (inherits(model, "lmerMod")) {
      
      cat("\nChecking residuals for mixed Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      plot(sim_res)  # Base R plot, should work fine
      
      readline(prompt = "Press Enter to proceed to the next model: ")
      
    } else if (inherits(model, "lm")) {
      
      cat("\nChecking residuals for linear Model at Location:", i, "Treatment:", j, "\n")
      
      # Generate autoplot residual diagnostics
      sim_res <- autoplot(model)
      
      # Explicitly print the ggplot object
      print(sim_res)
      
      readline(prompt = "Press Enter to proceed to the next model: ")
    }
  }
}
#seems fine

#############################VARIANCE IN GROWTH###################################

 #These residuals will always be bad because it is variance

# # Loop through each location (1 to 3)
# for (i in seq_along(var_mod)) {
#   # Loop through each treatment (1 to 4)
#   for (j in seq_along(var_mod[[i]])) {
# 
#     model <- var_mod[[i]][[j]]  # Extract the model
# 
#     if (inherits(model, "lmerMod")) {
# 
#       cat("\nChecking residuals for mixed Model at Location:", i, "Treatment:", j, "\n")
# 
#       sim_res <- simulateResiduals(model)
#       plot(sim_res)  # Base R plot, should work fine
# 
#       readline(prompt = "Press Enter to proceed to the next model: ")
# 
#     } else if (inherits(model, "lm")) {
# 
#       cat("\nChecking residuals for linear Model at Location:", i, "Treatment:", j, "\n")
# 
#       # Generate autoplot residual diagnostics
#       sim_res <- autoplot(model)
# 
#       # Explicitly print the ggplot object
#       print(sim_res)
# 
#       readline(prompt = "Press Enter to proceed to the next model: ")
#     }
#   }
# }


#############################Number of Fruits###################################

# Loop through each location (1 to 3)
for (i in seq_along(nfr_mod)) {
  # Loop through each treatment (1 to 4)
  for (j in seq_along(nfr_mod[[i]])) {
    
    model <- nfr_mod[[i]][[j]]  # Extract the model
    
    if (inherits(model, "lmerMod")) {
      
      cat("\nChecking residuals for mixed Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      plot(sim_res)  # Base R plot, should work fine
      
      readline(prompt = "Press Enter to proceed to the next model: ")
      
    } else if (inherits(model, "lm")) {
      
      cat("\nChecking residuals for linear Model at Location:", i, "Treatment:", j, "\n")
      
      # Generate autoplot residual diagnostics
      sim_res <- autoplot(model)
      
      # Explicitly print the ggplot object
      print(sim_res)
      
      readline(prompt = "Press Enter to proceed to the next model: ")
    }
  }
} #overall fine


#############################Prob of fruiting####################################
# Loop through each location (1 to 3)
for (i in seq_along(pfr_mod)) {
  # Loop through each treatment (1 to 4)
  for (j in seq_along(pfr_mod[[i]])) {
    
    model <- pfr_mod[[i]][[j]]  # Extract the model
    
    # Check if it's a valid model (glmerMod or glm)
    if (inherits(model, "glmerMod")) {
      
      cat("\nChecking residuals for GLMER Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      sim_res <- recalculateResiduals(sim_res)  # Adjust for discrete data
      
      plot(sim_res)  # Display residual plot
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      
    } else if (inherits(model, "glm") && family(model)$family == "binomial") {
      
      cat("\nChecking residuals for BINOMIAL GLM Model at Location:", i, "Treatment:", j, "\n")
      
      sim_res <- simulateResiduals(model)
      sim_res <- recalculateResiduals(sim_res)  # Adjust for discrete data
      plot(sim_res)  # Display residual plot
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      
    } else {
      cat("\nSkipping model at Location:", i, "Treatment:", j, " (Not glm or glmer)\n")
      readline(prompt = "Press Enter to proceed to the next model: ")
      
      next  # Skip to the next iteration if it's not a glm or glmer
    }
  }
}#all good


################################################################################
# PREDICT VITAL RATES using mean coefficient estimates (and no perturbation)----

all_vital_rates <- array(data= NA, dim= c(3, 5, 3, 4), 
                         dimnames= list(c("sc1", "sc2", "sc3"),
                                        c("survival", "fruiting rate","p1,i", "p2, i", "p3, i"),
                                        levels, trmts))
# NB that the p1,i values are meant to reflect the probability of transitioning FROM size class 1 to size class i 

for (i in 1:length(levels)){
  mylevel <- levels[i]
  for (j in 1:length(trmts)){
    mytrmt <- trmts[j]
    
    # making the dataframe to use for prediction
    forpredict <- data.frame(matrix(data= binmids))
    names(forpredict) <- "BM"
    
    # now, actually doing the predictin'
    sur_i_j <- predict_custom(sur_mod[[i]][[j]], binmids)
    
    
    fr_i_j1 <- predict_custom(pfr_mod[[i]][[j]], binmids)
    
    fr_i_j2 <- predict_custom(nfr_mod[[i]][[j]], binmids)
    
    
    fr_i_j <- fr_i_j1*exp(fr_i_j2)*forpredict$BM# prob of fruiting times no.fruits per initial biomass, given it fruited, times intial biomass
    
    # now, predicting transitions among size classes...
    # which involve both mean and variance in growth
    #forpredict$BM <- exp(forpredict$logbm)
    # transition from sc 1 to sc 1, 2, 3
    
    
    
    #####MUTE WHEN TESTING VARIANCE
    my_var <- predict_custom(var_mod[[i]][[j]], binmids)
    
    my_mean <- predict_custom(gr_mod[[i]][[j]], binmids)
    
    forpredict_i_j <- forpredict[1, "BM"]
    #for sc1, what is the cumulative likelihood of entering next size class given a starting size?
    
    
    if ((my_var[1])< 0) {
      # cat("replacing negative variance for level", i, "and trmt", j)
      predictedvar <- 0.00000001 # replace values that are negative with small #
    } else {
      predictedvar <- my_var[1]
    }
    
    growcdf <- pnorm(binedges, 
                     mean=my_mean[1]+forpredict_i_j,
                     sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
    
    grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    trans_from_sc1_i_j <- grows} else {trans_from_sc1_i_j <- NA} 
    
    
    #transition from sc 2 to sc 1, 2, 3
    forpredict_i_j <- forpredict[2,"BM" ]
    
    #####MUTE WHEN TESTING VARIANCE
    if (my_var[2] < 0) {
      # cat("replacing negative variance for level", i, "and trmt", j)
      predictedvar <- 0.00000001 # replace values that are negative with small #
    }else {
      predictedvar <- my_var[2]
    }
    
    growcdf <- pnorm(binedges, 
                     mean=my_mean[2]+forpredict_i_j,
                     sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
    
    grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
    
    if(sum(grows)>0){grows <- grows/sum(grows)
    trans_from_sc2_i_j <- grows} else {trans_from_sc2_i_j <- NA} 
    
    # transition from sc 3 to sc 1, 2, 3
    forpredict_i_j <- forpredict[3, "BM" ]
    
    #####MUTE WHEN TESTING VARIANCE
    
    if (my_var[3]< 0) {
      #  cat("replacing negative variance for level", i, "and trmt", j)
      predictedvar <- 0.00000001 # replace values that are negative with small #
    }else {
      predictedvar <- my_var[3]
    }
    
    growcdf <- pnorm(binedges, 
                     mean=my_mean[3]+forpredict_i_j,
                     sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
    
    grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
    
    if(sum(grows)>0){grows <- grows/sum(grows)
    trans_from_sc3_i_j <- grows} else {trans_from_sc3_i_j <- NA} 
    
    # let's store them, i guess in one array
    all_vital_rates[,1, i, j] <- sur_i_j # survival
    all_vital_rates[,2,i, j] <- fr_i_j # fruit
    all_vital_rates[,3,i, j] <- trans_from_sc1_i_j # transition from size class 1 to each of the three size classes
    all_vital_rates[,4,i, j] <- trans_from_sc2_i_j # and so on
    all_vital_rates[,5,i, j] <- trans_from_sc3_i_j
    
  }  }



transtoS<-array(data= all_vital_rates[1,3:5,,], dim= c(3, 1, 3, 4), 
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoS",
                               levels, trmts))

transtoM<-array(data= all_vital_rates[2,3:5,,], dim= c(3, 1, 3, 4), 
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoM",
                               levels, trmts))
transtoL<-array(data= all_vital_rates[3,3:5,,], dim= c(3, 1, 3, 4), 
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoL",
                               levels, trmts))

gr_alt <- array(data= c(transtoS, transtoM, transtoL), dim=c(3,3,3,4),
                dimnames=list(c("sc1", "sc2", "sc3"), c("transtoS", "transtoM", "transtoL"), 
                              levels, trmts))

gr_alt<-abind(transtoS, transtoM, transtoL, along=2)

#add everything back to original array
full_vrs<-abind(all_vital_rates, gr_alt, along=2)
#remove old gr structure
full_vrs_edit <-full_vrs[,-c(3:5),,]
#rename to fit existing structures
full_vrs_rn <- array(data=full_vrs_edit, dim=c(3,5,3,4),
                     dimnames=list(c("S", "M", "L"), c("sr", "fr", "transtoS", "transtoM", "transtoL"), 
                                   levels, trmts))



#     Adding in germination rates 
germ<-read.csv(file="RawData/22-4-25_levelspecific-germ.csv")

#add row for fruiting size class to full array
fruitna<-array(data=NA, dim=c(1,5,3,4),
               dimnames=list("F", c("sr", "fr", "transtoS", "transtoM", "transtoL"), 
                             levels, trmts))
full_vrs_rn1 <-abind(fruitna,full_vrs_rn,  along=1)


#add column for germ rate
germna<-array(data=NA, dim=c(4,1,3,4),
              dimnames=list(c("F", "S", "M", "L"), "germ",
                            levels, trmts))

full_vrs_2 <- abind(full_vrs_rn1, germna, along=2)

#add germination rates--should be different across N/C/S, but same across everything else
for(i in 1:3){
  mylevel <- levels[i]
  for (j in 1:4){
    mytrmt <- trmts[j]
    full_vrs_2["F", "germ", i, j] <- germ$germ[which(germ$level==mylevel)]
  }
}
full_vrs_allsites <- full_vrs_2
full_vrs_averaged <- apply(full_vrs_allsites, c(1,2,4), mean) # dimensions here are size class 1-4 ("F" "S" "M" "L"),
# vital rate 1-6 ("sr"       "fr"       "transtoS" "transtoM" "transtoL" "germ"    ), 
# and treatment 1-4 ("CONT" "LMH"  "MEGA" "MESO")



#############################################################################################################
# predict vital rate values when coefficients of vital rate functions are sampled from the var-cov matrix----
##########################################################################################################

rand_iter<-10000

all_vital_rates <- array(data= NA, dim= c(3, 5, 3, 4, rand_iter),
                         dimnames= list(c("sc1", "sc2", "sc3"),
                                        c("survival", "fruiting rate","p1,i", "p2, i", "p3, i"),
                                        levels, trmts, 1:rand_iter))

for (z in 1:rand_iter){
  myiterationofvrs<-z
  
  for (i in 1:length(levels)){
    mylevel <- levels[i]
    for (j in 1:length(trmts)){
      mytrmt <- trmts[j]
      
      #make copies of mod for replacing coefficients
      sur_mod1<-sur_mod[[i]][[j]]
      gr_mod1<-gr_mod[[i]][[j]]
      var_mod1<-var_mod[[i]][[j]]
      pfr_mod1<-pfr_mod[[i]][[j]]
      nfr_mod1<-nfr_mod[[i]][[j]]
      
      
      sur_mass <- MASS::mvrnorm(mu = summary(sur_mod[[i]][[j]])$coefficients[, "Estimate"],
                                Sigma = vcov(sur_mod[[i]][[j]])) # this command works for both glmers and glm
      
      gr_mass <- MASS::mvrnorm(mu = summary(gr_mod[[i]][[j]])$coefficients[, "Estimate"],
                               Sigma = vcov(gr_mod[[i]][[j]]))
      
      var_mass <- MASS::mvrnorm(mu = summary(var_mod[[i]][[j]])$coefficients[, "Estimate"],
                                Sigma = vcov(var_mod[[i]][[j]]))
      
      pfr_mass <- MASS::mvrnorm(mu = summary(pfr_mod[[i]][[j]])$coefficients[, "Estimate"],
                                Sigma = vcov(pfr_mod[[i]][[j]]))
      
      nfr_mass <- MASS::mvrnorm(mu = summary(nfr_mod[[i]][[j]])$coefficients[, "Estimate"],
                                Sigma = vcov(nfr_mod[[i]][[j]]))
      
      if (class(sur_mod1)[1]== "glm") { sur_mod1$coefficients <- sur_mass } else sur_mod1@beta <- sur_mass
      if (class(gr_mod1)[1]== "lm") { gr_mod1$coefficients <- gr_mass } else gr_mod1@beta <- gr_mass
      if (class(var_mod1)[1]== "lm") { var_mod1$coefficients <- var_mass } else var_mod1@beta <- var_mass
      if (class(pfr_mod1)[1]== "glm") { pfr_mod1$coefficients <- pfr_mass } else pfr_mod1@beta <- pfr_mass
      if (class(nfr_mod1)[1]== "lm") { nfr_mod1$coefficients <- nfr_mass } else nfr_mod1@beta <- nfr_mass
      
      
      
      # making the dataframe to use for prediction
      forpredict <- data.frame(matrix(data= (binmids)))
      names(forpredict) <- "BM"
      forpredict$level <- factor(mylevel)
      forpredict$trmtUHURU <- factor(mytrmt)
      
      # now, actually doing the predictin'
      sur_i_j <- predict_custom(sur_mod1, binmids)
      
      
      fr_i_j1 <- predict_custom(pfr_mod1, binmids)
      fr_i_j2 <- predict_custom(nfr_mod1, binmids)
      
      
      fr_i_j <- fr_i_j1*exp(fr_i_j2)*forpredict$BM# prob of fruiting times no.fruits per initial biomass, given it fruited, times intial biomass
      
      # now, predicting transitions among size classes...
      # which involve both mean and variance in growth
      #forpredict$BM <- exp(forpredict$logbm)
      # transition from sc 1 to sc 1, 2, 3
      
      
      
      #####MUTE WHEN TESTING VARIANCE
      my_var <- predict_custom(var_mod1, binmids)
      
      my_mean <- predict_custom(gr_mod1, binmids)
      
      forpredict_i_j <- forpredict[1, ]
      #for sc1, what is the cumulative likelihood of entering next size class given a starting size?
      
      
      if ((my_var[1])< 0) {
        #  cat("replacing negative variance for level", i, "and trmt", j)
        predictedvar <- 0.00000001 # replace values that are negative with small #
      } else {
        predictedvar <- my_var[1]
      }
      
      growcdf <- pnorm(binedges, 
                       mean=my_mean[1]+forpredict_i_j$BM,
                       sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
      
      grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
      if(sum(grows)>0){grows <- grows/sum(grows)
      trans_from_sc1_i_j <- grows} else {trans_from_sc1_i_j <- NA} 
      
      
      #transition from sc 2 to sc 1, 2, 3
      forpredict_i_j <- forpredict[2, ]
      
      #####MUTE WHEN TESTING VARIANCE
      if (my_var[2] < 0) {
        # cat("replacing negative variance for level", i, "and trmt", j)
        predictedvar <- 0.00000001 # replace values that are negative with small #
      }else {
        predictedvar <- my_var[2]
      }
      
      growcdf <- pnorm(binedges, 
                       mean=my_mean[2]+forpredict_i_j$BM,
                       sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
      
      grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
      
      if(sum(grows)>0){grows <- grows/sum(grows)
      trans_from_sc2_i_j <- grows} else {trans_from_sc2_i_j <- NA} 
      
      # transition from sc 3 to sc 1, 2, 3
      forpredict_i_j <- forpredict[3, ]
      
      #####MUTE WHEN TESTING VARIANCE
      
      if (my_var[3]< 0) {
        # cat("replacing negative variance for level", i, "and trmt", j)
        predictedvar <- 0.00000001 # replace values that are negative with small #
      }else {
        predictedvar <- my_var[3]
      }
      
      growcdf <- pnorm(binedges, 
                       mean=my_mean[3]+forpredict_i_j$BM,
                       sd=sqrt(predictedvar)) # NB this treatment of var/ sd is a bit different than the other 3 code
      
      grows <- growcdf[2:length(binedges)]-growcdf[1:(length(binedges)-1)]
      
      if(sum(grows)>0){grows <- grows/sum(grows)
      trans_from_sc3_i_j <- grows} else {trans_from_sc3_i_j <- NA} 
      
      # let's store them, i guess in one array
      all_vital_rates[,1, i, j,z] <- sur_i_j # survival
      all_vital_rates[,2,i, j,z] <- fr_i_j # fruit
      all_vital_rates[,3,i, j,z] <- trans_from_sc1_i_j # transition from size class 1 to each of the three size classes
      all_vital_rates[,4,i, j,z] <- trans_from_sc2_i_j # and so on
      all_vital_rates[,5,i, j,z] <- trans_from_sc3_i_j
      
    }  }
}

#reformatting randomly sampled growth rates to be more matrix friendly

transtoS<-array(data= all_vital_rates[1,3:5,,,], dim= c(3, 1, 3, 4, rand_iter),
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoS",
                               levels, trmts))

transtoM<-array(data= all_vital_rates[2,3:5,,,], dim= c(3, 1, 3, 4,rand_iter),
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoM",
                               levels, trmts))
transtoL<-array(data= all_vital_rates[3,3:5,,,], dim= c(3, 1, 3, 4, rand_iter),
                dimnames= list(c("sc1", "sc2", "sc3"), "transtoL",
                               levels, trmts))

gr_alt <- array(data= NA, dim=c(3,3,3,4, rand_iter),
                dimnames=list(c("sc1", "sc2", "sc3"), c("transtoS", "transtoM", "transtoL"),
                              levels, trmts, 1:rand_iter))


for (i in 1:rand_iter){
  for (j in 1:4){
    for (k in 1:3){
      transtoS_iteration_herb_lvl<-transtoS[,,k,j,i]
      transtoM_iteration_herb_lvl<-transtoM[,,k,j,i]
      transtoL_iteration_herb_lvl<-transtoL[,,k,j,i]
      
      yoink<- abind(transtoS_iteration_herb_lvl, transtoM_iteration_herb_lvl, transtoL_iteration_herb_lvl)
      
      gr_alt[,,k,j,i]<-yoink
      
    }
  }
}

#add everything back to original array
full_vrs<-abind(all_vital_rates, gr_alt, along=2)
#remove old gr structure
full_vrs_edit <-full_vrs[,-c(3:5),,,]
#rename to fit existing structures
full_vrs_rn <- array(data=full_vrs_edit, dim=c(3,5,3,4, rand_iter),
                     dimnames=list(c("S", "M", "L"), c("sr", "fr", "transtoS", "transtoM", "transtoL"),
                                   levels, trmts, 1:rand_iter))

# Adding in germination rates
germ<-read.csv(file="RawData/22-4-25_levelspecific-germ.csv")

#add row for fruiting size class to full array
fruitna<-array(data=NA, dim=c(1,5,3,4, rand_iter),
               dimnames=list("F", c("sr", "fr", "transtoS", "transtoM", "transtoL"),
                             levels, trmts, 1:rand_iter))
full_vrs_rn1 <-abind(fruitna,full_vrs_rn,  along=1)

#add column for germ rate
germna<-array(data=NA, dim=c(4,1,3,4, rand_iter),
              dimnames=list(c("F", "S", "M", "L"), "germ",
                            levels, trmts, 1:rand_iter))

full_vrs_2 <- abind(full_vrs_rn1, germna, along=2)

#add germination rates--should be different across N/C/S, but same across everything else
for(z in 1:rand_iter){
  for(i in 1:3){
    mylevel <- levels[i]
    for (j in 1:4){
      mytrmt <- trmts[j]
      full_vrs_2["F", "germ", i, j,z] <- germ$germ[which(germ$level==mylevel)]
    }
  }
}

#check to make sure iterating worked
hist(full_vrs_2[3,1,1,1,])
my_quantile <- function(input){if (is.null(input)) {output <- NA} else {quantile (input, c(0.5/2, 0.5, (1-0.05)/2), na.rm=TRUE)}} #custom function that returns NA if the vector is empty
vrs_random_CIs <- apply(full_vrs_2, c(1, 2, 3, 4), my_quantile) # this command gives you the median and CI's on each vital rate
save(vrs_random_CIs, file= "Data/3_mean_CIs_on_vital_rates.RData")
full_vrs_random <- full_vrs_2

# perturb rates for subsequent senstivity analyses (NB: you want to perturb the mean--across levels-- rate)----

perturbs<-c('high_sr1', 'low_sr1',
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

sizes <-c("FR", "S", "M", "L")

rates <- c("sr", "fr", "transtoS", "transtoM", "transtoL", "germ")

perturbarray <- array(data=NA, dim=c(4,6,4,32), 
                      dimnames=list(sizes, rates,trmts, perturbs))

perturbarray2<- abind(
  full_vrs_averaged, # takes the across-level mean vital rate
  perturbarray, along=4)

perturbs2<-c('ambient',
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


perturbarray3 <- array(data=perturbarray2, dim=c(4,6,4,33), 
                       dimnames=list(sizes, rates, trmts, perturbs2))


#now that we have an array with all the dimensions we need, 
#perturb each vital rate in said array 
#do this only for the non-germ stuff because of array structure
perturb_no <-c(2,8,14,20,26)


for(y in 1:4){
  ###survival###
  #small
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[2,1,y] ==0){
    yote[2,1,y]<-0.01
    yeet[2,1,y]<-0.01
  } else if (full_vrs_averaged[2,1,y] ==1){
    yote[2,1,y]<-0.99
    yeet[2,1,y]<-0.99
  } else {
    yote[2,1,y]<-1.05*(full_vrs_averaged[2,1,y])
    yeet[2,1,y]<-0.95*(full_vrs_averaged[2,1,y]) }
  
  perturbarray3[,,y,2]<-yote[,,y] #high perturb 
  perturbarray3[,,y,3]<-yeet[,,y] #low perturb 
  ###med###
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[3,1,y] ==0){
    yote[3,1,y]<-0.01
    yeet[3,1,y]<-0.01
  } else if (full_vrs_averaged[3,1,y] ==1){
    yote[3,1,y]<-0.99
    yeet[3,1,y]<-0.99
  } else {
    yote[3,1,y] <-1.05*(full_vrs_averaged[3,1,y]) 
    yeet[3,1,y] <-0.95*(full_vrs_averaged[3,1,y]) }
  
  perturbarray3[,,y,4]<-yote[,,y] #high perturb 
  perturbarray3[,,y,5]<-yeet[,,y]#low perturb 
  
  ###large###
  
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[4,1,y] ==0){
    yote[4,1,y]<-0.01
    yeet[4,1,y]<-0.01
  } else if (full_vrs_averaged[4,1,y] ==1){
    yote[4,1,y]<-0.99
    yeet[4,1,y]<-0.99
  } else {
    yote[4,1,y]<- 1.05*(full_vrs_averaged[4,1, y])
    yeet[4,1,y]<-0.95*(full_vrs_averaged[4,1,y])}
  
  perturbarray3[,,y,6]<- yote[,,y]#high perturb 
  perturbarray3[,,y,7]<- yeet[,,y]#low perturb 
  
  
  ###fr###
  #small
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[2,2,y] ==0){
    yote[2,2,y]<-0.01
    yeet[2,2,y]<-0.01
  } else {
    yote[2,2,y]<-1.05*(full_vrs_averaged[2,2,y])
    yeet[2,2,y]<-0.95*(full_vrs_averaged[2,2,y]) }
  
  perturbarray3[,,y,8]<-yote[,,y]#high perturb 
  perturbarray3[,,y,9]<-yeet[,,y]#low perturb 
  
  #med
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  
  if(full_vrs_averaged[3,2,y] ==0){
    yote[3,2,y]<-0.01
    yeet[3,2,y]<-0.01
  } else {
    yote[3,2,y]<-1.05*(full_vrs_averaged[3,2,y])
    yeet[3,2,y]<-0.95*(full_vrs_averaged[3,2,y]) }
  
  perturbarray3[,,y,10]<-yote[,,y] #high perturb 
  perturbarray3[,,y,11]<-yeet[,,y]#low perturb 
  
  #large
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[4,2,y] ==0){
    yote[4,2,y]<-0.01
    yeet[4,2,y]<-0.01
  } else {
    yote[4,2,y]<-1.05*(full_vrs_averaged[4,2,y])
    yeet[4,2,y]<-0.95*(full_vrs_averaged[4,2,y]) }
  
  perturbarray3[,,y,12]<-yote[,,y]#high perturb 
  perturbarray3[,,y,13]<-yeet[,,y]#low perturb 
  
  ###transtoS###
  
  #small
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[2,3,y] ==0){
    yote[2,3,y]<-0.01
    yeet[2,3,y]<-0.01
  } else {
    yote[2,3,y]<-1.05*(full_vrs_averaged[2,3,y])
    yeet[2,3,y]<-0.95*(full_vrs_averaged[2,3,y]) }
  
  perturbarray3[,,y,14]<-yote[,,y] #high perturb 
  perturbarray3[,,y,15]<-yeet[,,y]#low perturb 
  
  #med
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[3,3,y] ==0){
    yote[3,3,y]<-0.01
    yeet[3,3,y]<-0.01
  } else {
    yote[3,3,y]<-1.05*(full_vrs_averaged[3,3,y])
    yeet[3,3,y]<-0.95*(full_vrs_averaged[3,3,y]) }
  
  perturbarray3[,,y,16]<-yote[,,y] #high perturb
  perturbarray3[,,y,17]<-yeet[,,y] #low perturb 
  
  #large
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[4,3,y] ==0){
    yote[4,3,y]<-0.01
    yeet[4,3,y]<-0.01
  } else {
    yote[4,3,y]<-1.05*(full_vrs_averaged[4,3,y])
    yeet[4,3,y]<-0.95*(full_vrs_averaged[4,3,y]) }
  
  perturbarray3[,,y,18]<-yote[,,y]#high perturb 
  perturbarray3[,,y,19]<-yeet[,,y] #low perturb 
  
  ###transtoM###
  
  #small
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged 
  
  
  if(full_vrs_averaged[2,4,y] ==0){
    yote[2,4,y]<-0.01
    yeet[2,4,y]<-0.01
  } else {
    yote[2,4,y]<-1.05*(full_vrs_averaged[2,4,y])
    yeet[2,4,y]<-0.95*(full_vrs_averaged[2,4,y]) }
  
  
  perturbarray3[,,y,20]<-yote[,,y] #high perturb 
  perturbarray3[,,y,21]<-yeet[,,y]#low perturb 
  
  #med
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[3,4,y] ==0){
    yote[3,4,y]<-0.01
    yeet[3,4,y]<-0.01
  } else {
    yote[3,4,y]<-1.05*(full_vrs_averaged[3,4,y])
    yeet[3,4,y]<-0.95*(full_vrs_averaged[3,4,y]) }
  
  perturbarray3[,,y,22]<-yote[,,y] #high perturb 
  perturbarray3[,,y,23]<-yeet[,,y] #low perturb 
  
  #large
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[4,4,y] ==0){
    yote[4,4,y]<-0.01
    yeet[4,4,y]<-0.01
  } else {
    yote[4,4,y]<-1.05*(full_vrs_averaged[4,4,y])
    yeet[4,4,y]<-0.95*(full_vrs_averaged[4,4,y]) }
  
  perturbarray3[,,y,24]<-yote[,,y]#high perturb 
  perturbarray3[,,y,25]<-yeet[,,y] #low perturb 
  
  ###transtoL###
  
  #small
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[2,5,y] ==0){
    yote[2,5,y]<-0.01
    yeet[2,5,y]<-0.01
  } else {
    yote[2,5,y]<-1.05*(full_vrs_averaged[2,5,y])
    yeet[2,5,y]<-0.95*(full_vrs_averaged[2,5,y]) }
  
  perturbarray3[,,y,26]<-yote[,,y] #high perturb 
  perturbarray3[,,y,27]<-yeet[,,y]#low perturb 
  
  #med
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[3,5,y] ==0){
    yote[3,5,y]<-0.01
    yeet[3,5,y]<-0.01
  } else {
    yote[3,5,y]<-1.05*(full_vrs_averaged[3,5,y])
    yeet[3,5,y]<-0.95*(full_vrs_averaged[3,5,y]) }
  
  perturbarray3[,,y,28]<-yote[,,y] #high perturb 
  perturbarray3[,,y,29]<-yeet[,,y] #low perturb 
  
  #large
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[4,5,y] ==0){
    yote[4,5,y]<-0.01
    yeet[4,5,y]<-0.01
  } else {
    yote[4,5,y]<-1.05*(full_vrs_averaged[4,5,y])
    yeet[4,5,y]<-0.95*(full_vrs_averaged[4,5,y]) }
  
  perturbarray3[,,y,30]<-yote[,,y] #high perturb 
  perturbarray3[,,y,31]<-yeet[,,y] #low perturb 
  
  
  ###germination###
  yote<-full_vrs_averaged
  yeet<-full_vrs_averaged
  
  if(full_vrs_averaged[1,6,y] ==0){
    yote[1,6,y]<-0.01
    yeet[1,6,y]<-0.01
  } else {
    yote[1,6,y]<-1.05*(full_vrs_averaged[1,6,y])
    yeet[1,6,y]<-0.95*(full_vrs_averaged[1,6,y]) }
  
  perturbarray3[,,y,32]<-yote[,,y]
  perturbarray3[,,y,33]<-yeet[,,y]
  
}
full_vrs_perturb <- perturbarray3

# store the outputs---
save(sur_mod, gr_mod, var_mod, pfr_mod, nfr_mod, 
     full_vrs_allsites, 
     full_vrs_averaged,
     full_vrs_perturb, 
     full_vrs_random,
     file="Data/3_vital_rate_functions_and_vital_rate_values.RData") 


###############################################################################
                            #Table of VRs#
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


df_wide <- df %>%
  pivot_wider(names_from = c(Treatment, Site), values_from = Value, names_sep = "_")



# Create a header dataframe
header_df <- data.frame(
  col_keys = colnames(df_wide),  # Ensure col_keys match df_wide
  Treatment = c("", "","", "Control", "", "", "LMH", "","", "MEGA", "","", "MESO", ""),  # Top-level header
  Site = c("ID", "N", "C", "S", "N", "C", "S")  # Second-level header
)

# Create the flextable with grouped headers
ft <- flextable(df_wide) %>%
  set_header_df(mapping = header_df, key = "col_keys") %>%
  merge_h(part = "header") %>%
  theme_vanilla() %>%  # Optional theme
  autofit() %>%        # Adjust column width
  bold(j = 1:2, bold = TRUE) %>%  # Bold the first two columns
  set_header_labels("Size Class" = "Size Class", "Rate" = "Rate",
                    "CONT_N" = "Arid",
                    "CONT_C" = "Intermediate",
                    "CONT_S" = "Mesic",
                    "LMH_N" = "Arid",
                    "LMH_C" = "Intermediate",
                    "LMH_S" = "Mesic",
                    "MEGA_N" = "Arid",
                    "MEGA_C" = "Intermediate",
                    "MEGA_S" = "Mesic",
                    "MESO_N" = "Arid",
                    "MESO_C" = "Intermediate",
                    "MESO_S" = "Mesic")  # Adjust labels if needed

# Export to Word

ft  # Display flextable


# Export to Word document
doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_add_par("\n")  # Adds spacing after the table

print(doc, target = "vr_table_output.docx")

