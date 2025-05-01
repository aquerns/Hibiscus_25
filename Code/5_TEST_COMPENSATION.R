############################################################################
########                                                         ###########
########      9-14-23: Test for demographic compensation         ###########
########                                                         ###########
############################################################################
#Written by: Aleah Querns
# last modified: A Louthan 241220

#####################          TASKS                    #####################

#1: Get significant negative/positive correlations across N/C/S within herbivory/no herbivory/mega/meso
#2: Calculate percentage of negative/positive correlations out of total possible correlations (REALIZED PERCENT)
#3: Randomization test! Over 10,000 iterations, randomly reassign LTRE values across N/C/s
#4: For each iteration, calculate percentage of negative correlations + percentage of positive correlations
#5: Create a null distribution from iterations
#6: Compare realized percentage of neg/pos to null distributions

###Load libraries------
rm(list = ls())

  library(popbio)
  library(tidyverse)
  library(lsr)
  library(reshape2)
 # library(truncnorm)
  library(gridExtra)
  library(data.table)
library(plotly)

# loading in data 
load("Data/5_LTRE_contributions.RData")

#Number of randomizations for use in building null distributions
rando=10000# No. of reps in randomization

rates<- c("sr", "fr", "transtoS", "transtoM", "transtoL", "germ")
sizes <- c("FR", "S", "M", "L")
levels <- c("C", "N", "S")
trmts <- c("CONT" ,          "LMH",          "MEGA",          "MESO")


# Get observed number of significant neg/pos correlations & calculate percentage of significant neg/pos correlations ####

# first, positive correlations
dfstore_pos <- dfstore_neg <- array(data=NA, dim=c(4, 4,6,4,6,2), dimnames= list(trmts, sizes,rates,
                                                              sizes,rates,
                                                                c("P", "corr")))
countpos <-  countneg <- array(data=NA, dim=c(4, 4,6,4,6), dimnames= list(trmts, sizes,rates,
                                                         sizes,rates))

for (ii in 1:length(trmts)){ # runs through trmts
  LTRE_contributions_ii <- LTRE_contributions[,ii,] # getting LTRE contributions of treatment ii
  # then, reorganizing this array into a dataframe that works with AQ's code
  LTRE_contributions_ii_mod <- array(NA, dim=c(4, 6, 3), dimnames = list(sizes, rates, levels) ) 
  # NB: sizes <- c("FR", "S", "M", "L"); 
  # rates<- c("sr", "fr", "transtoS", "transtoM", "transtoL", "germ"); 
  # levels <- c("C", "N", "S")
  
    LTRE_contributions_ii_mod["FR", "germ",] <- LTRE_contributions_ii[,"germ"]
    LTRE_contributions_ii_mod["S", "sr",] <- LTRE_contributions_ii[,"sr1"]
    LTRE_contributions_ii_mod["S", "fr",] <- LTRE_contributions_ii[,"fr1"]
    LTRE_contributions_ii_mod["S", "transtoS",] <- LTRE_contributions_ii[,'S-S']
    LTRE_contributions_ii_mod["S", "transtoM",] <-LTRE_contributions_ii[,'S-M']
    LTRE_contributions_ii_mod["S", "transtoL",] <- LTRE_contributions_ii[,'S-L']
    LTRE_contributions_ii_mod["M", "sr",] <- LTRE_contributions_ii[,"sr2"]
    LTRE_contributions_ii_mod["M", "fr",] <- LTRE_contributions_ii[,"fr2"]
    LTRE_contributions_ii_mod["M", "transtoS",] <-LTRE_contributions_ii[,'M-S']
    LTRE_contributions_ii_mod["M", "transtoM",] <-LTRE_contributions_ii[,'M-M']
    LTRE_contributions_ii_mod["M", "transtoL",] <-LTRE_contributions_ii[,'M-L']
    LTRE_contributions_ii_mod["L", "sr",] <- LTRE_contributions_ii[,"sr3"]
    LTRE_contributions_ii_mod["L", "fr",] <- LTRE_contributions_ii[,"fr3"]
    LTRE_contributions_ii_mod["L", "transtoS",] <- LTRE_contributions_ii[,'L-S']
    LTRE_contributions_ii_mod["L", "transtoM",] <- LTRE_contributions_ii[,'L-M']
    LTRE_contributions_ii_mod["L", "transtoL",] <-LTRE_contributions_ii[,'L-L']
    # NB that the ORDER of levels is C-N-S here, rather than N-C-S. 
  # ALSO, same as LTRE_H, LTRE_contributions_i_mod now contains NA's when the vital rate DN apply
    
    #holder dataframe to store p-value and correlation value for a single pairwise comparison of vital rates
holder_NA<- data.frame(pval=NA, corval=NA) 


#note
#each array holding LTRE values is 4x6x3. 4 size classes, 6 vital rates, 3 sites (N/C/S)
#what I've opted to do is go spot by spot grabbing up all vital rates across North/Central/South in the first two dimensions [sizes, rates,], 
#and then I get correlations with every other vital rate across North/Central/South in the first two dimensions
#essentially it is pairwise comparisons of correlations between each vital rate across the range

for (i in sizes){
  for (k in rates){
    for (j in sizes){
      for (z in rates){
        holder_ikjz_pos <-  holder_ikjz_neg <- holder_NA
        p1<-LTRE_contributions_ii_mod[i, k,] #Vital rate 1 across N/C/S
        p2<-LTRE_contributions_ii_mod[j, z,] #Vital rate 2 across N/C/S
        
        if(i==j & k==z){ #we don't want to compare same vital rates: put NA for pval and corval
          
          holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
          holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
          
        } else if(any(is.na(p1)) | any(is.na(p2))) {#we don't want to run the cortest if there are NAs in vital rate 1
          
          holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
          holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
          
      }  else if (sd(p1)==0 | sd(p2)==0) {#If the standard deviation of LTRE values of either vr1 and vr2 is zero, add tiny random variation
          meow<-  cat("excluding because eiher size", i, " vital rate", k, "and size", j, "vital rate", z, "have no variation")
          print(meow)
          
          holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
          holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
          
        } 
 else {#if no NAs, no sd=0, and not comparing same vital rates, run cor-test
          
          
          cordf.pos<-cor.test(p1,p2, method="pearson", alternative= "greater")# run cor.test: pearson, one sided
          cordf.neg <-cor.test(p1,p2, method="pearson", alternative= "less")
          holder_ikjz_pos[1,1]<-cordf.pos[["p.value"]]
          holder_ikjz_pos[1,2]<- cordf.pos[["estimate"]]
          holder_ikjz_neg[1,1]<-cordf.neg[["p.value"]]
          holder_ikjz_neg[1,2]<- cordf.neg[["estimate"]]
          
        } 
        
        dfstore_pos[ii,i,k,j,z,]<- as.matrix(holder_ikjz_pos) #put the pval and corval between vr1 and vr2 in the proper place of comparison
        dfstore_neg[ii,i,k,j,z,]<- as.matrix(holder_ikjz_neg) #put the pval and corval between vr1 and vr2 in the proper place of comparison
        
        #COUNT POSITIVE CORRELATIONS
        if(is.na(dfstore_pos[ii, i,k,j,z,"P"])) {
          
          countpos[ii, i,k,j,z] <-NA #if P is NA, do not count
          
        } else if(
          
          !is.na(dfstore_pos[ii, i,k,j,z,"P"]) & #P not NA
          
          dfstore_pos[ii, i,k,j,z,"P"] <= 0.05 & #P is significant
          
          dfstore_pos[ii, i,k,j,z,"corr"]>0){ #corr value is positive
          
          countpos[ii, i,k,j,z]<-1 #count a pos correlation
          
        } else {
          
          countpos[ii, i,k,j,z] <-0
          
        } #if P is not NA but doesn't follow the other requirements, count as zero
        
        #COUNT NEGATIVE CORRELATIONS
        if(is.na(dfstore_neg[ii, i,k,j,z,"P"])) {
          
          countneg[ii,i,k,j,z] <- NA#if P is NA, do not count
          
        } else if(
          
          !is.na(dfstore_neg[ii,i,k,j,z,"P"]) &#P is NOT NA
          
          dfstore_neg[ii,i,k,j,z,"P"] <= 0.05 &#P is significant
          
          dfstore_neg[ii,i,k,j,z,"corr"]<0){ #correlation is negative
          
          countneg[ii,i,k,j,z]<-1 #count a negative correlation
          
        } else {
          countneg[ii,i,k,j,z]<- 0
        }
      }}}}

}

hist(dfstore_pos[1,,,,,"P"]) #what do the P's look like? (trmt 1, pos corr)

sum(!is.na(dfstore_pos[1, ,,,,]))/2# total # of spots which aren't NA (# of comparisons) (trmt 1, pos corr)

sum(dfstore_pos[1,,,,,"P"]<= 0.05, na.rm=TRUE) #total # of Ps < 0.05 --is there variation across runs? good to check because if so, code is broken


# observed proportion of pos & neg corr's
proppos <- propneg<- rep(NA, length.out= length(trmts))
for (ii in 1:length(trmts)){
proppos[ii] <- sum(countpos[ii,,,,], na.rm=T)/sum(!is.na(countpos[ii,,,,]))#Proportion positive corr out of total corr's that aren't NA
propneg[ii] <- sum(countneg[ii,,,,], na.rm=T)/sum(!is.na(countneg[ii,,,,]))} #Proportion neg corr out of total corr's that aren't NA

############################################################################
# generating null distributions using a randomization procedure ------
############################################################################

dfstore_pos_random <- dfstore_neg_random  <- array(data=NA, dim=c(4, 4,6,4,6,2, rando), dimnames= list(trmts, sizes,rates,
                                                                                                       sizes,rates,
                                                                                                       c("P", "corr"), 1:rando))
countpos_random <-  countneg_random <- array(data=NA, dim=c(4, 4,6,4,6, rando), dimnames= list(trmts, sizes,rates,
                                                                                            sizes,rates, 1:rando))

for (rep in 1:rando){
for (ii in 1:length(trmts)){ # runs through trmts
  LTRE_contributions_ii <- LTRE_contributions[,ii,] # getting LTRE contributions of treatment ii
  # then, reorganizing this array into a dataframe that works with AQ's code
  LTRE_contributions_ii_mod <- array(NA, dim=c(4, 6, 3), dimnames = list(sizes, rates, levels) ) 
  # NB: sizes <- c("FR", "S", "M", "L"); 
  # rates<- c("sr", "fr", "transtoS", "transtoM", "transtoL", "germ"); 
  # levels <- c("C", "N", "S")
  
  LTRE_contributions_ii_mod["FR", "germ",] <- LTRE_contributions_ii[,"germ"]
  LTRE_contributions_ii_mod["S", "sr",] <- LTRE_contributions_ii[,"sr1"]
  LTRE_contributions_ii_mod["S", "fr",] <- LTRE_contributions_ii[,"fr1"]
  LTRE_contributions_ii_mod["S", "transtoS",] <- LTRE_contributions_ii[,'S-S']
  LTRE_contributions_ii_mod["S", "transtoM",] <-LTRE_contributions_ii[,'S-M']
  LTRE_contributions_ii_mod["S", "transtoL",] <- LTRE_contributions_ii[,'S-L']
  LTRE_contributions_ii_mod["M", "sr",] <- LTRE_contributions_ii[,"sr2"]
  LTRE_contributions_ii_mod["M", "fr",] <- LTRE_contributions_ii[,"fr2"]
  LTRE_contributions_ii_mod["M", "transtoS",] <-LTRE_contributions_ii[,'M-S']
  LTRE_contributions_ii_mod["M", "transtoM",] <-LTRE_contributions_ii[,'M-M']
  LTRE_contributions_ii_mod["M", "transtoL",] <-LTRE_contributions_ii[,'M-L']
  LTRE_contributions_ii_mod["L", "sr",] <- LTRE_contributions_ii[,"sr3"]
  LTRE_contributions_ii_mod["L", "fr",] <- LTRE_contributions_ii[,"fr3"]
  LTRE_contributions_ii_mod["L", "transtoS",] <- LTRE_contributions_ii[,'L-S']
  LTRE_contributions_ii_mod["L", "transtoM",] <- LTRE_contributions_ii[,'L-M']
  LTRE_contributions_ii_mod["L", "transtoL",] <-LTRE_contributions_ii[,'L-L']
  # NB that the ORDER of levels is C-N-S here, rather than N-C-S. 
  # ALSO, same as LTRE_H, LTRE_contributions_i_mod now contains NA's when the vital rate DN apply
  
  #holder dataframe to store p-value and correlation value for a single pairwise comparison of vital rates
  holder_NA<- data.frame(pval=NA, corval=NA) 
  
  
  #note
  #each array holding LTRE values is 4x6x3. 4 size classes, 6 vital rates, 3 sites (N/C/S)
  #what I've opted to do is go spot by spot grabbing up all vital rates across North/Central/South in the first two dimensions [sizes, rates,], 
  #and then I get correlations with every other vital rate across North/Central/South in the first two dimensions
  #essentially it is pairwise comparisons of correlations between each vital rate across the range
  
  for (i in sizes){
    for (k in rates){
      for (j in sizes){
        for (z in rates){
          holder_ikjz_pos <-  holder_ikjz_neg <- holder_NA
          pa<-LTRE_contributions_ii_mod[i, k,] #Vital rate 1 across N/C/S
          pb<-LTRE_contributions_ii_mod[j, z,] #Vital rate 2 across N/C/S
          
          p1<- sample(pa)
          p2 <- sample(pb)
          if(i==j & k==z){ #we don't want to compare same vital rates: put NA for pval and corval
            
            holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
            holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
            
          } else if(any(is.na(p1)) | any(is.na(p2))) {#we don't want to run the cortest if there are NAs in vital rate 1
            
            holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
            holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
            
          }  else if (sd(p1)==0 | sd(p2)==0) {#If the standard deviation of LTRE values of either vr1 and vr2 is zero, add tiny random variation
            #meow<-  cat("excluding because eiher size", i, " vital rate", k, "and size", j, "vital rate", z, "have no variation")
            #print(meow)
            
            holder_ikjz_pos<-holder_ikjz_pos %>% transmute(pval=NA, corval=NA) 
            holder_ikjz_neg<-holder_ikjz_neg %>% transmute(pval=NA, corval=NA) 
            
          } 
          else {#if no NAs, no sd=0, and not comparing same vital rates, run cor-test
            
            
            cordf.pos<-cor.test(p1,p2, method="pearson", alternative= "greater")# run cor.test: pearson, one sided
            cordf.neg <-cor.test(p1,p2, method="pearson", alternative= "less")
            holder_ikjz_pos[1,1]<-cordf.pos[["p.value"]]
            holder_ikjz_pos[1,2]<- cordf.pos[["estimate"]]
            holder_ikjz_neg[1,1]<-cordf.neg[["p.value"]]
            holder_ikjz_neg[1,2]<- cordf.neg[["estimate"]]
            
          } 
          
          dfstore_pos_random[ii,i,k,j,z,, rep]<- as.matrix(holder_ikjz_pos) #put the pval and corval between vr1 and vr2 in the proper place of comparison
          dfstore_neg_random[ii,i,k,j,z,, rep]<- as.matrix(holder_ikjz_neg) #put the pval and corval between vr1 and vr2 in the proper place of comparison
          
          #COUNT POSITIVE CORRELATIONS
          if(is.na(dfstore_pos_random[ii, i,k,j,z,"P", rep])) {
            
            countpos_random[ii, i,k,j,z, rep] <-NA #if P is NA, do not count
            
          } else if(
            
            !is.na(dfstore_pos_random[ii, i,k,j,z,"P", rep]) & #P not NA
            
            dfstore_pos_random[ii, i,k,j,z,"P", rep] <= 0.05 & #P is significant
            
            dfstore_pos_random[ii, i,k,j,z,"corr", rep]>0){ #corr value is positive
            
            countpos_random[ii, i,k,j,z, rep]<-1 #count a pos correlation
            
          } else {
            
            countpos_random[ii, i,k,j,z, rep] <-0
            
          } #if P is not NA but doesn't follow the other requirements, count as zero
          
          #COUNT NEGATIVE CORRELATIONS
          if(is.na(dfstore_neg_random[ii, i,k,j,z,"P", rep])) {
            
            countneg_random[ii,i,k,j,z, rep] <- NA#if P is NA, do not count
            
          } else if(
            
            !is.na(dfstore_neg_random[ii,i,k,j,z,"P", rep]) &#P is NOT NA
            
            dfstore_neg_random[ii,i,k,j,z,"P", rep] <= 0.05 &#P is significant
            
            dfstore_neg_random[ii,i,k,j,z,"corr", rep]<0){ #correlation is negative
            
            countneg_random[ii,i,k,j,z, rep]<-1 #count a negative correlation
            
          } else {
            countneg_random[ii,i,k,j,z, rep]<- 0
          }
        }}}}
  
}


  print(rep)} #ends rep loop

save(list = ls(), file = xzfile("Data/6_test_for_comp_COMPRESS.RData.xz", compression = 9))


#############################################################################
############################################################################

          #Validation process--is there evidence for compensation?#
              # Restart here after initial randomization process

###########################################################################
###########################################################################
 load("Data/6_test_for_comp_COMPRESS.RData.xz")

proppos_random <- propneg_random<- matrix(NA, nrow= length(trmts), ncol= (rando))
upper_pos <- upper_neg <- rep(NA, length(trmts))

for (ii in 1:length(trmts)){
  
  for (rep in 1:rando){
        proppos_random[ii, rep] <- sum(countpos_random[ii,,,,, rep], na.rm=T)/sum(!is.na(countpos_random[ii,,,,, rep]))#Proportion positive corr out of total corr's that aren't NA
        propneg_random[ii, rep] <- sum(countneg_random[ii,,,,, rep], na.rm=T)/sum(!is.na(countneg_random[ii,,,,, rep]))#Proportion neg corr out of total corr's that aren't NA
  if( is.na(proppos_random[ii, rep] )){stop("pos is na")}   }
  upper_pos[ii] <- quantile(proppos_random[ii,], c(0.05)) # Aleah-- should this not be 0.05/2 ??? --- No, we want the lower 5th percentile, not the lower 2.5 percentile
      upper_neg[ii] <- quantile(propneg_random[ii,], c(0.95)) 
      } 

# checking evidence for compensation ------

for (ii in 1:length(trmts)){

# is there evidence via the proportion of positive being less than you expect by chance?
positive_result <- ifelse (proppos[ii] <= upper_pos[ii],"YES", "NO")# Aleah-- I changed this in these 2 lines to a less than or equal to sign, as I think that's right-- can you confirm? #Allison--looking at this, he says "as high as or higher" for negative, but "the same methodology" for positive, so I guess this is right?
#is there evidence via the proportion of negative corrs being more than you expect by chance?
negative_result <- ifelse(propneg[ii] >= upper_neg[ii], "YES","NO")
print(trmts[ii])
print(positive_result)
print(negative_result)
}

save(list = ls(), file = xzfile("Data/6_test_for_comp_COMPRESS.RData.xz", compression = 9))

# ##################################################################################
# ###                           GRAPHING                                     ######
# #################################################################################

load("Data/6_test_for_comp_COMPRESS.RData.xz")

# ##### Grazed (CONT)
# 

props_Hdf_pos<- as.data.frame(proppos_random[1,])
colnames(props_Hdf_pos) <- c("pos")
upper_pos_H<-upper_pos[1]
proppos_H<- proppos[1]

# #positive
p1<-ggplot()+
  geom_histogram(data=props_Hdf_pos, aes(x=pos, y=..count../sum(..count..)), fill="lightcoral", bins=10) +
  # geom_histogram(data=props_Hdf %>% filter ( pos>upper_pos_H), aes(x=pos, y=..count../sum(..count..)), fill="lightcoral", bins=10) +
   geom_vline(xintercept=upper_pos_H, linetype="dashed", color="firebrick", size=1.5)+ #5th perc of pos corr
  geom_vline(xintercept=proppos_H, linetype="solid", color="black", size=1.5)+ #actual percentile of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

#negative

props_Hdf_neg<- as.data.frame(propneg_random[1,])
colnames(props_Hdf_neg) <- c("neg")
upper_neg_H<-upper_neg[1]
propneg_H<- propneg[1]

p2<-ggplot(data=props_Hdf_neg)+
  geom_histogram(aes(x=neg, y=..count../sum(..count..)), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=upper_neg_H, linetype="dashed", color="royalblue3", size=1.5)+ #95th perc of neg corr
  geom_vline(xintercept=propneg_H, linetype="solid", color="black", size=1.5)+ #actual perc of neg corr
  xlab("Proportion of Negative Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

#####Ungrazed (LMH)

props_NHdf_pos<- as.data.frame(proppos_random[2,])
colnames(props_NHdf_pos) <- c("pos")
upper_pos_NH<-upper_pos[2]
proppos_NH<- proppos[2]

#positive
p3<-ggplot(data=props_NHdf_pos)+
  geom_histogram(aes(x=pos, y=..count../sum(..count..)), fill="lightcoral", bins=10) +
  geom_vline(xintercept=upper_pos_NH, linetype="dashed", color="firebrick", size=1.5)+#5th perc of pos corr
  geom_vline(xintercept=proppos_NH, linetype="solid", color="black", size=1.5)+#actual perc of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)
#negative
props_NHdf_neg<- as.data.frame(propneg_random[2,])
colnames(props_NHdf_neg) <- c("neg")
upper_neg_NH<-upper_neg[2]
propneg_NH<- propneg[2]

p4<-ggplot(data=props_NHdf_neg)+
  geom_histogram(aes(x=neg, y=..count../sum(..count..)), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=upper_neg_NH, linetype="dashed",  color="royalblue3",  size=1.5)+#95th perc of neg corr
  geom_vline(xintercept=propneg_NH, linetype="solid", color="black", size=1.5)+#actual perc of neg corr
  xlab("Proportion of Negative Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

####MEGA
props_MEGAdf_pos<- as.data.frame(proppos_random[3,])
colnames(props_MEGAdf_pos) <- c("pos")
upper_pos_MEGA<-upper_pos[3]
proppos_MEGA<- proppos[3]


#positive
p5<-ggplot(data=props_MEGAdf_pos)+
  geom_histogram(aes(x=pos, y=..count../sum(..count..)), fill="lightcoral", bins=10) +
  geom_vline(xintercept=upper_pos_MEGA, linetype="dashed", color="firebrick", size=1.5)+#5th perc of pos corr
  geom_vline(xintercept=proppos_MEGA, linetype="solid", color="black", size=1.5)+#actual percentile of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

props_MEGAdf_neg<- as.data.frame(propneg_random[3,])
colnames(props_MEGAdf_neg) <- c("neg")
upper_neg_MEGA<-upper_neg[3]
propneg_MEGA<- propneg[3]

#negative
p6<-ggplot(data=props_MEGAdf_neg)+
  theme_bw(base_size=10)+
  theme(axis.text.x  = element_text(hjust = 0.7))+
  geom_histogram(aes(x=neg, y=..count../sum(..count..)), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=upper_neg_MEGA, linetype="dashed", color="royalblue3", size=1.5)+#95th perc of neg corr
  geom_vline(xintercept=propneg_MEGA, linetype="solid", color="black", size=1.5)+#actual perc of neg corr
  xlab("Proportion of Negative Correlations")+
  ylab("Frequency")

# ###MESO
props_MESOdf_pos<- as.data.frame(proppos_random[4,])
colnames(props_MESOdf_pos) <- c("pos")
upper_pos_MESO<-upper_pos[4]
proppos_MESO<- proppos[4]

#positive
p7<-ggplot(data=props_MESOdf_pos)+
  geom_histogram(aes(x=pos, y=..count../sum(..count..)), fill="lightcoral", bins=10) +
  geom_vline(xintercept=upper_pos_MESO, linetype="dashed", color="firebrick", size=1.5)+#5th perc of pos corr
  geom_vline(xintercept=proppos_MESO, linetype="solid", color="black", size=1.5)+#actual percentile of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

props_MESOdf_neg<- as.data.frame(propneg_random[4,])
colnames(props_MESOdf_neg) <- c("neg")
upper_neg_MESO<-upper_neg[4]
propneg_MESO<- propneg[4]

#negative
p8<-ggplot(data=props_MESOdf_neg)+
  geom_histogram(aes(x=neg, y=..count../sum(..count..)), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=upper_neg_MESO, linetype="dashed", color="royalblue3", size=1.5)+#95th perc of neg corr
  geom_vline(xintercept=propneg_MESO, linetype="solid", color="black", size=1.5)+#actual perc of neg corr
  xlab("Proportion of Negative Correlations")+
  ylab("Frequency")+
  theme_bw(base_size=10)

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow=4)
grid.arrange(p3,p4,p7,p8,p5,p6, p1,p2, nrow=4) #actual order for ppr
#no herbs, meso, mega, all herbs

jpeg(filename = "25-2-20_null_dist.jpg", width = 7.5, height = 7.5, units = "in", res = 300)
grid.arrange(p3,p4,p7,p8,p5,p6, p1,p2, nrow=4) #actual order for ppr
dev.off()

#no lmh, meso, mega, control

#ONLY presence and absence of herbivory


props_Hdf_pos$Herbtrt <- "All Herbivores"

props_NHdf_pos$Herbtrt <- "No Herbivores"

props_Hdf <-cbind(props_Hdf_pos,props_Hdf_neg )
props_NHdf <- cbind(props_NHdf_pos, props_NHdf_neg )

allcorrs <- rbind(props_Hdf, props_NHdf)


# ################HERBS
#positive
p1<-ggplot()+
  geom_histogram(data=props_Hdf %>% filter ( pos<=upper_pos_H), aes(x=pos), fill="firebrick", bins=10) +
  geom_histogram(data=props_Hdf %>% filter ( pos>upper_pos_H), aes(x=pos), fill="lightcoral", bins=10) +
  geom_vline(xintercept=proppos_H, linetype="dashed", color="black", size=1.5)+ #actual percentile of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Count")+
  theme_bw(base_size=16)

#negative
p2<-ggplot()+
  geom_histogram(data=props_Hdf %>% filter ( neg>=upper_neg_H), aes(x=neg), fill="royalblue3", bins=10) +
  geom_histogram(data=props_Hdf %>% filter ( neg<upper_neg_H), aes(x=neg), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=propneg_H, linetype="dashed", color="black", size=1.5)+ #actual percentile of pos corr
  xlab("Proportion of Negative Correlations")+
  ylab("Count")+
  theme_bw(base_size=16)
# 
# 
# ################No HERBS
#positive
p3<-ggplot()+
  geom_histogram(data=props_NHdf %>% filter ( pos<=upper_pos_NH), aes(x=pos), fill="firebrick", bins=10) +
  geom_histogram(data=props_NHdf %>% filter ( pos>upper_pos_NH), aes(x=pos), fill="lightcoral", bins=10) +
  geom_vline(xintercept=proppos_NH, linetype="dashed", color="black", size=1.5)+ #actual percentile of pos corr
  xlab("Proportion of Positive Correlations")+
  ylab("Count")+
  theme_bw(base_size=16)

#negative
p4<-ggplot()+
  geom_histogram(data=props_NHdf %>% filter ( neg>=upper_neg_NH), aes(x=neg), fill="royalblue3", bins=10) +
  geom_histogram(data=props_NHdf %>% filter ( neg<upper_neg_NH), aes(x=neg), fill="lightsteelblue", bins=10) +
  geom_vline(xintercept=propneg_NH, linetype="dashed", color="black", size=1.5)+ #actual percentile of pos corr
  xlab("Proportion of Negative Correlations")+
  ylab("Count")+
  theme_bw(base_size=16)

# save 1200*650 or 12 * 6.5
grid.arrange(p1,p2,p3,p4, nrow=2)
