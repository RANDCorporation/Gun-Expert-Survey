########################################################################################
### CompareModelsPublic.r ###

## last update: 04/10/21
## authors:     Andrew R.Morral
##
## description: Compares models of gun policy preferences with and without interactions by expert group
## code updated for PUF file November 2021 by Robin Beckman

##############################################################

library(tidyr)
library(reshape2)
library(patchwork)
library(polycor)
library(MASS)
library(ggplot2)

########################################################
##### Functions
########################################################

### Winsorize predictors ####
ends = c(95,105)
winsor<-function(u){
  u[which(u>ends[2])]<-ends[2]
  u[which(u<ends[1])]<-ends[1]
  return(u)
}

### Compare predicted with actual policy ratings on Q13 ###
compare.dist = function(u){   
  matchglm = ifelse(dat$Q13==u,1,0) # Prediction matches true value
  nearglm = ifelse(abs(as.numeric(dat$Q13)-as.numeric(u))<2,1,0) # Prediction within 1 scale point
  tot.match = sum(matchglm)/length(dat$cat2) #Mean perfect match
  cat0.match= sum(matchglm[dat$cat2==0])/sum(dat$cat2 == 0) #Mean perfect match "permissive"
  cat1.match = sum(matchglm[dat$cat2==1])/sum(dat$cat2 == 1) #Mean perfect match "restrictive"
  tot.near = sum(nearglm)/length(dat$cat2) #Mean within 1 scale point
  cat0.near = sum(nearglm[dat$cat2==0])/sum(dat$cat2 == 0) #Mean within 1 "permissive"
  cat1.near = sum(nearglm[dat$cat2==1])/sum(dat$cat2 == 1) #Mean within 1 "restrictive"
  mle = cbind.data.frame(Total = c(tot.match, tot.near), "Permissive" = c(cat0.match, cat0.near), "Restrictive"= c(cat1.match, cat1.near))
  mle$match = c("Exact match", "Within one")
  mle = mle[,c(4,1:3)]
  return(mle)
}

### Return mean and interquartile range for policy ratings and predicted ratings ###
summarize = function(u,true.rat=FALSE, set){
  temp = data.frame(matrix(u,nrow=length(u)))
  temp$class = dat$cat2
  temp$policy = dat$policy
  # Calculate mean rating by policy and expert group 
  summ1 = tapply(as.numeric(temp[,1]), temp[,c("class","policy")], mean)
  # Create separate data frames for permissive and restrictive groups
  ret0 = data.frame(summ1[1,])
  ret0$class = "Permissive"
  ret1 = data.frame(summ1[2,])
  ret1$class = "Restrictive"
  ret0$policy = ret1$policy = colnames(summ1)
  names(ret1)[1]=names(ret0)[1]="Rating"
  # If true ratings, calculate inter-quartile ranges by policy and expert group
  if (true.rat==TRUE){
    summ1lo = tapply(as.numeric(temp[,1]), temp[,c("class","policy")], quantile,probs=c(.25))
    summ1hi = tapply(as.numeric(temp[,1]), temp[,c("class","policy")], quantile,probs=c(.75))
    ret0$LCI50 = summ1lo[1,] ; ret1$LCI50 = summ1lo[2,]
    ret0$UCI50 = summ1hi[1,] ; ret1$UCI50 = summ1hi[2,]
  }
  # Stack data.frames for permissive and restrictive groups
  ret2 = rbind.data.frame(ret0,ret1)
  ret2$set = set
  return(ret2)
}

########################################################
##### Assemble Data 
########################################################

cant.miss = 6 #Exclude policies for which there are 6 or more imputed predictors

x.wide<-read.csv("gun_policy_puf.csv", stringsAsFactors = FALSE)


### Index of policy effects (predictors) ###

preds=""
for (i in 1:19) preds = c(preds, (3+(i-1)*11):(13+(i-1)*11))
preds = as.numeric(preds[-1])

### Index of q13 missing value indicators
q13missIndicator = intersect(grep("Q13",names(x.wide)),grep("ms_",names(x.wide)))

### Index of variables listing number of policy effects rated per policy 
outcomes_ansInd = grep("outcomes_ans_", names(x.wide))

### Reduce data frame to just model variables ###
x.wide.small = x.wide[,c(1,2,preds,q13missIndicator, outcomes_ansInd)]


#### Convert wide data frame to long ###########
x.long = gather( x.wide.small, variable, effect, OP1Q1:NP5Q13, 
                 ms_OP1Q13:ms_NP5Q13,outcomes_ans_OP1:outcomes_ans_NP5)  

### Create index for start of policy effect number in outcome name ###
locationQ = regexpr("Q",x.long[,3])
for (i in names(attributes(locationQ))) attr(locationQ,i) = NULL
locationP = regexpr("P",x.long[,3])
for (i in names(attributes(locationP))) attr(locationP,i) = NULL

### Construct policy and outcome variables ###
x.long$policy = ifelse(locationQ != -1, 
                       substr(x.long[,3],locationP-1,locationQ-1),
                       substr(x.long[,3],locationP-1,nchar(x.long[,3])))
x.long$outcome = ifelse(locationQ != -1,
                       substr(x.long[,3],locationQ, nchar(x.long[,3])),
                       "Effects Rated")
x.long$outcome = ifelse(substr(x.long[,3],1,2)=="ms",
                        paste0(x.long$outcome,"miss"),
                        x.long$outcome)
                        
### Create data frame with one row per policy and respondent ###
dat = dcast(x.long, participant_ID + policy +cat2 ~ outcome , value.var="effect", fun.aggregate=sum)
dat = dat[, c(1,2,3,5,10:16,6,7,8,9,4)] # Reorder variables

### Winsorize quantitative IVs and normalize all IVs ###
dat[,c(4:11)]<-data.frame(apply(dat[,c(4:11)],2,winsor))
dat[,4:13] = scale(dat[4:13])

### Make interaction terms ###
dat$cat2 = dat$cat2-1
clus = dat$cat2 - mean(dat$cat2)
dat[,paste0("int",1:10)] = apply(dat[,4:13],2, function(u){
  clus * (u - mean(u))})

#### Drop cases with missing Q13 or with fewer than cant.miss predictions ####
dat = dat[ ! (dat$Q13miss==1 | dat[,"Effects Rated"]<= 10-cant.miss),]


########################################################
##### Model Policy Favorability Ratings, Q13 
########################################################

### Model with interactions ###
dat$Q13 = as.factor(dat$Q13)
ol1 = polr(Q13 ~ ., dat[,c(4:14,16:25)])
summary(ol1)
pred.glm1 = predict(ol1, dat[,c(4:14,16:25)])
compare.dist(pred.glm1)

### Model without interactions ###
ol2 = polr(Q13 ~ ., dat[,c(4:14)])
summary(ol2)
pred.glm2 = predict(ol2, dat[,c(4:14)])
compare.dist(pred.glm2)

### Check correlations between actual and predicted ratings ####

dat1 = cbind(dat$Q13,pred.glm1)
dat2 = cbind(dat$Q13,pred.glm2)

## Model with Interactions ##
# Full sample correlations
polychor(dat1[,1],dat1[,2], ML=TRUE) 
# Permissive group correlations
polychor(dat1[dat$cat2==0,1],dat1[dat$cat2==0,2], ML=TRUE) 
# Restrictive group correlations
polychor(dat1[dat$cat2==1,1],dat1[dat$cat2==1,2], ML=TRUE)

## Model without interactions ##
# Full sample correlations
polychor(dat2[,1],dat2[,2], ML=TRUE) 
# Permissive group correlations
polychor(dat2[dat$cat2==0,1],dat2[dat$cat2==0,2], ML=TRUE)
# Restrictive group correlations
polychor(dat2[dat$cat2==1,1],dat2[dat$cat2==1,2], ML=TRUE)

########################################################
##### Visualize Model Predictions 
########################################################

### Figure labels ###
plist <- c("1. Universal background checks",
           "2. Ban on sale of 'assault weapons'",
           "3. Stand your ground law",
           "4. Mental health prohibitions",
           "5. Report lost or stolen firearms",
           "6. License to purchase required",
           "7. Report firearms sales",
           "8. Child access-prevention law",
           "9. Surrender of firearms by prohibited possessors",
           "10. Taxes on firearms and ammunition",
           "11. Minimum age requirements",
           "12. Permitless carry",
           "13. Ten-day waiting period to purchase",
           "14. Elimination of gun-free zones",
           "15. Extreme-risk protection order",
           "16. Domestic violence restraining order",
           "17. Arming school personnel",
           "18. Gun purchase limits",
           "19. Prosecution of prohibitted possessors")


polist = c("OP1","OP2", "OP3", "OP4", "OP5", "OP6", "OP8", "OP9","OP11", "OP12","OP13","OP14",
           "OP15","OP16","NP1","NP2","NP3","NP4","NP5")

with.ints = summarize(pred.glm1,set="With interactions")
without.ints = summarize(pred.glm2,set="Without interactions")
true.rats =  summarize(dat$Q13,true.rat = TRUE, set="True ratings")
rats = rbind.data.frame(with.ints, without.ints, true.rats[,c(1:3,6)])
rats$LCI50[rats$set=="True ratings"] = true.rats[,4]
rats$UCI50[rats$set=="True ratings"] = true.rats[,5]
rats$policy=factor(rats$policy, levels = polist[19:1], labels=plist[19:1])
rats$set = factor(rats$set, levels = c("Without interactions", "With interactions", "True ratings"))

perm.set = rats[rats$class=="Permissive",]
g1 = ggplot(perm.set, aes(x=Rating, y=policy, color=set))+
  geom_pointrange(aes(xmin=LCI50, xmax=UCI50), position = position_dodge(width = 0.5), size=.4)+
  theme(aspect.ratio=8/3, legend.position="none",
        plot.margin =  margin(rep(0,4)),
        axis.text.y = element_text(size = 12)) +
  ggtitle("Permissive class") +
  scale_color_manual(values=c("goldenrod","darkred","olivedrab4")) +
  scale_x_continuous(limits = c(1,5)) + 
  guides(col = guide_legend(reverse = TRUE)) + 
  xlab("Mean rating") + ylab("")

perm.set = rats[rats$class=="Restrictive",]
g2 = ggplot(perm.set, aes(x=Rating, y=policy, color=set))+
  geom_pointrange(aes(xmin=LCI50, xmax=UCI50), position = position_dodge(width = 0.5), size=.4)+
  theme(aspect.ratio=8/3,  axis.text.y = element_blank(), legend.position="right", legend.title=element_blank(),
        plot.margin = margin(rep(0,4)), legend.text = element_text(size=10))+
  ggtitle("Restrictive class") +
  scale_color_manual(values=c("goldenrod","darkred","olivedrab4")) +
  scale_x_continuous(limits = c(1,5)) + 
  guides(col = guide_legend(reverse = TRUE)) + 
  xlab("Mean rating") + ylab("")
g3 = g1|g2   #Make sure plot viewer is wide enough
g3
ggsave("ratings_puf.png")

