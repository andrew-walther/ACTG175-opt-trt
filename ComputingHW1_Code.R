#BIOS 740 Computing HW 1 Code
#Andrew Walther
#9/15/2021

##Data & Packages
set.seed(420)
#for ML method implementation
library(DynTxRegime)
library(randomForest)
#data manipulation/visualization
library(tidyverse)
#summary statistics
library(mosaic)
#cross-validation
library(caret)
#AIDS Clinical Trials Group Study 175 data (ACTG175)
library(speff2trial)
data(ACTG175)

#select observations in only arms 1,2 & eliminate unused factors as directed (included 'arms')
dataCD496 <- ACTG175 %>% filter(!is.na(cd496), arms %in% 1:2) %>% 
        select(-pidnum,-str2, -offtrt, -cd420, -r, -cd820, -cens,
               -days, -treat, -oprior, -z30, -zprior, -preanti)
##########################################################################
##Preprocessing & TEST/TRAIN Sets
#remove observations where 'cd496 = 0' (drops 1 observation)
dataCD496 <- dataCD496[dataCD496$cd496 != 0, ]
#add column for log transform of CD4 cell count (check if this is necessary)
data_log_cd496 <- dataCD496 %>% mutate(log_cd496 = log(cd496))

#Create 10 folds of the dataset and add to TEST/TRAIN lists
folds <- createFolds(dataCD496$cd496, 10)
TEST_DATA <- list()
TRAIN_DATA <- list()
for(i in 1:10){
        TEST_DATA[[i]] <- dataCD496[folds[[i]],]
        TRAIN_DATA[[i]] <- dataCD496[-folds[[i]],]}
##########################################################################
##Exploratory Analysis / Feature Selection
#summary stats for CD490 & log(CD4)
favstats(dataCD496$cd496)
favstats(data_log_cd496$log_cd496)

#hist of CD4 Cell counts, possible right skew in the data
dataCD496 %>% ggplot() + geom_histogram(aes(x=cd496),binwidth = 50)
+labs(title = "Histogram of CD4 cell counts",x = "CD4 Cell Count (cells/mm^3)")
#hist of log(CD4 cell), now there's left skew in the data (check for normality)
data_log_cd496 %>% ggplot() + geom_histogram(aes(x=log(cd496)),binwidth = 0.1)
+labs(title = "Histogram of log CD4 cell counts",x = "CD4 Cell Count (cells/mm^3)")

#K-S test for normality - H0: data follows Normal distribution (p>0.05 so confirm CD496 is normal)
ks.test(dataCD496$cd496, "pnorm", mean=mean(dataCD496$cd496), sd=sd(dataCD496$cd496))

#Feature Selection: Train RPart Model to compute variable importance to decide what to keep
rPartMod <- train(cd496 ~ ., data=dataCD496, method="rpart")
rpartImp <- varImp(rPartMod)
print(rpartImp)
##########################################################################
##Random Forest
#Initialize lists to hold mean predictions
RF_PREDS1 <- list()
RF_PREDS2 <- list()
#loop over 10 data folds with Random forest method
for(i in 1:10){
        #filter just observations in arm 1
        TRAIN_DATA1 <- TRAIN_DATA[[i]][TRAIN_DATA[[i]]$arms==1,]
        #filter just observations in arm 2
        TRAIN_DATA2 <- TRAIN_DATA[[i]][TRAIN_DATA[[i]]$arms==2,]
        #train RF model on just arm 1 observations
        rf.obj1 <- randomForest(cd496 ~ ., data = TRAIN_DATA1, 
                                importance = TRUE, proximity = TRUE)
        #train RF model on just arm 2 observations
        rf.obj2 <- randomForest(cd496 ~ ., data = TRAIN_DATA2, 
                                importance = TRUE, proximity = TRUE)
        #predict on RF models with full set of test data
        RF_PREDS1[[i]] <- mean(predict(rf.obj1, newdata = TEST_DATA[[i]]))
        RF_PREDS2[[i]] <- mean(predict(rf.obj2, newdata = TEST_DATA[[i]]))}
#compute mean & sd for each fold prediction mean (larger mean will serve as the value estimate)
mean.pred1 <- mean(as.numeric(RF_PREDS1))
mean.pred2 <- mean(as.numeric(RF_PREDS2))
se.pred1 <- sd(as.numeric(RF_PREDS1))
se.pred2 <- sd(as.numeric(RF_PREDS2))
##########################################################################
##Augmented Inverse Probability Weighted Estimator
aipw.estimates <- list()
opt.Tx <- list()
for(i in 1:10){
        #specify propensity model (moPropen)
        moPropen <- buildModelObj(model = ~ 1, solver.method = 'glm', 
                                  solver.args = list('family'='binomial'), 
                                  predict.method = 'predict.glm', 
                                  predict.args = list(type='response'))
        #specify outcome model (moMain)
        moMain <- buildModelObj(model = ~age+wtkg+hemo+homo+drugs+
                                        karnof+race+gender+strat+symptom+cd40+cd80
                                , solver.method = 'lm')
        #specify regimes
        data <- TRAIN_DATA[[i]]
        regime <- function(eta1, data) {
                tst <- {data$cd496 > eta1}
                rec <- rep(1, nrow(x=data))
                rec[!tst] <- 2
                return( rec )}
        #inverse probability weighted estimator object
        fit.AIPW <- optimalSeq(moPropen = moPropen, moMain = moMain, 
                               regimes = regime, data = data, 
                               response = data$cd496, txName = 'arms',
                               Domains = cbind(1,1062),
                               pop.size = 100, starting.values = 1)
        #value function estimates & test data treatment assignments
        aipw.estimates[[i]] <- as.numeric(estimator(fit.AIPW))
        opt.Tx[[i]]<- optTx(fit.AIPW, newdata = TEST_DATA[[i]])}
#compute mean & sd of estimated value functions (mean = 355.4116, se = 2.482627)
mean.aipw <- mean(as.numeric(aipw.estimates))
se.aipw <- sd(as.numeric(aipw.estimates))
#count patients assigned to trt 1 or 2 and compute majority proportion
trt1.total <- 0
trt2.total <- 0
for(i in 1:10){
        trt1 <- count(opt.Tx[[i]]$optimalTx == 1)
        trt2 <- count(opt.Tx[[i]]$optimalTx == 2)
        trt1.total <- trt1.total + trt1
        trt2.total <- trt2.total + trt2}
print(trt1.total)
print(trt2.total)
print(trt2.total/(trt1.total+trt2.total))
##########################################################################
##Residual Weighted Learning
#Initialize lists to hold value estimates & optimal treatment outcomes
rwl.estimates <- list()
opt.Tx <- list()
for(i in 1:10){
        #specify propensity model w/ all parameters (moPropen)
        moPropen <- buildModelObj(model = ~age+wtkg+hemo+homo+drugs+
                                          karnof+race+gender+strat+symptom+cd40+cd80,
                                  solver.method = 'glm',
                                  solver.args = list('family'='binomial'),
                                  predict.method = 'predict.glm',
                                  predict.args = list(type='response'))
        #specify outcome model w/ all parameters (moMain)
        moMain <- buildModelObj(model = ~age+wtkg+hemo+homo+drugs+
                                        karnof+race+gender+strat+symptom+
                                        cd40+cd80, solver.method = 'lm')
        #residual weighted learning object (cvfolds = 1 for efficiency & include all covariates)
        fit.rwl <- rwl(moPropen = moPropen, moMain = moMain, data = TRAIN_DATA[[i]], 
                       reward = TRAIN_DATA[[i]]$cd496, txName = 'arms', 
                       regime = ~age+wtkg+hemo+homo+drugs+
                               karnof+race+gender+strat+symptom+cd40+cd80, cvFolds = 1,
                       responseType = 'continuous')
        #value function estimates
        rwl.estimates[[i]] <- as.numeric(estimator(fit.rwl))
        #optimal treatment assignments from test data
        opt.Tx[[i]]<- optTx(fit.rwl, newdata = TEST_DATA[[i]])}
#compute mean & sd of estimated value functions (mean = 369.8081, se = 12.60558)
mean.rwl <- mean(as.numeric(rwl.estimates))
se.rwl <- sd(as.numeric(rwl.estimates))
#count patients assigned to trt 1 or 2 and compute majority proportion
trt1.total <- 0
trt2.total <- 0
for(i in 1:10){
        trt1 <- count(opt.Tx[[i]]$optimalTx == 1)
        trt2 <- count(opt.Tx[[i]]$optimalTx == 2)
        trt1.total <- trt1.total + trt1
        trt2.total <- trt2.total + trt2}
print(trt1.total)
print(trt2.total)
print(trt2.total/(trt1.total+trt2.total))