## This code is the implementation of the data analysis in the following paper: 
## "Predicting Corticosteroid-Free Biologic Remission with Vedolizumab in Crohn's Disease" (2017)
## Authors: Akbar K. Waljee, Boang Liu, Kay Sauder, Ji Zhu, Shail M. Govani, Ryan W. Stidham, Peter D.R. Higgins
## Author of the code: Boang Liu
## Date: April 17 2017

library("randomForest")
library("pROC")

# Read in the data. y contains outcome. x_bsl contains baseline variables. 
# x_cs contains week 6 labs and other baseline variables. x_lgt contains x_cs and longitudinal variables.

n = nrow(y) #n is the sample size
nTree = 1000 #number of trees for random forest

# Split the data into training and testing for the 50 replications
train_all = matrix(0,50,floor(n*0.7))
test_all = matrix(0,50,n-floor(n*0.7))
set.seed(5)
for(m in 1:50){
  train_all[m,] = sample(n,floor(n*0.7)) 
  test_all[m,] = (1:n)[-train_all[m,]] 
}

##########################################################################
# Week6 model
var_lgt = c("USUBJID","CALPRO_Week0","CRP_DIFF","PC_DIFF") #CALPRO_Week0: baseline FCP, CRP_DIFF: slope of CRP, PC_DIFF: slope of VDZ level
set.seed(35)
seeds = sample(2000,50)
auc_cs = rep(0,50)
for(m in 1:50){
  print(m)
  train = train_all[m,]
  test = test_all[m,]
  
  var_cs = c(names(x_cs)[-1],var_lgt[-1])
  
  set.seed(seeds[m])
  rf_cs = randomForest(x=x_lgt[train,var_cs],y=as.factor(y[train,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
  rf_cs_pred = predict(rf_cs,newdata=x_lgt[test,var_cs],type="prob")[,2]
  rf_cs_roc = roc(y[test,"OUTCOME"],rf_cs_pred)
  auc_cs[m] = rf_cs_roc$auc
}
auc_cs
mean(auc_cs) #average AuROC over 50 replications

ind_best = which.min(abs(auc_cs-mean(auc_cs))) ##the split that has closest AuROC to the average
train = train_all[ind_best,]
test = test_all[ind_best,]
var_cs = c(names(x_cs)[-1],var_lgt[-1])
set.seed(seeds[ind_best])
rf_cs = randomForest(x=x_lgt[train,var_cs],y=as.factor(y[train,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
rf_cs_pred = predict(rf_cs,newdata=x_lgt[test,var_cs],type="prob")[,2]
rf_cs_roc = roc(y[test,"OUTCOME"],rf_cs_pred,ci=T)
rf_cs_roc #ROC

rf_cs_roc$direction
co_cs = coords(rf_cs_roc,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","1-npv","tp","tn","fp","fn"))
co_cs #get the best cutoff and other quantities

# Variable importance
var_cs = c(names(x_cs)[-1],var_lgt[-1])
set.seed(20)
rf_cs = randomForest(x=x_lgt[,var_cs],y=as.factor(y[,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
rf_cs_imp = rf_cs$importance[,4]*100/sum(rf_cs$importance[,4])
rf_cs_imp = sort(rf_cs_imp) #variable importance scores

# Partial dependence plots
x = x_lgt[,var_cs]
imp_names = names(sort(rf_cs_imp, decreasing = TRUE))
impdata = x[imp_names]
length_v = rep(0, length(rf_cs_imp))
uniqueV = list()
for(i in 1:length(length_v)){
  uniqueV[[i]] = unique(impdata[,i])
  length_v[i] = length(uniqueV[[i]])
  if(length_v[i]>50){
    uniqueV[[i]] = quantile(impdata[,i], seq(0.01, 0.99, 0.02))
    length_v[i] = length(uniqueV[[i]])
  }
}

est = list()
for(i in 1:length(length_v)){
  print(i)
  est[[i]] = rep(0, length_v[i])
  for(j in 1:length_v[i]){
    newdata_j = x
    newdata_j[imp_names[i]] = uniqueV[[i]][j]
    est[[i]][j] = mean(predict(rf_cs, newdata_j, type = 'prob')[,2])
    gc()
  }
}

# given i, make partial dependence plot for the ith predictor
plot(uniqueV[[i]], est[[i]], xlab = imp_names[i], ylab = "Success Probability")

#################################################################################
# Baseline model
set.seed(10)
seeds = sample(2000,50)
auc_bsl = rep(0,50)
for(m in 1:50){
  print(m)
  train = train_all[m,]
  test = test_all[m,]
  
  var_bsl = names(x_bsl)[-1]
  
  set.seed(seeds[m])
  rf_bsl = randomForest(x=x_bsl[train,var_bsl],y=as.factor(y[train,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
  rf_bsl_pred = predict(rf_bsl,newdata=x_bsl[test,var_bsl],type="prob")[,2]
  rf_bsl_roc = roc(y[test,"OUTCOME"],rf_bsl_pred)
  auc_bsl[m] = rf_bsl_roc$auc
}
auc_bsl
mean(auc_bsl) #average AuROC over 50 replications 
sd(auc_bsl)

# use the split that is selected by the week 6 model
train = train_all[ind_best,]
test = test_all[ind_best,]
var_bsl = names(x_bsl)[-1]
set.seed(10)
rf_bsl = randomForest(x=x_bsl[train,var_bsl],y=as.factor(y[train,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
rf_bsl_pred = predict(rf_bsl,newdata=x_bsl[test,var_bsl],type="prob")[,2]
rf_bsl_roc = roc(y[test,"OUTCOME"],rf_bsl_pred,ci=T)
rf_bsl_roc #ROC

rf_bsl_roc$direction # "<" 
co_bsl = coords(rf_bsl_roc,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","1-npv","tp","tn","fp","fn"))
co_bsl #get the best cutoff and other quantities

# Variable importance
set.seed(10)
rf_bsl = randomForest(x=x_bsl[,-1],y=as.factor(y[,"OUTCOME"]),ntree=nTree,importance=T,replace=F)
rf_bsl_imp = rf_bsl$importance[,4]*100/sum(rf_bsl$importance[,4])
rf_bsl_imp = sort(rf_bsl_imp) #variable importance scores

# Partial dependence plots
x = x_bsl[,var_bsl]
imp_names = names(sort(rf_bsl_imp, decreasing = TRUE))
impdata = x[imp_names]
length_v = rep(0, length(rf_bsl_imp))
uniqueV = list()
for(i in 1:length(length_v)){
  uniqueV[[i]] = unique(impdata[,i])
  length_v[i] = length(uniqueV[[i]])
  if(length_v[i]>50){
    uniqueV[[i]] = quantile(impdata[,i], seq(0.01, 0.99, 0.02))
    length_v[i] = length(uniqueV[[i]])
  }
}

est = list()
for(i in 1:length(length_v)){
  print(i)
  est[[i]] = rep(0, length_v[i])
  for(j in 1:length_v[i]){
    newdata_j = x
    newdata_j[imp_names[i]] = uniqueV[[i]][j]
    est[[i]][j] = mean(predict(rf_bsl, newdata_j, type = 'prob')[,2])
    gc()
  }
}

# given i, make partial dependence plot for the ith predictor
plot(uniqueV[[i]], est[[i]], xlab = imp_names[i], ylab = "Success Probability")

###############################################################################
# Simpler model: (HGB*ALB*VDZ)/(CRP*Weight) at week 6
x = rep(NA,nrow(x_lgt))
for(i in 1:length(x)){
  if(x_lgt$CRP_Week6[i]*x_lgt$WEIGHT[i]!=0){
    x[i] = (x_lgt$HGB_Week6[i]*x_lgt$ALB_Week6[i]*x_lgt$PC_Week6[i])/(x_lgt$CRP_Week6[i]*x_lgt$WEIGHT[i])
  }
}
scs_roc = roc(y[,"OUTCOME"],x,ci=T)
scs_roc #ROC

scs_roc$direction #"<"
co_scs = coords(scs_roc,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","1-npv","tp","tn","fp","fn"))
co_scs #get the best cutoff and other quantities

#############################################################################