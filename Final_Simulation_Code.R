#### Setup ####

### Load in packages, functions, and data
# Install SSVS package from Github
devtools::install_github("mahmoud-mfahmy/SSVSforPsych")

# Packages
library(plyr)
library(SSVSforPsych)
library(data.table)
library(glmnet)
library(bayestestR)
library(gtools)
library(MASS)


#### Adaptive LASSO Stuff ----

## Mylars function (dependency of adalasso)
mylars<-function (X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE) 
{
  x<-X
  n<-length(y)
  all.folds <- split(sample(1:n),rep(1:k,length=n))
  
  if (use.Gram==TRUE){
    type="covariance"
  }
  if (use.Gram==FALSE){
    type="naive"
  }
  globalfit<-glmnet(x,y,family="gaussian",standardize=normalize,type.gaussian=type,intercept=intercept)
  lambda<-globalfit$lambda
  residmat <- matrix(0, length(lambda), k)
  for (i in seq(k)) {
    omit <- all.folds[[i]]
    fit <- glmnet(x[-omit, ,drop=FALSE], y[-omit],type.gaussian=type,standardize=normalize,family="gaussian",intercept=intercept)
    fit <- predict(fit, newx=x[omit, , drop = FALSE], type = "response", 
                   s = lambda)
    if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
  }
  cv <- apply(residmat, 1, mean)
  cv.lasso<-min(cv)
  cv.error <- sqrt(apply(residmat, 1, var)/k)
  lambda.opt<-lambda[which.min(cv)]
  coefficients=predict(globalfit,type="coefficients",s=lambda.opt)
  inter=coefficients[1]
  coefficients=coefficients[-1]
  names(coefficients)=1:ncol(X)
  object <- list(lambda=lambda,cv=cv,lambda.opt=lambda.opt,cv.lasso=cv.lasso,intercept=inter,coefficients=coefficients)
  invisible(object)
}

##### Adaptive LASSO function 
adalasso<-function(X, y,k=10,use.Gram=TRUE,both=TRUE,intercept=TRUE){
  colnames(X)=1:ncol(X)
  n<-length(y)
  cv.adalasso<-NULL
  globalfit<-mylars(X,y,k=k,use.Gram=use.Gram,normalize=TRUE,intercept=intercept)
  coefficients.lasso=globalfit$coefficients
  intercept.lasso=globalfit$intercept
  cv.lasso<-globalfit$cv.lasso
  lambda<-globalfit$lambda
  lambda.lasso<-globalfit$lambda.opt
  coefficients.adalasso=NULL
  lambda.adalasso<-intercept.adalasso<-NULL
  if (use.Gram==TRUE){
    type="covariance"
  }
  if (use.Gram==FALSE){
    type="naive"
  }
  if (both==TRUE){ 
    # cross-validation for adaptive lasso
    all.folds <- split(sample(1:n),rep(1:k,length=n))
    residmat <- matrix(0, length(lambda), k)
    
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      Xtrain<-X[-omit,,drop=FALSE]
      ytrain<-y[-omit]
      Xtest<-X[omit,,drop=FALSE]
      ytest<-y[omit]
      my.lars<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram,intercept=intercept)
      coef.lasso<-my.lars$coefficients
      weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
      #cat(paste("-- non-zero weights ",length(weights),"\n"))
      if (length(weights)==0){
        residmat[,i]<-mean((mean(ytrain)-ytest)^2)
      }
      if (length(weights)==1){
        residmat[,i]=mean((ytest -my.lars$intercept - Xtest%*%coef.lasso)^2)
      }
      if (length(weights)>1){
        XXtrain <- Xtrain[ , names(weights), drop=FALSE]
        XXtest<-Xtest[ , names(weights), drop=FALSE]
        XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
        XXtest<-scale(XXtest, center=FALSE, scale=weights)
        #cat(paste("ncol of XXtrain: ",ncol(XXtrain),"\n"))
        fit<-glmnet(XXtrain,ytrain,type.gaussian=type,standardize=FALSE,intercept=intercept)
        pred<-predict(fit, newx=XXtest, type = "response",s = lambda)
        if (length(omit) == 1){
          pred <- matrix(pred, nrow = 1)
        }
        residmat[, i] <- apply((ytest - pred)^2, 2, mean)
      }
    }
    cv <- apply(residmat, 1, mean)
    cv.adalasso<-min(cv)
    weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
    coefficients.adalasso<-rep(0,ncol(X))
    names(coefficients.adalasso)<-1:ncol(X)
    if (length(weights)>0){
      XX <- X[ , names(weights), drop=FALSE]
      if ( length(weights)==1 )  XX <- XX/weights        
      else  XX <- scale(XX, center=FALSE, scale=weights)
      if (length(weights)<=1){
        intercept.adalasso=intercept.lasso 
        coefficients.adalasso<-coefficients.lasso
        lambda.adalasso=0
      }
      else{
        fit<-glmnet(XX,y,type.gaussian=type,standardize=FALSE,intercept=intercept)
        lambda.adalasso<-lambda[which.min(cv)]
        coefficients=predict(fit,type="coefficients",s=lambda.adalasso)
        intercept.adalasso<-coefficients[1]
        coefficients.adalasso[names(weights)]<-coefficients[-1]/weights
      }
    }
  }
  return(list(cv.lasso=cv.lasso,lambda.lasso=lambda.lasso,cv.adalasso=cv.adalasso,lambda.adalasso=lambda.adalasso,intercept.lasso=intercept.lasso, intercept.adalasso=intercept.adalasso, coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso))
}


# Make results reproducible ----
set.seed(13513548)

### Prepare data

# read in the data
bigData <- read.csv("sim_p50_t10_corr0.8_first10.csv")


# split data into manageable dataframes
myArray <- with(bigData, 
                split(bigData, 
                      f = do.call(paste, 
                                  bigData[,c("r","N")])))

# subset dataframe
subset <- myArray[1]

#### For loop ####

# Identify the predictor columns we'll use (For p= 50)
predCols <- paste0("X",1:50)

# Create an empty data frame to save values
finalValues <- NULL

### Run the for loop
for (dataset in 1:length(subset)){
  
  # Scale all variables in the dataset
  temp <- scale(subset[[dataset]][,c("y",predCols)])
  subset[[dataset]] <- as.data.frame(cbind(temp,subset[[dataset]][,c("N","r")]))
  
  #### Adaptive LASSO ####
  
  adaResults <- NULL
  
  # Adaptive lasso with 10-fold CV
  adalassoTemp <- adalasso(X = as.matrix(scale(subset[[dataset]][,predCols])),
                           y = scale(subset[[dataset]][,"y"]))
  adalassoTemp$coefficients.adalasso
  
  # Save the coefficients
  adalassoTempResults <- as.data.frame(adalassoTemp$coefficients.adalasso)
  
  # Save the last column
  adalassoSelected <- as.data.frame(adalassoTempResults[,ncol(adalassoTempResults)])
  
  # Add a column to name the variables
  adalassoSelectedNew <- as.data.frame(cbind(adalassoSelected[,1], paste0("X",1:50)))
  
  # Name all columns
  colnames(adalassoSelectedNew) <- c("coefficients", "Preds")
  
  # Make preds numeric
  adalassoSelectedNew[,1] <- as.numeric(as.character(adalassoSelectedNew[,1]))
  
  # Save the nonzero
  number.ada.nonzero <- abs(adalassoSelectedNew$coefficients)  != 0
  save.ada.nonzero <- c(which(number.ada.nonzero))
  save.ada.zero <- c(which(!number.ada.nonzero))
  
  # Assign a 1 to all nonzero.vars
  n.ada.nonzero <- cbind(paste0("X",save.ada.nonzero),
                         rep(1,length(paste0("X",save.ada.nonzero))))
  n.ada.nonzero <- as.data.frame(n.ada.nonzero)
  n.ada.nonzero[,2] <- as.numeric(n.ada.nonzero[,2])
  colnames(n.ada.nonzero) <- c("Preds", "Selected")
  
  # Assign a 0 to all zero.vars
  n.ada.zero <- cbind(paste0("X",save.ada.zero),
                      rep(0,length(paste0("X",save.ada.zero))))
  n.ada.zero <- as.data.frame(n.ada.zero)
  n.ada.zero[,2] <- as.numeric(n.ada.zero[,2])
  colnames(n.ada.zero) <- c("Preds", "Selected")
  if (nrow(n.ada.zero)>0){
    n.ada.zero[,2] <- 0
  }
  # Save results
  adaResults <- rbind(n.ada.nonzero,n.ada.zero)
  
  # Order variable names alphabetically
  adaResults <- adaResults[order(adaResults$Preds),]
  
  # Arrange columns so x1 -> x50 are numerically ordered
  adaResults <- adaResults[mixedorder(as.character(adaResults$Preds)),]
  
  # Save results to a datarame
  adaResults <- cbind(adaResults,
                      adalassoSelectedNew$coefficients)
  rownames(adaResults) <- adaResults$Preds
  

  
  
  #### SSVS ####
  
  
  # Run the SSVS analysis
  ssvs.results <- SSVSforPsych::SSVS(x = scale(subset[[dataset]][,predCols]),
                       y = scale(subset[[dataset]][,"y"]))
  
  # Create a dataframe with all post-burn-in beta balues for saving
  temp.beta.frame <- as.data.frame(ssvs.results[["beta"]])
  
  ### Save the MIP values 
  inc_prob <- as.data.frame(apply(ssvs.results$beta!=0,2,mean))
  inc_prob$var<-rownames(inc_prob)
  names(inc_prob)<-c("MIP","Variable_Name")
  
  ### Save the average betas (including the zero values)
  # Loop
  average.beta <- NULL
  lower.credibility <- NULL
  upper.credibility <- NULL
  for (m in names(temp.beta.frame)){
    average.beta[m] <- mean(temp.beta.frame[,m])
    # 95% credibility interval lower
    lower.credibility[m] <- ci(temp.beta.frame[,m], method = "HDI",ci = .95)[[2]]
    # 95% credibility interval upper
    upper.credibility[m] <- ci(temp.beta.frame[,m], method = "HDI",ci = .95)[[3]]
  }
  
  ### Save the median beta values (including the zero values)
  # Loop
  median.beta <- NULL
  for (m in names(temp.beta.frame)){
    # Obtain mean
    median.beta[m] <- median(temp.beta.frame[,m])
  }
  
  ### Save the average non-zero betas
  # Make the zero values into NAs
  temp.beta.frame.nonzero <- temp.beta.frame
  is.na(temp.beta.frame.nonzero) <- temp.beta.frame.nonzero==0
  
  # Loop 
  average.nonzero.beta <- NULL
  for (m in names(temp.beta.frame.nonzero)){
    # Obtain mean
    average.nonzero.beta[m] <- mean(temp.beta.frame.nonzero[,m], na.rm = TRUE)
  }
  
  ### Save the proportion of runs that produced non-zero values
  number.nonzero <- colSums(temp.beta.frame != 0)
  proportion.nonzero <- number.nonzero/nrow(temp.beta.frame)
  
  #### LASSO min ####
  
  ## Run the lASSO
  # Set the number of folds
  foldid <- sample(rep(1:10, 
                       length.out = length(subset[[dataset]][,"y"])))
  
  # Run the LASSO
  cvlasso = cv.glmnet(x = as.matrix(subset[[dataset]][,predCols]),
                      y = as.matrix(subset[[dataset]][,"y"]), 
                      foldid = foldid,
                      alpha=1)
  
  # Extract the coefficients produced by the solution that minimizes lambda
  lassocoef.min = coef(cvlasso, s="lambda.min")
  
  ## Create a dataframe of the results (75 columns, 1 row)
  lassoCoefs <- as.data.frame(as.matrix(lassocoef.min))[-1,]
  temp <- as.data.frame(lassoCoefs)
  lassoCoefsTranspose <- t(temp)
  colnames(lassoCoefsTranspose) <- colnames(subset[[dataset]][,predCols])
  
  ## Save the proportion of runs that produced non-zero values in the lASSO. This should be a single value (1 or 0) since the LASSO is based on just a single run
  number.nonzero.lasso <- colSums(lassoCoefsTranspose != 0)
  proportion.nonzero.lasso <- number.nonzero.lasso/nrow(lassoCoefsTranspose)
  
  
  #### LASSO 1se ####
  
  # Running the lasso is identical as above
  
  # Extract the coefficients produced by 1se lambda
  lassocoef.1se = coef(cvlasso, s="lambda.1se")
  
  ## Create a dataframe of the results (75 columns, 1 row)
  lassoCoefs.1se <- as.data.frame(as.matrix(lassocoef.1se))[-1,]
  temp <- as.data.frame(lassoCoefs.1se)
  lassoCoefs.1seTranspose <- t(temp)
  colnames(lassoCoefs.1seTranspose) <- colnames(subset[[dataset]][,predCols])
  
  ## Save the proportion of runs that produced non-zero values in the lASSO. This should be a single value (1 or 0) since the LASSO is based on just a single run
  number.nonzero.lasso.1se <- colSums(lassoCoefs.1seTranspose != 0)
  proportion.nonzero.lasso.1se <- number.nonzero.lasso.1se/nrow(lassoCoefs.1seTranspose)
  
  #### Correlation with p < 0.1 ####
  corr_results <- matrix(nrow=50,ncol=4)
  colnames(corr_results) <- c('Pred','Coefficient','P-value','Selected')
  
  x = subset[[dataset]][,predCols]
  y = subset[[dataset]][,"y"]
  for (col in 1:50) {
    cor_ <- cor.test(x[,col],y)
    corr_results[col,1] <- paste0('X',col)
    corr_results[col,2] <- cor_$estimate
    corr_results[col,3] <- cor_$p.value
    corr_results[col,4] <- as.integer(cor_$p.value < 0.1)
  }
  corr_results <- as.data.frame(corr_results)
  rownames(corr_results) <- corr_results$Pred
  
  #### Correlation with p < 0.05 ####
  corr_results_05 <- matrix(nrow=50,ncol=4)
  colnames(corr_results_05) <- c('Pred','Coefficient','P-value','Selected')
  
  x = subset[[dataset]][,predCols]
  y = subset[[dataset]][,"y"]
  for (col in 1:50) {
    cor_ <- cor.test(x[,col],y)
    corr_results_05[col,1] <- paste0('X',col)
    corr_results_05[col,2] <- cor_$estimate
    corr_results_05[col,3] <- cor_$p.value
    corr_results_05[col,4] <- as.integer(cor_$p.value < 0.05)
  }
  corr_results_05 <- as.data.frame(corr_results_05)
  rownames(corr_results_05) <- corr_results_05$Pred
  
  
  # Transpose the SSVS values and change rownames
  inc_tran = as.data.frame(t(inc_prob[,1]))
  names(inc_tran) = rownames(inc_prob)
  beta_tran <- as.data.frame(t(average.beta))
  beta_tran <- sapply(beta_tran,as.numeric)
  lower_tran <- as.data.frame(t(lower.credibility))
  upper_tran <- as.data.frame(t(upper.credibility))
  average_tran <- as.data.frame(t(average.nonzero.beta))
  median_tran <- as.data.frame(t(median.beta))
  prop_tran <- as.data.frame(t(proportion.nonzero))
  
  #### Save values to a data frame ####
  loopValues <- rbind(
    
    ## Adaptive LASSO values
    
    as.data.frame(t(adaResults))[2,],
    as.data.frame(t(adaResults))[3,],
    
    
    ## SSVS values
    as.data.frame(inc_tran),
    beta_tran,
    as.data.frame(lower_tran),
    as.data.frame(upper_tran),
    as.data.frame(average_tran),
    as.data.frame(median_tran),
    as.data.frame(prop_tran),
    
    ## LASSO min values
    as.data.frame(t(proportion.nonzero.lasso)),
    as.data.frame(lassoCoefsTranspose),
    
    ## LASSO 1se values
    as.data.frame(t(proportion.nonzero.lasso.1se)),
    as.data.frame(lassoCoefs.1seTranspose),
  
    ## Correlation p < 0.1 values
    as.data.frame(t(corr_results))[2,],
    as.data.frame(t(corr_results))[3,],
    as.data.frame(t(corr_results))[4,],
    
    ## Correlation p < 0.05 values
    as.data.frame(t(corr_results_05))[2,],
    as.data.frame(t(corr_results_05))[3,],
    as.data.frame(t(corr_results_05))[4,]
    
  )
  
  # Denote which measures were taken
  loopValues$measures <- c(
    # Adaptive Lasso values
    "Adaptive Lasso selected",
    "Adaptive Lasso coefficient",
    
    
    # SSVS values
    "SSVS MIP",
    "SSVS average Beta",
    "SSVS average Beta lower credibility interval (HDI)",
    "SSVS average Beta upper credibility interval (HDI)",
    "SSVS average nonzero Beta",
    "SSVS median Beta",
    "SSVS proportion nonzero variables",
    
    # LASSO min values
    "LASSO min selected",
    "LASSO min coefficient",
    
    # LASSO 1se values
    "LASSO 1se selected",
    "LASSO 1se coefficient",

    # Correlation p < 0.1
    "Correlation 0.1 coefficient",
    "Correlation 0.1 p-value",
    "Correlation 0.1 selected",
    
    # Correlation p < 0.05
    "Correlation 0.05 coefficient",
    "Correlation 0.05 p-value",
    "Correlation 0.05 selected"
    
  )
  
  ### Model-level values
  
  # Save the number of predictors >0.5 that were selected for SSVS
  loopValues$numSvsPreds <- c(sum(inc_prob[,1]>.5), 
                              rep(NA,(length(loopValues$measures)-1)))
  
  # Denote the rep number
  loopValues$repNumber <- c(dataset, 
                            rep(NA,(length(loopValues$measures)-1)))
  
  # Save the N value for the dataset
  loopValues$N <- c(unique(subset[[dataset]][,"N"]), 
                    rep(NA,(length(loopValues$measures)-1)))
  
  # Save the r value for the dataset
  loopValues$r <- c(unique(subset[[dataset]][,"r"]), 
                    rep(NA,(length(loopValues$measures)-1)))
  
  ## Add the loop to the established dataframe
  finalValues <- rbind(finalValues, loopValues)
  
}

#### After the loop is finished ####

# Name the columns 
colnames(finalValues) <- c(colnames(subset[[dataset]][,predCols]), 
                           "Measures",
                           "SSVSforPsych model size", 
                           "Dataset",
                           "N",
                           "r")

# Reorder the columns
finalValues2 <- finalValues[c(53,54,55,52,51,1:50)]


