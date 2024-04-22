## Load the required libraries 
if(1){
  library(boot)
  library(caret)
  library(glmnet)
  library(BayesLogit)
  library(invgamma)
  library(MASS)
  library(statmod)
  library(truncnorm)
  library(coda)
  library(broom)
  library(tibble)
  library(ggimage)
  library(rsvg)
  library(pROC)
}

## Read in the data 
## The data is arranged with the following columns
## Col. 1 : Event Longtitude
## Col. 2 : Event Latitude
## Col. 3 : Station Longitude
## Col. 4 : Station Latitude
## Col. 5 : Recorded Efficiency (0: blocked, 1: inefficient, 2: efficient) 

data1=read.table("lg_test.txt")
colnames(data1)[1]="S_Long"
colnames(data1)[2]="S_Lat"
colnames(data1)[3]="R_Long"
colnames(data1)[4]="R_Lat"
colnames(data1)[5]="efficiency"
one=which(data1$efficiency==1)
zero=which(data1$efficiency==0)
two=which(data1$efficiency==2)

## Check the distribution of the recorded efficiency
## table(data1$efficiency)

## Filter the data based on the application
## Case 1: inefficient paths are removed 
## Case 2: inefficient paths are treated as blocked
## Case 3: inefficient paths are treated as efficient

data2=data1[c(two,zero),]
table(data2$efficiency)

#data2=data1[c(one,zero,two),]
#data2$efficiency[data2$efficiency==1]=0
#table(data2$efficiency)

#data2=data1[c(one,zero,two),]
#data2$efficiency[data2$efficiency==1]=2
#table(data2$efficiency)

## Split into training and validation sets ###
## Training set is obtained by randomly sampling 75% of the data 
## Remaining 25% is used as the validation set to determine the decision rule. 

train=sample(nrow(data2),nrow(data2)*0.75)
train.data=data2[train,]
val.data=data2[-train,]

## Define the model components 

dist=train.data[,c(6:ncol(train.data))] 
X=dist*pi*1/3.5
X=as.matrix(X)

Y=as.factor(train.data$efficiency)
y=model.matrix(~Y)
y=y[,-1]
Y=as.matrix(y)

## Use cross-validation to determine the tuning parameter #################
## Set alpha=1 in cv.glmnet() to choose Lasso.

mod=train.data[,-c(1:4)]
xo=model.matrix(efficiency~.,data=mod)
xo=xo[,-1]
yo=train.data$efficiency
cv=cv.glmnet(xo,yo,alpha=1)
lambda=cv$lambda.min

###################### Define the sampler ##############################

updateSig2 <- function(a,b,Y,X,Beta,Tau){
  N <- length(Y)
  p <- length(Beta)
  D <- solve(diag(Tau))
  a.post <- (N - 1 + p)*0.5 + a
  b.post <- 0.5*(t(Y - X%*%Beta)%*%(Y - X%*%Beta) + t(Beta)%*%D%*%Beta) + b
  return(rinvgamma(1,a.post,rate=b.post))
} 
updateTau <- function(Lambda,Sig2,Beta){
  dispersion <- Lambda^2
  mu <- sqrt(Lambda^2 * Sig2 / Beta^2)
  return(1/rinvgauss(length(mu),mean = mu, dispersion = dispersion))
}

## A prior distribution can also be used to estimate the tuning parameter, lambda
#updateLambda <- function(p,r,Delta,Tau){
#  a.post <- p + r
#  b.post <- sum(Tau^2)*0.5 + Delta
#  return(sqrt(rgamma(1,a.post,b.post)))
#}

logit.R <- function(y, X, samp=4000, burn=1000, verbose=1000)
{
  X = as.matrix(X);
  y = as.numeric(y)
  
  p = ncol(X)
  N = nrow(X)
  n=rep(1, length(y))
  m0=rep(0, ncol(X))
  P0=matrix(0, nrow=ncol(X), ncol=ncol(X))
  
  kappa = (y-1/2)*n
  
  Z = colSums(X*kappa) + P0 %*% m0;
  
  w = rep(0,N)
  beta = rep(0, p)
  
  output <- list(w = matrix(nrow=samp, ncol=N), beta = matrix(nrow=samp, ncol=p))
  
  for ( j in 1:(samp+burn) )
  {
    
    ## draw w
    psi = drop(X%*%beta)
    w = rpg.devroye(N, n, psi);
    
    tau=updateTau(lambda,1,beta)
    D_tau=diag(tau)
    
    ## draw beta - Joint Sample.
    PP = t(X) %*% (X * w) + solve(D_tau);
    S = chol2inv(chol(PP));
    m = S %*% (t(X)%*%alpha+solve(D_tau)%*%m0)
    beta = m + t(chol(S)) %*% rtruncnorm(p,a=0);
    
    # Save the post burnin samples
    if (j>burn) {
      output$w[j-burn,] <- w
      output$beta[j-burn,] <- beta
    }
    
    if (j %% verbose == 0) { print(paste("LogitPG: Iteration", j)); }
  }
  
  ## Produce the output 
  output$"y" = y;
  output$"X" = X;
  output$"n" = n;
  
  output
}


####################### Run the sampler ########################

fit=logit.R(y,X, samp=4000, burn=1000)

betas=fit$beta
betas=as.matrix(betas)

####################### Determine the decision rule using the validation set ###########################
## Define a grid between 0 and 1
dec_grid=seq(0,1,by=0.025)

pred=matrix(0,nrow(val.data),length(dec_grid))
AUC=rep(NA,length(dec_grid))
Acc=rep(NA,length(dec_grid))

dist4=val.data[,c(6:ncol(val.data))]
X4=as.matrix(dist4)
X4=X4*pi*1/3.5

g4=X4%*%t(betas)
h4=inv.logit(g4)

prob=rowMeans(h4)

for(i in 1:length(dec_grid))
{
  
  for(j in 1:nrow(X4))
  {
    if(prob[j]==0){prob[j]=prob[j]+1e-6}
    else if(prob[j]==1) {prob[j]=prob[j]-1e-6}
    if(prob[j]>dec_grid[i])
    {pred[j,i]=2}
    else {pred[j,i]=0}
  }
  if(i==1 | i==length(dec_grid))
  {
    AUC[i]=0
    Acc[i]=0
  }
  else
  {
    AUC[i]=auc(roc(pred[,i],val.data$efficiency))
    Acc[i]=confusionMatrix(as.factor(pred[,i]),as.factor(val.data$efficiency))$overall[1]
  }
  
}

results=data.frame(dec_grid,AUC,Acc)
results=na.omit(results)

max(results$Acc)

dr=dec_grid[which(Acc==max(results$Acc))]  ## classification rule corresponding to the highest accuracy

##################### Predictions ############################
## Produce the final classification table and compute the AUC under the ROC

dist3=data2[,c(6:ncol(data2))]
X3=as.matrix(dist3)
X3=X3*pi*1/3.5

g3=X3%*%t(betas)
h3=inv.logit(g3)
y_post3=rowMeans(h3)

A3=rep(0,nrow(X3))
for(i in 1:nrow(X3))
{
  if(y_post3[i]>dr)
  {A3[i]=2}
  else {A3[i]=0}
}

table(A3)
confusionMatrix(as.factor(A3),as.factor(data2$efficiency))
auc(roc(A3,data2$efficiency))

