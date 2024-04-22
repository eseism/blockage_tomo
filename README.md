# Bayesian Lasso-Logit model

This Rscript can be used to produce the to reproduce the main results in the journal article ''Predicting Lg Blockage in the Middle East using a Bayesian Lasso Logistic Regression Model".

The model is a Bayesian Logistic regression model with LASSO penalty and uses polya-gamma data augmentation strategies to handle non-conjugacy. 


## Authors

- Saikat Nandy [@SaikatNandy](https://github.com/SaikatNandy)


## Installation


```bash
#install.packages("devtools",dependencies=TRUE) #install devtools from CRAN
#library(devtools) #load devtools
libraries=c("boot","caret","BayesLogit","invgamma","MASS","statmod",
                   "coda","broom","tibble","rsvg","pROC","ggimage")
install.packages(libraries,dependencies=TRUE)
sapply(libraries,require,character=TRUE)

```

    
## Data
The data set contains the Lg efficiency data to reproduce the main results in the journal article ''Predicting Lg Blockage in the Middle East using a Bayesian Lasso Logistic Regression Model".

Column information of the data: \
1 Event Longitude \
2 Event Latitude\
3 Station Longitude\
4 Station Latitude\
5 Recorded Efficiency (0: blocked, 1: inefficient, 2: efficient) 

```bash
data1=read.table("lg_test.txt")
colnames(data1)[1]="S_Long"
colnames(data1)[2]="S_Lat"
colnames(data1)[3]="R_Long"
colnames(data1)[4]="R_Lat"
colnames(data1)[5]="efficiency"
one=which(data1$efficiency==1)
zero=which(data1$efficiency==0)
two=which(data1$efficiency==2)
```
## Usage/Examples

Filter the data depending on how the inefficient rays are handled \
Case 1: inefficient paths are removed \
Case 2: inefficient paths are treated as blocked\
Case 3: inefficient paths are treated as efficient


```javascript
data2=data1[c(two,zero),]
table(data2$efficiency)

#data2=data1[c(one,zero,two),]
#data2$efficiency[data2$efficiency==1]=0
#table(data2$efficiency)

#data2=data1[c(one,zero,two),]
#data2$efficiency[data2$efficiency==1]=2
#table(data2$efficiency)
```

Split the data into training and validation sets. Training set is obtained by randomly sampling 75% of the data and the remaining 25% is used as the validation set to determine the decision rule. 

```javascript
train=sample(nrow(data2),nrow(data2)*0.75)
train.data=data2[train,]
val.data=data2[-train,]
```

Prepare the response dummy variables and the subdistance matrix.
```javascript
dist=train.data[,c(6:ncol(train.data))] 
X=dist*pi*1/3.5
X=as.matrix(X)

Y=as.factor(train.data$efficiency)
y=model.matrix(~Y)
y=y[,-1]
Y=as.matrix(y)
```

Use cross-validation to determine the tuning parameter . Set ```alpha=1 ``` to choose LASSO penalty.

```javascript
mod=train.data[,-c(1:4)]
xo=model.matrix(efficiency~.,data=mod)
xo=xo[,-1]
yo=train.data$efficiency
cv=cv.glmnet(xo,yo,alpha=1)
lambda=cv$lambda.min
```

Run the MCMC algorithm using the ```logit.R()``` function defined in the Rscript. ```5000``` samples are drawn and the first ```1000``` are discarded as burnin.  

```javascript
fit=logit.R( y, X, samp=4000, burn=1000)

betas=fit$beta
betas=as.matrix(betas)
```

Determine the decision rule using the validation set. Further details about this can be found in Hui et al. (2023).

```javascript
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

```

Produce the final predictions, confusion matrix, and compute the AUROC

```javascript
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

# y_post3 is the vector of estimated probability of blockage.
# A3 is vector of classified ray efficiencies.

table(A3)
confusionMatrix(as.factor(A3),as.factor(data2$efficiency))
auc(roc(A3,data2$efficiency))
```


## References


Hongjun Hui, Saikat Nandy, Scott H. Holan, Jingjing Pan, Duyi Li, Eric A. Sandvol; A Bayesian Lasso Logistic Regression Model for Predicting the Probability of Regional Seismic Phase Observation Using Sn in the Middle East and East Asia as Examples. Bulletin of the Seismological Society of America 2023; 113 (2): 562–576. 

Polson, Nicholas G., James G. Scott, and Jesse Windle. “Bayesian Inference for Logistic Models Using Pólya–Gamma Latent Variables.” Journal of the American Statistical Association 108, no. 504 (2013): 1339–49.
