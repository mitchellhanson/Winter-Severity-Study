rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(splines)
library(xtable)
library(stargazer)
allweatherdata <- read.csv("/Users/xuanruizhang/Desktop/voltage/Thesis/Winter-And-Pop.csv")
myset <- setdiff(1:nrow(allweatherdata), (14224:14288))
allweatherdata <- allweatherdata[myset,]
allweatherdata$Location = paste0(allweatherdata$State, allweatherdata$City)
allweatherdata$maxwinter <- rep(0,nrow(allweatherdata))
allweatherdata$AdjAWSSI <- rep(0,nrow(allweatherdata))
allweatherdata$SmallAWSSI <- allweatherdata$AWSSI/1000
for (i in 1:nrow(allweatherdata)){
  temp <- allweatherdata$Location[i]
  allweatherdata$maxwinter[i] <-  max(allweatherdata[allweatherdata$Location == temp,]$AWSSI)
  allweatherdata$AdjAWSSI[i] <- allweatherdata$AWSSI[i]/allweatherdata$maxwinter[i]
}

#Average AdjAWSSI
avgadja <- (aggregate(allweatherdata$AdjAWSSI, list(allweatherdata$Season), FUN=mean))[1:69,]
avgg <- (aggregate(allweatherdata$AWSSI, list(allweatherdata$Season), FUN=mean))[1:69,]
avlg<- avgg <- (aggregate(log(allweatherdata$AWSSI + 1), list(allweatherdata$Season), FUN=mean))[1:69,]


#Population weighted AdjAWSSI
totalpop <- aggregate(allweatherdata$X2010Pop,list(allweatherdata$Season), FUN=sum )[1:69,]
paww <-aggregate((allweatherdata$X2010Pop) * allweatherdata$SmallAWSSI  ,list(allweatherdata$Season), FUN=sum )[1:69,]
paww2 <-aggregate((allweatherdata$X2010Pop) * allweatherdata$AdjAWSSI  ,list(allweatherdata$Season), FUN=sum )[1:69,]
paww $adjx <- paww$x / totalpop$x
paww2 $adjx <- paww2$x / totalpop$x

##Autocorrelation Utilities
#J-th Order Covariance
calc_autocov <- function(j,X){
  m <- mean(X)
  count <- 0
  for (i in (j+1):length(X)){
    count <- count + (X[i]-m)*(X[i-j]-m)
  }
  return(count/length(X))
}

#Ljung-Box
LBS <- function(p,X,rou){
  count <- 0
  n <- length(X)
  for (j in 1:p){
    count <- count + rou[j]*rou[j]/(n-j)
  }
  return(count*n*(n+2))
}

#Ljung-Box P-value
calc_lbp <- function(mydata,rou,lth){
  rou_l <- rou[-1]
  LBSS<- rep(0,lth)
  LBP <- rep(0,lth)
  for (t in 1:lth){
    LBSS[t] <- LBS(t,mydata,rou_l)
    LBP[t] <- 1-pchisq(LBSS[t],t)
  }
  return(LBP)
}

##Import GDP Data
dfgdp <- read.csv("/Users/xuanruizhang/Desktop/voltage/Thesis/gdp.csv")[1:289,]
#GDP Used (2000-2019)
newgdp = dfgdp$GDP [212:289]
covars_gdp = rep(0,21)
for (i in 1:21){
  covars_gdp[i] = calc_autocov(i-1,newgdp)
}
corr_gdp <- covars_gdp/covars_gdp[1]
LB <- rep(0,20)
for (i in 1:20){
  LB[i] = LBS(i,newgdp,corr_gdp)
}
LBP <- calc_lbp(newgdp, corr_gdp, length(corr_gdp) - 1)
mydf <- cbind.data.frame(0:20, corr_gdp) 
colnames(mydf) <- c('Index', 'Corr')
#GDP Corrrelogram
plot1 <- ggplot() + geom_line(aes(x = mydf$Index , y=mydf$Corr)) +
  labs(x = "Period", y = "Autocorrelation Coefficients", title = "Correlogram With 20 Lags") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=20,  family="CMU Serif"))  

ggsave(filename = "plot1.png", plot = plot1, width = 8, height = 5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

#output table
cto <- cbind.data.frame(corr_gdp[2:21], LB,LBP)
xtable(cto)

##Setup of OLS
q2gdp <- rep(0,69)
q3gdp <- rep(0,69)
q4gdp <- rep(0,69)
q1gdp <- rep(0,69) 

for (i in seq(16,288,4)){
  q1gdp[(i/4)-3] <-  dfgdp$GDP[i] 
  q4gdp[(i/4)-3] <-  dfgdp$GDP[i-1] 
  q3gdp[(i/4)-3] <-  dfgdp$GDP[i-2] 
  q2gdp[(i/4)-3] <-  dfgdp$GDP[i-3] 
}
#Run OLS
#Adj
OLS1 <- lm(q1gdp[51:69]~avgadja $x[51:69] + q4gdp[51:69]  + q3gdp[51:69])
#
OLS2 <- lm(q1gdp[51:69]~avgg $x[51:69] + q4gdp[51:69]  + q3gdp[51:69])
#Pa
OLS3 <- lm(q1gdp[51:69]~paww $adjx[51:69] + q4gdp[51:69]  + q3gdp[51:69])
#Adj,Pa
OLS4 <- lm(q1gdp[51:69]~paww2 $adjx[51:69] + q4gdp[51:69]  + q3gdp[51:69])
stargazer(OLS1,OLS2, title="Results", align=TRUE)
stargazer(OLS3,OLS4, title="Results", align=TRUE)

##Local Linear Estimator with gaussian kernel
#Gaussian Kernel

#Gaussian Kernel
my_gaussian_kernel <- function(sigma,mu,x){
  res <- exp(-0.5 * ((x-mu)/sigma)^2)
  return ((1/sqrt(2*pi)) * res)
}

#Local Linear Estimator
myll <- function(bw, X, Y, L){
  estimator <- matrix(0,2,length(L))
  for (i in 1:length(L)){
    alldta <- cbind(X,Y)
    sel <- c()
    for(j in 1:length(X)){
 #     if ((abs(L[i] - X[j]) < bw) && (abs(L[i] - X[j])>0) ){
        sel <- rbind(sel,alldta[j,])
#      }
    }
    weight <- rep(0,nrow(sel))
    for (j in 1:nrow(sel)){
      weight[j] = my_gaussian_kernel(1,0,(sel[j,1] - L[i])/bw)
    }
    W <- diag(weight)
    Z <- cbind(rep(1,nrow(sel)),(sel[,1] - rep(L[i], nrow(sel))))
    ll_coefs <- solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W %*% sel[,2] 
    estimator[,i] <- ll_coefs
  }
  return (estimator)
}

#llcv,leave one out CV, returns optimal bandwidth
myll_cv <- function(X,Y,d,s,e){
  # d is density
  min_mse <- 2147483647
  optimal_bw <- 2147483647
  for (i in seq(s,e, by = d)){
    r <- myll(i,X,Y,X)
    residues <- r[1,] - Y
    mse <- t(residues) %*% residues
    if (mse < min_mse){
      min_mse = mse
      optimal_bw = i
    }
  }
  return(optimal_bw)
}

myll_gauss_sigma <-function(X,Y,bw,L){
  res <- rep(0, length(L))
  #esti <- myll(bw,X,Y,X)[1,] - Y
  #e2 <- esti ^ 2
  Rk <- 1/(2*sqrt(pi))
  for (i in 1:length(L)){
    E <- c(L[i])
    #esti <- cbind(rep(1,length(X)),X) %*% myll(bw,X,Y,E) - Y
    esti <- myll(bw,X,Y,X)[1,] - Y
    e2 <- esti^2
    upper <- 0
    lower <- 0
    for (j in 1:length(X)){
      upper <- upper + my_gaussian_kernel(1,0,(X[j] - L[i])/bw) * e2[j]
      lower <- lower + my_gaussian_kernel(1,0,(X[j] - L[i])/bw)
    }
    mid <-  Rk * (upper/lower)
    res[i] <- mid/lower
  }
  return(sqrt(res))
}

Y <- q1gdp[51:69]
X <- avgadja$x[51:69]
G <- seq(0.3,0.65, by=0.01)
opt_bandwidth  <- myll_cv(X,Y,0.01,0.06,0.9)
LL <- myll(opt_bandwidth,X,Y,G)[1,] 
ll_95cl <- 1.96 * myll_gauss_sigma(X,Y,opt_bandwidth,G)
ll95_low  <- LL - ll_95cl
ll95_high <- LL + ll_95cl

alldata = cbind.data.frame(X,Y)
output = cbind.data.frame(G, ll95_low, ll95_high,LL)
plot2 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output,mapping=aes(x=G,ymin=ll95_low,ymax=ll95_high),fill = 'rosybrown',alpha=0.5)+
  geom_line(data=output, mapping=aes(x=G,y=LL),size=1.5)+
  labs(x = "Average Adjusted AWSSI", y = "Q1 GDP Growth", title = "Local Linear Regression") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot2.png", plot = plot2, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

X <- avgg$x[51:69]
G <- seq(5.3,6.2 ,by=0.01)
opt_bandwidth  <- myll_cv(X,Y,0.01,0.1,2)
LL <- myll(opt_bandwidth,X,Y,G)[1,]
ll_95cl <- 1.96 * myll_gauss_sigma(X,Y,opt_bandwidth,G)
ll95_low  <- LL - ll_95cl
ll95_high <- LL + ll_95cl
alldata = cbind.data.frame(X,Y)
output = cbind.data.frame(G, ll95_low, ll95_high,LL)
plot3 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output,mapping=aes(x=G,ymin=ll95_low,ymax=ll95_high),fill = 'rosybrown',alpha=0.5)+
  geom_line(data=output, mapping=aes(x=G,y=LL),size=1.5)+
  labs(x = "Average AWSSI", y = "Q1 GDP Growth", title = "Local Linear Regression") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot3.png", plot = plot3, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

X <- paww$adjx[51:69]
G <- seq(0.15,0.5 ,by=0.01)
opt_bandwidth  <- myll_cv(X,Y,0.01,0.06,0.9)
LL <- myll(opt_bandwidth,X,Y,G)[1,]
ll_95cl <- 1.96 * myll_gauss_sigma(X,Y,opt_bandwidth,G)
ll95_low  <- LL - ll_95cl
ll95_high <- LL + ll_95cl
alldata = cbind.data.frame(X,Y)
output = cbind.data.frame(G, ll95_low, ll95_high,LL)
plot4 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output,mapping=aes(x=G,ymin=ll95_low,ymax=ll95_high),fill = 'rosybrown',alpha=0.5)+
  geom_line(data=output, mapping=aes(x=G,y=LL),size=1.5)+
  labs(x = "Average AWSSI (Pop Adjusted)", y = "Q1 GDP Growth", title = "Local Linear Regression") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot4.png", plot = plot4, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")


X <- paww2$adjx[51:69]
G <- seq(0.15,0.7 ,by=0.01)
opt_bandwidth  <- myll_cv(X,Y,0.01,0.06,0.9)
LL <- myll(opt_bandwidth,X,Y,G)[1,]
ll_95cl <- 1.96 * myll_gauss_sigma(X,Y,opt_bandwidth,G)
ll95_low  <- LL - ll_95cl
ll95_high <- LL + ll_95cl
alldata = cbind.data.frame(X,Y)
output = cbind.data.frame(G, ll95_low, ll95_high,LL)
plot5 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output,mapping=aes(x=G,ymin=ll95_low,ymax=ll95_high),fill = 'rosybrown',alpha=0.5)+
  geom_line(data=output, mapping=aes(x=G,y=LL),size=1.5)+
  labs(x = "Average AWSSI (Pop Adjusted)", y = "Q1 GDP Growth", title = "Local Linear Regression") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot5.png", plot = plot5, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")


##B - Spline Estimation
my_spline<- function(p, X, Y){
  quant_list <- seq(1/(p+1),1,length.out = p+1)[-(p+1)]
  K <-  unname(quantile(X,quant_list))
  bsp <- lm(Y~bs(X,knots=K))
  return (bsp)
}

my_bscv <- function(X, Y){
  pgrid <- seq(1,20, length.out=20)
  least_mse <- 2147483647
  p_optimal <- 1
  for (j in 1:length  (pgrid)){
    p_val <-  pgrid[j]
    mse <- 0
    for (i in 1:length(X)){
      x_data <- X[i]
      y_data <- Y[i]
      x_set <- X[-i]
      y_set <- Y[-i]
      bsp_model <- my_spline(p_val, x_set, y_set)
      loo_bspline <- predict(bsp_model,data.frame(X=c(x_data)))
      mse <- mse + (y_data-loo_bspline)^2
    }
    if (mse <  least_mse){
      least_mse <- mse
      p_optimal <- p_val
    }
  }
  return   (p_optimal)
}

#BS Estimator
X <- avgadja$x[51:69]
G <- seq(0.3,0.65, by=0.01)
optimal_bs <- my_bscv(X,Y)
BS <- unname(predict(my_spline(optimal_bs,X,Y),data.frame(X=G),interval="confidence"))
alldata = cbind.data.frame(X,Y)
output2 <- cbind.data.frame(G,BS)
plot6 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output2,mapping=aes(x=G,ymin=BS[,2],ymax=BS[,3]),fill = 'olivedrab',alpha=0.5)+
  geom_line(data=output2, mapping=aes(x=G,y=BS[,1],color="BS Estimator"),size=1.5)+
  labs(x = "Average Adjusted AWSSI", y = "v", color = "Method", title = "B-Spline Estimation") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot6.png", plot = plot6, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

X <- avgg$x[51:69]
G <- seq(5.3,6.2 ,by=0.01)
optimal_bs <- my_bscv(X,Y)
BS <- unname(predict(my_spline(optimal_bs,X,Y),data.frame(X=G),interval="confidence"))
alldata = cbind.data.frame(X,Y)
output2 <- cbind.data.frame(G,BS)
plot7 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output2,mapping=aes(x=G,ymin=BS[,2],ymax=BS[,3]),fill = 'olivedrab',alpha=0.5)+
  geom_line(data=output2, mapping=aes(x=G,y=BS[,1],color="BS Estimator"),size=1.5)+
  labs(x = "Average AWSSI", y = "Q1 GDP Growth", color = "Method", title = "B-Spline Estimation") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot7.png", plot = plot7, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

X <- paww$adjx[51:69]
G <- seq(0.15,0.5 ,by=0.01)
optimal_bs <- my_bscv(X,Y)
BS <- unname(predict(my_spline(optimal_bs,X,Y),data.frame(X=G),interval="confidence"))
alldata = cbind.data.frame(X,Y)
output2 <- cbind.data.frame(G,BS)
plot8 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output2,mapping=aes(x=G,ymin=BS[,2],ymax=BS[,3]),fill = 'olivedrab',alpha=0.5)+
  geom_line(data=output2, mapping=aes(x=G,y=BS[,1],color="BS Estimator"),size=1.5)+
  labs(x = "Average AWSSI (Pop. Adjusted)", y = "Q1 GDP Growth", color = "Method", title = "B-Spline Estimation") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot8.png", plot = plot8, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

X <- paww2$adjx[51:69]
G <- seq(0.15,0.68 ,by=0.01)
optimal_bs <- my_bscv(X,Y)
BS <- unname(predict(my_spline(optimal_bs,X,Y),data.frame(X=G),interval="confidence"))
alldata = cbind.data.frame(X,Y)
output2 <- cbind.data.frame(G,BS)
plot9 <- ggplot() + geom_point(data=alldata, mapping=aes(x=X,y=Y))+
  geom_ribbon(data=output2,mapping=aes(x=G,ymin=BS[,2],ymax=BS[,3]),fill = 'olivedrab',alpha=0.5)+
  geom_line(data=output2, mapping=aes(x=G,y=BS[,1],color="BS Estimator"),size=1.5)+
  labs(x = "Average Adjusted AWSSI (Pop. Adjusted)", y = "Q1 GDP Growth", color = "Method", title = "B-Spline Estimation") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=22,  family="CMU Serif")) 
ggsave(filename = "plot9.png", plot = plot9, width = 12, height = 7.5, path = "/Users/xuanruizhang/Desktop/voltage/Thesis")

