rm(list = ls())
library(MatchIt)
library(haven)
library(dplyr)
library(lmtest)

## Import data
data <- read_dta("data.dta")

# outcome variable: CEO compensation (tot_compen)
# unit variable: cik 
# time variable :fyear
# treatment group:Stayers (Stayers = 1)
# control group: Switchers (Stayers = 0)
# covariates: stock price (stockprice), net income (net_income), total assets (total_assets)

## Generate dummy variables for treatment, post-intervention, and the interaction between the two 
data$trt <- data$Stayers

data$post <- 0 
data$post[which(data$fyear>=2009)] <- 1
data <- data[complete.cases(data), ]  ## Remove NA

data$did <- data$post*data$Stayers ## treatment * Post

## outcome variable and covariates 
data$x1 <- as.numeric((data$stockprice))
data$x2 <- as.numeric((data$net_income))
data$x3 <- as.numeric((data$total_assets))
data$y <- data$tot_compen

data <- subset(data, select=c("cik","fyear", "trt","y", "x1", "x2", "x3", "post", "did"))

## Generate dummy variables for each post-intervention period
data$td0 <- 0 ; data$td0[data$fyear==2009]<-1;  
data$td1 <- 0 ; data$td1[data$fyear==2010]<-1;
data$td2 <- 0 ; data$td1[data$fyear==2011]<-1;

###### Conventional DiD 
# DiD estimation (model 1)
m1.qmle <- glm(y ~ post + trt + did + x1 + x2 + x3, data = data, family = quasipoisson(link="log"))

# Robust standard error
coeftest(m1.qmle, vcov = vcovHC(m1.qmle, type="HC1"))

# Mean Squared Error 
m1.qmle.mse <- mean((m1.qmle[["residuals"]])^2)

###### Conventional hybrid approach
## Stage 1: Matching on covariates
## using the year 2008 for illustrative purposes
data.pre <- subset(data, data$fyear==2008, select = c("cik", "x1","x2","x3", "trt"))
m.out1 <- matchit(trt ~ 0 + x1 + x2 + x3, data = data.pre, caliper = 0.2)
matched_df1 <- match.data(m.out1)
id.m1 <- unique(matched_df1$cik)

# Matched Sample
dat.match1 <- subset(data, cik %in% id.m1)

## Stage 2: DiD estimation (model 2)
m2.qmle <- glm(y ~ post + trt + did + x1 + x2 + x3 , data = dat.match1, family = quasipoisson(link="log"))

# Robust standard error
coeftest(m2.qmle, vcov = vcovHC(m2.qmle, type="HC1"))

# Mean Squared Error 
m2.qmle.mse <- mean((m2.qmle[["residuals"]])^2)

###### GM hybrid approach
# Step for coefficients
coefX <- function(yp,Xp,iXpXp) {
  betaX <- iXpXp %*% crossprod(Xp,yp)
  betaX
}
# Step for unit fixed-effects  
coefCit <- function(ym,Xm, yi,Xi,yt,Xt,betaX) {
  coefi <- yi - Xi%*%betaX
  coeff <- list(coefi)
  names(coeff) <- c("coefi")
  coeff
}
# Estimation of parameters
coefsbai <- function(y,X, N, T) {
  
  N <- N
  T <- T
  
  numi <- seq(1:N) %x% rep(1,T)
  numt <- rep(1,N) %x% seq(1:T)
  ym <- apply(y,2,mean)
  yi <- as.matrix(tapply(y,numi,mean),N,1)
  yp <- y-as.matrix(yi %x% rep(1,T))
  
  Xm <- t(apply(X,2,mean))
  Xi <- matrix(0,N,ncol(X))
  for (k in 1:ncol(X)) 
  { Xi[,k] <- tapply(X[,k],numi,mean) }
  Xp   <- X-Xi%x%rep(1,T)
  iXpXp <- solve(crossprod(Xp))
  
  coefX <- coefX(yp,Xp,iXpXp)
  cf<- coefCit(ym,Xm, yi,Xi,yt,Xt,coefX)
  coeff <- list(coefX,cf$coefi)
  names(coeff) <- c("coefX","coefi")
  coeff
}

## Stage 1: Extract unit fixed effects and matching on covariates, outcomes, and the estimated unit fixed effects 
# Using Pre-intervention data (should be balanced)
data.pre <- subset(data, fyear < 2009)

# Make a balanced panel 
balanced.id <- data.pre %>% dplyr::count(cik)
balanced.id <- subset(balanced.id, n==2)
balanced.id <- unique(balanced.id$cik)
data.bpre <- subset(data.pre, cik %in% balanced.id) # Keep units that were observed during all pre-intervention periods

# Number of units and number of pre-intervention periods
N <- length(unique(data.bpre$cik))
T0 <- 2

# Estimate the CEO fixed effects
beta <- coefsbai(as.matrix(data.bpre$y),as.matrix(cbind(data.bpre$x1, data.bpre$x2, data.bpre$x3)),N,T0)
res.ai <- data.frame(unique(data.bpre$cik), beta$coefi)
names(res.ai) <- c("cik", "mu")

data.bpre <- data.bpre %>% 
  left_join(res.ai, by = "cik") 
data <- data %>% 
  left_join(res.ai, by = "cik") 

## Matching on covariates, outcomes, and the estimated unit fixed effects 
## using the year 2008 for illustrative purposes
data.pre2 <- subset(data.bpre, data.bpre$fyear==2008,  select = c("cik", "x1", "x2", "x3", "trt", "y", "mu"))
m.out2 <- matchit(trt ~ 0 + x1 + x2 + x3 + y + mu, data = data.pre2, caliper = 0.2)
matched_df2 <- match.data(m.out2)
id.m2 <- unique(matched_df2$cik)
dat.match2 <- subset(data, cik %in% id.m2)

## Stage 2: DiD estimation (model 3)
m3.qmle <- glm(y ~ post + trt + did + x1 + x2 + x3
               + mu + 
                 mu*td0 + mu*td1 + mu*td2, data = dat.match2, family=quasipoisson(link='log'))

# Robust standard error
coeftest(m3.qmle, vcov = vcovHC(m3.qmle, type="HC1"))

# Mean Squared Error 
m3.qmle.mse <- mean((m3.qmle[["residuals"]])^2)



