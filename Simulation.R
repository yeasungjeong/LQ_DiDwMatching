rm(list = ls())

library(dplyr)
library(MatchIt)
library(plm)
library(tidyverse)
library(phtt)
library(Synth)
library("gsynth")
library(synthdid)

####### Prepare for the GM method 
# Step for coefficients
coefX <- function(yp,Xp,iXpXp) {
  betaX <- iXpXp %*% crossprod(Xp,yp)
  betaX
}
# Step for CEO fixed-effects  
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

# Number of Iterations
k <- 1000 

# Number of CEOs (N)
n <- 200  

# Total time period (T)
total.T <- c(5) 

## Data generating process
res.m1 <- res.sd1 <- NULL;
res.m2 <- res.sd2 <- NULL;

for (s in 1:2) {
  for (a in 1:length(total.T)) {
    max.time <- total.T[a]
    trt.time <- max.time - 2  # Time of intervention
    res1 <- res2 <- NULL;
    for (i in 1:k) {
    tryCatch({
      set.seed(i)
      print(i)
      
      ############  Data generating process ############  
      N <- n 
      T <- max.time
      T0 <- trt.time - 1
      
      dat <- expand.grid(id = 1:n, tp = 1:T) %>%
        arrange(id,tp) %>% group_by(id) %>%
        mutate(
          # Covariates
          x=rnorm(T, mean = 0, sd = 1),  
          # CEO fixed-effects
          a = rnorm(1, -1, 1),           
          x.m = ifelse(tp==T0, x, -100),
          x.m = max(x.m),
          # Coefficient of the covariate in the propensity score model
          beta.ps = 0.2,
        ) %>%
        ungroup()
      
      # Error term
      dat$err <- rnorm(N*T, 0, 0.3)
      
      if (s == 1) {
        # Absence of time-variant USD (Scenario 1)
        dat$ft <- 0
        
        # Generate pre-intervention outcomes
        dat$y.pre <- 0.3*dat$x + dat$a + dat$err
        
        # Generate propensity scores and assign treatment
        dat <- dat %>% arrange(id,tp) %>% group_by(id) %>%
          mutate(
            y.m = ifelse(tp==T0, y.pre, -100), 
            y.m = max(y.m),
            # Functional form for the propensity score model
            logit.treat = beta.ps*x.m + a, 
            p.treat =  exp(logit.treat) / (exp(logit.treat) + 1),
            # assigning the treatment group
            trt1 = (rbinom(1, 1, p.treat)), 
          ) %>%
          ungroup()
        
        # Generate post-intervention dummies
        dat$trt <- dat$trt1
        dat$td0<-0; dat$td0[dat$tp==(T0+1)]<-1;
        dat$td1<-0; dat$td1[dat$tp==(T0+2)]<-1;
        dat$td2<-0; dat$td2[dat$tp==(T0+3)]<-1;

        # Generate post-intervention outcomes
        beta <- 0.3
        dat <- dat %>% mutate(y = ifelse(tp>T0, trt + beta*x + a + err, y.pre)) 
        
        # Generate a treatment dummy
        dat$p <- 0
        dat$p[dat$tp>T0] <- 1 
        dat$did <- dat$p*dat$trt
        
        # Use only pre-intervention data for extracting unit fixed-effects
        dat.m <- subset(dat, dat$p==0) 
      } else {
        # Presence of time-variant USD (Scenario 1)
        dat$ft <- 3*sin(2*pi*dat$tp/T)
        dat$ft[which(dat$tp<=T0)] <- 0
        
        # Generate pre-intervention outcomes
        dat$y.pre <- 0.3*dat$x + dat$a + dat$a*dat$ft + dat$err
        
        # Generate propensity scores and assign treatment
        dat <- dat %>% arrange(id,tp) %>% group_by(id) %>%
          mutate(
            y.m = ifelse(tp==T0, y.pre, -100), 
            y.m = max(y.m),
            # Functional form for the propensity score model
            logit.treat = beta.ps*x.m + a, 
            p.treat =  exp(logit.treat) / (exp(logit.treat) + 1),
            # assigning the treatment group
            trt1 = (rbinom(1, 1, p.treat)), 
          ) %>%
          ungroup()
        
        # Generate post-intervention dummies
        dat$trt <- dat$trt1
        dat$td0<-0; dat$td0[dat$tp==(T0+1)]<-1;
        dat$td1<-0; dat$td1[dat$tp==(T0+2)]<-1;
        dat$td2<-0; dat$td2[dat$tp==(T0+3)]<-1;

        
        # Generate post-intervention outcomes
        beta <- 0.3
        dat <- dat %>% mutate(y = ifelse(tp>T0, trt + beta*x + a + ft*a + err, y.pre)) 
        
        # Generate a treatment dummy
        dat$p <- 0
        dat$p[dat$tp>T0] <- 1 
        dat$did <- dat$p*dat$trt
        
        # Use only pre-intervention data for extracting unit fixed-effects
        dat.m <- subset(dat, dat$p==0) 
      }
      
      
      ##### 1. DiD 
      out1 <- plm(y~did + x, data = dat, 
                  index = c("id", "tp"), effect = "twoways")
      
      # Obtain bias (true treatment effect = 1)
      bias.m1 <- (out1$coefficients[1]-1)
      
      ##### 2. Conventional hybrid model 
      #### Matching on covariates only 
      m.out1 <- matchit(trt ~ 0 + x.m, data = dat.m, caliper = 0.2)
      
      # Matched sample 
      matched_df1 <- match.data(m.out1)
      matched_df1 <- subset(matched_df1, select=-c(weights, subclass))
      
      # Keep only units that are matched   
      id.m1 <- unique(matched_df1$id)
      dat.match1 <- subset(dat, id %in% id.m1)
      
      # DiD estimator
      out2 <- lm(y ~ p + trt + did + x, data = dat.match1) 
      
      # Obtain bias (true treatment effect = 1)
      TE.match2 <- as.matrix(out2$coefficients[4])  # estimated treatment effect
      bias.m2 <- as.numeric(TE.match2 - 1) # bias
      
      ##### 3. The GM hybrid approach
      ### First stage: Extract unit fixed-effects 
      beta <- coefsbai(as.matrix(dat.m$y),as.matrix(dat.m$x),N,T0)
      id <- c(1:n)
      res.ai <- data.frame(id, beta$coefi)
      
      ## Add the estimated unit fixed-effects to the pre-intervention data 
      dat.m <- dat.m %>% 
        left_join(res.ai, by = "id") 
      
      m.out3 <- matchit(trt ~ 0 + x.m + y.m + beta.coefi, data = dat.m, caliper = 0.2)
      matched_df3 <- match.data(m.out3)
      matched_df3 <- subset(matched_df3, select=-c(weights, subclass))
      id.m3 <- unique(matched_df3$id)
      dat.match3 <- subset(dat, id %in% id.m3)
      dat.match3 <- dat.match3 %>% 
        left_join(res.ai, by = "id") 
      
      ### Second stage: DiD with the interaction terms
      out3 <- lm(y ~  p + trt + did + x + beta.coefi +
                   beta.coefi*td0 + beta.coefi*td1
                 , data = dat.match3) 
      summary(out3)
      TE.match3 <- as.matrix(out3$coefficients[4]) 
      bias.m3 <- as.numeric(TE.match3 - 1)

      ##### 4. IFE
      T.ife <- length(unique(dat$tp))
      N.ife <- length(unique(dat$id))
      Y <- matrix(dat$y, T.ife, N.ife)
      did <- matrix(dat$did,  T.ife, N.ife)
      x <- matrix(dat$x, T.ife, N.ife)
      TE.4 <- Eup(Y ~ -1 + did + x, 
                        additive.effects = "twoways")
      bias.m4 <- (as.numeric(TE.4$slope.para[1])-1)
      
      ########## SCM 
      dat.scm <- as.data.frame(dat)
      scm.trt <- subset(dat.scm, trt==1, select=c("id", "tp", "y", "x"))
      scm.trt2 <- scm.trt %>% group_by(tp) %>% mutate(mean.y = mean(y))
      scm.trt2 <- subset(scm.trt2, select=c("tp", "mean.y", "x"))
      scm.trt2 <- scm.trt2 %>% distinct(tp, mean.y, .keep_all = T)
      names(scm.trt2) <- c("tp", "y", "x")
      scm.trt2$id <- 0 
      
      dat.scm <- subset(dat.scm, trt == 0, select=c("id", "tp", "y", "x"))
      dat.scm <- rbind(scm.trt2, dat.scm)
      
      dat.scm$cms_id <- as.numeric(dat.scm$id)
      dat.scm$year <- as.numeric(dat.scm$tp)
      dat.scm$cms_cha <- as.character(dat.scm$cms_id)
      dat.scm <- data.frame(dat.scm)
      control.id <- unique(dat.scm$id[which(dat.scm$id!=0)])
      
        dataprep.out <-
          dataprep(dat.scm,
                   predictors = c("x"),
                   dependent     = "y",
                   unit.variable = "cms_id",
                   time.variable = "year",
                   unit.names.variable = "cms_cha",
                   treatment.identifier  = 0,
                   controls.identifier   = control.id,
                   time.predictors.prior = c(1:T0),
                   time.optimize.ssr     = c(1:T0),
                   time.plot             = c(1:T)
          )
        
        # Run synth
        synth.out <- synth(dataprep.out)
        gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot%*%synth.out$solution.w)
        te.scm <- mean(gaps[trt.time:T])
        bias.m5 <- (as.numeric(te.scm)-1)
      
      ########## GSCM 
      
      if (a == 2) {
          gsynth.out <- gsynth(y ~ did + x, data = dat,
                               index = c("id","tp"), 
                               CV = TRUE, r = c(0,1), se = TRUE, 
                               inference = "parametric", nboots = 1000,
                               parallel = TRUE,  min.T0 = 3)
          bias.m6 <- (as.numeric(gsynth.out$est.avg[1])-1)} 
      
      if (a == 3) {
        gsynth.out <- gsynth(y ~ did + x, data = dat,
                             index = c("id","tp"), 
                             CV = TRUE, r = c(0,3), se = TRUE, 
                             inference = "parametric", nboots = 1000,
                             parallel = TRUE,  min.T0 = 5)
        bias.m6 <- (as.numeric(gsynth.out$est.avg[1])-1)} 
      
      if (a > 3) {
                                 
        gsynth.out <- gsynth(y ~ did + x, data = dat,
                             index = c("id","tp"),
                             CV = TRUE, se = TRUE, 
                             inference = "parametric", nboots = 1000,
                             parallel = TRUE)
        bias.m6 <- (as.numeric(gsynth.out$est.avg[1])-1)}
        
      ########## SDID 
      dat <- as.data.frame(dat)
      panel.sim <- subset(dat, select=c("id", "tp", "y", "did"))
      setup = panel.matrices(panel.sim)
      dat <- dat %>% arrange(tp, id)
      N.sdid <- length(unique(dat$id))
      T.sdid <- length(unique(dat$tp))
      cov.sdid <-  array(dat$x, dim =c(N.sdid, T.sdid, 1))
      tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0, cov.sdid)
      bias.m7 <- (as.numeric(tau.hat[1]-1))
      
      if (a == 1 ) {
      res1 <-rbind(res1, data.frame(cbind(scenario = as.numeric(s),
                                       bias.did = bias.m1, 
                                       bias.conv_hybrid = bias.m2, 
                                       bias.proposed = bias.m3,
                                       bias.ife = bias.m4, 
                                       bias.scm = bias.m5, 
                                       bias.sdid = bias.m7)))} else {
      res2 <-rbind(res2, data.frame(cbind(scenario = as.numeric(s),
                                          bias.did = bias.m1, 
                                          bias.conv_hybrid = bias.m2, 
                                          bias.proposed = bias.m3,
                                          bias.ife = bias.m4, 
                                          bias.scm = bias.m5, 
                                          bias.gscm = bias.m6,
                                          bias.sdid = bias.m7)))}
      
    }, error=function(e){})
  }
    
    if (a == 1) {
    m.bias1 <-cbind(scenario=as.numeric(s),
                   T0 = T0,
                   t(round(apply(res1, 2, mean)[2:7],3)))
    
    sd.bias1 <-cbind(scenario=as.numeric(s),
                   T0=T0,
                   t(round(apply(res1, 2, sd)[2:7],3) ))
    
    res.m1<-rbind(res.m1, m.bias1);
    res.sd1<-rbind(res.sd1, sd.bias1)} else {
    m.bias2<-cbind(scenario=as.numeric(s),
                     T0 = T0,
                     t(round(apply(res2, 2, mean)[2:8],3)))
      
    sd.bias2 <-cbind(scenario=as.numeric(s),
                      T0=T0,
                      t(round(apply(res2, 2, sd)[2:8],3) ))
      
     res.m2<-rbind(res.m2, m.bias2);
     res.sd2<-rbind(res.sd2, sd.bias2)}
  
  }
}

res.m1; res.sd1
res.m2; res.sd2

library(foreign)
write.dta(data.frame(res.m), "C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_mean1.dta")
write.dta(data.frame(res.sd), "C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_sd1.dta")

###############################################################################################

res <- NULL; 
k <- 1000
n <- 200  # sample size for each iteration
#total.T <- c(5, 6) # Time period
total.T <- c(8, 13, 18, 23) # Time period

## Data generating process
res.m <- res.sd <- NULL;
for (s in 1:2) {
  for (a in 1:length(total.T)) {
    max.time <- total.T[a]
    trt.time <- max.time - 2 # time point for the intervention
    res <- NULL;
    for (i in 1:k) {
      tryCatch({
        set.seed(i)
        print(i)
        
        ############  Data generating process ############  
        N <- n 
        T <- max.time
        T0 <- trt.time - 1
        
        dat <- expand.grid(id = 1:n, tp = 1:T) %>%
          arrange(id,tp) %>% group_by(id) %>%
          mutate(
            # Covariates
            x=rnorm(T, mean = 0, sd = 1),  
            # Unit fixed-effects
            a = rnorm(1, -1, 1),           
            x.m = ifelse(tp==T0, x, -100),
            x.m = max(x.m),
            # Coefficient of the covariate in the propensity score model
            beta.ps = 0.2,
          ) %>%
          ungroup()
        
        # Error term
        dat$err <- rnorm(N*T, 0, 0.3)
        
        if (s == 1) {
          # Absence of time-variant USD (Scenario 1)
          dat$ft <- 0
          
          # Generate pre-intervention outcomes
          dat$y.pre <- 0.3*dat$x + dat$a + dat$err
          
          # Generate propensity scores and assign treatment
          dat <- dat %>% arrange(id,tp) %>% group_by(id) %>%
            mutate(
              y.m = ifelse(tp==T0, y.pre, -100), 
              y.m = max(y.m),
              err.ps = rnorm(1, 0, 0.5),
              # Functional form for the propensity score model
              logit.treat = beta.ps*x.m + a, 
              p.treat =  exp(logit.treat) / (exp(logit.treat) + 1),
              # assigning the treatment group
              trt1 = (rbinom(1, 1, p.treat)), 
            ) %>%
            ungroup()
          
          # Generate post-intervention dummies
          dat$trt <- dat$trt1
          dat$td0<-0; dat$td0[dat$tp==(T0+1)]<-1;
          dat$td1<-0; dat$td1[dat$tp==(T0+2)]<-1;
          dat$td2<-0; dat$td2[dat$tp==(T0+3)]<-1;
          
          # Generate post-intervention outcomes
          beta <- 0.3
          dat <- dat %>% mutate(y = ifelse(tp>T0, trt + beta*x + a + err, y.pre)) 
          
          # Generate a treatment dummy
          dat$p <- 0
          dat$p[dat$tp>T0] <- 1 
          dat$did <- dat$p*dat$trt
          
          # Use only pre-intervention data for extracting unit fixed-effects
          dat.m <- subset(dat, dat$p==0) 
        } else {
          # Presence of time-variant USD (Scenario 1)
          dat$ft <- 3*sin(2*pi*dat$tp/T)
          dat$ft[which(dat$tp<=T0)] <- 0
          
          # Generate pre-intervention outcomes
          dat$y.pre <- 0.3*dat$x + dat$a + dat$a*dat$ft + dat$err
          
          # Generate propensity scores and assign treatment
          dat <- dat %>% arrange(id,tp) %>% group_by(id) %>%
            mutate(
              y.m = ifelse(tp==T0, y.pre, -100), 
              y.m = max(y.m),
              err.ps = rnorm(1, 0, 0.5),
              # Functional form for the propensity score model
              logit.treat = beta.ps*x.m + a, 
              p.treat =  exp(logit.treat) / (exp(logit.treat) + 1),
              # assigning the treatment group
              trt1 = (rbinom(1, 1, p.treat)), 
            ) %>%
            ungroup()
          
          # Generate post-intervention dummies
          dat$trt <- dat$trt1
          dat$td0<-0; dat$td0[dat$tp==(T0+1)]<-1;
          dat$td1<-0; dat$td1[dat$tp==(T0+2)]<-1;
          dat$td2<-0; dat$td2[dat$tp==(T0+3)]<-1;
          
          
          # Generate post-intervention outcomes
          beta <- 0.3
          dat <- dat %>% mutate(y = ifelse(tp>T0, trt + beta*x + a + ft*a + err, y.pre)) 
          
          # Generate a treatment dummy
          dat$p <- 0
          dat$p[dat$tp>T0] <- 1 
          dat$did <- dat$p*dat$trt
          
          # Use only pre-intervention data for extracting unit fixed-effects
          dat.m <- subset(dat, dat$p==0) 
        }
        
        ##### 1. DiD 
        out1 <- plm(y~did + x, data = dat, 
                    index = c("id", "tp"), effect = "twoways")
        
        # Obtain bias (true treatment effect = 1)
        bias.m1 <- (out1$coefficients[1]-1)
        
        ##### 2. Conventional hybrid model 
        #### Matching on covariates only 
        m.out1 <- matchit(trt ~ 0 + x.m, data = dat.m, caliper = 0.2)
        
        # Matched sample 
        matched_df1 <- match.data(m.out1)
        matched_df1 <- subset(matched_df1, select=-c(weights, subclass))
        
        # Keep only units that are matched   
        id.m1 <- unique(matched_df1$id)
        dat.match1 <- subset(dat, id %in% id.m1)
        
        # DiD estimator
        out2 <- lm(y ~ p + trt + did + x, data = dat.match1) 
        
        # Obtain bias (true treatment effect = 1)
        TE.match2 <- as.matrix(out2$coefficients[4])  # estimated treatment effect
        bias.m2 <- as.numeric(TE.match2 - 1) # bias
        
        ##### 3. The GM hybrid approach
        ### First stage: Extract unit fixed-effects 
        N <- n 
        T0 <- trt.time - 1
        
        beta <- coefsbai(as.matrix(dat.m$y),as.matrix(dat.m$x),N,T0)
        id <- c(1:n)
        res.ai <- data.frame(id, beta$coefi)
        
        ## Add the estimated unit fixed-effects to the pre-intervention data 
        dat.m <- dat.m %>% 
          left_join(res.ai, by = "id") 
        
        m.out3 <- matchit(trt ~ 0 + x.m + y.m + beta.coefi, data = dat.m, caliper = 0.2)
        matched_df3 <- match.data(m.out3)
        matched_df3 <- subset(matched_df3, select=-c(weights, subclass))
        id.m3 <- unique(matched_df3$id)
        dat.match3 <- subset(dat, id %in% id.m3)
        dat.match3 <- dat.match3 %>% 
          left_join(res.ai, by = "id") 
        
        ### Second stage: DiD with the interaction terms
        out3 <- lm(y ~  p + trt + did + x + beta.coefi +
                     beta.coefi*td0 + beta.coefi*td1 + beta.coefi*td2
                   , data = dat.match3) 
        summary(out3)
        TE.match3 <- as.matrix(out3$coefficients[4]) 
        bias.m3 <- as.numeric(TE.match3 - 1)
        
        ##### 4. IFE
        T.ife <- length(unique(dat$tp))
        N.ife <- length(unique(dat$id))
        Y <- matrix(dat$y, T.ife, N.ife)
        did <- matrix(dat$did,  T.ife, N.ife)
        x <- matrix(dat$x, T.ife, N.ife)
        TE.4 <- Eup(Y ~ -1 + did + x, 
                    additive.effects = "twoways")
        bias.m4 <- (as.numeric(TE.4$slope.para[1])-1)
        
        ########## SCM 
        dat.scm <- as.data.frame(dat)
        scm.trt <- subset(dat.scm, trt==1, select=c("id", "tp", "y", "x"))
        scm.trt2 <- scm.trt %>% group_by(tp) %>% mutate(mean.y = mean(y))
        scm.trt2 <- subset(scm.trt2, select=c("tp", "mean.y", "x"))
        scm.trt2 <- scm.trt2 %>% distinct(tp, mean.y, .keep_all = T)
        names(scm.trt2) <- c("tp", "y", "x")
        scm.trt2$id <- 0 
        
        dat.scm <- subset(dat.scm, trt == 0, select=c("id", "tp", "y", "x"))
        dat.scm <- rbind(scm.trt2, dat.scm)
        
        dat.scm$cms_id <- as.numeric(dat.scm$id)
        dat.scm$year <- as.numeric(dat.scm$tp)
        dat.scm$cms_cha <- as.character(dat.scm$cms_id)
        dat.scm <- data.frame(dat.scm)
        control.id <- unique(dat.scm$id[which(dat.scm$id!=0)])
        
        dataprep.out <-
          dataprep(dat.scm,
                   predictors = c("x"),
                   dependent     = "y",
                   unit.variable = "cms_id",
                   time.variable = "year",
                   unit.names.variable = "cms_cha",
                   treatment.identifier  = 0,
                   controls.identifier   = control.id,
                   time.predictors.prior = c(1:T0),
                   time.optimize.ssr     = c(1:T0),
                   time.plot             = c(1:T)
          )
        
        # Run synth
        synth.out <- synth(dataprep.out)
        gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot%*%synth.out$solution.w)
        te.scm <- mean(gaps[trt.time:T])
        bias.m5 <- (as.numeric(te.scm)-1)
        
        
        ########## GSCM 
        if (a == 1) {
          gsynth.out <- gsynth(y ~ did + x, data = dat,
                               index = c("id","tp"), 
                               CV = TRUE, r = c(0,3), se = TRUE, 
                               inference = "parametric", nboots = 1000,
                               parallel = TRUE,  min.T0 = 5)} else {
                                 
          gsynth.out <- gsynth(y ~ did + x, data = dat,
                               index = c("id","tp"),
                               CV = TRUE, se = TRUE, 
                               inference = "parametric", nboots = 1000,
                               parallel = TRUE)}
        
          bias.m6 <- (as.numeric(gsynth.out$est.avg[1])-1)
        
        ########## SDID 
        library(synthdid)
        
          dat <- as.data.frame(dat)
          panel.sim <- subset(dat, select=c("id", "tp", "y", "did"))
          setup = panel.matrices(panel.sim)
          dat <- dat %>% arrange(tp, id)
          N.sdid <- length(unique(dat$id))
          T.sdid <- length(unique(dat$tp))
          cov.sdid <-  array(dat$x, dim =c(N.sdid, T.sdid, 1))
          tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0, cov.sdid)
          bias.m7 <- (as.numeric(tau.hat[1]-1))
        
        res <-rbind(res, data.frame(cbind(scenario = as.numeric(s),
                                          bias.did = bias.m1, 
                                          bias.conv_hybrid = bias.m2, 
                                          bias.proposed = bias.m3,
                                          bias.ife = bias.m4, 
                                          bias.scm = bias.m5, 
                                          bias.gscm = bias.m6,
                                          bias.sdid = bias.m7)))
      }, error=function(e){})
    }
    
    m.bias <-cbind(scenario=as.numeric(s),
                   T0 = T0,
                   t(round(apply(res, 2, mean)[2:8],3)))
    
    sd.bias <-cbind(scenario=as.numeric(s),
                    T0=T0,
                    t(round(apply(res, 2, sd)[2:8],3) ))
    
    res.m<-rbind(res.m, m.bias);
    res.sd<-rbind(res.sd, sd.bias);
  }
}

res.m; res.sd

library(foreign)
write.dta(data.frame(res.m), "C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_mean2.dta")
write.dta(data.frame(res.sd), "C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_sd2.dta")



dd <- read.dta("C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_mean1.dta")
dd2 <- read.dta("C:/Users/yj547242/Dropbox/[1] Projects/2023/DiD-LQ/RR/Responses/result_sd1.dta")


