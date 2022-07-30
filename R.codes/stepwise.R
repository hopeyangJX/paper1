data1<-read.csv("mydataB2n-.csv")
library(jsonlite)
install.packages("jsonlite")
library(carData)
library(lme4)
library(Matrix)
library(lmerTest)
fix(data1)
fit0<-lmer(LnosZ~LDOC+LNO3+LCEC+LB+LE+LR+(1|block)+(1|Year), data=data1)
vif(fit0)
fit<-lm(LnosZ~LDOC+LNO3+LCEC+LB+LE+LR, data=data1)
summary(fit0)
summary(fit)
AIC(fit0,fit)
install.packages("MuMIn")
install.packages("sjPlot")
library(sjPlot)
?MuMIn
library(MuMIn)
r.squaredGLMM(fit0)
tab_model(fit0)
vif(fit)

data1<-read.csv("Fungidata.csv")
head(data1)
fit3<-lm(PATrichness~DOC+NH4+NO3+AP+CEC+pH+Biomass+Richness+Shannon+Eveness,data = data1)
vif(fit3)
#library(MASS)
fit3<-lm(PATrichness~DOC+NH4+NO3+pH+Biomass+Richness+Eveness,data = data1)
vif(fit3)
stepAIC(fit3,direction = "backward")
fit<-lm( PATrichness~DOC+NH4+NO3+pH+Biomass+Richness+Eveness,data = data1)
summary(fit)
relweights <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  barplot(t(import),names.arg=lbls,
          ylab="% of R-Square",
          xlab="Predictor Variables",
          main="Relative Importance of Predictor Variables",
          sub=paste("R-Square=", round(rsquare, digits=3)),
          ...)
  return(import)
}##计算预测变量的相对权重
relweights(fit, col="lightgrey")
summary(fit)
**
