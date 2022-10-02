### Analysis of tumor data
### Taking tumor volume values for 0,1,...,8 weeks (the complete cases only)
### Performing a linear mixed model based analysis

# calling the data:
dat1 <- read.csv("tumor_data_mixed_model.csv",header = TRUE)
dat1 <- dat1[,-1]
head(dat1)
#attach(dat1)

library(lme4)
mod0 <- lmer(tumor_vol~block*treat+as.factor(week)+(1|rat_ID)+(1|rat_ID:tumor_nos),data = dat1)
summary(mod0)

#variance of rat_ID is zero, so, let's ignore rat effect:
mod1 <- lmer(tumor_vol~block*treat+as.factor(week)+(1|rat_ID:tumor_nos),data = dat1)
summary(mod1)

#what if we ignore tumor effect:
mod2 <- lmer(tumor_vol~block*treat+as.factor(week)+(1|rat_ID),data = dat1)
summary(mod2)

#which model amongst mod0, mod1 and mod2:
anova(mod0,mod1,mod2) #mod1 has minimum AIC & BIC

#diagnostic plots:
fit <- fitted(mod1)
par(mfrow=c(2,2))
plot(dat1$tumor_vol,fit,type = "p",pch=20,xlab="tumor volume",ylab="fitted volume")
abline(a=0,b=1,col=2,lwd=2)

qqnorm(resid(mod1),main = "residuals",pch=20)
qqnorm(ranef(mod1)[[1]][,1],main = "tumor effect (random)",pch=20)

plot(fit,resid(mod1),xlab="fitted",ylab="residuals",pch=20)
abline(h=0,col=2,lwd=2)
par(mfrow=c(1,1))


#gamma regression using glm
glm1 <- glm((tumor_vol+0.001)~block*treat+as.factor(week),data = dat1,
            family = Gamma(link = "log"))
summary(glm1)

glm2 <- glm((tumor_vol+0.001)~block+treat+as.factor(week),data = dat1,
            family = Gamma(link = "log"))
summary(glm2)

glm3 <- glm((tumor_vol+0.0001)~block*treat+as.factor(week),data = dat1,
            family = Gamma(link = "inverse"))
summary(glm3)

#Replacing 0 in the data by 1e-3
dat2=dat1
dat2$tumor_vol[dat1$tumor_vol==0] = 1e-3

#gamma regression using glm
glm1 <- glm(tumor_vol~block+treat+as.factor(week),data = dat2,
            family = Gamma(link = "log"))
summary(glm1)
par(mfrow=c(2,2)); plot(glm1); par(mfrow=c(1,1))

#mixed effect gamma regression:
library(lme4)
glm0 <- glmer((tumor_vol+0.001)~block*treat+as.factor(week)+(1|rat_ID)+(1|rat_ID:tumor_nos),
              data = dat1,family = Gamma(link = "log"))
summary(glm0)
#replacing 0 by 1e-3
glm00 <- glmer(tumor_vol~block*treat+as.factor(week)+(1|rat_ID)+(1|rat_ID:tumor_nos),
              data = dat2,family = Gamma(link = "log"))
summary(glm00)

glm1 <- glmer((tumor_vol+0.001)~block*treat+as.factor(week)+(1|rat_ID:tumor_nos),
              data = dat1,family = Gamma(link = "log"))
summary(glm1)

glm2 <- glmer((tumor_vol+0.001)~block*treat+as.factor(week)+(1|rat_ID),
              data = dat1,family = Gamma(link = "log"))
summary(glm2)





######################################################################
############ Zero inflated mixed effects gamma regression ############
######################################################################
library(glmmTMB) #library(glmmADMB)

# calling the data:
dat1 <- read.csv("tumor_data_mixed_model.csv",header = TRUE)
dat1 <- dat1[,-1]
head(dat1)

par(mfrow=c(2,2))
hist(dat1$tumor_vol[dat1$block=="L" & dat1$treat=="AdLib"],breaks = 20,
     freq = FALSE,xlab="tumor volume",main="Lean group - AdLib diet")
hist(dat1$tumor_vol[dat1$block=="L" & dat1$treat=="WM"],breaks = 20,
     freq = FALSE,xlab="tumor volume",main="Lean group - WM diet")
hist(dat1$tumor_vol[dat1$block=="OB" & dat1$treat=="AdLib"],breaks = 20,
     freq = FALSE,xlab="tumor volume",main="Obese group - AdLib diet")
hist(dat1$tumor_vol[dat1$block=="OB" & dat1$treat=="WM"],breaks = 20,
     freq = FALSE,xlab="tumor volume",main="Obese group - WM diet")
par(mfrow=c(1,1))

#dat1$rat_ID=as.factor(dat1$rat_ID)
#dat1$tumor_nos=as.factor(dat1$tumor_nos)
#dat1$week=as.factor(dat1$week)
#glm0 <- glmmadmb(tumor_vol~block*treat+as.factor(week)+(1|rat_ID)+(1|rat_ID:tumor_nos),data = dat1,family = "gamma",zeroInflation = TRUE)

glm0 <- glmmTMB(tumor_vol~block*treat+as.factor(week),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
summary(glm0)

glm1 <- glmmTMB(tumor_vol~block*treat+as.factor(week),
                data = dat1,ziformula = ~block*treat+as.factor(week),
                family = ziGamma(link = "log"))
summary(glm1)

glm2 <- glmmTMB(tumor_vol~block*treat+as.factor(week),
                data = dat1,ziformula = ~block*treat,
                family = ziGamma(link = "log"))
summary(glm2)
#par(mfrow=c(2,2)); plot(glm0); par(mfrow=c(1,1))

glm0 <- glmmTMB(tumor_vol~block+treat+as.factor(week)+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
glm1 <- glmmTMB(tumor_vol~block+treat+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
glm2 <- glmmTMB(tumor_vol~as.factor(week)+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
#-2*(as.vector(logLik(glm1))-as.vector(logLik(glm0)))
#qchisq(0.05,df.residual(glm1)-df.residual(glm0),lower.tail = FALSE)
pchisq(as.numeric(-2*(as.vector(logLik(glm1))-as.vector(logLik(glm0)))),
       df.residual(glm1)-df.residual(glm0),lower.tail = FALSE)
pchisq(as.numeric(-2*(as.vector(logLik(glm2))-as.vector(logLik(glm0)))),
       df.residual(glm2)-df.residual(glm0),lower.tail = FALSE)


#Quadratic time effect:
glm0 <- glmmTMB(tumor_vol~block*treat+week+I(week^2)+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
glm1 <- glmmTMB(tumor_vol~block*treat+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
glm2 <- glmmTMB(tumor_vol~week+I(week^2)+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
#-2*(as.vector(logLik(glm1))-as.vector(logLik(glm0)))
#qchisq(0.05,df.residual(glm1)-df.residual(glm0),lower.tail = FALSE)
pchisq(as.numeric(-2*(as.vector(logLik(glm1))-as.vector(logLik(glm0)))),
       df.residual(glm1)-df.residual(glm0),lower.tail = FALSE)
pchisq(as.numeric(-2*(as.vector(logLik(glm2))-as.vector(logLik(glm0)))),
       df.residual(glm2)-df.residual(glm0),lower.tail = FALSE)


#fitting ZIgamma for 8 seperate weeks:
glm.w0 <- glmmTMB(tumor_vol~block*treat,data=subset(dat1,week==0),
                  ziformula = ~1,family = ziGamma(link = "log"))
summary(glm.w0)

#plotting p-values:
dp=dz=integer()
for(i in 0:8){
  glm.w <- glmmTMB(tumor_vol~block*treat+(1|rat_ID),data=subset(dat1,week==i),
                   ziformula = ~1,family = ziGamma(link = "log"))
  dp=rbind(dp,summary(glm.w)$coefficients$cond[,4])
  dz=rbind(dz,summary(glm.w)$coefficients$cond[,3])
}
dp <- data.frame(dp); names(dp) <- c("Intercept","blockOB","treatWM","blockOB.treatWM")
dz <- data.frame(dz); names(dz) <- c("Intercept","blockOB","treatWM","blockOB.treatWM")

plot(0:8,dp$treatWM,type="b",pch=20, lwd = 2, cex = 1.5,
     main="Statistical significance of WM over week",
     xlab="Week",ylab="P-value")
abline(h=c(0.005,0.05),lty=2,col=2, lwd = 2)
text(c(0,0), c(0.013, 0.058), c('P=0.005', 'P=0.05'), pos = 4,
     cex = 1.2, col = 'blue')

plot(0:8,dz$treatWM,type="b",pch=20,main="Z-values over week",
     xlab="week",ylab="treatment : WM Group")


#rat table:
# for(wk in 0:8){
#   print(paste('Wk', wk))
#   d.wk = dat1[dat1$week==wk,]
#   print(sort(table(d.wk$tumor_nos)))
# }


# tracking WM-week interaction effects:
dat1[,8]=as.factor(paste(dat1$treat,dat1$week,sep="-"))
names(dat1)=c(names(dat1)[-8],"new.treat")
head(dat1)
glm0 <- glmmTMB(tumor_vol~block+new.treat+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
summary(glm0)

mat=model.matrix(glm0)
mat1=mat[,1:10]
g1=glmmTMB(dat1$tumor_vol ~0+mat+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
g2=glmmTMB(dat1$tumor_vol ~0+mat1+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
as.numeric(-2*(logLik(g2)-logLik(g1))) #log-likelihood ratio value
qchisq(0.005,9,lower.tail = FALSE)
pchisq(as.numeric(-2*(logLik(g2)-logLik(g1))),9,lower.tail = FALSE) #significant p-value
##### interaction effect due to week and diet is significant!!!
mat2=mat[,-(3:11)]
g3=glmmTMB(dat1$tumor_vol ~0+mat2+(1|rat_ID/tumor_nos),
                data = dat1,ziformula = ~1,family = ziGamma(link = "log"))
logLik(g3)
logLik(g2)




#detach(dat1)