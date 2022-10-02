### Analysis of tumor data
### Taking tumor volume values for 0,1,...,8 weeks (the complete cases only)
### Performing a panel linear model based analysis

# calling the data:
dat1 <- read.csv("tumor_data_panel.csv",header = TRUE)
dat1 <- dat1[,-1]
head(dat1)
boxplot(dat1$tumor_vol~dat1$block)
boxplot(dat1$tumor_vol~dat1$treat)
boxplot(dat1$tumor_vol~dat1$week)
combo1=paste(dat1$block,"+",dat1$treat)
boxplot(dat1$tumor_vol~combo1)

par(mfrow=c(2,1))
boxplot(y[x2=="AdLib"]~wk[x2=="AdLib"],ylim=range(y))
boxplot(y[x2=="WM"]~wk[x2=="WM"],ylim=range(y))
par(mfrow=c(1,1))

library(plm)
data1=pdata.frame(dat1,index=c("ID","week"))
id=data1$ID
wk=data1$week
y=data1$tumor_vol
x1=data1$block
x2=data1$treat

library(lme4)
lme0 <- lmer(y~x1+x2+factor(wk)+(1|id),data=data1)
summary(lme0)

pool=plm(y~x1+x2,data=data1,index=c("ID","week"),model="pooling",effect="twoways")

#Test for panel effects
#H0: No panel effect
plmtest(pool,effect="twoway",type="ghm")		#Test for panel effect
plmtest(pool,effect="twoway",type="kw")

#Testing for time effect
plmtest(pool,effect="time",type="bp") #Significant time effects

#Test for individual effect
plmtest(pool,effect="individual",type="bp") #Significant individual effects

fix.t=plm(y~x1+x2,data=data1,effect="time",model="within",random.method=NULL)
ran.t=plm(y~x1+x2,data=data1,effect="time",model="random",random.method="walhus")
phtest(fix.t,ran.t) #Fail to reject the null => random effect models (for time effects)

ran.tw=plm(y~x1+x2,data=data1,effect="twoways",model="random",random.method="walhus",index=c("ID","week"))
summary(ran.tw) #both individual and time effects are random

fix.tw=plm(y~x1+x2+factor(wk),data=data1,effect="individual",model="random",random.method="walhus",index=c("ID","week"))
summary(fix.tw) #random individual effects and fixed time effects

# MODEL:
# random=plm(y~x1+x2,data=data1,index=c("ID","week"),model="random",random.method="walhus")
# random.time=plm(y~x1+x2+factor(week),data=data1,index=c("ID","week"),model="random",random.method="walhus")
# pFtest(random.time,random) #Null is rejected => use time-fixed
# summary(random.time)


data1.new = data1
data1.new$tumor_vol[data1$tumor_vol==0] = 1e-4
id=data1.new[,1]
wk=data1.new[,2]
y=log(data1.new[,3])
x1=data1.new[,4]
x2=data1.new[,5]

par(mfrow=c(2,2))
hist(y[x1=="L" & x2=="AdLib"],xlab="",main="L+AdLib",breaks=20)
hist(y[x1=="L" & x2=="WM"],xlab="",main="L+WM",breaks=20)
hist(y[x1=="OB" & x2=="AdLib"],xlab="",main="OB+AdLib",breaks=20)
hist(y[x1=="OB" & x2=="WM"],xlab="",main="OB+WM",breaks=20)
par(mfrow=c(1,1))

random.time.new=plm(y~x1+x2+factor(wk),data=data1.new,index=c("ID","week"),model="random",random.method="walhus")
summary(random.time.new)

library(lme4)
lme1 <- lmer(y~x1+x2+factor(wk)+(1|id),data=data1.new)
summary(lme1)


data.new2 = data1
data.new2 = data.new2[-which(data1$tumor_vol==0),]
id=data.new2[,1]
wk=data.new2[,2]
y=log(data.new2[,3])
x1=data.new2[,4]
x2=data.new2[,5]

par(mfrow=c(2,2))
hist(y[x1=="L" & x2=="AdLib"],xlab="",main="L+AdLib",breaks=20)
hist(y[x1=="L" & x2=="WM"],xlab="",main="L+WM",breaks=20)
hist(y[x1=="OB" & x2=="AdLib"],xlab="",main="OB+AdLib",breaks=20)
hist(y[x1=="OB" & x2=="WM"],xlab="",main="OB+WM",breaks=20)
par(mfrow=c(1,1))

random.new=plm(y~x1+x2+factor(wk),data=data.new2,index=c("ID","week"),model="random",random.method="walhus")
summary(random.new)

library(lme4)
lme2 <- lmer(y~x1+x2+factor(wk)+(1|id),data=data.new2)
summary(lme2)




