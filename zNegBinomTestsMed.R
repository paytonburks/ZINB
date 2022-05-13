####	Demand for Medical Care by the Elderly
#### Forward selection using Wald Test and LRT (Zero-inflated Poisson Regression Model)

library(MASS)
library("pscl")
library("stats4")
library("splines")
library("VGAM")
library(foreign)


#### Data
data<- read.dta("racd06data1healthcare.dta")
dt <- data[, c(1, 6:8, 13, 15, 16, 11, 18)]
X=dt[,2:9]
y=dt[,1]

library(olsrr)
ols_step_forward()


############    Wald Statisitcs - chi sq w 8df, p = 0.05 ~> 15.507  #################

######### Wald Statistic for beta1=0   #####
fit1= vglm(y ~ X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1)[4]
B=vcov(fit1)
w1=beta1hat^2/B[4,4]

######### Wald Statistic for beta3=0   #####
fit3= vglm(y ~ X[,3], family=zinegbinomial(zero = 1))
beta1hat=coef(fit3)[4]
B=vcov(fit3)
w3=beta1hat^2/B[4,4]

######### Wald Statistic for beta5=0   #####
fit5= vglm(y ~ X[,5], family=zinegbinomial(zero = 1))
beta1hat=coef(fit5)[4]
B=vcov(fit5)
w5=beta1hat^2/B[4,4]

######### Wald Statistic for beta6=0   #####
fit6= vglm(y ~ X[,6], family=zinegbinomial(zero = 1))
beta1hat=coef(fit6)[4]
B=vcov(fit6)
w6=beta1hat^2/B[4,4]


######### Wald Statistic for beta7=0   #####
fit7= vglm(y ~ X[,7], family=zinegbinomial(zero = 1))
beta1hat=coef(fit7)[4]
B=vcov(fit7)
w7=beta1hat^2/B[4,4]

c(w1,w3,w5,w6,w7)

### x1 added  ###


####  To test whether beta3=0 or beta5=0 or beta6=0 or beta7=0 when x1 already in the model ####

######### Wald Statistic for beta3=0 when x1 in the model   #####
fit1_3= vglm(y ~ X[,3]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_3)[4]
B=vcov(fit1_3)
w1_3=beta1hat^2/B[4,4]

######### Wald Statistic for beta5=0 when x1 in the model   #####
fit1_5= vglm(y ~ X[,5]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_5)[4]
B=vcov(fit1_5)
w1_5=beta1hat^2/B[4,4]

######### Wald Statistic for beta6=0 when x1 in the model   #####
fit1_6= vglm(y ~ X[,6]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_6)[4]
B=vcov(fit1_6)
w1_6=beta1hat^2/B[4,4]

######### Wald Statistic for beta7=0 when x1 in the model   #####
fit1_7= vglm(y ~ X[,7]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_7)[4]
B=vcov(fit1_7)
w1_7=beta1hat^2/B[4,4]


c(w1_3,w1_5, w1_6, w1_7)

### x3 added ###


####  To test whether beta5=0 or beta6=0 or beta7=0 when x1 and x3 is already in the model ####

######### Wald Statistic for beta5=0 when x1 and x3 in the model   #####
fit1_3_5= vglm(y ~ X[,5]+X[,3]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_3_5)[4]
B=vcov(fit1_3_5)
w1_3_5=beta1hat^2/B[4,4]

######### Wald Statistic for beta6=0 when x1 and x3 in the model   #####
fit1_3_6= vglm(y ~ X[,6]+X[,3]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_3_6)[4]
B=vcov(fit1_3_6)
w1_3_6=beta1hat^2/B[4,4]

######### Wald Statistic for beta7=0 when x1 and x3 in the model   #####
fit1_3_7= vglm(y ~ X[,7]+X[,3]+X[,1], family=zinegbinomial(zero = 1))
beta1hat=coef(fit1_3_7)[4]
B=vcov(fit1_3_7)
w1_3_7=beta1hat^2/B[4,4]

c(w1_3_5, w1_3_6, w1_3_7)

### none added, everything < chi sq. score ###


############Likelihood Ratio Statistics - chi sq w 8df, p = 0.05 ~> 15.507 #####

###   LRT Statistics for beta1=0   #####
fit0=vglm(y ~ 1, zinegbinomial(zero = 1))
fit1= vglm(y ~ X[,1], family=zinegbinomial(zero = 1))
lrt1=2*(logLik(fit1)-logLik(fit0))
#lrt1=lrtest(fit0, fit1) 

###   LRT Statistics for beta3=0   #####
fit0=vglm(y ~ 1, zinegbinomial(zero = 1))
fit3= vglm(y ~ X[,3], family=zinegbinomial(zero = 1))
lrt3=2*(logLik(fit3)-logLik(fit0))

###   LRT Statistics for beta5=0   #####
fit0=vglm(y ~ 1, zinegbinomial(zero = 1))
fit5= vglm(y ~ X[,5], family=zinegbinomial(zero = 1))
lrt5=2*(logLik(fit5)-logLik(fit0))


###   LRT Statistics for beta6=0   #####
fit0=vglm(y ~ 1, zinegbinomial(zero = 1))
fit6= vglm(y ~ X[,6], family=zinegbinomial(zero = 1))
lrt6=2*(logLik(fit6)-logLik(fit0))

###   LRT Statistics for beta7=0   #####
fit0=vglm(y ~ 1, zinegbinomial(zero = 1))
fit7= vglm(y ~ X[,7], family=zinegbinomial(zero = 1))
lrt7=2*(logLik(fit7)-logLik(fit0))


c(lrt1,lrt3,lrt5,lrt6,lrt7)

### x1 added ###


####  To test whether beta3=0 or beta5=0 or beta7=0 when x1 is already in the model ####

#########   LRT Statistic for beta3=0 when x1 in the model   #####
fit1=vglm(y ~ X[,1], zinegbinomial(zero = 1))
fit1_3= vglm(y ~ X[,1]+X[,3], family=zinegbinomial(zero = 1))
lrt1_3=2*(logLik(fit1_3)-logLik(fit1))

#########   LRT Statistic for beta5=0 when x1 in the model   #####
fit1=vglm(y ~ X[,1], zinegbinomial(zero = 1))
fit1_5= vglm(y ~ X[,1]+X[,5], family=zinegbinomial(zero = 1))
lrt1_5=2*(logLik(fit1_5)-logLik(fit1))

#########   LRT Statistic for beta5=0 when x1 in the model   #####
fit1=vglm(y ~ X[,1], zinegbinomial(zero = 1))
fit1_6= vglm(y ~ X[,1]+X[,6], family=zinegbinomial(zero = 1))
lrt1_6=2*(logLik(fit1_6)-logLik(fit1))

#########   LRT Statistic for beta5=0 when x1 in the model   #####
fit1=vglm(y ~ X[,1], zinegbinomial(zero = 1))
fit1_7= vglm(y ~ X[,1]+X[,7], family=zinegbinomial(zero = 1))
lrt1_7=2*(logLik(fit1_7)-logLik(fit1))

c(lrt1_3,lrt1_5, lrt1_6,lrt1_7)

### x3 added ###


####  To test whether beta5=0 or beta6=0 or beta7=0 when x1 and x3 is already in the model ####

######### LRT Statistic for beta5=0 when x1 and x3 in the model   #####
fit1_3=vglm(y ~ X[,1]+X[,3], zinegbinomial(zero = 1))
fit1_3_5= vglm(y ~ X[,1]+X[,3]+X[,5], family=zinegbinomial(zero = 1))
lrt1_3_5=2*(logLik(fit1_3_5)-logLik(fit1_3))

######### LRT Statistic for beta6=0 when x1 and x3 in the model   #####
fit1_3=vglm(y ~ X[,1]+X[,3], zinegbinomial(zero = 1))
fit1_3_6= vglm(y ~ X[,1]+X[,3]+X[,6], family=zinegbinomial(zero = 1))
lrt1_3_6=2*(logLik(fit1_3_6)-logLik(fit1_3))

######### LRT Statistic for beta7=0 when x1 and x3 in the model   #####
fit1_3=vglm(y ~ X[,1]+X[,3], zinegbinomial(zero = 1))
fit1_3_7= vglm(y ~ X[,1]+X[,3]+X[,7], family=zinegbinomial(zero = 1))
lrt1_3_7=2*(logLik(fit1_3_7)-logLik(fit1_3))

c(lrt1_3_5, lrt1_3_6, lrt1_3_7)

### x5 added ###

####  To test whether beta6=0 or beta7=0 when x1, x3 and x5 is already in the model ####
######### LRT Statistic for beta6=0 when x1, x3 and x5 in the model   #####
fit1_3_5=vglm(y ~ X[,1]+X[,3]+X[,5], zinegbinomial(zero = 1))
fit1_3_5_6= vglm(y ~ X[,1]+X[,3]+X[,5]+X[,6], family=zinegbinomial(zero = 1))
lrt1_3_5_6=2*(logLik(fit1_3_5_6)-logLik(fit1_3_5))


######### LRT Statistic for beta7=0 when x1, x3 and x5 in the model   #####
fit1_3_5=vglm(y ~ X[,1]+X[,3]+X[,5], zinegbinomial(zero = 1))
fit1_3_5_7= vglm(y ~ X[,1]+X[,3]+X[,5]+X[,7], family=zinegbinomial(zero = 1))
lrt1_3_5_7=2*(logLik(fit1_3_5_7)-logLik(fit1_3_5))

c(lrt1_3_5_6, lrt1_3_5_7)

### x6 added ###

####  To test whether beta7=0 when x1, x3, x5 and x6 is already in the model ####
######### LRT Statistic for beta7=0 when x1, x3, x5 and x6 in the model   #####
fit1_3_5_6=vglm(y ~ X[,1]+X[,3]+X[,5]+X[,6], zinegbinomial(zero = 1))
fit1_3_5_6_7= vglm(y ~ X[,1]+X[,3]+X[,5]+X[,6]+X[,7], family=zinegbinomial(zero = 1))
lrt1_3_5_6_7=2*(logLik(fit1_3_5_6_7)-logLik(fit1_3_5_6))

lrt1_3_5_6_7

### x7 not added ###


####### BIC Tests #######

X=data[,7:22]
y=data[,1]

# 1st step
bc=c()
for(i in 1:16)
{
  d=BIC( vglm(y ~ X[,i], family=zinegbinomial(zero = 1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x3 added - 24568.07 ###

# 2nd step
X1=X[,-3]
bc=c()
for(i in 1:15)
{
  d=BIC( vglm(y ~ X[,3]+X1[,i], family=zinegbinomial(zero = 1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x15 added BIC = 24506.03 ###

# 3rd step
X2=X1[,-14]
bc=c()
for(i in 1:14)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X2[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x2 added BIC = 24473.25 ###


# 4th step
X3=X2[,-2]
bc=c()
for(i in 1:13)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X3[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x12 added BIC = 24454.38 ###


# 5th step
X4=X3[,-10]
bc=c()
for(i in 1:12)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X[,12]+X4[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x16 added BIC = 24435.50 ###


# 6th step
X5=X4[,-12]
bc=c()
for(i in 1:11)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X[,12]+X[,16]+X5[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x1 added BIC = 24419.46 ###


# 7th step
X6=X5[,-1]
bc=c()
for(i in 1:10)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X[,12]+X[,16]+X[,1]+X6[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x10 added BIC = 24414.96 


# 8th step
X7=X6[,-7]
bc=c()
for(i in 1:9)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X[,12]+X[,16]+X[,1]+X[,10]+X7[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)


### x11 added BIC = 24413.47


# 9th step
X8=X7[,-7]
bc=c()
for(i in 1:8)
{
  d=BIC( vglm(y ~ X[,3]+X[,15]+X[,2]+X[,12]+X[,16]+X[,1]+X[,10]+X[,11]+X8[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

# STOP, don't add any more - BIC = 24416.82


###### AIC Tests ######

X=dt[,2:9]
y=data[,1]

# 1st step
bc=c()
for(i in 1:8)
{
  d=AIC( vglm(y ~ X[,i], family=zinegbinomial(zero = 1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x1 added AIC = 24751.58 ###

# 2nd step
X1=X[,-1]
bc=c()
for(i in 1:7)
{
  d=AIC( vglm(y ~ X[,1]+X1[,i], family=zinegbinomial(zero = 1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x3 added AIC = 24669.25 ###

# 3rd step
X2=X1[,-2]
bc=c()
for(i in 1:6)
{
  d=AIC( vglm(y ~ X[,1]+X[,3]+X2[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

# x8 added AIC = 24615.98 ###

# 4th step
X3=X2[,-6]
bc=c()
for(i in 1:5)
{
  d=AIC( vglm(y ~ X[,1]+X[,3]+X[,8]+X3[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x2 added AIC = 24567.01 ###

# 5th step
X4=X3[,-1]
bc=c()
for(i in 1:4)
{
  d=AIC( vglm(y ~ X[,1]+X[,2]+X[,3]+X[,8]+X4[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x5 added AIC = 24541.39 ###

# 6th step
X5=X4[,-2]
bc=c()
for(i in 1:3)
{
  d=AIC( vglm(y ~ X[,1]+X[,2]+X[,3]+X[,5]+X[,8]+X5[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x6 added AIC = 24519.88 ###

# 7th step
X6=X5[,-2]
bc=c()
for(i in 1:2)
{
  d=AIC( vglm(y ~ X[,1]+X[,2]+X[,3]+X[,5]+X[,6]+X[,8]+X6[,i], family=zinegbinomial(zero=1)) )
  bc=c(bc,d)
}
bc
min(bc)

### x7 added AIC = 24515.61 ###
