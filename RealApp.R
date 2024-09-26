setwd("D:/Joint works/Touqeer Files/Discrete modeling/code/Github_code")
source("Functions.R")

############################################
#########Upheld complaints data-------------
##############################################
data<-read.csv("Automobile_Insurance_Company_Complaint_Rankings__Beginning_2009.csv")


x<- data$Upheld.Complaints
#print(x)
x<- na.omit(x)
plot(table(x1),main = paste("Insurance complaints of New York city"),  xlab = paste("Number of upheld complaints"),ylab = paste("Frequency"))
##degpd1--------------------------
start.time <- Sys.time()
degpd1<-fit.model(x, model =1, init=c(sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_degpd1<- round(end.time - start.time,2)
value<-3642.467
paste(round(c(degpd1$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")
##KS test
from.pdist.to.stepfun = function(pdist,a,b){
  stopifnot(is.finite(pdist(b)) & is.finite(pdist(a)))
  x = a:b
  y = c(0,sapply(a:(b-1),pdist),1)
  return(stepfun(x,y,f=0)) # f=0: right-continuous
}
stepFun.test = from.pdist.to.stepfun(pdist=function(x) (x+1)/6,a=0,b=2)
stopifnot(stepFun.test(-0.001)==0  & stepFun.test(0)>0 & stepFun.test(1)<1 & stepFun.test(2)==1)

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pdiscegpd(x, kappa = degpd1$fit$mle[1], sigma = degpd1$fit$mle[2], xi = degpd1$fit$mle[3], type = 1),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))

##degpd2--------------------------
start.time <- Sys.time()
degpd2<-fit.model(x, model =2, init=c(sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_degpd2<- round(end.time - start.time,2)
value<-3642.309
paste(round(c(degpd2$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pdiscegpd(x, kappa = degpd2$fit$mle[1], sigma = degpd2$fit$mle[2], xi = degpd2$fit$mle[3], type = 2),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))
##degpd3--------------------------
start.time <- Sys.time()
degpd3<-fit.model(x, model =3, init=c(sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_degpd3<- round(end.time - start.time,2)
value<-3642.442
paste(round(c(degpd3$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pdiscegpd(x, kappa = degpd3$fit$mle[1], sigma = degpd3$fit$mle[2], xi = degpd3$fit$mle[3], type = 3),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))
##
#Return Level plots--------DEGPD
pd<- rdiscegpd(length(x), kappa = degpd1$fit$mle[1], sigma = degpd1$fit$mle[2], xi=degpd1$fit$mle[3], type = 1)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd2$fit$mle[1], sigma = degpd2$fit$mle[2], xi=degpd2$fit$mle[3], type = 2)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd3$fit$mle[1], sigma = degpd3$fit$mle[2], xi=degpd3$fit$mle[3], type = 3)
returnlevelplot(x, pd)


############################################
#########Doctors visit data-------------
##############################################
library(zic)
data("docvisits")
x<- docvisits$docvisits
x<- na.omit(x)
plot(table(x), main = paste("Doctors visits to the hospital"),
     xlab = paste("Number of doctors visits"),ylab = paste("Frequency"),)
##zidegpd1--------------------------
start.time <- Sys.time()
zidegpd1<-fit.model(x, model =1, init=c(0.2,sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("zidiscegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_zidegpd1<- round(end.time - start.time,2)
value<-3865.037
paste(round(c(zidegpd1$fit$mle,aic(value,4),bic(value,4, length(x)),value),2),collapse = " & ")

###
stepFun.test = from.pdist.to.stepfun(pdist=function(x) (x+1)/6,a=0,b=2)
stopifnot(stepFun.test(-0.001)==0  & stepFun.test(0)>0 & stepFun.test(1)<1 & stepFun.test(2)==1)

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pzidiscegpd(x,pi=zidegpd1$fit$mle[1], kappa = zidegpd1$fit$mle[2], sigma = zidegpd1$fit$mle[3], xi = zidegpd1$fit$mle[4], type = 1),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))
##zidegpd2--------------------------
start.time <- Sys.time()
zidegpd2<-fit.model(x, model =2, init=c(0.2,sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("zidiscegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_zidegpd2<- round(end.time - start.time,2)
value<-3864.801
paste(round(c(zidegpd2$fit$mle,aic(value,4),bic(value,4, length(x)),value),2),collapse = " & ")

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pzidiscegpd(x,pi=zidegpd2$fit$mle[1], kappa = zidegpd2$fit$mle[2], sigma = zidegpd2$fit$mle[3], xi = zidegpd2$fit$mle[4], type = 2),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))

##zidegpd3--------------------------
start.time <- Sys.time()
zidegpd3<-fit.model(x, model =3, init=c(0.2,sqrt(sqrt(mean(x))),sqrt(sqrt(sd(x))),0.2), family=c("zidiscegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
end.time <- Sys.time()
time.taken_zidegpd3<- round(end.time - start.time,2)
value<-3864.437
paste(round(c(zidegpd3$fit$mle,aic(value,4),bic(value,4, length(x)),value),2),collapse = " & ")

##
B.bootstrap = 200
v = rep(NA, B.bootstrap)  # Initialize a vector to store the p-values
for (j in 1:B.bootstrap) {
  # Define the step function for the discrete EGPD fitted model
  stepFun = from.pdist.to.stepfun(
    pdist = function(x) pzidiscegpd(x,pi=zidegpd3$fit$mle[1], kappa = zidegpd3$fit$mle[2], sigma = zidegpd3$fit$mle[3], xi = zidegpd3$fit$mle[4], type = 3),
    a = 0, 
    b = max(x) + 1
  )
  
  # Perform the KS test with bootstrapping
  ks_result = dgof::ks.test(
    x = x, 
    y = stepFun, 
    alternative = "two.sided", 
    simulate.p.value = TRUE, 
    B = 200
  )
  
  # Store the p-value from the KS test
  v[j] = ks_result$p.value
}

# Output the vector of p-values from the bootstrapped KS tests
print(mean(v))
##
#Return Level plots-------ZIDEGPD
pd<- rzidiscegpd(length(x),  pi = zidegpd1$fit$mle[1],kappa = zidegpd1$fit$mle[2], sigma = zidegpd1$fit$mle[3], xi=zidegpd1$fit$mle[4], type = 1)
returnlevelplot(x, pd)
pd<- rzidiscegpd(length(x),pi = zidegpd2$fit$mle[1],kappa = zidegpd2$fit$mle[2], sigma = zidegpd2$fit$mle[3], xi=zidegpd2$fit$mle[4], type = 2)
returnlevelplot(x, pd)
pd<- rzidiscegpd(length(x),pi = zidegpd3$fit$mle[1],kappa = zidegpd3$fit$mle[2], sigma = zidegpd3$fit$mle[3], xi=zidegpd3$fit$mle[4], type = 3)
returnlevelplot(x, pd)


############################################
#########Gaming and batting offensces data-------------
##############################################
data<- read.csv("Gaming and batting offensces.csv")
d<-data$Gaming_and_bating_ofense
plot(table(d),main = paste("Betting and gaming offences at\n New South Wales
Australia"),
     xlab = paste("Number of offences"),ylab = paste("Frequency"))

##Threshold at 10% quntile
u=floor(quantile(d, 0.1)) #threshold at 10% qunatile
u
x = d[d>=u]-u
t = table(x)
hist(x, main = paste("Betting and gaming offences at\n New South Wales
Australia"),
     xlab = paste("Number of offences"),ylab = paste("Frequency"), col = "lightblue", border = "black", breaks = 200)

degpd1<-fit.model(x, model =1, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =T, R = 100, ncpus = 1, plots =TRUE)
value<-1123.932
paste(round(c(degpd1$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")
degpd2<-fit.model(x, model =2, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
value<-1124.184
paste(round(c(degpd2$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

degpd3<-fit.model(x, model =3, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)

value<-1123.979
paste(round(c(degpd3$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

##Return Level plots
#Return Level plots--------DEGPD
pd<- rdiscegpd(length(x), kappa = degpd1$fit$mle[1], sigma = degpd1$fit$mle[2], xi=degpd1$fit$mle[3], type = 1)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd2$fit$mle[1], sigma = degpd2$fit$mle[2], xi=degpd2$fit$mle[3], type = 2)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd3$fit$mle[1], sigma = degpd3$fit$mle[2], xi=degpd3$fit$mle[3], type = 3)
returnlevelplot(x, pd)


##Threshold at 20% quntile
u=floor(quantile(d, 0.2)) #threshold at 20% qunatile
u
x = d[d>=u]-u
t = table(x)
hist(x, main = paste("Betting and gaming offences at\n New South Wales
Australia"),
     xlab = paste("Number of offences"),ylab = paste("Frequency"), col = "lightblue", border = "black", breaks = 200)

degpd1<-fit.model(x, model =1, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =T, R = 100, ncpus = 1, plots =TRUE)
value<-939.45
paste(round(c(degpd1$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")
degpd2<-fit.model(x, model =2, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)
value<-939.7068
paste(round(c(degpd2$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

degpd3<-fit.model(x, model =3, init=c(0.5,0.9,0.2), family=c("discegpd"), confint =TRUE, R = 1000, ncpus = 1, plots = TRUE)

value<-939.6411
paste(round(c(degpd3$fit$mle,aic(value,3),bic(value,3, length(x)),value),2),collapse = " & ")

##Return Level plots
#Return Level plots--------DEGPD
pd<- rdiscegpd(length(x), kappa = degpd1$fit$mle[1], sigma = degpd1$fit$mle[2], xi=degpd1$fit$mle[3], type = 1)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd2$fit$mle[1], sigma = degpd2$fit$mle[2], xi=degpd2$fit$mle[3], type = 2)
returnlevelplot(x, pd)
pd<- rdiscegpd(length(x), kappa = degpd3$fit$mle[1], sigma = degpd3$fit$mle[2], xi=degpd3$fit$mle[3], type = 3)
returnlevelplot(x, pd)