library(fGarch)

library(goftest)

library(timeSeries)

library(VineCopula)

library(KScorrect)

library(stats)

library(ADGofTest)

library("readxl")

library("tseries")

library(readxl)
Data <- read_excel("Data.xlsx")

# Question 2a

port = 1/2*Data$`S&PLogReturns` + 1/2*Data$DAXLogReturns

sigma2_SP = var(Data$`S&PLogReturns`)
sigma2_DAX = var(Data$DAXLogReturns)
sigma_SP = sqrt(sigma2_SP)
sigma_DAX = sqrt(sigma2_DAX)
r_SPDAX = cor(Data$DAXLogReturns,Data$`S&PLogReturns`)
mu_SP = mean(Data$`S&PLogReturns`)
mu_DAX = mean(Data$DAXLogReturns)
r_SPDAX

sigma2_port = 1/4*sigma2_SP + 1/4*sigma2_DAX + 1/2*r_SPDAX*sigma_DAX*sigma_SP
sigma_port = sqrt(sigma2_port)
mu_port = 1/2*mu_SP + 1/2*mu_DAX

norminv_95 = qnorm(0.95,mean=0,sd=1)
norminv_99 = qnorm(0.99,mean=0,sd=1)

VaR_95 = norminv_95*sigma_port - mu_port #VaR 95
VaR_99 = norminv_99*sigma_port - mu_port #VaR 99
VaR_95
VaR_99

VaR_95NR = exp(VaR_95)-1
VaR_99NR = exp(VaR_99)-1
VaR_95NR
VaR_99NR

# Question 2b

# We form a time series plot to visually observe the data over time

time <- seq(2000,2018,18/989)
plot(time,Data$`S&PLogReturns`,type='l') # Time series plot
plot(time,Data$DAXLogReturns,type='l') # Time series plot

# Clearly we see some evidence of autocorrelation and clustered volatility

SPLR = Data$`S&PLogReturns`
DAXLR = Data$DAXLogReturns

# We carry out the Jarque-Bera test to check if the log returns of the data are normally distributed

jarque.bera.test(SPLR)

# Clearly we cannot use the normal distribution, as the p-values do not pass the test at the 0.05 sig level

# Now, we plot the autocorrelation functions of the data and the squared data, to test for autocorrelation and GARCH effects

acf(SPLR)

# The ACF plot shows significant lag at 1, though the p-value is very slightly above 0.05 

pacf(SPLR)

# The PACF plot shows a significant lag at 1, so we think we should proceed with an AR(1) model

acf(SPLR^2)

# Also, we see GARCH effects, so we must incorporate this into our model

# We try various AR(1)-GARCH models, with different conditional distributions, and check which of these has the best information criteria

modelSPnorm = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "norm")

modelSPstd = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "std")

modelSPged = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "ged")

modelSPsstd = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "sstd")

modelSPsnorm = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "snorm")

modelSPsged = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "sged")

ic_SP = rbind(modelSPnorm@fit$ics, modelSPstd@fit$ics, modelSPged@fit$ics, modelSPsstd@fit$ics, modelSPsnorm@fit$ics, modelSPsged@fit$ics)
row.names(ic_SP) <- c('norm','std', 'ged', 'sstd', 'snorm', 'sged')
ic_SP

# We see that 'sstd', or Skew Student-t distribution has the lowest AIC/BIC, so we proceed with this choice

# Now, we test various GARCH(p,q) models with 'sstd' as the conditional distribution, to check their information criteria

modelSP11 = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "sstd")
modelSP12 = garchFit(formula= ~arma(1,0) + garch(1,2),data=SPLR,trace = F, cond.dist = "sstd")
modelSP13 = garchFit(formula= ~arma(1,0) + garch(1,3),data=SPLR,trace = F, cond.dist = "sstd")
modelSP21 = garchFit(formula= ~arma(1,0) + garch(2,1),data=SPLR,trace = F, cond.dist = "sstd")
modelSP22 = garchFit(formula= ~arma(1,0) + garch(2,2),data=SPLR,trace = F, cond.dist = "sstd")
modelSP23 = garchFit(formula= ~arma(1,0) + garch(2,3),data=SPLR,trace = F, cond.dist = "sstd")
modelSP31 = garchFit(formula= ~arma(1,0) + garch(3,1),data=SPLR,trace = F, cond.dist = "sstd")
modelSP32 = garchFit(formula= ~arma(1,0) + garch(3,2),data=SPLR,trace = F, cond.dist = "sstd")
modelSP33 = garchFit(formula= ~arma(1,0) + garch(3,3),data=SPLR,trace = F, cond.dist = "sstd")

ic_pqSP = rbind(modelSP11@fit$ics, modelSP12@fit$ics, modelSP13@fit$ics, modelSP21@fit$ics, modelSP22@fit$ics, modelSP23@fit$ics, modelSP31@fit$ics, modelSP32@fit$ics, modelSP33@fit$ics)
row.names(ic_pqSP) <- c('1,1','1,2','1,3','2,1','2,2','2,3','3,1','3,2','3,3')
ic_pqSP

# We find that AR(1)-GARCH(1,1) is the best fit for the data, so we proceed with AR(1)-GARCH(1,1) with sstd distribution

modelSP11 = garchFit(formula= ~arma(1,0) + garch(1,1),data=SPLR,trace = F, cond.dist = "sstd")

# Now, just to check, we extract the residuals, and check for autocorrelation in squared residuals

res_SP11 <- residuals(modelSP11, standardize = TRUE)
acf(res_SP11^2, main = "ACF of squared residuals of AR(1)-GARCH(1,1) model for S&P 500 Log Returns")

# We have significant lag = 1, so it turns out our model was not adequate.

# Instead, we try the GARCH(p,q) model with next best information criteria, namely GARCH(1,2)

modelSP12 = garchFit(formula= ~arma(1,0) + garch(1,2),data=SPLR,trace = F, cond.dist = "sstd")

res_SP12 <- residuals(modelSP12, standardize = TRUE)
acf(res_SP12^2, main = "ACF of squared residuals of AR(1)-GARCH(1,2) model for S&P 500 Log Returns")

# We still have significant lags. Perhaps we need to try a different model?
# Come to think of it, the lag significance at lag 1 in the original ACF/PACF was quite small, so it may just be randomness which caused this
# Furthermore, we did not see any sign of a geometrically decaying ACF, which would generally suggest true autocorrelation

# What happens if we try AR(0)-GARCH(1,2) instead?

modelSP12 = garchFit(formula= ~arma(0,0) + garch(1,2),data=SPLR,trace = F, cond.dist = "sstd")

res_SP12 <- residuals(modelSP12, standardize = TRUE)
acf(res_SP12^2, main = "ACF of squared residuals of AR(0)-GARCH(1,2) model for S&P 500 Log Returns")

# We no longer have significant lags, so we proceed with this model
# We carry out the Ljung-Box test just to be sure

Box.test(res_SP12, lag=10, type = "Ljung-Box", fitdf = 0)

# We pass the Ljung-Box test, so we proceed

# Now, to apply the PIT to the data. We do this by first plotting a histogram of the residuals, to see what may fit

hist(res_SP12,breaks=100)

# This seems like it has a slight skew;

# We use 'sstdFit' to fit an appropriate skew Student's t distribution

SSTDFITSP <- sstdFit(res_SP12)

# Now, we apply the probability integral transform

u_SP12 <- psstd(res_SP12, mean = SSTDFITSP$estimate["mean"], sd = SSTDFITSP$estimate["sd"], nu = SSTDFITSP$estimate["nu"], xi= SSTDFITSP$estimate["xi"])
hist(u_SP12)

# This looks uniform, but we apply the KS test just to be sure

KStest_SP12 <- LcKS(u_SP12, cdf = "punif")
KStest_SP12$p.value

# We pass the KS test, so we finalise this model

# To summarise, we found that the best fit was a GARCH(1,2) model with conditional distribution Skewed t
# We apply skewed t to the residuals, using SstdFit, and found a uniform distribution, implying that our model choice was accurate




# DAX Returns model fitting

# We carry out the Jarque-Bera test to check if the log returns of the data are normally distributed

jarque.bera.test(DAXLR)

# Clearly we cannot use the normal distribution, as the p-value does not pass the test at the 0.05 sig level

# Now, we plot the autocorrelation functions of the data and the squared data, to test for autocorrelation and GARCH effects

acf(DAXLR, main = "Autocorrelation function of DAX log-returns")

# The ACF plot shows significant lag at lag 3
# However, note that the ACF plot does not show geometrical decay, so it does not imply clear autocorrelation

# Now, we plot the PACF to double check

pacf(DAXLR, main = "Partial autocorrelation function of DAX log-returns")

# Indeed there is still significant lag at lag 3

# This suggests we could proceed with AR(3)

acf(DAXLR^2)

# Also, we see GARCH effects, so we must incorporate this into our model

# We try various AR(3)-GARCH models, with different conditional distributions, and check which of these has the best information criteria

modelDAXnorm = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "norm")

modelDAXstd = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "std")

modelDAXged = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "ged")

modelDAXsstd = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sstd")

modelDAXsnorm = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "snorm")

modelDAXsged = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sged")

ic_DAX = rbind(modelDAXnorm@fit$ics, modelDAXstd@fit$ics, modelDAXged@fit$ics, modelDAXsstd@fit$ics, modelDAXsnorm@fit$ics, modelDAXsged@fit$ics)
row.names(ic_DAX) <- c('norm','std', 'ged', 'sstd', 'snorm', 'sged')
ic_DAX

# We see that 'sstd', or Skew Student-t distribution has the lowest AIC/BIC, so we proceed with this choice

# Now, we test various GARCH(p,q) models with 'sstd' as the conditional distribution, to check their information criteria

modelDAX11 = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX12 = garchFit(formula= ~arma(3,0) + garch(1,2),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX13 = garchFit(formula= ~arma(3,0) + garch(1,3),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX21 = garchFit(formula= ~arma(3,0) + garch(2,1),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX22 = garchFit(formula= ~arma(3,0) + garch(2,2),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX23 = garchFit(formula= ~arma(3,0) + garch(2,3),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX31 = garchFit(formula= ~arma(3,0) + garch(3,1),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX32 = garchFit(formula= ~arma(3,0) + garch(3,2),data=DAXLR,trace = F, cond.dist = "sstd")
modelDAX33 = garchFit(formula= ~arma(3,0) + garch(3,3),data=DAXLR,trace = F, cond.dist = "sstd")

ic_pqDAX = rbind(modelDAX11@fit$ics, modelDAX12@fit$ics, modelDAX13@fit$ics, modelDAX21@fit$ics, modelDAX22@fit$ics, modelDAX23@fit$ics, modelDAX31@fit$ics, modelDAX32@fit$ics, modelDAX33@fit$ics)
row.names(ic_pqDAX) <- c('1,1','1,2','1,3','2,1','2,2','2,3','3,1','3,2','3,3')
ic_pqDAX

# We find that GARCH(1,1) is the best fit for the data, so we proceed with GARCH(1,1) with sstd distribution

modelDAX11 = garchFit(formula= ~arma(3,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sstd")

# Now, just to check, we extract the residuals, and check for autocorrelation in squared residuals

res_DAX11 <- residuals(modelDAX11, standardize = TRUE)
acf(res_DAX11^2)

# It seems that we do not have significant lags, so we could proceed with this model

# We carry out the Ljung-Box test just to be sure

Box.test(res_DAX11, lag=10, type = "Ljung-Box", fitdf = 3)

# Unfortunately, we do not pass the Ljung-Box test

# We try other GARCH models

res_DAX12 <- residuals(modelDAX12, standardize = TRUE)
Box.test(res_DAX12, lag=10, type = "Ljung-Box", fitdf = 3)

res_DAX13 <- residuals(modelDAX13, standardize = TRUE)
Box.test(res_DAX13, lag=10, type = "Ljung-Box", fitdf = 3)

res_DAX21 <- residuals(modelDAX21, standardize = TRUE)
Box.test(res_DAX21, lag=10, type = "Ljung-Box", fitdf = 3)

res_DAX22 <- residuals(modelDAX22, standardize = TRUE)
Box.test(res_DAX22, lag=10, type = "Ljung-Box", fitdf = 3)

res_DAX23 <- residuals(modelDAX23, standardize = TRUE)
Box.test(res_DAX23, lag=10, type = "Ljung-Box", fitdf = 3)

# Perhaps we need to try a new model? Come to think of it, the lag at 3 was not overwhelmingly significant
# Also, we did not see a geometrically decaying ACF plot, so it does not imply clear autocorrelation
# This implies that the lag at 3 may have just been a freak scenario, since the model is clearly not an appropriate fit

# Let us try AR(0), to see if it works

modelDAX11 = garchFit(formula= ~arma(0,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sstd")
res_DAX11 <- residuals(modelDAX11, standardize = TRUE)
acf(res_DAX11^2)

# Still, with AR(0) we do not see significant lags, but we must try the Ljung-Box test

# Let us try Ljung-Box test to be sure

Box.test(res_DAX11, lag=10, type = "Ljung-Box", fitdf = 0)

# We do pass the Ljung-Box test, which implies that our choice was accurate

# We continue. We note that our model choice was a guess, so we may have to further tweak the model along the way

# We look to apply PIT to the data, and then test for uniformity, as this would give us a clue as to whether it is appropriate

# Now, to apply the PIT to the data. We do this by first plotting a histogram of the residuals, to see what may fit

hist(res_DAX11,breaks=100)

# Now, we try sstdFit

SSTDFITDAX <- sstdFit(res_DAX11)

# We can apply probability integral transform

u_DAX11 <- psstd(res_DAX11, mean = SSTDFITDAX$estimate["mean"], sd = SSTDFITDAX$estimate["sd"], nu = SSTDFITDAX$estimate["nu"], xi= SSTDFITDAX$estimate["xi"])
hist(u_DAX11)

# This looks uniform, but we apply the KS test to be sure

KStest_DAX11 <- LcKS(u_DAX11, cdf = "punif")
KStest_DAX11$p.value

# We pass the KS test, so we can continue

# To summarise, we found that the best fit was a GARCH(1,1) model with conditional distribution Skewed Student t
# We apply skewed Student t using sstdFit to the residuals, and found a uniform distribution, implying that our model choice was accurate

# We save our final models in a clear format, to work with copula effectively

modelSP = garchFit(formula= ~arma(0,0) + garch(1,2),data=SPLR,trace = F, cond.dist = "sstd")
res_SP <- residuals(modelSP, standardize = TRUE)
u_SP <- psstd(res_SP, mean = SSTDFITSP$estimate["mean"], sd = SSTDFITSP$estimate["sd"], nu = SSTDFITSP$estimate["nu"], xi= SSTDFITSP$estimate["xi"])

modelDAX = garchFit(formula= ~arma(0,0) + garch(1,1),data=DAXLR,trace = F, cond.dist = "sstd")
res_DAX <- residuals(modelDAX, standardize = TRUE)
u_DAX <- psstd(res_DAX11, mean = SSTDFITDAX$estimate["mean"], sd = SSTDFITDAX$estimate["sd"], nu = SSTDFITDAX$estimate["nu"], xi= SSTDFITDAX$estimate["xi"])

# We now proceed to copula modelling

# Just to get an idea of the copula, we plot the uniform residuals against each other

plot(u_DAX, u_SP)

# This roughly resembles the T-copula, due to the tail dependence

# This makes sense, since financial returns generally have lower tail dependence, which makes the t-copula good for modelling financial returns

# We apply BiCopSelect to fit a copula

copmodel = BiCopSelect(u_DAX, u_SP, familyset = NA, selectioncrit="AIC", indeptest = TRUE, level = 0.05, se=TRUE)
copmodel

# As suspected we obtain a t-copula 

# Now, we proceed to simulating the copula

N=1000
cop_sim = BiCopSim(N, family = copmodel$family, copmodel$par, copmodel$par2)
resDAX_sim <- qsstd(cop_sim[,1], mean = SSTDFITDAX$estimate["mean"], sd = SSTDFITDAX$estimate["sd"], nu = SSTDFITDAX$estimate["nu"], xi= SSTDFITDAX$estimate["xi"])
resSP_sim = qsstd(cop_sim[,2], mean = SSTDFITSP$estimate["mean"], sd = SSTDFITSP$estimate["sd"], nu = SSTDFITSP$estimate["nu"], xi= SSTDFITSP$estimate["xi"])

# We have our simulated residuals. Now we must reintroduce the GARCH effects that we found, to find the value at risk.

modelDAX_omega = modelDAX@fit$par[2]
modelDAX_alpha1 = modelDAX@fit$par[3]
modelDAX_lastres = modelDAX@residuals[990]
modelDAX_lastres2 = modelDAX_lastres^2
modelDAX_beta1 = modelDAX@fit$par[4]
modelDAX_lastsigma = modelDAX@sigma.t[990]
modelDAX_lastsigma2 = modelDAX_lastsigma^2

DAXsim_sigma991sq = modelDAX_omega +modelDAX_alpha1*modelDAX_lastres2 + modelDAX_beta1*modelDAX_lastsigma2
DAXsim_sigma991 = sqrt(DAXsim_sigma991sq)

modelDAX_mu = modelDAX@fit$par[1]
DAXLR_sim = modelDAX_mu + DAXsim_sigma991*resDAX_sim

modelSP_omega = modelSP@fit$par[2]
modelSP_alpha1 = modelSP@fit$par[3]
modelSP_lastres = modelSP@residuals[990]
modelSP_lastres2 = modelSP_lastres^2
modelSP_beta1 = modelSP@fit$par[4]
modelSP_lastsigma = modelSP@sigma.t[990]
modelSP_lastsigma2 = modelSP_lastsigma^2
modelSP_beta2 = modelSP@fit$par[5]
modelSP_seclastsigma = modelSP@sigma.t[989]
modelSP_seclastsigma2 = modelSP_seclastsigma^2

SPsim_sigma991sq = modelSP_omega + modelSP_alpha1*modelSP_lastres2 + modelSP_beta1*modelSP_lastsigma2 + modelSP_beta2*modelSP_seclastsigma2
SPsim_sigma991 = sqrt(SPsim_sigma991sq)

modelSP_mu = modelSP@fit$par[1]
SPLR_sim = modelSP_mu + SPsim_sigma991*resSP_sim

# Now, we have found simulated values for the S&P500 and DAX log returns

# All that is left to do is to find the VaR

plot(DAXLR_sim, SPLR_sim)
port_sim = 1/2*DAXLR_sim + 1/2*SPLR_sim
var_sim_95 = quantile(port_sim, 0.05)*-1
var_sim_99 = quantile(port_sim, 0.01)*-1

var_sim_99
var_sim_95

# Notice that these are higher than the VaR that we found using the parametric approach.
# This is because the parametric approach does not account for lower tail dependence, but the t-copula that we used does

