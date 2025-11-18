# Portfolio Optimization with Minimization of CVaR, where the returns distribution is fitted with MNTS distribution. 


library(foreach)
library(doParallel)
library(xtable)


library(rugarch)
library(temStaR)
library(quadprog)
library(quantmod)
library(functional)
library(nloptr)
library(pracma)
library(NlcOptim)
library(evmix)
library(spatstat)
library(Matrix)
library(ggplot2)

setwd("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization")

source("readStockData.R")
source("func_findOptport.R")




begindate <- as.Date("2020-01-01")
enddate<-as.Date("2025-01-01")
intervalDateNumber <- 22


# Dow Jones tickers
yTickers<-c("AAPL","AMGN","AXP","BA","CAT","CRM","CSCO","CVX","DIS","DOW","GS","HD","HON","IBM","JPM","KO","MCD","MMM","MRK","MSFT","NKE","PG","TRV","UNH","V","VZ","WBA","WMT","INTC","JNJ")
numofelements<-length(yTickers)
ret <- readdata(yTickers, numofelements, begindate, enddate, intervalDateNumber)


# Fitting the data to an mnts distribution 
options(warn=-1)
st <- fitmnts(returndata = ret, n = numofelements)
options(warn=0)
st

# Using the defined function by professor to find the optimal portfolio with mnts's distributed returns  

findOptPort_MNTS <- function(stMNTS, eta, benchmarkReturn){
  numofelements <- stMNTS$ndim
  w0 <- rep(1/numofelements, numofelements)
  optport <- fmincon(x0 = w0, 
                     fn = functional::Curry(portfolioCVaRmnts, 
                                            strPMNTS = stMNTS,
                                            eta = eta),
                     A = -matrix(stMNTS$mu, nrow = 1, ncol = numofelements), 
                     b = -benchmarkReturn,
                     Aeq = matrix(data = 1, nrow = 1, ncol = numofelements), 
                     beq = 1,
                     lb = rep(0, numofelements), 
                     ub = rep(1, numofelements))
  
  ExpectedReturn <- sum(optport$par*stMNTS$mu)
  returnvalue <- list(weight=optport$par, 
                      reward = ExpectedReturn, 
                      risk = optport$value)
  return( returnvalue )
}

# Benchmark return

mu0<-0.018

retv<-findOptPort_MNTS(stMNTS=st,eta=0.01,benchmarkReturn = mu0)
retv

# This is my single portfolio that is obtained through the portfolio optimization process

# We can plot the efficient frontier by generating portfolios on a list of target returns
target_returns <- seq(min(st$mu),max(st$mu),length.out=20)

reward_list <- numeric(length(target_returns))
risk_list   <- numeric(length(target_returns))

for (i in seq_along(target_returns)) {
  mu0 <- target_returns[i]
  
  # run optimization for this target return
  optport <- findOptPort_MNTS(
    stMNTS = st,
    eta = 0.01,
    benchmarkReturn = mu0
  )
  
  # save reward and risk
  reward_list[i] <- optport$reward
  risk_list[i]   <- optport$risk
}

# combine into data frame
frontier_df <- data.frame(
  Return = reward_list,
  CVaR   = risk_list
)

frontier_df

# Plotting the efficient frontier
ggplot(frontier_df, aes(x = CVaR, y = Return)) +
  geom_line(color = "red", size = 1.2) +
  geom_point(color = "blue") +
  labs(title = "CVaR Minimized-Efficient Frontier (MNTS)",
       x = "CVaR (tail risk)",
       y = "Expected Return") +
  theme_minimal()

# Doing the portfolio optimization for distribution that is not distributed as MNTS


er <- colMeans(ret)
covMtx <- cov(ret)

empiricalCVaR<-function(w,returns,eta=0.01){
  port_ret<-as.numeric(returns %*% w)
  q<-quantile(port_ret,probs=eta,type=7)
  cvar<-mean(port_ret[port_ret<=q])
  return(-cvar)
}

findOptPort_empirical <- function(returns, eta, benchmarkReturn){
  numofelements <- ncol(returns)
  mu <- colMeans(returns)   # historical mean returns
  w0 <- rep(1/numofelements, numofelements)
  
  optport <- fmincon(
    x0 = w0,
    fn = functional::Curry(empiricalCVaR, returns = returns, eta = eta),
    A  = -matrix(mu, nrow = 1, ncol = numofelements), 
    b  = -benchmarkReturn,
    Aeq = matrix(1, nrow = 1, ncol = numofelements), 
    beq = 1,
    lb = rep(0, numofelements), 
    ub = rep(1, numofelements)
  )
  
  ExpectedReturn <- sum(optport$par * mu)
  returnvalue <- list(
    weight = optport$par,
    reward = ExpectedReturn,
    risk   = optport$value
  )
  return(returnvalue)
}


target_returns_empherical <- seq(min(er),max(er),length.out=20)

reward_list_empherical <- numeric(length(target_returns_empherical))
risk_list_empherical   <- numeric(length(target_returns_empherical))

findOptPort_empirical(
  returns = ret,
  eta = 0.01,
  benchmarkReturn = target_returns_empherical[17]
)

bad_idx <- c(2, 6, 9, 15, 17, 18, 19, 20)

target_returns_clean <- target_returns_empherical[-bad_idx]

target_returns_clean


for (i in seq_along(target_returns_clean)) {
  mu0 <- target_returns_clean[i]
  
  # run optimization for this target return
  optport <- findOptPort_empirical(
    returns = ret,
    eta = 0.01,
    benchmarkReturn = mu0
  )
  
  # save reward and risk
  reward_list_empherical[i] <- optport$reward
  risk_list_empherical[i]   <- optport$risk
}

# combine into data frame
frontier_df_empherical <- data.frame(
  Return = reward_list_empherical,
  CVaR   = risk_list_empherical
)

frontier_df_empherical


ggplot(frontier_df_empherical, aes(x = CVaR, y = Return)) +
  geom_line(color = "red", size = 1.2) +
  geom_point(color = "blue") +
  labs(title = "CVaR Minimized-Efficient Frontier (NOT MNTS)",
       x = "CVaR (tail risk)",
       y = "Expected Return") +
  theme_minimal()
frontier_df_empherical[1:12,]

plot(frontier_df_empherical[1:12,]
$CVaR, frontier_df_empherical[1:12,]
$Return,
     type = "b", col = "blue", pch = 16,
     xlab = "CVaR", ylab = "Expected Return",
     main = "Efficient Frontiers (Empirical vs MNTS)")

lines(frontier_df$CVaR, frontier_df$Return,
      type = "b", col = "red", pch = 17)

legend("bottomright", legend = c("Empirical", "MNTS"),
       col = c("blue", "red"), pch = c(16, 17), lty = 1)

# Let's try a sharpe based portfolio optimization

rf_annual <- 0.0418
rf_monthly_log <- log(1 + rf_annual) / 12

portfolioSharpe <-function(w,mu,covMtx,rf=rf_monthly_log){
  port_ret<-sum(w*mu)
  port_vol<-sqrt(t(w) %*% covMtx %*% w)
  sharpe <- (port_ret - rf) / port_vol
  return(-sharpe)
}
optSharpe <- fmincon(
  x0 = rep(1/length(st$mu), length(st$mu)),
  fn = functional::Curry(portfolioSharpe, mu = st$mu, covMtx = st$CovMtx, rf = rf_monthly_log),
  Aeq = matrix(1, nrow = 1, ncol = length(st$mu)),
  beq = 1,
  lb = rep(0, length(st$mu)),
  ub = rep(1, length(st$mu))
)
weights <- optSharpe$par
weights[abs(weights) < 1e-6] <- 0
# Expected portfolio return
port_ret <- sum(weights * st$mu)

# Portfolio volatility
port_vol <- sqrt(t(weights) %*% st$CovMtx %*% weights)

# Sharpe ratio
sharpe <- (port_ret - rf_monthly_log) / port_vol

c(Return = port_ret, Volatility = port_vol, Sharpe = sharpe)


# Doing the ARMA GARCH NTS Fit

source("miscellaneous_tools.R")
source("fcts_mntsARMAGARCH_paramest.R")

yTickers <- c("^GSPC", "^IXIC", "^RUT") #,"^DJI"
numofelements <- length(yTickers)

begindate <- as.Date("2019-1-1")
enddate <- as.Date("2021-12-31")

ret <- readdata(yTickers, numofelements, begindate, enddate, 1 )

# Modifying the function before using it

paramest_mntsarmagarch <- function(ret, numofelements, yTickers=NULL){
  muvec <- matrix(nrow = numofelements, ncol = 1)
  sigvec <- matrix(nrow = numofelements, ncol = 1)
  stdret <- matrix(nrow = dim(ret)[1], ncol = numofelements)
  tgarchparam <- matrix(nrow = numofelements, ncol = 7)
  rownames(tgarchparam) <- yTickers
  
  tspec <- rugarch::ugarchspec(
    mean.model = list(armaOrder = c(1,1)),
    variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
    distribution.model = "std"
  )
  
  for(n in 1:numofelements) {
    tf <- rugarch::ugarchfit(spec = tspec, data = ret[,n])
    a <- coef(tf)
    
    if (n == 1) {
      colnames(tgarchparam) <- names(a)   # set column names from first asset
    }
    
    tgarchparam[n,] <- a
    v1 <- tf@fit$sigma
    stdret[,n] <- tf@fit$residuals / v1   # standardized residuals
    tforc <- rugarch::ugarchforecast(tf, n.ahead=1)
    muvec[n] <- tforc@forecast$seriesFor[1]
    sigvec[n] <- tforc@forecast$sigmaFor[1]
  }
  
  options(warn=-1)
  st <- temStaR::fitmnts(returndata = stdret, n = numofelements, stdflag = TRUE)
  st$mu <- as.numeric(muvec)
  st$sigma <- as.numeric(sigvec)
  st$CovMtx <- diag(as.numeric(sigvec)) %*% st$CovMtx %*% diag(as.numeric(sigvec))
  options(warn=0)  
  
  colnames(st$Rho) <- yTickers
  rownames(st$Rho) <- yTickers
  colnames(st$CovMtx) <- yTickers
  rownames(st$CovMtx) <- yTickers
  names(st$mu) <- yTickers
  names(st$sigma) <- yTickers
  names(st$beta) <- yTickers
  
  return(list(st = st, tagparam = tgarchparam))
}


marketparam <- paramest_mntsarmagarch(ret,numofelements,yTickers)

stMNTS <- indStdNTS2StdMNTS(marketparam$st$alpha,
                            marketparam$st$theta, 
                            marketparam$st$beta,
                            marketparam$st$mu,
                            marketparam$st$sigma,
                            marketparam$st$CovMtx, 
                            yTickers = yTickers)

stMNTS
