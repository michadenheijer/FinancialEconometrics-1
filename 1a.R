rm(list=ls())    
# Jochem installeer eerst deze
# install.packages("here")
setwd(here::here())
library("tseries") # For Jarque Bera

dax = read.table("DAX.txt",header=FALSE)
dax_log_returns = diff(log(dax[,1]) * 100)


# 1

plot(dax[,1], type="l", xlab="t", ylab="Dax Index Value")
plot(dax_log_returns, type="l", xlab="t", ylab="Dax Index Single Day Return (%)")
dax_log_squared = dax_log_returns**2
plot(dax_log_squared, type="l")
plot(acf(dax_log_returns,main="")[1:35])
plot(acf(dax_log_squared,main="")[1:35])

# 3
# This file contains the average log-likelihood function. 
# The function below takes the data, labeled "x", and the parameter vector, 
# labeled "par", as input and gives the average log-likelihood, labeled "llik", 
# as output. 

llik_fun_GARCH <- function(par,x, plot_filter=FALSE){
  
  n <- length(x)
  
  #set paramter values from the input par using link functions for restrictions
  omega <- exp(par[1])                   #exp() to ensure omega>0
  alpha_1 <- exp(par[2])/(1+exp(par[2]))   #logistic()=exp()/(exp()) for 0<alpha<1
  alpha_2 <- exp(par[3])/(1+exp(par[3]))   #logistic()=exp()/(exp()) for 0<alpha<1
  beta_1 <- exp(par[4])/(1+exp(par[4]))    #logistic()=exp()/(exp()) for 0<beta<1
  beta_2 <- exp(par[5])/(1+exp(par[5]))    #logistic()=exp()/(exp()) for 0<beta<1
  
  ## Filter Volatility
  sig2 <- rep(0,n)
  sig2[1] <- var(x) #initialize volatility at unconditional variance
  sig2[2] <- var(x)
  
  for(t in 3:n){
    sig2[t] <- omega + alpha_1*x[t-1]^2  + alpha_2*x[t-2]^2 + beta_1*sig2[t-1] + beta_2*sig2[t-2] 
  }
  
  if (plot_filter){
    plot(sig2, type="l", ylab="cond. variance", xlab="t")
    
    u <- x/sqrt(sig2)
    
    # homosked. test
    plot(acf(u^2, main="")[1:35])
    
    # norm. test
    print(jarque.bera.test(u))
  }
  
  ## Calculate Log Likelihood Values
  
  #construct sequence of log-likelihood contributions
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  
  llik <- mean(l)  # obtain the average log-likelihood
  
  return(llik) # return the average log-likelihood as output
}

a <- 0.3
b <- 0.3  

omega <- var(dax_log_returns)*(1-a-b) 
par_ini <- c(log(omega),log(a/(1-a)), log(a/(1-a)), log(b/(1-b)), log(b/(1-b)))
est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,dax_log_returns), method = "BFGS")

omega_hat <- exp(est$par[1])
alpha_1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
alpha_2_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta_1_hat <- exp(est$par[4])/(1+exp(est$par[4]))
beta_2_hat <- exp(est$par[5])/(1+exp(est$par[5]))

llik_fun_GARCH(est$par, dax_log_returns, plot_filter = TRUE)

print(paste("Omega_hat", omega_hat))
print(paste("alpha_1_hat", alpha_1_hat))
print(paste("alpha_2_hat", alpha_2_hat))
print(paste("beta_1_hat", beta_1_hat))
print(paste("beta_2_hat", beta_2_hat))

########### 4
# Inlcuded in prev.

########## 5
# Ja zeer naar taakje weer want nu mag ik die hele functie weer herschrijven

llik_fun_GARCH_p_q <- function(p, q, par, x, plot_filter=FALSE, print_aic_bic=FALSE){
  
  n <- length(x)
  
  #set paramter values from the input par using link functions for restrictions
  omega <- exp(par[1])                   #exp() to ensure omega>0
  alpha = exp(par[2:(2+p-1)])/(1+exp(par[2:(2+p-1)]))
  beta = exp(par[(2+p):(2+p+q-1)])/(1+exp(par[(2+p):(2+p+q-1)]))
  
  ## Filter Volatility
  sig2 <- rep(var(x),n) #initialize volatility at unconditional variance
  
  for(t in (max(c(p,q)) + 1):n){
    sig2[t] <- omega
    
    # Is deze code niet cool gevectorized, ja
    # Heb ik zin om dit nu te gaan oplossen in R, nee
    for(i in 1:p){
      sig2[t] <- sig2[t] + alpha[i]*x[t-i]^2
    }
    
    for(j in 1:q){
      sig2[t] <- sig2[t] + beta[j]*sig2[t-j]
    }
  }
  
  if (plot_filter){
    plot(sig2, type="l", ylab="cond. variance", xlab="t")
  }
  
  ## Calculate Log Likelihood Values
  
  #construct sequence of log-likelihood contributions
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  
  llik <- mean(l)  # obtain the average log-likelihood
  
  if (print_aic_bic){
    llik_val <- -llik*n # Deze code komt van week 2
    
    # Echt heel dom maar de formule stond natuurlijk
    # verkeerd in de code van de opdracht
    k <- p + q + 1
    aic <- 2*k-2*llik_val
    bic <- log(n)*k-2*llik_val
    
    cat("The AIC is:")
    print(aic)
    
    cat("The BIC is:")
    print(bic)
  }
  
  return(llik) # return the average log-likelihood as output
}

est_garch <- function(p, q, x){
  a <- 0.2
  b <- 0.2  
  
  par_ini <- rep(0, 1 + p + q)
  par_ini[1] <- var(dax_log_returns)*(1-a-b)
  par_ini[2:(2+p-1)] <- log(a/(1-a))
  par_ini[(2+p):(2+p+q-1)] <- log(b/(1-b))

  est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_p_q(p, q, par, dax_log_returns), method = "BFGS")
  omega_hat <- exp(est$par[1])
  alpha_hat <- exp(est$par[2:(2+p-1)])/(1+exp(est$par[2:(2+p-1)]))
  beta_hat <- exp(est$par[(2+p):(2+p+q-1)])/(1+exp(est$par[(2+p):(2+p+q-1)]))
  
  print(paste("p:", p))
  print(paste("q:", q))
  print("omega_hat")
  print(omega_hat)
  print("alpha_hat")
  print(alpha_hat)
  print("beta_hat")
  print(beta_hat)
  
  llik_fun_GARCH_p_q(p, q, est$par, x, print_aic_bic = TRUE)
}

# Het wordt lekker naar geprint, dus succes daarmee
for (i in 1:4){
  for (j in 1:4) {
    est_garch(i, j, dax_log_returns)
  }
}


###### 7
llik_fun_GARCH_mod <- function(par,x){
  
  n <- length(x)
  
  #set paramter values from the input par using link functions for restrictions
  omega <- exp(par[1])                   #exp() to ensure omega>0
  alpha <- exp(par[2])/(1+exp(par[2]))   #logistic()=exp()/(exp()) for 0<alpha<1
  beta <- exp(par[3])/(1+exp(par[3]))    #logistic()=exp()/(exp()) for 0<beta<1
  gamma <- par[4]                        # Natuurlijk staat er niet gegeven wat gamma kan zijn
                                         # dus ja zal wel alles zijn.
  
  ## Filter Volatility
  sig2 <- rep(0,n)
  sig2[1] <- var(x) #initialize volatility at unconditional variance
  
  for(t in 2:n){
    sig2[t] <- omega + alpha*x[t-1]^2 + beta*sig2[t-1] + gamma * (x[t-1] < 0) * x[t-1]^2
  }
  
  ## Calculate Log Likelihood Values
  
  #construct sequence of log-likelihood contributions
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  
  llik <- mean(l)  # obtain the average log-likelihood
  
  return(llik) # return the average log-likelihood as output
}

par_ini <- c(log(omega),log(a/(1-a)), log(b/(1-b)), 0)
est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_mod(par,dax_log_returns), method = "BFGS")


# Nu moeten we nog vergelijken, dit kan natuurlijk met de AIC
# maar ja hebben we daar zin in? Ja geen idee, dit is dus lekker een probleem voor jou
# maar ik heb het wel geprogrammeerd, zodat je kan kiezen.
llik_fun_GARCH_p_q(1, 1, c(log(omega), log(a/(1-a)), log(b/(1-b))), dax_log_returns, print_aic_bic = TRUE)

n <- length(dax_log_returns)
llik_val <- -est$value*n # Deze code komt van week 2

# Echt heel dom maar de formule stond natuurlijk
# verkeerd in de code van de opdracht
# Volgens mij zit hier een rekenfout, maar ja kan wel
k <- 4
print(paste("AIC", 2*k-2*llik_val)) # AIC
print(paste("BIC", log(n)*k-2*llik_val)) # BIC

# Het kan ook zijn dat je natuurlijk gewoon de likelyhood wilt printen
# Die kan je gewoon vinden in deze code natuurlijk

# Hier zijn de ge-estimate waarden:
omega <- exp(est$par[1])                   #exp() to ensure omega>0
alpha <- exp(est$par[2])/(1+exp(est$par[2]))   #logistic()=exp()/(exp()) for 0<alpha<1
beta <- exp(est$par[3])/(1+exp(est$par[3]))   #logistic()=exp()/(exp()) for 0<alpha<1
gamma <- est$par[4]

print(paste("omega", omega))
print(paste("alpha", alpha))
print(paste("beta", beta))
print(paste("gamma", gamma))
