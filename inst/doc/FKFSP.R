## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
##The packages 'FKF', 'stats' and 'NFCP' are required for this Vignette:
library(FKF.SP)
library(FKF)
library(stats)
library(NFCP)

## -----------------------------------------------------------------------------
# Set constants:
## Length of series
n <- 10000

## AR parameters
AR <- c(ar1 = 0.6, ar2 = 0.2, ma1 = -0.2, sigma = sqrt(0.2))

# Generate observations:
set.seed(1)
a <- stats::arima.sim(model = list(ar = AR[c("ar1", "ar2")], ma = AR["ma1"]), n = n,
            innov = rnorm(n) * AR["sigma"])

## -----------------------------------------------------------------------------
arma21ss <- function(ar1, ar2, ma1, sigma) {
Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
Zt <- matrix(c(1, 0), ncol = 2)
ct <- matrix(0)
dt <- matrix(0, nrow = 2)
GGt <- matrix(0)
H <- matrix(c(1, ma1), nrow = 2) * sigma
HHt <- H %*% t(H)
a0 <- c(0, 0)
## Diffuse assumption
P0 <- matrix(1e6, nrow = 2, ncol = 2)
return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
            HHt = HHt))}

## -----------------------------------------------------------------------------
# The objective function passed to 'optim'
objective <- function(theta, yt, SP) {
param <- arma21ss(theta["ar1"], theta["ar2"], theta["ma1"], theta["sigma"])
# Kalman filtering through the 'fkf.SP' function:
if(SP){
 ans <- - fkf.SP(a0 = param$a0, P0 = param$P0, dt = param$dt, ct = param$ct, 
               Tt = param$Tt, Zt = param$Zt, HHt = param$HHt, GGt = param$GGt, 
               yt = yt)
 }
# Kalman filtering through the 'fkf' function:
 else{
 ans <- - fkf(a0 = param$a0, P0 = param$P0, dt = param$dt, ct = param$ct, Tt = param$Tt,
            Zt = param$Zt, HHt = param$HHt, GGt = param$GGt, yt = yt)$logLik
   
 }
 return(ans)
}
##Optim minimizes functions by default, so the negative is returned


## -----------------------------------------------------------------------------
#This test estimates parameters through 'optim'.
#Please run the complete chunk for a fair comparison:

#Initial values:
theta <- c(ar = c(0, 0), ma1 = 0, sigma = 1)

###MLE through the 'fkf' function:
start <- Sys.time()
set.seed(1)
FKF_estimation <- optim(theta, objective, yt = rbind(a), hessian = TRUE, SP = F)
FKF_runtime <- Sys.time() - start

###MLE through the 'fkf.SP' function:
start <- Sys.time()
set.seed(1)
FKF.SP_estimation <- optim(theta, objective, yt = rbind(a), hessian = TRUE, SP = T)
FKF.SP_runtime <- Sys.time() - start


## -----------------------------------------------------------------------------
print(rbind(FKF.SP = FKF.SP_estimation$par, FKF = FKF_estimation$par))

## -----------------------------------------------------------------------------
print(c(FKF.SP = FKF.SP_estimation$counts[1], FKF = FKF_estimation$counts[1]))

## -----------------------------------------------------------------------------
print(c(FKF.SP = FKF.SP_runtime, FKF = FKF_runtime))

## -----------------------------------------------------------------------------
## Transition equation:
## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)
## Measurement equation:
## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)

##Complete Nile Data - no NA's
y_complete <- y_incomplete <- Nile
##Incomplete Nile Data - two NA's are present:
y_incomplete[c(3, 10)] <- NA


## Set constant parameters:
dt <- ct <- matrix(0)
Zt <- Tt <- matrix(1)
a0 <- y_incomplete[1]   # Estimation of the first year flow
P0 <- matrix(100)     # Variance of 'a0'

## Parameter estimation - maximum likelihood estimation:
Nile_MLE <- function(yt, SP){
##Unknown parameters initial estimates:
GGt <- HHt <- var(yt, na.rm = TRUE) * .5
set.seed(1)
# Kalman filtering through the 'fkf.SP' function:
if(SP){
  return(optim(c(HHt = HHt, GGt = GGt),
        fn = function(par, ...)
             -fkf.SP(HHt = matrix(par[1]), GGt = matrix(par[2]), ...),
             yt = rbind(yt), a0 = a0, P0 = P0, dt = dt, ct = ct,
             Zt = Zt, Tt = Tt))
} else {
# Kalman filtering through the 'fkf' function:
  return(optim(c(HHt = HHt, GGt = GGt),
        fn = function(par, ...)
             -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
             yt = rbind(yt), a0 = a0, P0 = P0, dt = dt, ct = ct,
             Zt = Zt, Tt = Tt))
}}

## -----------------------------------------------------------------------------
fkf.SP_MLE_complete <- Nile_MLE(y_complete, SP = T)
fkf_MLE_complete <- Nile_MLE(y_complete, SP = F)

## -----------------------------------------------------------------------------
print(fkf.SP_MLE_complete[1:3])

## -----------------------------------------------------------------------------
print(fkf_MLE_complete[1:3])

## -----------------------------------------------------------------------------
fkf.SP_MLE_incomplete <- Nile_MLE(y_incomplete, SP = T)
fkf_MLE_incomplete <- Nile_MLE(y_incomplete, SP = F)

## -----------------------------------------------------------------------------
print(fkf.SP_MLE_incomplete[1:3])

## -----------------------------------------------------------------------------
print(fkf_MLE_incomplete[1:3])

## -----------------------------------------------------------------------------
#Number of NA values:
NA_values <- length(which(is.na(y_incomplete)))

print( 0.5 * NA_values * log(2 * pi))

## -----------------------------------------------------------------------------
#This test uses estimated parameters of complete data. 
#Please run the complete chunk for a fair comparison:

#'fkf' 
set.seed(1)
start <- Sys.time()
for(i in 1:1e4) fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fkf_MLE_complete$par[1]),
                    GGt = matrix(fkf_MLE_complete$par[2]), yt = rbind(y_complete))
FKF_runtime <- Sys.time() - start

#'fkf.SP'
set.seed(1)
start = Sys.time()
for(i in 1:1e4) fkf.SP(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fkf.SP_MLE_complete$par[1]),
                       GGt = matrix(fkf.SP_MLE_complete$par[2]), yt = rbind(y_complete))
fkf.SP_runtime <- Sys.time() - start

print(c(FKF.SP = fkf.SP_runtime, FKF = FKF_runtime))


## -----------------------------------------------------------------------------
#This test estimates parameters 10 times through 'optim'.
#Please run the complete chunk for a fair comparison:

## Transition equation:
## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)
## Measurement equation:
## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)

## tree-ring widths in dimensionless units
y <- treering

## Set constant parameters:
dt <- ct <- matrix(0)
Zt <- Tt <- matrix(1)
a0 <- y[1]            # Estimation of the first width
P0 <- matrix(100)     # Variance of 'a0'

##Time comparison - Estimate parameters 10 times:

###MLE through the 'fkf' function:
start = Sys.time()
set.seed(1)
for(i in 1:10)  fit_fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                     GGt = var(y, na.rm = TRUE) * .5),
                   fn = function(par, ...)
                     -fkf(HHt = array(par[1],c(1,1,1)), GGt = array(par[2],c(1,1,1)), ...)$logLik,
                   yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                   Zt = Zt, Tt = Tt)

run_time_FKF = Sys.time() - start

###MLE through the 'fkf.SP' function:
start = Sys.time()
set.seed(1)
for(i in 1:10)  fit_fkf.SP <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                        GGt = var(y, na.rm = TRUE) * .5),
                      fn = function(par, ...)
                        -fkf.SP(HHt = array(par[1],c(1,1,1)), GGt = matrix(par[2]), ...),
                      yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                      Zt = Zt, Tt = Tt)
run_time_FKF.SP = Sys.time() - start

print(c(fkf.SP = run_time_FKF.SP, fkf = run_time_FKF))

## Filter tree ring data with estimated parameters using 'fkf':
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = array(fit_fkf$par[1],c(1,1,1)),
               GGt = array(fit_fkf$par[2],c(1,1,1)), yt = rbind(y))


## -----------------------------------------------------------------------------

yt = t(log(NFCP::SS_oil$contracts)) # quoted log futures prices
delta_t <- NFCP::SS_oil$dt # Discrete time step
##time to maturity of quoted futures contracts:
TTM <- t(NFCP::SS_oil$contract_maturities)

a0 <- yt[1,1]     # initial estimate
P0 <- matrix(100) # Variance of 'a0'

## GBM Function
gbm_mle <- function(theta, SP){

ct <- theta["alpha_rn"] * TTM
dt <- (theta["alpha"] - 0.5 * theta["sigma"]^2) * delta_t
Zt <- matrix(1, nrow(yt))
HHt <- matrix(theta["sigma"]^2 * delta_t)
Tt <- matrix(1)

##'fkf.SP' requires a vector of the diagonal elements of the variances of the measurement error 
if(SP){
GGt = rep(theta["ME_1"]^2, nrow(yt))
} else {
##'fkf' instead requires a matrix of the elements of the variances of the measurement error 
GGt = diag(theta["ME_1"]^2, nrow(yt))
}

##'fkf.SP' returns only the log-likelihood numeric value, whilst 'fkf' returns a list of filtered values
logLik = ifelse(SP,
                - fkf.SP(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = yt),
                - fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = yt)$logLik
                )
return(logLik)
}


## -----------------------------------------------------------------------------
#This test estimates parameters through 'optim'.
#Please run the complete chunk for a fair comparison:

#Initial estimates
gbm_par <- c(alpha = 0, alpha_rn = 0.01, sigma = 0.1, ME_1 = 0.05)

###MLE through the 'fkf.SP' function:
set.seed(1)
start = Sys.time()
fkf.SP.gbm = optim(par = gbm_par, fn = gbm_mle, SP = T)
fkf.SP_runtime <- Sys.time() - start

###MLE through the 'fkf' function:
set.seed(1)
start = Sys.time()
fkf.gbm = optim(par = gbm_par, fn = gbm_mle, SP = F)
fkf_runtime <- Sys.time() - start


## -----------------------------------------------------------------------------
print(rbind(FKF.SP = - fkf.SP.gbm$value, FKF = - fkf.gbm$value))

## -----------------------------------------------------------------------------
print(rbind(FKF.SP = fkf.SP.gbm$par, FKF = fkf.gbm$par))

## -----------------------------------------------------------------------------
print(c(FKF.SP = fkf.SP.gbm$counts[1], FKF = fkf.gbm$counts[1]))

## -----------------------------------------------------------------------------
print(c(FKF.SP = fkf.SP_runtime, FKF = fkf_runtime))

