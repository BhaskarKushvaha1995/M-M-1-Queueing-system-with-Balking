# Load necessary libraries
if (!require("MASS")) install.packages("MASS", dependencies=TRUE)
if (!require("zipfR")) install.packages("zipfR", dependencies=TRUE)
if (!require("hypergeo")) install.packages("hypergeo", dependencies=TRUE)
if (!require("gsl")) install.packages("gsl", dependencies=TRUE)
if (!require("pracma")) install.packages("pracma", dependencies=TRUE)
if (!require("expint")) install.packages("expint", dependencies=TRUE)
if (!require("cascsim")) install.packages("cascsim", dependencies=TRUE)
if (!require("mltools")) install.packages("mltools", dependencies=TRUE)

# Attach the libraries
library(MASS)
library(zipfR)
library(hypergeo)
library(gsl)
library(pracma)
library(expint)
library(cascsim)
library(mltools)

# Define the simulation function
gerRo <- function(tau, rho) {
  return(rpois(tau, rho))
}

# Maximum likelihood estimate for rho
MLERho<-function(samp) {
n<-length(samp)
y<-sum(samp)
return(y/n)} 
tau <- c(10, 20, 50, 100, 200)
rho <- c(0.20, 0.5, 0.8, 5)
r.MLERho<-TabMtCa(tau,rho,MLERho)
r.MLERho
#write.csv(r.MLERho, "MMLERho.csv")


# Function to calculate mean and mean squared error (MSE)
MtCa <- function(tau, rho, fEst, ...) {
  rep <- 10000
  est <- numeric(rep)
  samp <- numeric(tau)
  for (i in 1:rep) {
    samp <- gerRo(tau, rho)
    est[i] <- fEst(samp, ...)
  }
  return(c(mean(est), mean((rho - est)^2)))  # Calculate MSE
}

# Function to create a table of results
TabMtCa <- function(tau, rho, fEst, ...) {
  res <- matrix(nrow = length(rho), ncol = 2 * length(tau))
  for (i in 1:length(rho)) {
    for (j in 1:length(tau)) {
      res[i, (2 * j - 1):(2 * j)] <- MtCa(tau[j], rho[i], fEst, ...)
    }
  }
  return(res)
}

# Bayesian estimate rho ~ Gamma
GammaRhoselfA <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  return((1 + a + y) / (b + n))  # Correct calculation
}

# Define parameters
tau <- c(10, 20, 50, 100, 200)
rho <- c(0.20,0.5,0.8, 5)
a <- 0.4
b <- 0.4

# Calculate results
r.GammaRhoselfA <- TabMtCa(tau, rho, GammaRhoselfA, a, b)
r.GammaRhoselfA
#write.csv(r.GammaRhoselfA, "GammaRhoselfA.csv")

GammaRhoselfB <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  return((1 + a + y) / (b + n))  # Correct calculation
}

# Define parameters
a <- 0.4
b <- 1

# Calculate results
r.GammaRhoselfB <- TabMtCa(tau, rho, GammaRhoselfB, a, b)
r.GammaRhoselfB
#write.csv(r.GammaRhoselfB, "GammaRhoselfB.csv")

GammaRhoselfC <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  return((1 + a + y) / (b + n))  # Correct calculation
}

# Define parameters
a <- 1
b <- 0.4

# Calculate results
r.GammaRhoselfC <- TabMtCa(tau, rho, GammaRhoselfC, a, b)
r.GammaRhoselfC 
#write.csv(r.GammaRhoselfC, "GammaRhoselfC.csv")

GammaRhoselfD <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  return((1 + a + y) / (b + n))  # Correct calculation
}

# Define parameters
a <- 1
b <- 1

# Calculate results
r.GammaRhoselfD <- TabMtCa(tau, rho, GammaRhoselfD, a, b)
r.GammaRhoselfD
#write.csv(r.GammaRhoselfD, "GammaRhoselfD.csv")

GammaRHOPRECA <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  
  # Use lgamma for more numerically stable computation
  numerator <- lgamma(2 + a + y)
  denominator <- 2 * log(b + n) + lgamma(a + y)
  
  # Calculate and handle potential numerical issues
  value <- numerator - denominator
  if (value < -700) return(0)  # Avoid exp() underflow
  
  return(sqrt(exp(value)))
}

GammaRHOPRECB <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  
  numerator <- lgamma(2 + a + y)
  denominator <- 2 * log(b + n) + lgamma(a + y)
  
  value <- numerator - denominator
  if (value < -700) return(0)  # Avoid exp() underflow
  
  return(sqrt(exp(value)))
}

GammaRHOPRECC <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  
  numerator <- lgamma(2 + a + y)
  denominator <- 2 * log(b + n) + lgamma(a + y)
  
  value <- numerator - denominator
  if (value < -700) return(0)  # Avoid exp() underflow
  
  return(sqrt(exp(value)))
}

GammaRHOPRECD <- function(samp, a, b) {
  n <- length(samp)
  y <- sum(samp)
  if (y == 0) return(0)  # Handle zero-sum case
  
  numerator <- lgamma(2 + a + y)
  denominator <- 2 * log(b + n) + lgamma(a + y)
  
  value <- numerator - denominator
  if (value < -700) return(0)  # Avoid exp() underflow
  
  return(sqrt(exp(value)))
}
write.csv(r.GammaRHOPRECD, "GammaRHOPRECD.csv")
# Calculate results with corrected functions
r.GammaRHOPRECA <- TabMtCa(tau, rho, GammaRHOPRECA, a, b)
r.GammaRHOPRECA
#write.csv(r.GammaRHOPRECA, "GammaRHOPRECA.csv")

r.GammaRHOPRECB <- TabMtCa(tau, rho, GammaRHOPRECB, a, b)
r.GammaRHOPRECB
#write.csv(r.GammaRHOPRECB, "GammaRHOPRECB.csv")

r.GammaRHOPRECC <- TabMtCa(tau, rho, GammaRHOPRECC, a, b)
r.GammaRHOPRECC
#write.csv(r.GammaRHOPRECC, "GammaRHOPRECC.csv")

r.GammaRHOPRECD <- TabMtCa(tau, rho, GammaRHOPRECD, a, b)
r.GammaRHOPRECD
#write.csv(r.GammaRHOPRECD, "GammaRHOPRECD.csv")






# Corrected BetaRhoself Functions with internal calculation of y and n
BetaRhoselfA <- function(samp, a, b) { 
  y <- sum(samp)
  n <- length(samp)
  num <- (a + y) * hyperg_U(a + b, b - y, n * log(exp(1))) 
  denom <- n * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1)) 
  result <- num / denom 
  return(result) 
}

BetaRhoselfB <- function(samp, a, b) {
  y <- sum(samp)
  n <- length(samp)
  num <- (a + y) * hyperg_U(a + b, b - y, n * log(exp(1))) 
  denom <- n * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1)) 
  result <- num / denom 
  return(result) 
}

BetaRhoselfC <- function(samp, a, b) {
  y <- sum(samp)
  n <- length(samp)
  num <- (a + y) * hyperg_U(a + b, b - y, n * log(exp(1))) 
  denom <- n * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1)) 
  result <- num / denom 
  return(result) 
}

BetaRhoselfD <- function(samp, a, b) { 
  y <- sum(samp)
  n <- length(samp)
  num <- (a + y) * hyperg_U(a + b, b - y, n * log(exp(1))) 
  denom <- n * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1)) 
  result <- num / denom 
  return(result) 
}

# Generate and save results for BetaRhoself functions
tau <- c(10, 20, 50, 100, 200) 
rho <- c(0.20, 0.5, 0.8, 5) 

# Running simulations and saving results
results_list_self <- list()

# Example usage with revised BetaRhoself functions
results_A <- TabMtCa(tau, rho, BetaRhoselfA, a = 2, b = 2)
results_B <- TabMtCa(tau, rho, BetaRhoselfB, a = 2, b = 1)
results_C <- TabMtCa(tau, rho, BetaRhoselfC, a = 1, b = 2)
results_D <- TabMtCa(tau, rho, BetaRhoselfD, a = 1, b = 1)

# Add results to the list
results_list_self$BetaRhoselfA <- results_A
results_list_self$BetaRhoselfB <- results_B
results_list_self$BetaRhoselfC <- results_C
results_list_self$BetaRhoselfD <- results_D

# Print the results
print(results_list_self)


# Corrected BetaRhoself Functions with internal calculation of y and n
BetaRhoPLFA <- function(samp, a, b) { 
  y <- sum(samp)
  n <- length(samp)
  numerator <- (a + y) * (1 + a + y) * hyperg_U(a + b, -1 + b - y, n * log(exp(1))) 
denominator <- n^2 * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1))^2 
result <- sqrt(numerator / denominator) 
 return(result) 
}

BetaRhoPLFB <- function(samp, a, b) {
  y <- sum(samp)
  n <- length(samp)
   y <- sum(samp)
  n <- length(samp)
  numerator <- (a + y) * (1 + a + y) * hyperg_U(a + b, -1 + b - y, n * log(exp(1))) 
denominator <- n^2 * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1))^2 
result <- sqrt(numerator / denominator) 
  return(result) 
}

BetaRhoPLFC <- function(samp, a, b) {
   y <- sum(samp)
  n <- length(samp)
  numerator <- (a + y) * (1 + a + y) * hyperg_U(a + b, -1 + b - y, n * log(exp(1))) 
denominator <- n^2 * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1))^2 
result <- sqrt(numerator / denominator)
  return(result) 
}

BetaRhoPLFD <- function(samp, a, b) { 
   y <- sum(samp)
  n <- length(samp)
  numerator <- (a + y) * (1 + a + y) * hyperg_U(a + b, -1 + b - y, n * log(exp(1))) 
denominator <- n^2 * hyperg_U(a + b, 1 + b - y, n * log(exp(1))) * log(exp(1))^2 
result <- sqrt(numerator / denominator)
  return(result) 
}

# Generate and save results for BetaRhoself functions
tau <- c(10, 20, 50, 100, 200) 
rho <- c(0.20, 0.5, 0.8, 5) 

# Running simulations and saving results
results_list_PLF <- list()

# Example usage with revised BetaRhoself functions
results_A <- TabMtCa(tau, rho, BetaRhoPLFA, a = 2, b = 2)
results_B <- TabMtCa(tau, rho, BetaRhoPLFB, a = 2, b = 1)
results_C <- TabMtCa(tau, rho, BetaRhoPLFC, a = 1, b = 2)
results_D <- TabMtCa(tau, rho, BetaRhoPLFD, a = 1, b = 1)

# Add results to the list
results_list_PLF$BetaRhoPLFA <- results_A
results_list_PLF$BetaRhoPLFB <- results_B
results_list_PLF$BetaRhoPLFC <- results_C
results_list_PLF$BetaRhoPLFD <- results_D

# Print the results
print(results_list_PLF)


# Revised Jeffreys prior estimate for rho
JefRhoself <- function(samp) {
  n <- length(samp)
  y <- sum(samp)
  
  # Logarithmic form for stability
  log_numerator <- lgamma(y + 3 / 2)
  log_denominator <- log(n) + lgamma(y + 1 / 2)
  
  # Calculate the ratio in the exponential scale for stability
  est <- exp(log_numerator - log_denominator)
  
  if (is.nan(est) || is.infinite(est)) {
    return(NA)
  }
  
  return(est)
}

# Revised Jeffreys prior precision estimate for rho
JefRhoprec <- function(samp) {
  n <- length(samp)
  y <- sum(samp)
  
  # Logarithmic form for stability
  log_numerator <- lgamma(y + 5 / 2)
  log_denominator <- 2 * log(n) + lgamma(y + 1 / 2)
  
  # Calculate the ratio in the exponential scale for stability
  est <- exp(0.5 * (log_numerator - log_denominator))
  
  if (is.nan(est) || is.infinite(est) || est < 0) {
    return(NA)
  }
  
  return(est)
}

# Define parameters
tau <- c(10, 20, 50, 100, 200)
rho <- c(0.20, 0.5, 0.8, 5)

# Run calculations with the improved functions
r.JefRhoself <- TabMtCa(tau, rho, JefRhoself)
print(r.JefRhoself)

r.JefRhoprec <- TabMtCa(tau, rho, JefRhoprec)
print(r.JefRhoprec)

# Optionally save to CSV
 #write.csv(r.JefRhoself, "JefRhoself.csv")
# write.csv(r.JefRhoprec, "JefRhoprec.csv")







# Load required libraries
library(gsl)       # For gamma functions
library(knitr)     # For displaying tables

# Define parameters for the prior distributions
a <- 1       # Parameter for Gamma and Beta second kind prior
b <- 1       # Parameter for Gamma prior
d <- 0.1     # Parameter for Beta second kind prior
m_values <- 0:4 # Values of m for which we compute the predictive distribution
rho_values <- c(0.2, 0.5, 0.8, 5) # Different rate parameters for Poisson distribution
sample_sizes <- c(10) # Different sample sizes

# Initialize a data frame to store the results
results <- data.frame(
  Rho = numeric(),
  Sample_Size = integer(),
  m = integer(),
  P1 = numeric(),
  P2 = numeric(),
  P3 = numeric(),
  B12 = numeric(),
  B13 = numeric(),
  B23 = numeric()
)

# Predictive distributions
# 1. Predictive distribution P1 for Beta second kind prior
predictive_beta2 <- function(m, a, b, y, n) {
numerator <- gamma(a + m + y) * hyperg_U(a + m + y, 1 - b + m + y, (1 + n) * log(exp(1))) 
denominator <- factorial(m) * gamma(a + y) * hyperg_U(a + y, 1 - b + y, n * log(exp(1))) 
return(numerator / denominator)
}

# 2. Predictive distribution P2 for Gamma prior
predictive_gamma <- function(m, a, b, y, n) {
numerator <- gamma(a + m + y) * ((b + n) * log(exp(1)))^(a + y) * ((1 + b + n) * log(exp(1)))^(-a - m - y) 
denominator <- factorial(m) * gamma(a + y) 
return(numerator / denominator)
}

# 3. Predictive distribution P3 for Jeffreys prior
predictive_jeffreys <- function(m, y, n) {
numerator <- gamma(1/2 + m + y) * (n * log(exp(1)))^(1/2 + y) * ((1 + n) * log(exp(1)))^(-(1/2) - m - y) 
denominator <- factorial(m) * gamma(1/2 + y) 
return(numerator / denominator)
}

# Loop over different values of rho and sample sizes
for (rho in rho_values) {
  for (n in sample_sizes) {
    # Simulate Poisson data
    sample_data <- rpois(n, rho)
    y <- sum(sample_data)     # Sufficient statistic for Poisson data (sum of counts)
    n<- length(sample_data)
# Loop over values of m
    for (m in m_values) {
      # Predictive distributions with checks
      P1 <- predictive_beta2(m, a, b, y, n)
      P2 <- predictive_gamma(m, a, b, y, n)
      P3 <- predictive_jeffreys(m, y, n)
      
      # Bayes Factors with NaN and Inf handling
      B12 <- ifelse(is.nan(P1 / P2) | is.infinite(P1 / P2), NA, P1 / P2)
      B13 <- ifelse(is.nan(P1 / P3) | is.infinite(P1 / P3), NA, P1 / P3)
      B23 <- ifelse(is.nan(P2 / P3) | is.infinite(P2 / P3), NA, P2 / P3)
      
      # Append the results to the data frame
      results <- rbind(results, data.frame(
        Rho = rho,              # Current value of rho
        Sample_Size = n,        # Current sample size
        m = m,                  # Current value of m
        P1 = P1,
        P2 = P2,
        P3 = P3,
        B12 = B12,
        B13 = B13,
        B23 = B23
      ))
    }
  }
}

# Display results for rho = 0.2 for different values of m
print(kable(subset(results, Rho == 0.2), caption = "Predictive Distributions and Bayes Factors for rho = 0.2"))

# Display results for rho = 0.5 for different values of m
print(kable(subset(results, Rho == 0.5), caption = "Predictive Distributions and Bayes Factors for rho = 0.5"))

# Display results for rho = 0.8 for different values of m
print(kable(subset(results, Rho == 0.8), caption = "Predictive Distributions and Bayes Factors for rho = 0.8"))

# Display results for rho = 5 for different values of m
print(kable(subset(results, Rho == 5), caption = "Predictive Distributions and Bayes Factors for rho = 5"))

# Save the results to a CSV file
#write.csv(results, file = "predictive_distributions_and_bayes_factors_all.csv", row.names = FALSE)

# Print a confirmation message
print("All results saved to 'predictive_distributions_and_bayes_factors_all.csv'")


