
################################################################################
# Antoniak (?) glm
################################################################################

#------------------------------------------------------------ 
# Part 1 - define the functions
#------------------------------------------------------------ 

# Internal function to simulate from a Dirichlet process
rDP <- function(alpha, size) {
  n_seq <- 0:(size - 1)
  sum(rbinom(size, 1, alpha / (alpha + n_seq)))
}
rDP <- Vectorize(FUN = rDP, vectorize.args = c("alpha", "size"))

# Useful function that is present in the glm.fit...
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Function to make the link for the AntoniakGLM
make.link.Antoniak <- function() {
  linkfun <- function(mu) {
    # Get the size parameter from the global environment
    get("size", envir = parent.env(environment()))
    g_inv <- function(mu_target, s){
      if(mu_target <= 1){
        sol <- 1e-8
      } else if (mu_target >= s) {
        sol <- 1e8
      } else {
        sol <- uniroot(function(x) x * (digamma(x + s) - digamma(x)) - mu_target, c(1e-12, 1e12))$root
      }
      log(sol)
    }
    g_inv <- Vectorize(g_inv, vectorize.args = c("mu_target", "s")) 
    g_inv(mu, size)
  }
  
  linkinv <- function(eta, size.) {
    exp_eta <- pmax(exp(eta), .Machine$double.eps)
    pmax(exp_eta * (digamma(exp_eta + size.) - digamma(exp_eta)), .Machine$double.eps)
  }
  
  mu.eta <- function(eta, size.) {
    exp_eta <- pmax(exp(eta), .Machine$double.eps)
    pmax(exp_eta * (digamma(exp_eta + size.) - digamma(exp_eta)) + 
           exp_eta^2 * (trigamma(exp_eta + size.) - trigamma(exp_eta)), .Machine$double.eps)
  }
  
  valideta <- function(eta, size.) TRUE
  
  environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                 valideta = valideta, name = "rarefaction"), class = "link-glm")
}


antoniak <- function() {
  
  family <- "Antoniak"
  # Make the accumulation curve link for the Antoniak glm
  stats <- make.link.Antoniak()
  linktemp <- "rarefaction"
  
  # Variance function
  variance <- function(mu, size.){
    exp_eta <- exp(stats$linkfun(mu, size.))
    mu + exp_eta^2 * (trigamma(exp_eta + size.) - trigamma(exp_eta))
  }
  
  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)
  
  # residual deviance
  dev.resids <- function(y, mu, wt, size.) {
    alpha_y <- exp(stats$linkfun(y, size.))
    alpha_mu <- exp(stats$linkfun(mu, size.))
    wt * (y * log(alpha_y/alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) + 
            lgamma(alpha_mu + size.) - lgamma(alpha_mu))
  }
  
  # aic
  aic <- function(y, n, mu, wt, dev, size.) {
    alpha <- exp(stats$linkfun(mu, size.))
    - 2 * sum(y * log(alpha) - lgamma(alpha + size.) + lgamma(alpha) * wt)
  } 
  
  # Initialization function
  initialize <- expression({
    if (ncol(y) != 2) stop("must specify cbind(n, y) as dependent variable for the 'antoniak' family")
    #n <- y[, 1]
    #y <- y[, 2]
    #if (any(y < 1 | y > n)) stop("values not in range for the 'antoniak' family")
    #n <- rep.int(1, nobs)
    name_size <- colnames(y)[1]
    name_y <- colnames(y)[2]
    size <- y[, 1]
    mustart <- y[, 2]
    y <- y[, 2]
  })
  
  # Simulation function
  simfun <- function(object, nsim, size.) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    alpha <- exp(stats$linkfun(ftd))
    sapply(1:nsim, function(i) rDP(alpha = alpha, size.))
  }
  
  # Return the structure for thr glm type
  structure(list(family = family, link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, 
                 dispersion = 1), class = "family")
}


#------------------------------------------------------------ 
# Part 2 - try it!
#------------------------------------------------------------ 

#------------------------------------------------------- Simulated data
# Run the experiment to replicate data
simulate_data <- function(N, p, lambda){
  X <- cbind(1, matrix(rnorm(N * (p - 1), sd = 10), nrow = N, ncol = p - 1))
  beta <- c(2, rnorm(p-1, sd = 0.1))
  alphas <- exp(X %*% beta)
  n <- rpois(N, lambda = lambda)
  y <- rDP(alpha = alphas, size = n)
  return(list(n = n, y = y, beta = beta, X = X, df = data.frame(n = n, y = y, X)))
}


# Try the function
set.seed(10)
data <- simulate_data(N = 1200, p = 10, lambda = 1000)

fit <- glm(cbind(n, y) ~ . -1,
           family = antoniak, 
           data = data$df[1:1000, ],      # <------ estimate using the first 1000
           method = "glm.fit")   # <------ Specify this fitting method
# Print
print(fit)

# Summary
summary(fit)
plot(fit$coefficients, data$beta)
plot(predict(fit, type = "response"), data$y[1:1000])

# plot
plot(fit, size = data$n) # <---- Still working on it! it is very boring. The quick solution is to have the parameter size when calling the function

# Prediction (in-sample)
predict(fit)
predict(fit, se.fit = TRUE)

predict(fit, type = "response")
predict(fit, type = "response", se.fit = TRUE)

# Prediction (out-of-sample)
predict(fit, newdata = data$df[1000:1200, ], se.fit = TRUE)
predict(fit, newdata = data$df[1000:1200, ], type = "response", se.fit = TRUE)

# Check if we have retrieved parameters
plot(data$df[1000:1200, "y"], predict(fit, newdata = data$df[1000:1200, ], type = "response"))
plot(fit$coefficients, data$beta)

#------------------------------------------------------- Barro Colorado Island dataset
# Original function using tommi's code
fisher_NR <- function(X, y, n, tol = 1e-16, beta_start = NULL, maxiter = 10000) {
  
  g_inv <- function(mu_target, size){
    if(mu_target <= 1){
      return(1e-7)
    } else if (mu_target >=size) {
      return(1e8)
    } else {
      uniroot(function(x) x * (digamma(x + size) - digamma(x)) - mu_target, c(1e-10, 1e8))$root
    }
  }
  g_inv <- Vectorize(g_inv, vectorize.args = c("mu_target","size"))
  
  loglik <- numeric(maxiter)
  
  # Initialization (If null, implicitely initialized at beta=0)
  
  # Initialization
  mu <- y #c(X %*% solve(crossprod(X), crossprod(X, y)))
  alpha <- g_inv(mu, n)
  eta <- log(alpha)
  w <- mu + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
  z <- eta + (y - mu) / w
  
  # First value of the likelihood
  loglik[1] <- sum(y * eta - lgamma(alpha + n) + lgamma(alpha))
  
  # Iterative procedure
  for (t in 2:maxiter) {
    beta <- solve(qr(crossprod(X * w, X)), crossprod(X * w, z))
    eta <- c(X %*% beta)
    alpha <- exp(eta)
    mu <- alpha * (digamma(alpha + n) - digamma(alpha))
    w <- mu + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
    z <- eta + (y - mu) / w
    
    loglik[t] <- sum(y * eta - lgamma(alpha + n) + lgamma(alpha))
    if (loglik[t] - loglik[t - 1] < tol) {
      return(list(beta = beta, Convergence = cbind(Iteration = (1:t) - 1, Loglikelihood = loglik[1:t])))
    }
  }
  stop("The algorithm has not reached convergence")
}

library(vegan)
data("BCI")
data("BCI.env")
y <- apply(BCI, 1, function(x) sum(x>0))
n <- rowSums(BCI)
X <- model.matrix(Precipitation ~ EnvHet + Habitat, data = BCI.env)
c(beta_NR <- fisher_NR(X = X, y = y, n = n, maxiter = 25, tol = 1e-10)$beta)

dataBCI <- data.frame(y = y, n = n, BCI.env)
fitBCI <- glm(cbind(n, y) ~ EnvHet + Habitat + Stream, family = antoniak(),
              data = dataBCI, method = "glm.fit.antoniak")

summary(fitBCI)

plotAntoniak <- function(fit){
  size <<- fit$size
  plot(fit)
  rm(size, envir = globalenv())
}

cbind("NR_tommi" = c(beta_NR), "glm routine" = fit$coefficients)

plot(predict(fitBCI, type = "response"), dataBCI$y)
abline(a = 0, b = 1, col = "red")



set.seed(10)
data <- simulate_data(N = 1200, p = 10, lambda = 1e6)

fit <- glm(cbind(n, y) ~ . -1,
           family = antoniak, 
           data = data$df[1:1000, ],      # <------ estimate using the first 1000
           method = "glm.fit.antoniak")   #


summary(fit)

fitPois <- glm(y ~ . -1, family = "poisson",
               offset = -log(data$df[1:1000, ]$n), 
               data = data$df[1:1000, -1])
summary(fitPois)
round(fitPois$coefficients, 5)





