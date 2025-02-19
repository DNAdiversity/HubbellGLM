# We test the package with some simulations

# Load packages
library(HubbellGLM)
library(tidyverse)

# Source useful functions
source("additional_functions.R")

#------------------------------------------------------------
# Part 1 - test the package on the Barro Colorado island
#------------------------------------------------------------
# We will use the BarroColorado dataset
train_index <- sort(sample(1:nrow(BarroColorado), size = 40))
y <- BarroColorado$y[train_index]
n <- BarroColorado$n[train_index]
X <- model.matrix(cbind(n, y) ~ ., BarroColorado[train_index, ])

# Test the functions with sigma = 0

#----------- Tommi's formulation (canonical link)
out_tommi <- fisher_NR(X = X, y = y, n = n)
#----------- Ale's formulation (polyseries link)
out_ale <- Hubbell_NR(X = X, y = y, n = n, sigma = 0)
#----------- HubbellGLM package
out_hub <- HubbellGLM(cbind(n, y) ~ ., family = hubbell(sigma = 0), data = BarroColorado[train_index, ])

# Check if coefficients are the same
cbind("tommi" = out_tommi$beta[, 1],
      "ale" =  out_ale$beta[, 1],
      "glm" = coef(out_hub))

# Check if logliks are the same
out_tommi$loglik
out_ale$loglik
logLik(out_hub) # <------- this is different because it uses the Stirling numbers

#------------ Predict with HubbellGLM
predict(out_hub)
predict(out_hub, type = "response")
plot(predict(out_hub, type = "response"), y)
# Predict and also adding standard errors
predict(out_hub, type = "response", se.fit = TRUE)
# Out of sample prediction
predict(out_hub, type = "response", newdata = BarroColorado[-train_index, ])

#--------------
# Run the method for a grid of sigma values, and pick the one
# with the highest AIC (on full dataset)
sigma_all <- seq(-3, 0.99, 0.05)
aics <- numeric(length(sigma_all))
for(i in 1:length(sigma_all)){
  print(sigma_all[i])
  # Run HubbellGLM
  out <- HubbellGLM(cbind(n, y) ~ .,
                    family = hubbell(sigma = sigma_all[i]),
                    data = BarroColorado)
  aics[i] <- AIC(out)
}
plot(sigma_all, aics)
best_sigma <- sigma_all[which.min(aics)]

# The lowest AIC appears to be happening with sigma = -1.2
out_best <- HubbellGLM(cbind(n, y) ~ .,
                  family = hubbell(sigma = best_sigma),
                  data = BarroColorado)
family(out_best)
summary(out_best)
plot(predict(out_best, type = "response"), BarroColorado$y)

plot(density(BarroColorado$y))
lines(density(predict(out_best, type = "response")))

#---- Test now HubbellGLM vs linear model
X <- model.matrix(cbind(n, y) ~ ., BarroColorado)
inverse_alpha <- HubbellGLM:::inv_polyseries(mu_target = BarroColorado$y,
                            size = BarroColorado$n, sigma = best_sigma)
linear_model <- lm(log(inverse_alpha) ~ X - 1)
cbind("linar_model" = coef(linear_model), "glm" = coef(out_best))
# Coefficients are very similar in this case, as we expected

#------------------------------------------------------------
# Part 2 - test the package on simulated data
#------------------------------------------------------------
# We simulate data from the polyseries link, and check if the
# AIC retrieves the correct model
set.seed(10)
sigma_true <- -1
data <- simulate_data(N = 1000, p = 10,
                      lambda = 1000,
                      sigma = sigma_true)

sigma_all <- seq(-2, 0.99, 0.05)
aics <- sapply(sigma_all, function(sigma) {
  print(sigma)
  out <- HubbellGLM(cbind(n, y) ~ . - 1,
                    family = hubbell(sigma = sigma),
                    data = data$df)
  AIC(out)
})

plot(sigma_all, aics)
abline(v = sigma_true)
abline(v = sigma_all[which.min(aics)], col = "red")

sigma_all[which.min(aics)]

out_true <- HubbellGLM(cbind(n, y) ~ . - 1,
                  family = hubbell(sigma = sigma_true),
                  data = data$df)
out_best <- HubbellGLM(cbind(n, y) ~ . - 1,
                       family = hubbell(sigma =  sigma_all[which.min(aics)]),
                       data = data$df)

cbind("true_beta" =  data$beta,
      "true_sigma" = coef(out_true),
      "best_sigma" = coef(out_best))


sigma_all[which.min(aics)]

# Try now with Hubbell_NR
logliks <- sapply(sigma_all, function(sigma) {
  print(sigma)
  Hubbell_NR(data$X, data$y, data$n, sigma = sigma)$loglik
})

plot(sigma_all, logliks)
abline(v = sigma_true)
abline(v = sigma_all[which.min(aics)], col = "red")


