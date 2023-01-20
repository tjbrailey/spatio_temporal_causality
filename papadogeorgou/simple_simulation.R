# Date: Oct 2, 2019
# Description: Simple simulations where the treatment is generated from a homogeneous
# poisson process with constant intensity and the outcome is generated from a homogeneous
# poisson process with intensity the number of observed treatment points.
#
# Note: Consider an intervention that assigns treatments as HPP(h). Then, the
# estimand is equal to h|B|.
#
# Update date: May 31, 2020.
# Update description: Include commenting.

rm(list = ls())

library(spatstat)
set.seed(123)

# Number of time points over which data are measured.
time_points <- 500

# Size of square spatial window.
window_length <- 3
# True value for the intensity: observed treatments are generated from HPP(true_h)
# with expected number of points true_h * window_length ^ 2
true_h <- 1
# Values of counterfactual intensity for which causal quantities will be calculated.
counter_h <- seq(0.7, 1.3, by = 0.05)

window <- owin(c(0, window_length), c(0, window_length))
# The set B I will consider will be the whole window.
B_window <- window
# Since I assume that Y is generated from HPP with intensity the number of treatment locations,
# The true value of the estimands are:
true_values <- counter_h * window_length ^ 2


# ----- GENERATING THE DATA ------- #

# Treatments and outcomes:
Wt <- lapply(1 : time_points, function(x) rpoispp(lambda = true_h, win = window))
Yt <- lapply(1 : time_points, function(x) rpoispp(lambda = Wt[[x]]$n / (window_length ^ 2), win = window))

# Average number of treatment and outcome points
mean(sapply(Wt, function(x) x$n))
# should be approximately true_h * window_length ^ 2
true_h * window_length ^ 2

# Turning into hyperframes to be used in the spatstat R package:
Wt <- hyperframe(Points = Wt)
Yt <- hyperframe(Points = Yt)

plot(Wt[1 : 4])  # Plotting first 4 treatments
plot(Yt[1 : 4])  # Plotting first 4 outcomes

# Smoothing the outcome process with pre-specified bandwidth 0.5.
Yt$Dens <- with(Yt, density(Points, adjust = 0.5, edge = TRUE, diggle = TRUE))
plot(Yt$Dens[1 : 4])

# Fitting the constant propensity score model:
ps_data <- hyperframe(Points = Wt$Points)
ps_fit <- mppm(Points ~ 1, data = ps_data)
fit_coefs <- summary(ps_fit)$coef
true_coefs <- log(true_h)



# ----- Calculating the estimator ------ #

# The estimates at each time point (should not be used for inference)
estimates <- array(NA, dim = c(2, time_points, length(counter_h)))
dimnames(estimates) <- list(ps_val = c('true', 'estimated'),
                            t =  1 : time_points,
                            counter_h = counter_h)

# Denominator: propensity score at the observed treatments based on fitted and true ps values:
fit_denoms <- sapply(Wt$Points, function(x) dpois(x$n, lambda = exp(fit_coefs) * (window_length ^ 2)))
true_denoms <- sapply(Wt$Points, function(x) dpois(x$n, lambda = exp(true_coefs) * (window_length ^ 2)))
# Integrated smoothed outcome surface over the set B:
integrals <- with(Yt, integral(Dens, domain = B_window))
# Numerator: density of the observed treatment under the counterfactual scenarios:
numerators <- sapply(counter_h, function(hh) sapply(Wt$Points, function(x) dpois(x$n, lambda = hh * (window_length ^ 2))))

for (hh in 1 : length(counter_h)) {
  estimates[1, , hh] <- numerators[, hh] / true_denoms * integrals
  estimates[2, , hh] <- numerators[, hh] / fit_denoms * integrals
}

# Averaging the estimates over time:
est <- apply(estimates, c(1, 3), mean)
# Variance upper bound:
cons_var <- apply(estimates ^ 2, c(1, 3), mean) / time_points

result <- abind::abind(est, cons_var, est - 1.96 * sqrt(cons_var),
                       est + 1.96 * sqrt(cons_var), along = 3)
dimnames(result)[[3]] <- c('est', 'cons_var', 'LB', 'UB')


# Plot estimates against true values:
plot(counter_h, true_values, ylim = range(c(true_values, result[, , 1])))
points(counter_h, result[1, , 1], col = 'blue')
points(counter_h, result[1, , 1], col = 'red', cex = 0.5, pch = 16)
legend('topleft', pch = c(1, 1, 16), col = c('black', 'blue', 'red'),
       legend = c('truth', 'true ps', 'est ps'))


which_h <- 5
result[1, which_h, c(1, 3, 4)]
true_values[which_h]

