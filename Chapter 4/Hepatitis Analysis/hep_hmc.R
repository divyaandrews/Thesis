library(rstan)
library(loo)
library(scoringRules)
library(ggplot2)
library(tidyr)
library(dplyr)
library(VGAM)
library(LaplacesDemon)
data_hep<- read.csv("/Users/apple/Documents/Paper3/R files/hepb.csv")
data_1 <- data_hep$data
data_1 <- data_1[1:110]
#Tlen <- length(data_1)
#df <- data.frame(1:110,data_1)
Y <- data_1
N<-length(data_1)
#noge-ingarch
model <- stan_model('/Users/apple/Documents/Paper3/R files/hep_noge.stan')
options(mc.cores=4)
fit <- sampling(model,list(N=N,Y=Y),algorithm="NUTS", iter = 2000,warmup=1000, chains=4,thin=1)
s<- summary(fit)
est<- unlist(s[1]$summary[c(1,5,6,7,8),1])
sd<-unlist(s[1]$summary[c(1,5,6,7,8),2])
rhat<- unlist(s[1]$summary[c(1,5,6,7,8),10])
neff<- unlist(s[1]$summary[c(1,5,6,7,8),9])
log_lik_4 <- extract_log_lik(fit,
                             parameter_name = "log_lik",
                             merge_chains = FALSE)
col_means <- apply(log_lik_4[, , 1], 2, mean)
ebic_ng <- -2*(sum(col_means))/4 + 4*log(110)

posterior <- rstan::extract(fit)
dnoge <- function(x, theta, phi) {
  delta <- as.numeric(x == 0)
  delta * phi + (1 - delta) * (1 - phi) * (1 - theta)^(x - 1) * theta
}
pnoge <- function(x, theta, phi) {
  sapply(x, function(xi) {
    if (xi >= 0) {
      sum(dnoge(0:xi, theta, phi))
    } else {
      0
    }
  })
}


phi_0 <- mean(posterior$phi)
# posterior$lambda: G x T matrix of posterior λ_t for each draw g
lambda_samples <- posterior$cmean  # dimensions G x T
G <- dim(lambda_samples)[1]
T <- dim(lambda_samples)[2]

log_cpo_vec <- numeric(T)

for (t in 1:T) {
  liks <- dnoge(Y[t], theta = (1-phi_0)/lambda_samples[, t], phi=phi_0)  # vector of length G
  cpo_t <- 1 / mean(1 / liks)
  log_cpo_vec[t] <- log(cpo_t)
}

# Final log-CPO score
log_CPO <- sum(log_cpo_vec)

max_y <- max(Y) + 10  # support for CRPS computation
crps_values <- numeric(T)

for (t in 1:T) {
  obs <- Y[t]
  crps_t <- numeric(G)
  
  for (g in 1:G) {
    lambda_t <- lambda_samples[g, t]
    support <- 0:max_y
    F_vals <- pnoge(support,theta = (1-phi_0)/lambda_t, phi= phi_0)
    
    # CRPS for discrete distribution at x = obs
    term1 <- sum(F_vals[1:(obs+1)]^2)
    term2 <- sum((F_vals[(obs+2):(max_y+1)] - 1)^2)
    crps_t[g] <- term1 + term2
  }
  
  crps_values[t] <- mean(crps_t)
}

# Average CRPS over time
CRPS_avg <- mean(crps_values)

lambda_110_samp <- posterior$cmean[,110]
alpha0_samp <- posterior$alpha0
alpha1_samp <- posterior$alpha1
beta1_samp <- posterior$beta1

H <- 27       # number of forecast steps
G <- 1000      # number of posterior samples
support <- 0:1000  # Support of the predictive distribution
# True values (test set)
Y_test <- data_hep$data[111:137]

# Preallocate
pred_mean <- numeric(H)
pred_median <- numeric(H)
pred_lower <- numeric(H)
pred_upper <- numeric(H)
pred_matrix <- matrix(NA, nrow = H, ncol = G)  # Store sampled values (optional)

# Initial values
x_prev <- data_1[110]
lambda_prev <- lambda_110_samp  # Vector of length G

for (h in 1:H) {
  # Step 1: Compute lambda_{t+1}^{(g)} for all G posterior draws
  lambda_next <- alpha0_samp + alpha1_samp * x_prev + beta1_samp * lambda_prev  # Vector of length G
  
  # Step 2: Construct predictive mixture distribution
  # For each k in support, average the probabilities from all G Poisson(lambda_next)
  pred_probs <- sapply(support, function(k) {
    mean(dnoge(k, theta =(1-phi_0)/lambda_next, phi =phi_0))  # Mixture probability at support point k
  })
  
  # Ensure probabilities sum to 1 (numerical stability)
  pred_probs <- pred_probs / sum(pred_probs)
  
  # Step 3: Compute predictive summaries
  pred_mean[h]   <- sum(support * pred_probs)  # E[Y]
  cdf <- cumsum(pred_probs)
  pred_median[h] <- support[which.min(abs(cdf - 0.5))]  # Find closest to median
  pred_lower[h]  <- support[which(cdf >= 0.025)[1]]
  pred_upper[h]  <- support[which(cdf >= 0.975)[1]]
  
  # (Optional) Simulate samples from the predictive distribution for storage
  pred_matrix[h, ] <- sample(support, G, replace = TRUE, prob = pred_probs)
  
  # Step 4: Update for next step
  x_prev <- Y_test[h]  # Plug in predictive mean
  lambda_prev <- lambda_next  # Reuse same lambda_next values
}



# RMSE
rmse_mean <- sqrt(mean((Y_test - pred_mean)^2)); rmse_mean
rmse_median <- sqrt(mean((Y_test - pred_median)^2))

# MAE
mae_mean <- mean(abs(Y_test - pred_mean)); 
mae_median <- mean(abs(Y_test - pred_median)); mae_median

# Create a data frame
plot_data <- data.frame(
  Time = 111:137,
  Observed = Y_test,
  PredMean = pred_mean,
  PredMedian = pred_median,
  Lower = pred_lower,
  Upper = pred_upper
)
colnames(plot_data) <- c("Time", "Observed", "Predictive Mean", "Predictive Median", "Lower", "Upper")

# Convert to long format
plot_long <- plot_data %>%
  pivot_longer(cols = c("Observed", "Predictive Mean", "Predictive Median"),
               names_to = "Type", values_to = "Value")

# Create the plot
plot1<-ggplot(plot_long, aes(x = Time, y = Value, color = Type, shape = Type)) +
  geom_point(size = 2, fill = "white", stroke = 1) +  # Stroke for border, fill="white" for unfilled shapes
  
  # Add ribbon for 95% credible interval
  geom_ribbon(data = plot_data,
              aes(x = Time, ymin = Lower, ymax = Upper),
              inherit.aes = FALSE,
              fill = "red", alpha = 0.2) +
  
  scale_color_manual(values = c(
    "Observed" = "black",
    "Predictive Mean" = "green",
    "Predictive Median" = "darkblue"
  )) +
  scale_shape_manual(values = c(
    "Observed" = 16,    # filled circle
    "Predictive Mean" = 1,     # open circle
    "Predictive Median" = 2    # open triangle
  )) +
  ylim(0, 30)+
  labs(
    x = "Week",
    y = "Hep-B Cases",
    title = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
   legend.title = element_blank(),
    panel.border = element_blank(),        # Remove the box
  panel.grid = element_blank(),          # Remove grid lines
    axis.line = element_line(color = "black"),  # Draw x and y axis lines
    axis.ticks = element_line(color = "black")  # Keep axis ticks
  )
