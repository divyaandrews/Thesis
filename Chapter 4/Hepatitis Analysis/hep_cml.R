library(readr)
library(forecast)
library("cowplot")
library("miscTools")
library("maxLik")
library(ggplot2)
library(ggfortify)
library(dplyr)
library(gridExtra)
library(grid)

  d <- data.frame(x = x,
                  y = cumsum(y) / sum(y),
                  upper = 1 / xm * x + crit,
                  lower = 1 / xm * x - crit)
  p <- ggplot2::ggplot(data = d, mapping = ggplot2::aes_string(x = 'x', y = 'y')) +
    geom_line(colour = colour, linetype = linetype) +
    ggplot2::scale_x_continuous(name = 'Frequency', limits = c(0, xm)) +
    ggplot2::scale_y_continuous(name = 'Cum. Periodogram', limits = c(0, 1))
  p <- plot_confint(p = p, data = d, conf.int = conf.int,
                    conf.int.colour = conf.int.colour,
                    conf.int.linetype = conf.int.linetype,
                    conf.int.fill = conf.int.fill, conf.int.alpha = conf.int.alpha)
  p
}



#Import data
data_hep<- read.csv("C:\\Users\\USER\\Downloads\\hepb.csv")
data_1 <- data_hep$data
data_1 <- data_1[1:110]
Tlen <- length(data_1)
df <- data.frame(1:110,data_1)
plot1<-ggplot(df, aes(x = 1:110, y = data_1)) +
  geom_line(color="darkmagenta")+
  geom_point(size=0.2)+
  scale_y_continuous(breaks = seq(0, 14, by = 2))+
  scale_x_continuous(breaks = seq(0, 110, by = 10))+
  labs(x = "Week", y = "Hep-B Cases", title = "") 

acf_values <- acf(data_1, plot = FALSE)
pacf_values <- pacf(data_1, plot = FALSE)
# Create a data frame including lag 0
acf_df <- data.frame(
  Lag = acf_values$lag[1:22],
  ACF = acf_values$acf[1:22]
)
# Create a data frame including lag 0
pacf_df <- data.frame(
  Lag = pacf_values$lag,
  PACF = pacf_values$acf
)

# Plot using ggplot2
plot2<-ggAcf(
  data_1,
  lag.max = 20,
  type = c("correlation", "covariance", "partial"),
  plot = TRUE,
  na.action = na.contiguous,
  demean = TRUE,
)+ ggtitle("")

plot3<-ggplot(pacf_df, aes(x = Lag, y = PACF)) +
  geom_bar(stat = "identity", position = "identity", width = 0.05, fill = "black") +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(data_1)), linetype = "dashed", color = "blue") +
  labs(title = "", x = "Lag", y = "PACF") +
  scale_x_continuous(breaks = seq(0, 21, by = 2))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  theme_minimal()+
  theme(
    panel.grid = element_blank(),                   # Removes grid lines
    #panel.border = element_rect(color = "black", fill = NA, linetype = "solid"),  # Adds outline box
    axis.ticks = element_line(linetype = "dashed")  # Dashed lines at y-axis ticks
  )

# Combining 3 plots in one window
ggdraw() +
  draw_plot(plot1,  x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(plot2, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(plot3, x = .5, y = 0, width = .5, height = .5)

hist(data_1-0.5)
maxval<-max(data_1)
#Absolute frequencies
absfreq <- tabulate(data_1+1) #+1 to include 0
plot(0:maxval, absfreq/Tlen, type="h", xlab = "k", ylab = expression(paste("estimated P(X"[t],"=k)")), lwd=4, ylim=c(0,0.12))

#acf(data_1, plot=FALSE)[[1]][2:11]

rho1 <- acf(data_1, plot=FALSE)[[1]][2]
rho1 #0.4050859
rho2 <- acf(data_1, plot=FALSE)[[1]][3]
rho2 #0.3398476

#Observations' mean
barX <- mean(data_1)
barX

#Observations' variance
sX <- var(data_1)
sX
#6.435696

#Start with moment estimates:
#delta1 := alpha1+beta1
delta1mm <- rho2/rho1
delta1mm #0.6238413

alpha0mm <- barX*(1-delta1mm)
alpha0mm #1.313136

alpha1mm <- (1-delta1mm^2 - sqrt(1-delta1mm^2) * sqrt(1-delta1mm^2+4*rho1*(delta1mm-rho1)) )/(-2*(delta1mm-rho1))

beta1mm <- delta1mm-alpha1mm

dnoge <- function(x, theta, phi) {
  delta <- as.numeric(x == 0)
  delta * phi + (1 - delta) * (1 - phi) * (1 - theta)^(x - 1) * theta
}
pnoge <- function(x, theta, phi) {
  value <- 0
  if (x >= 0) {
    for (y in 0:x) {
      value <- value + dnoge(y, theta, phi)
    }
  }
  value
}
llnogepingarch11 <- function(par, data) {
  # par = (alpha0, beta1, alpha1, phi, lambda1)
  alpha0 <- par[1]
  beta1  <- par[2]
  alpha1 <- par[3]
  phi    <- par[4]
  cmean <- par[5] # initialization of λ₁
  
  T <- length(data)
  value <- 0
  
  for (t in 2:T) {
    cmean <- alpha0 + alpha1 * data[t-1] + beta1 * cmean
    theta <- (1 - phi) / cmean
    prob <- dnoge(data[t], theta, phi)
    value <- value - log(prob)
  }
  value
}





#Log-likelihood of NoGe-INGARCH(1,1) model:
llnogeingarch11mod <- function(par,data){
  #par is vector (beta0,beta1,alpha1,p,m1)
  llnogeingarch11(reparam(par),data)
}
estmlnoge <- suppressWarnings(optim(c(alpha0ml,alpha1ml+beta1ml,beta1ml/(alpha1ml+beta1ml),0.4,barX), llnogeingarch11mod, method="L-BFGS-B", lower=rep(1e-4,5), upper=c(2,0.9999,0.9999,0.9999,10), control=list(ndeps=rep(1e-4,5)), data=data_1))
estmlnoge$par
#0.2929597  0.9700822  0.8634678  0.3380541 13.0994283
reparam(estmlnoge$par)


# Parameter Estimates:
round(estmlnoge$par, 4)

estmlng <- suppressWarnings(optim(reparam(estmlnoge$par), llnogepingarch11, method="L-BFGS-B", lower=reparam(estmlnoge$par)-1e-3, upper=reparam(estmlnoge$par)+1e-3, control=list(ndeps=rep(1e-4,5)), data=data_1, hessian=TRUE))
estmlng$par

# Standard Errors:
se_noge <- sqrt(diag(solve(estmlng$hessian)))
round(se_noge, 4)

# Log-likelihood, AIC, BIC
neglik <- estmlng$value
Tlen <- length(data_1)
AIC_noge <- 2 * neglik + 2 * 4
BIC_noge <- 2 * neglik + log(Tlen) * 4

c(LogLik = -neglik, AIC = AIC_noge, BIC = BIC_noge)


alpha0mlng <- estmlng$par[[1]]
beta1mlng <- estmlng$par[[2]]
alpha1mlng <- estmlng$par[[3]]
pmlng <- estmlng$par[[4]]
ofiestng <- estmlng$hessian #inverse covariance
neglmaxng <- estmlng$value
estcovng <- solve(ofiestng)

#Estimates:
round(c(alpha0mlng,alpha1mlng,beta1mlng, pmlng),4)
#1.5811977 0.1786201 0.3806499 0.7238932

#Estimated standard errors:
c(sqrt(diag(estcovng)))
#0.84524914  0.29130282  0.10107850  0.09952219 22.53507260

#AIC and BIC:
AICng <- 2*neglmaxng+2*4
BICng <- 2*neglmaxng+log(Tlen)*4

#Pearson residuals
res <- c()
mt <- estmlng$par[[5]]

for(t in c(2:Tlen)){
  mt <- alpha0mlng+beta1mlng*mt+alpha1mlng*data_1[t-1]
  res <- c(res, (data_1[t]-mt)/sqrt(mt)*sqrt(pmlng))
}
acf_values <- acf(res, plot = FALSE)

# Create a data frame including lag 0
acf_df <- data.frame(
  Lag = acf_values$lag[1:11],
  ACF = acf_values$acf[1:11]
)

plot_10<-ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", position = "identity", width = 0.03, fill = "black") +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(res)), linetype = "dashed", color = "blue") +
  labs(title = "NoGe-INGARCH", x = "Lag", y = "ACF") +
  scale_x_continuous(breaks = seq(0, 10, by = 1))+
  theme_minimal()+
  theme(
    panel.grid = element_blank(),                   # Removes grid lines
    panel.border = element_rect(color = "black", fill = NA, linetype = "solid"),  # Adds outline box
    axis.ticks = element_line(linetype = "dashed")  # Dashed lines at y-axis ticks
  )
plot(res)
mean(res)
#-0.001765745
var(res)
#1.052954
