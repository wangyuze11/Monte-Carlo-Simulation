S_0 <- 100
r <- 0.05
sigma <- 0.2
K <- 100
T_total <- 1
t1 <- 0.5
t2 <- 0.75
F_1_afterdiv <- 97.53
F_2_afterdiv <- 93.76
# calculate the forward prices right before the dividends
F_1 <- S_0 * exp(r * t1)
F_2 <- F_1_afterdiv * exp(r * (t2 - t1))
# calculate the dividends' cash value
d1 <- F_1 - F_1_afterdiv
d2 <- F_2 - F_2_afterdiv
y1 <- F_1_afterdiv / F_1
y2 <- F_2_afterdiv / F_2
# calculate the dividends' ratios
d1.ratio <- d1 / F_1
d2.ratio <- d2 / F_2
F_T <- F_2_afterdiv * exp(r * (T_total - t2)) # forward price after 2nd dividend
d <- (log(F_T / K) - 0.5 * sigma^2 * T_total)/(sigma * sqrt(T_total))
p1 <- pnorm(d + sigma * sqrt(T_total))
p2 <- pnorm(d)
p <- exp(-r * T_total) * (F_T * p1 - K * p2) # theoretical value of the option
p
W <- matrix(rnorm(300000),nrow=3)
S_1 <- S_0 * exp((r - 0.5 * sigma^2) * t1 + sigma * sqrt(t1) * W[1,]) * y1
S_2 <- S_1 * exp((r - 0.5 * sigma^2) * (t2-t1) + sigma * sqrt(t2-t1) * W[2,]) * y2
S_T <- S_2 * exp((r - 0.5 * sigma^2) * (T_total-t2) + sigma * sqrt(T_total-t2) * W[3,])
payoff <- (S_T - K)
payoff[payoff<0] <- 0
X <- exp(-r * T_total) * payoff
d1.ratio # dividend 1's proportion
d2.ratio # dividend 2's proportion
mean(X) # estimated option price
W <- matrix(rnorm(300000), nrow=3)
S_1 <- S_0 * exp((r - 0.5 * sigma^2) * t1 + sigma * sqrt(t1) * W[1,]) - d1
S_2 <- S_1 * exp((r - 0.5 * sigma^2) * (t2-t1) + sigma * sqrt(t2-t1) * W[2,]) - d2
S_T <- S_2 * exp((r - 0.5 * sigma^2) * (T_total-t2) + sigma * sqrt(T_total-t2) * W[3,])
payoff <- (S_T - K)
payoff[payoff<0] <- 0
Y <- exp(-r * T_total) * payoff
mu <- mean(Y)
sd <- sqrt(var(Y))
d1 # cash value of dividend 1
d2 # cash value of dividend 2
mu # estimated option price
mu + c(-1,1) * qnorm(0.975) * sd / sqrt(100000) # confidence interval
#X Y assigned in previous chunks
lambda <- lm(Y~X)$coefficient[2]
Y_lambda <- Y - lambda * (X - p)
mu <- mean(Y_lambda)
sd <- sqrt(var(Y_lambda))
mu
mu + c(-1,1) * qnorm(0.975) * sd / sqrt(100000) # confidence interval

library(MASS)
T_total <- 1 #assign values
S1_0 <- 1
S2_0 <- 1
K <- 2
r <- 0.05
mu <- r
c <- c(0, 0.5, -0.5)
sigma <- sqrt(2 + 2 * c)
b <- sqrt(2 + 2 * c)
a <- 2 * r - 1
#compute the theoretical values for Y options
d <- (log(S1_0 * S2_0 / K) + a * T_total)/(sqrt(T_total) * b)
EY <- S1_0 * S2_0 * exp((r + c) * T_total) * pnorm(d + b * sqrt(T_total)) - K * exp(-r *T_total)
#simulate the Y and Z option prices
Y_price <- rep(0, times = 3)
Z_price <- rep(0, times = 3)
Z_lambda_price <- rep(0, times = 3)
rho <- rep(0, times = 3)
factor <- rep(0, times = 3)
sd <- rep(0, times = 3)
for (i in 1:3){
  corr <- c[i]
  SIGMA <- matrix(c(1,corr,corr,1), nrow = 2) #decompose the covariance matrix
  W <- mvrnorm(50000, mu = c(0, 0), Sigma = SIGMA)
  S1 <- S1_0 * exp((r - 0.5 * 1^2) * T_total + 1 * sqrt(T_total) * W[,1])
  S2 <- S2_0 * exp((r - 0.5 * 1^2) * T_total + 1 * sqrt(T_total) * W[,2])
  Y <- exp(-r * T_total) * (S1 * S2 - K)
  Y[Y<0] <- 0
  Z <- exp(-r * T_total) * (S1 + S2 - K)
  Z[Z<0] <- 0
  Z_sd <- sqrt(var(Z))
  Y_sd <- sqrt(var(Y))
  cov <- cov(Y, Z)
  rho[i] <- cov/(Z_sd * Y_sd)
  lambda <- lm(Z ~ Y)$coefficient[2]
  Z_lambda <- Z - lambda * (Y - EY[i])
  sd[i] <- sqrt(var(Z_lambda))
  Y_price[i] <- mean(Y)
  Z_lambda_price[i] <- mean(Z_lambda)
  Z_price[i] <- mean(Z)
}
factor <- round(1/(1-rho^2), digits = 4)
width <- qnorm(0.975) * sd / sqrt(50000)
ci_l <- Z_lambda_price - width
ci_u <- Z_lambda_price + width
row.name <- c("Expectation of Y", "sample mean of Y", "sample mean of Z",
              "sample mean of Lambda Z", "Reduction Factor",
              "95% CI lower bound", "95% CI upper bound")
table <- as.data.frame(rbind(EY, Y_price, Z_price, Z_lambda_price, factor, ci_l, ci_u))
table <- cbind(row.name, table)
names(table) <- c("Variable", "c=0", "c=0.5", "c=-0.5")
library(knitr)
kable(table, align = 'c', row.names = F, digits = 4)
factor
