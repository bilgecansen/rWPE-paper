
library(tseries)

# Equilibrium dynamics
x1 <- rick_grow(mu_r = 1, sd = 0.001, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 50)
plot(x1, type= "l")
acf(x1, lag = 1, pl = F)

x2 <- sample(x1)
plot(x2, type= "l")
acf(x2, lag = 1, pl = F)

x3 <- random_grow(x1)
plot(x3, type= "l")
acf(x3, lag = 1, pl = F)

# Chaos r = 3.0
x1 <- rick_grow(mu_r = 3, sd = 0.001, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 50)
plot(x1, type= "l")
acf(x1, lag = 1, pl = F)
calculate_pe(x1)

x2 <- sample(x1)
plot(x2, type= "l")
acf(x2, lag = 1, pl = F)
calculate_pe(x2)

x3 <- random_grow(x1)
plot(x3, type= "l")
acf(x3, lag = 1, pl = F)
calculate_pe(x3)

# Chaos r = 3.4
x1 <- rick_grow(mu_r = 3.4, sd = 0.001, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 50)
plot(x1, type= "l")
acf(x1, lag = 1, pl = F)
calculate_pe(x1)

x2 <- sample(x1)
plot(x2, type= "l")
acf(x2, lag = 1, pl = F)
calculate_pe(x2)

x3 <- random_grow(x1)
plot(x3, type= "l")
acf(x3, lag = 1, pl = F)
calculate_pe(x3)
