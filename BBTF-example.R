library(fdrtool)

## function------
quant=function(x){ quantile(x, prob=c(0.025,0.975)) }
source("BBTF.HS.function.R")
## --------------

# true upper boundary
n <- 100
x <- 1:n
th0 <- sqrt(x)/2

# mcmc
mc <- 10500
burn <- 500
th <- 5
thin <- seq(from=1, to=mc-burn, by=th)

# sigmoid approximation
eta <- 500

# order of TF
k <- 1

# data generate
set.seed(1)
sd <- 1
theta <- (1/sd)*sqrt(pi/2)
noise <- rhalfnorm(n, theta)
y <- th0 - noise

# fitting
fit <- BBTF.HS(y, x=x, mc=mc, burn=burn, eta=eta, k=k, upper=TRUE, shape="NI")
Est <- apply(fit[[1]][thin,], 2, mean)
CI <- apply(fit[[1]][thin,], 2, quant)


# plot
plot(th0, type="l", ylim = c(min(y), max(th0)))
points(y)
points(Est, col=2, type="l")
points(th0, type="l")
polygon(c(1:100,rev(1:100)), c(CI[1,],rev(CI[2,])), col=adjustcolor(2, alpha=0.4), border=NA)



