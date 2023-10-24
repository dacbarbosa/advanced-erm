## Chap 2 - Selection on Observable - Introduction to nonparametric regression

## This routine provides a smooth introduction to nonparametric regression
# - An example of empirical cdf that needs to be smoothed
# - Several examples of Kernel distribution
# - Regression methods: local average, local regression, kernel regression

## ---- Empirical CDF
set.seed=44935718
F10 <- ecdf(rnorm(100))
plot(function(x) pnorm(x), -3, 3)
lines(F10, verticals = TRUE, do.points = FALSE)

## ---- Kernel distributions
library(evmix)
x<-seq(-2,2,0.01) 
plot(function(x) kdepanechnikov(x),-2,2)
lines(x,kdgaussian(x))
lines(x,kduniform(x))
lines(x,kdtriangular(x))

## ---- plot Gaussian kernel function
xstar <- 0.3
plot(x, dnorm((x - xstar) / (1 / 60)), 
     t = "l", ylab = "K([x - 0.3] / h)", cex.lab = 1.5, 
     cex.axis = 1.25, ylim = c(0, 0.4))
for(k in 2:4){
  lines(x, dnorm((x - xstar) / (k / 60)) , lty = k)
}
legend("topright", legend = c("h = 1/60", "h = 2/60", "h = 3/60", "h = 4/60"),
       lty = 1:4, bty = "n")


## An example from http://users.stat.umn.edu/~helwig/notes/smooth-notes.html
## ---- Local averaging

# define function and data
set.seed(1)
n <- 101
x <- seq(0, 1, length.out = n)
fx <- sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# define x* and color for window
xstar <- 0.3
cols <- rgb(190/255, 190/255, 190/255, alpha = 0.5)

# set-up 2 x 2 subplot
par(mfrow = c(2,2))

# loop through spans (0.1, 0.2, 0.3, 0.4)
for(s in c(0.1, 0.2, 0.3, 0.4)){
  
  # plot data and true function
  plot(x, y, main = paste0("span = ", s), ylim = c(-2.5, 2.5),
       cex.lab = 1.5, cex.axis = 1.25)
  lines(x, fx, col = "blue", lwd = 2)
  
  # plot window
  window <- c(xstar - s / 2, xstar + s / 2)
  rect(window[1], -3, window[2], 3, col = cols)
  
  # define weights - exclude points that are outside the window for the calculation of the average
  w <- rep(0, n)
  w[x >= window[1] & x <= window[2]] <- 1
  
  # plot estimate
  ystar <- sum(y * w) / sum(w)
  points(xstar, ystar, pch = 17, col = "red", cex = 1)
  
  # add legend
  legend("topright", legend = c("data", "truth"),
         pch = c(1, NA), lty = c(NA, 1), col = c("black", "blue"), bty = "n")
  legend("bottomright", legend = c("estimate", "window"),
         pch = c(17, 15), col = c("red", "gray"), bty = "n")

}

## ---- Local linear Regression - uniform weights
# set-up 2 x 2 subplot
par(mfrow = c(2,2))

# loop through spans (0.1, 0.2, 0.3, 0.4)
for(s in c(0.1, 0.2, 0.3, 0.4)){
  
  # plot data and true function
  plot(x, y, main = paste0("span = ", s), ylim = c(-2.5, 2.5),
       cex.lab = 1.5, cex.axis = 1.25)
  lines(x, fx, col = "blue", lwd = 2)
  
  # plot window
  window <- c(xstar - s / 2, xstar + s / 2)
  rect(window[1], -3, window[2], 3, col = cols)
  
  # define weights - exclude points that are outside the window for the calculation of the average
  w <- rep(0, n)
  w[x >= window[1] & x <= window[2]] <- 1
  
  # plot estimate
  X.w <- sqrt(w) * cbind(1, x)
  y.w <- sqrt(w) * y
  beta <- solve(crossprod(X.w)) %*% crossprod(X.w, y.w)
  ystar <- as.numeric(cbind(1, xstar) %*% beta)
  points(xstar, ystar, pch = 17, col = "red", cex = 1)
  
  # add regression line
  abline(beta, lty = 3)
  
  # add legend
  legend("topright", legend = c("data", "truth"),
         pch = c(1, NA), lty = c(NA, 1), col = c("black", "blue"), bty = "n")
  legend("bottomright", legend = c("estimate", "window"),
         pch = c(17, 15), col = c("red", "gray"), bty = "n")
  
}

## ---- Local quadratic regression - with uniform weights

# set-up 2 x 2 subplot
par(mfrow = c(2,2))

# loop through spans (0.1, 0.2, 0.3, 0.4)
for(s in c(0.1, 0.2, 0.3, 0.4)){
  
  # plot data and true function
  plot(x, y, main = paste0("span = ", s), ylim = c(-2.5, 2.5),
       cex.lab = 1.5, cex.axis = 1.25)
  lines(x, fx, col = "blue", lwd = 2)
  
  # plot window
  window <- c(xstar - s / 2, xstar + s / 2)
  rect(window[1], -3, window[2], 3, col = cols)
  
  # define weights
  w <- rep(0, n)
  w[x >= window[1] & x <= window[2]] <- 1
  
  # plot estimate
  X <- cbind(1, x - 0.5, (x - 0.5)^2)
  X.w <- sqrt(w) * X
  y.w <- sqrt(w) * y
  beta <- solve(crossprod(X.w)) %*% crossprod(X.w, y.w)
  ystar <- as.numeric(cbind(1, xstar - 0.5, (xstar - 0.5)^2) %*% beta)
  points(xstar, ystar, pch = 17, col = "red", cex = 1)
  
  # add regression line
  lines(x, X %*% beta, lty = 3)
  
  # add legend
  legend("topright", legend = c("data", "truth"),
         pch = c(1, NA), lty = c(NA, 1), col = c("black", "blue"), bty = "n")
  legend("bottomright", legend = c("estimate", "window"),
         pch = c(17, 15), col = c("red", "gray"), bty = "n")
  
}

## ---- Local Kernel regression

# set-up 2 x 2 subplot
par(mfrow = c(2,2))

# loop through h = c(1, 2, 3, 4) / 60
for(h in c(1:4)/60){
  
  # plot data and true function
  plot(x, y, main = paste0("h = ", h * 60, "/60"), ylim = c(-2.5, 2.5),
       cex.lab = 1.5, cex.axis = 1.25)
  lines(x, fx, col = "blue", lwd = 2)
  
  # plot 99% window
  window <- c(xstar - 3 * h, xstar + 3 * h)
  rect(window[1], -3, window[2], 3, col = cols)
  
  # define weights
  w <- dnorm((x - xstar) / h) 
  w <- w / sum(w)
  
  # plot estimate
  ystar <- sum(y * w)
  points(xstar, ystar, pch = 17, col = "red", cex = 1)
  
  # add legend
  legend("topright", legend = c("data", "truth"),
         pch = c(1, NA), lty = c(NA, 1), col = c("black", "blue"), bty = "n")
  legend("bottomright", legend = c("estimate", "99% area"),
         pch = c(17, 15), col = c("red", "gray"), bty = "n")
  
}


# Exercise, change the weights to use a epaechnikov instead of gaussian. Draw the whole function on the interval [0.2;0.8]