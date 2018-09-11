
library(tidyverse)

euler = function(f, y0, a, b, h) {
  ts = seq(a, b, by = h)
  
  y = matrix(0, nrow = length(ts), ncol = length(y0))

  y[1,] = y0
  for (i in 1:(length(ts) - 1)) {
    y[i + 1,] = y[i,] + h*f(ts[i], y[i,])
  }
  
  m = cbind(ts, y)
  colnames(m) = c("t", if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0))
  m
}

fixedPoint = function(f, x0, tol = 1e-10, steps = 20) {
  while (TRUE) {
    xn = f(x0)
    if (near(xn, x0, tol) || steps == 0)
      break
    x0 = xn
    steps = steps - 1
  }
  xn
}

bwEuler = function(f, y0, a, b, h) {
  ts = seq(a, b, by = h)
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  
  y[1,] = y0
  for (i in 1:(length(ts) - 1)) {
    y[i + 1,] = fixedPoint(function(x) y[i,] + h*f(ts[i + 1], x), y[i,])
  }
  
  m = cbind(ts, y)
  colnames(m) = c("t", if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0))
  m
}

crankNicolson = function(f, y0, a, b, h) {
  ts = seq(a, b, by = h)
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  
  y[1,] = y0
  for (i in 1:(length(ts) - 1)) {
    y[i + 1,] = (y[i,] + h*f(ts[i], y[i,]) + fixedPoint(function(x) y[i,] + h*f(ts[i + 1], x), y[i,]))/2
  }
  
  m = cbind(ts, y)
  colnames(m) = c("t", if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0))
  m
}

huen = function(f, y0, a, b, h) {
  ts = seq(a, b, by = h)
  
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  
  y[1,] = y0
  for (i in 1:(length(ts) - 1)) {
    y[i + 1,] = y[i,] + h/2*(f(ts[i], y[i,]) + f(ts[i + 1], y[i,] + h*f(ts[i], y[i,])))
  }
  
  m = cbind(ts, y)
  colnames(m) = c("t", if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0))
  m
}
