
library(tidyverse)

euler = function(f, y0, ts) {
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  colnames(y) = if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0)
  y[1,] = y0
  
  for (i in 1:(length(ts) - 1)) {
    h = ts[i + 1] - ts[i]
    y[i + 1,] = y[i,] + h*f(ts[i], y[i,])
  }
  
  y
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

bwEuler = function(f, y0, ts) {
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  colnames(y) = if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0)
  y[1,] = y0
  
  for (i in 1:(length(ts) - 1)) {
    h = ts[i + 1] - ts[i]
    y[i + 1,] = fixedPoint(function(x) y[i,] + h*f(ts[i + 1], x), y[i,])
  }
  
  y
}

crankNicolson = function(f, y0, ts) {
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  colnames(y) = if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0)
  y[1,] = y0
  
  for (i in 1:(length(ts) - 1)) {
    h = ts[i + 1] - ts[i]
    y[i + 1,] = (y[i,] + h*f(ts[i], y[i,]) + fixedPoint(function(x) y[i,] + h*f(ts[i + 1], x), y[i,]))/2
  }
  
  y
}

huen = function(f, y0, ts) {
  y = matrix(0, nrow = length(ts), ncol = length(y0))
  colnames(y) = if (is.null(names(y0))) paste0("y", 1:length(y0)) else names(y0)
  y[1,] = y0
  
  for (i in 1:(length(ts) - 1)) {
    h = ts[i + 1] - ts[i]
    y[i + 1,] = y[i,] + h/2*(f(ts[i], y[i,]) + f(ts[i + 1], y[i,] + h*f(ts[i], y[i,])))
  }

  y
}

rk_mat = function(M, f, y0, a, b, h) {
  ts = seq(a, b, by = h)
  
  y = double(length(ts))
  y[1] = y0
  
  
}

M = matrix(c(0, 0, 0, 0, 0, 
  .5, .5, 0, 0, 0,
  .5, 0, .5, 0, 0,
  1, 0, 0, 1, 0,
  0, 1/6, 1/3, 1/3, 1/6), byrow = TRUE, nrow = 5)


