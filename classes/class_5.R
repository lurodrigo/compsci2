
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

tb = euler(function(t, u) (u - t)/(u + t), c("u" = 1), 0, 7, 0.1)
ggplot(tb %>% as.data.frame(), aes(x = t, y = u)) + geom_path()

sigma = 10
rho = 27
beta = 8/3
tb = euler(function(t, y) c(sigma*(y[2]- y[1]), 
                            y[1]*(rho - y[3]) - y[2], 
                            y[1]*y[2] - beta*y[3]), 
           c(1, 0, 0), 0, 30, 0.005) %>% as.data.frame()

ggplot(tb %>% as.data.frame(), aes(x = t, y = y1)) + geom_path()
ggplot(tb %>% as.data.frame(), aes(x = y2, y = y1)) + geom_path()

## aula 5

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

tb = euler(function(t, u) (u - t)/(u + t), c("u" = 1), 0, 7, 0.1)
ggplot(tb %>% as.data.frame(), aes(x = t, y = u)) + geom_path()

tb = bwEuler(function(t, u) (u - t)/(u + t), c("u" = 1), 0, 7, 0.1)
ggplot(tb %>% as.data.frame(), aes(x = t, y = u)) + geom_path()

tb = crankNicolson(function(t, u) (u - t)/(u + t), c("u" = 1), 0, 7, 0.1)
ggplot(tb %>% as.data.frame(), aes(x = t, y = u)) + geom_path()

tb = huen(function(t, u) (u - t)/(u + t), c("u" = 1), 0, 7, 0.1)
ggplot(tb %>% as.data.frame(), aes(x = t, y = u)) + geom_path()

####
arg = list(function(t, u) -5*u, c("u" = 1), 0, 1, 0.025)
tb = invoke_map(list(euler, bwEuler, crankNicolson, huen), 
                  list(arg, arg, arg, arg)) %>%
  map(as.data.frame) %>%
  bind_cols %>%
  select(t, mt_euler = u, mt_bwEuler = u1, mt_crankNicolson = u2, mt_huen = u3) %>%
  mutate(u = exp(-5*t)) %>%
  tail(1)

abs((tb %>% select(starts_with("mt")) %>% as.matrix) - tb$u[nrow(tb)]) 
