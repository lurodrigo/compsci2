
library(tidyverse)

# source('R/ivp_methods.R')

central_difference = function(f, ts, ya, yb, h) {
  n = length(ts) - 2
  
  fs_between = f(ts[2:(length(ts) - 1)])
  
  A = diag(rep(2, n))
  for (i in 1:(n-1)) {
    A[1 + i, i] = -1
    A[i, 1 + i] = -1
  }
  
  c(ya, solve(A, h*h*fs_between), yb)
}

## Example 1

f = function(t) sin(pi * t)
exact = function(t) { sin(pi * t) / (pi * pi) }

map_dfr(2^(2:7) - 1, function(N) {
  ts = seq(0, 1, length.out = N + 2)
  h = ts[2] - ts[1]
  numerical = central_difference(f, ts, 0, 0, h)
  exact_sol = exact(ts)
  
  half = (N + 3) / 2
  
  tibble(N = N, h = h, e = abs(exact_sol[half] - numerical[half]))
})

## Example 2

f = function(t) sin(pi * t)
exact = function(t) { sin(pi * t) / (pi * pi) }
