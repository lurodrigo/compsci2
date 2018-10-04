
source("classes/ivp_methods.R")

solveIVP = function(u0, du0, h = 0.001) {
  ys = euler(f = function(t, y) c(y[2], -exp(y[1] + 1)), 
             y0 = c("u" = u0, "u'" = du0), 
             ts = seq(0, 1, by = h))
  ys
}

secant = function(f, x0, x1, tol = 1e-6, maxsteps = 100) {
  i = 0
  while (TRUE) {
    i = i + 1
    x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))
    
    if (abs(x1 - x2) < tol || i > maxsteps)
      break
    
    x0 = x1
    x1 = x2
  }
  
  x2
}

solveBVP = function(u0, u1) {
  s = secant(function(du0) tail(solveIVP(u0, du0), 1)[1] - u1, 0, 1)
  solveIVP(u0, s)
}