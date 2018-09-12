
lsq = function(t, u) {
  n = length(u)
  w1 = (sum(u)*sum(t) - n*sum(u*t))/(sum(t)^2 - n*sum(t^2))
  w0 = (sum(u) - w1 * sum(t))/n
  c("intercept" = w0, "slope" = w1)
}

lsqm = function(A, b) {
  tA = t(A)
  solve(tA %*% A, tA %*% b)
}

lsqm_svd = function(A, b) {
  decomp = svd(A)
  x = diag(1/decomp$d) %*% t(decomp$u) %*% b
  decomp$v %*% x
}

tb = data.frame(t = 1:5, u = c(1, 2, 1.3, 3.75, 2.25))
lsqsqm(cbind(1, tb$t), tb$u)
lsqsqm(cbind(1, tb$t, tb$t^2), tb$u)
lsqsqm_svd(cbind(1, tb$t, tb$t^2), tb$u)
