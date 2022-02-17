############################################################################
############################################################################
###########                                                      ###########
###########                   Theory (cmb-lm)                    ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########           b.planterosejimenez@erasmusmc.nl           ###########
###########                                                      ###########
############################################################################
############################################################################


# Load libraries
library(MASS)
library(fastmatrix)
library(penalized)
library(ridge)
library(matrixcalc)

# Generate X
set.seed(1)
n = 100
m = 10
mu = rnorm(m) # Generate mean vector
Sigma = rWishart(n = 1, df = m, Sigma = diag(m))[,,1] # Generate Covariance matrix
X = mvrnorm(n = n, mu = mu, Sigma = Sigma) # Without intercept column
X_ = cbind(1, X) # With intercept column

# Generate y
THETA = rnorm(m+1)
y = X_ %*% THETA + rnorm(n, sd = 3)
X_y = cbind(X_, y) # Extended X matrix

# Reduced QR decomposition of X
dim(X_) # 100  11
QR = qr(X_)
R = qr.R(QR, complete = F); dim(R) # 11 11
Q = qr.Q(QR, complete = F); dim(Q) # 100  11

# OLS estimates
(theta_ols = coef(lm(y ~ X))) # lm function
solve(t(X_) %*% X_) %*% t(X_) %*% y # From Normal equation
ginv(X_) %*% y # From pseudoinverse
backsolve(R, t(Q) %*% y) # From QR
sweep.operator(t(X_y) %*% X_y, k = 1:(m+1))[1:(m+1),m+2] # From sweep operation

# Exclude columns operation
omega = c(2,4,7)
I = diag(m+1)[,-omega]
mean(X_ %*% I == X_[,-omega]) # 1
I %*% t(I) # Unitary matrix with zeros in Ajj for j in omega
t(I) %*% I # Unitary matrix of order m+1-|omega|
head(X_ %*% I %*% t(I)) # substitutes columns in omega by zeros

# Interpreting R
R[1,1]^2; nrow(X_) # n
R[1,]/R[1,1]; colMeans(X_) # mu
1/(R[1,1]^2 - 1) * (t(R) %*% R - R[1,] %*% t(R[1,])); cov(X_) # covariance matrix
chol((n-1)*cov(X_) + n*colMeans(X_) %*% t(colMeans(X_))); R # R from n, mu and covariance matrix; sign does not match because of sign multiplicity in Cholesky (but does not affect in any way)
mean(y); t(R[1,]/R[1,1]) %*% theta_ols # mean(y)

# Interpreting Q
round(colMeans(Q), 8); c(1/R[1,1], rep(0, m))
round(cov(Q), 8); diag(c(0, rep(1/(n-1), m)))

# Interpreting QTy
(Qty = t(Q) %*% y); R %*% theta_ols
Qty; (n-1)*ginv(t(R)) %*% sapply(1:ncol(X_), function(x) cov(X_[,x], y)) + c(R[1,1]*mean(y), rep(0, m))

# Extracting covariance between X_,j and y from R and theta_ols
sapply(1:ncol(X_), function(x) cov(X_[,x], y))
round(1/(n-1)*(t(R) %*% R %*% theta_ols - t(R) %*% c(R[1,] %*% theta_ols, rep(0,m))), 10)

# Cmb-lm estimates
omega = c(2,4,7); omega_c = (1:(m+1))[!(1:(m+1) %in% omega)]
(theta_omega = coef(lm(y ~ X[,-(omega-1)]))) # lm function; omega-1 as it does not include the intercept
ginv(R %*% diag(m+1)[,-omega]) %*% R %*% theta_ols # ginv
XtX = rbind(cbind(t(R) %*% R, t(R) %*% R %*% theta_ols), c(t(R) %*% R %*% theta_ols, t(y) %*% y))
diag(m+2)[-c(omega, m+2),] %*% sweep.operator(XtX, k = omega_c) %*% diag(m+2)[,-(1:(m+1))]# From sweep operation

# Cmb-lm standard errors
I = diag(m+1)[,-omega]
sigma2_epsilon = (t(y) %*% y - t(theta_omega) %*% t(I) %*% t(R) %*% R %*% I %*% theta_omega)/(n-(m+1)+length(omega))
sqrt(sum((y-X_[,-omega] %*% theta_omega)^2)/(n-(m+1)+length(omega)))
sigma(lm(y ~ X[,-(omega-1)])); sqrt(sigma2_epsilon) # Residual standard error
as.numeric(sigma2_epsilon) * ginv(t(I) %*% t(R) %*% R %*% I); unname(vcov(lm(y ~ X[,-(omega-1)]))) # squared standard errors
sigma2_epsilon*(n-(m+1)+length(omega)); diag(m+2)[-(1:(m+1)),] %*% sweep.operator(XtX, k = omega_c) %*% diag(m+2)[,-(1:(m+1))] # From sweep operation

# Confidence and prediction interval
x0 = X_[1,omega_c] # First sample as example
add = qt(p = 1-0.05/2, df = n-(m+1)+length(omega))*sqrt(sigma2_epsilon)*sqrt(t(x0) %*% ginv(t(I) %*% t(R) %*% R %*% I) %*% x0)
fit = x0 %*% theta_omega
c(fit = fit, lwr = fit - add, upr = fit + add); predict(lm(y ~ X[,-(omega-1)]), interval = 'confidence')[1,]
add = qt(p = 1-0.05/2, df = n-(m+1)+length(omega))*sqrt(sigma2_epsilon)*sqrt(t(x0) %*% ginv(t(I) %*% t(R) %*% R %*% I) %*% x0+1)
c(fit = fit, lwr = fit - add, upr = fit + add); predict(lm(y ~ X[,-(omega-1)]), interval = 'prediction')[1,]

# Ridge regression
I_ridge = diag(m+1); I_ridge[1,1] = 0
lambda = 100
(theta_ridge = coef(penalized(response = y, penalized = X, lambda1 = 0, lambda2 = lambda)))
coef(linearRidge(y ~ X, lambda = lambda, scaling = "none"))
coef(MASS::lm.ridge(y ~ X, lambda = lambda, scales = "F"))
solve(t(X_) %*% X_ + lambda*I_ridge) %*% t(X_) %*% y # Normal equation
R_l = chol(t(X_) %*% X_ + lambda*I_ridge)
XtX_l = rbind(cbind(t(R_l) %*% R_l, t(R) %*% R %*% theta_ols), c(t(R) %*% R %*% theta_ols, t(y) %*% y))
sweep.operator(XtX_l, k = 1:(m+1))[1:(m+1), m+2] # From sweep

# Cmb-ridge estimates
omega = c(2,4,7); omega_c = (1:(m+1))[!(1:(m+1) %in% omega)]; I = diag(m+1)[,-omega]
(theta_ridge_omega = coef(penalized(response = y, penalized = X[,-(omega-1)], lambda1 = 0, lambda2 = lambda)))
ginv(R_l %*% I) %*% R_l %*% theta_ridge
coef(linearRidge(y ~ X[,-(omega-1)], lambda = lambda, scaling = "none"))

# Degrees of freedom in Ridge
(df1 = n-(m+1)+length(omega))
H_cmb = X_[,-omega] %*% solve(t(X_[,-omega]) %*% X_[,-omega] + lambda*I_ridge[-omega,-omega]) %*% t(X_[,-omega])
(df2 = n-sum(diag(H_cmb))) # Effective degrees of freedom
n-sum(diag(t(I) %*% t(R) %*% R %*% I %*% ginv(t(I) %*% t(R_l) %*% R_l %*% I)))

# Cmb-ridge standard errors
(sigma2_epsilon = sum((y-X_[,-omega] %*% theta_ridge_omega)^2)/df2) # Residuals standard error
((t(y) %*% y - t(theta_ridge_omega) %*% t(I) %*% t(R_l) %*% R_l %*% I %*% theta_ridge_omega - lambda*norm(theta_ridge_omega[-1], '2')^2)/df2) # From cmb-lm
(sweep.operator(XtX_l, k = omega_c)[m+2, m+2] -lambda*norm(theta_ridge_omega[-1], '2')^2)/df2 # From sweep
# Standard error - classic formula
X_omega = X_[,-omega]; I_ridge_omega = diag(m+1-length(omega)); I_ridge_omega[1,1] = 0
(std.err2 = as.numeric(sigma2_epsilon) * solve(t(X_omega) %*% X_omega + lambda*I_ridge_omega) %*% t(X_omega) %*% X_omega %*% solve(t(X_omega) %*% X_omega + lambda*I_ridge_omega))
# Standard error - Cmb-lm
I2 = diag(m+1-length(omega))[,-1]
sigma2_epsilon * ginv(t(I) %*% t(R_l) %*% R_l %*% I) %*% (diag(m+1-length(omega)) - lambda*I2 %*% t(I2) %*% ginv(t(I) %*% t(R_l) %*% R_l %*% I))



####
# Mistakes on the ridge R-package!
# It penalizes the intercept and removes factor t(X_omega) %*% X_omega %*% solve(t(X_omega) %*% X_omega + lambda*diag(m+1-length(omega))))
vcov(linearRidge(y ~ X[,-(omega-1)], lambda = lambda, scaling = "none"))
sigma2_epsilon*df2/df1 * solve(t(X_omega) %*% X_omega + lambda*diag(m+1-length(omega)))
# Reported as an issue on Github https://github.com/SteffenMoritz/ridge/issues/15
####

# Confidence and prediction intervals
x0 = X_[1,omega_c] # First sample as example
add = qt(p = 1-0.05/2, df = df2)*sqrt(t(x0) %*% std.err2 %*% x0)
fit = x0 %*% theta_ridge_omega
c(fit = fit, lwr = fit - add, upr = fit + add)
add = qt(p = 1-0.05/2, df = df2)*sqrt(t(x0) %*% std.err2 %*% x0 + sigma2_epsilon)
c(fit = fit, lwr = fit - add, upr = fit + add)

# OLS Correlation formula
cor(X_ %*% theta_ols, y)^2
(t(theta_ols) %*% t(R) %*% R %*% theta_ols - n*mean(y)^2)/(t(y) %*% y -n*mean(y)^2)

# Ridge Correlation formula
cor(X_ %*% theta_ridge, y)^2
nom = (t(theta_ridge) %*% t(R_l) %*% R_l %*% theta_ridge - n*mean(y)^2)^2
denom = (sqrt(nom) - lambda*norm(theta_ridge[-1], "2")^2)*(t(y) %*% y -n*mean(y)^2)
nom/denom

# M is negative semi-definite
M = solve(t(R_l) %*% R_l) - solve(t(R) %*% R)
is.positive.semi.definite(round(-M,13))
mean(round(eigen(M)$values, 16) <= 1e-16)
-lambda * solve(t(R) %*% R) %*% I_ridge %*% solve(t(R_l) %*% R_l); M

# MSE
t(y-X_ %*% theta_ols) %*% (y-X_ %*% theta_ols)/n
(t(y) %*% y - t(theta_ols) %*% t(R) %*% R %*% theta_ols)/n

t(y-X_ %*% theta_ridge) %*% (y-X_ %*% theta_ridge)/n
(t(y) %*% y + t(theta_ridge) %*% t(R_l) %*% R_l %*% theta_ridge -2* t(theta_ols) %*% t(R) %*% R %*% theta_ols + 2*lambda*t(theta_ols[-1]) %*% theta_ridge[-1]- lambda*norm(theta_ridge[-1], "2")^2)/n
