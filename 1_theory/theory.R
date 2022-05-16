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

########################### Simulate data ###########################

# Generate X
set.seed(1)
n = 100
m = 10
mu = rnorm(m) # Generate mean vector
Sigma = rWishart(n = 1, df = m, Sigma = diag(m))[,,1] # Generate random Covariance matrix
X = mvrnorm(n = n, mu = mu, Sigma = Sigma) # Without intercept column
X_ = cbind(1, X) # With intercept column

# Generate y
THETA = rnorm(m+1)
y = X_ %*% THETA + rnorm(n, sd = 3)
X_y = cbind(X_, y) # Extended X matrix

########################### Exclude operation ###########################

# Exclude columns operation
omega = c(2,4,7)
I = diag(m+1)[,-omega]
X_[,-omega] == X_ %*% I # Equivalent to column exclusion
I %*% t(I) # Unitary matrix with zeros in Ajj for j in omega
t(I) %*% I # Unitary matrix of order m+1-|omega|
head(X_ %*% I %*% t(I)) # substitutes columns in omega by zeros


########################### Obtain cmb-lm parameters ###########################

# Reduced QR decomposition of X
dim(X_) # 100  11
QR = qr(X_)
R = qr.R(QR, complete = F); dim(R) # 11 11
Q = qr.Q(QR, complete = F); dim(Q) # 100  11

# OLS estimates; can be obtained from normal eq;, Planitz equ., QR or sweep operator.
(theta_ols = coef(lm(y ~ X))) # lm function
solve(t(X_) %*% X_) %*% t(X_) %*% y # From Normal equation
ginv(X_) %*% y # From pseudoinverse
backsolve(R, t(Q) %*% y) # From QR
sweep.operator(t(X_y) %*% X_y, k = 1:(m+1))[1:(m+1),m+2] # From sweep operator.

########################### Cmb-lm parameter interpretation ###########################

# Interpreting R
R[1,1]^2; nrow(X_) # n
R[1,]/R[1,1]; colMeans(X_) # mu
1/(R[1,1]^2 - 1) * (t(R) %*% R - R[1,] %*% t(R[1,])); cov(X_) # covariance matrix
chol((n-1)*cov(X_) + n*colMeans(X_) %*% t(colMeans(X_))); R # R from n, mu and covariance matrix; sign does not match because of sign multiplicity in Cholesky (but does not affect in anyway)


########################### Cmb-lm estimates ###########################

omega = c(2,4,7); omega_c = (1:(m+1))[!(1:(m+1) %in% omega)]

# From submodel
unname((theta_omega = coef(lm(y ~ X[,-(omega-1)])))) # lm function; omega-1 as it does not include the intercept

# From cmb-lm
ginv(R %*% diag(m+1)[,-omega]) %*% R %*% theta_ols # ginv

# From sweep operator
XtX = rbind(cbind(t(R) %*% R, t(R) %*% R %*% theta_ols), c(t(R) %*% R %*% theta_ols, t(y) %*% y))
diag(m+2)[-c(omega, m+2),] %*% sweep.operator(XtX, k = omega_c) %*% diag(m+2)[,-(1:(m+1))]# From sweep operation


########################### Cmb-lm std. errors - cmb-lm ###########################

I = diag(m+1)[,-omega]

# Variance of residuals (submodel)
sigma(lm(y ~ X[,-(omega-1)]))^2 # Residual standard error
sum((y-X_[,-omega] %*% theta_omega)^2)/(n-(m+1)+length(omega))

# Variance of residuals (cmb-lm)
(sigma2_epsilon = (t(y) %*% y - t(theta_omega) %*% t(I) %*% t(R) %*% R %*% I %*% theta_omega)/(n-(m+1)+length(omega)))

# Variance of residuals (sweep operator)
(sigma2_epsilon = sweep.operator(XtX, k = omega_c)[m+2, m+2]/(n-(m+1)+length(omega))) # From sweep operation

# Standard errors (submodel)
unname(vcov(lm(y ~ X[,-(omega-1)])))

# Standard errors (cmb-lm or sweep)
as.numeric(sigma2_epsilon) * ginv(t(I) %*% t(R) %*% R %*% I) # squared standard errors


########################### Confidence/Prediction intervals - cmb-lm ###########################

x0 = X_[1, omega_c] # First sample as example
dim(x0) = c(1, length(omega_c))
fit = x0 %*% theta_omega

### Confidence intervals
# Submodel
predict(lm(y ~ X[,-(omega-1)]), interval = 'confidence')[1,]

# Cmb-lm
add = qt(p = 1-0.05/2, df = n-(m+1)+length(omega))*sqrt(sigma2_epsilon)*sqrt(x0 %*% ginv(t(I) %*% t(R) %*% R %*% I) %*% t(x0))
c(fit = fit, lwr = fit - add, upr = fit + add)

### Prediction intervals
# Submodel
predict(lm(y ~ X[,-(omega-1)]), interval = 'prediction')[1,]

# Cmb-lm
add = qt(p = 1-0.05/2, df = n-(m+1)+length(omega))*sqrt(sigma2_epsilon)*sqrt(x0 %*% ginv(t(I) %*% t(R) %*% R %*% I) %*% t(x0)+1)
c(fit = fit, lwr = fit - add, upr = fit + add)


########################### Ridge regression ###########################

lambda = 100 # Regularization parameter

# From submodel
(theta_ridge = coef(penalized(response = y, penalized = X, lambda1 = 0, lambda2 = lambda)))
# Cannot compare with MASS::lm.ridge since it scales input matrix

# From closed-form solution
I_ridge = diag(m+1); I_ridge[1,1] = 0
solve(t(X_) %*% X_ + lambda*I_ridge) %*% t(X_) %*% y # Normal equation

########################### Cmb-Ridge regression ###########################

# From submodel
(theta_ridge_omega = coef(penalized(response = y, penalized = X[,-(omega-1)], lambda1 = 0, lambda2 = lambda)))

# Cmb-ridge estimates
omega = c(2,4,7); omega_c = (1:(m+1))[!(1:(m+1) %in% omega)]; I = diag(m+1)[,-omega]
R_l = chol(t(X_) %*% X_ + lambda*I_ridge)
ginv(R_l %*% I) %*% R_l %*% theta_ridge

# From sweep
XtX_l = rbind(cbind(t(R_l) %*% R_l, t(R_l) %*% R_l %*% theta_ridge), c(t(R_l) %*% R_l %*% theta_ridge, t(y) %*% y))
sweep.operator(XtX_l, k = omega_c)[omega_c, m+2] # From sweep

########################### Cmb-Ridge d.f. ###########################

# Degrees of freedom in Ridge
H_cmb = X_[,-omega] %*% solve(t(X_[,-omega]) %*% X_[,-omega] + lambda*I_ridge[-omega,-omega]) %*% t(X_[,-omega])
(df_ridge = n-sum(diag(H_cmb))) # Effective degrees of freedom
n-sum(diag(I %*% ginv(R_l %*% I) %*% R_l %*% solve(t(R_l) %*% R_l) %*% t(R) %*% R))

########################### Standard errors - cmb-ridge ###########################

# Variance residuals (submodel)
(sigma2_epsilon = sum((y-X_[,-omega] %*% theta_ridge_omega)^2)/df_ridge) # Residuals standard error

# Variance residuals (cmb-lm)
((t(y) %*% y - t(theta_ridge_omega) %*% t(I) %*% t(R_l) %*% R_l %*% I %*% theta_ridge_omega - lambda*norm(theta_ridge_omega[-1], '2')^2)/df_ridge) # From cmb-lm

# Variance residuals (sweep operator)
(sweep.operator(XtX_l, k = omega_c)[m+2, m+2] -lambda*norm(theta_ridge_omega[-1], '2')^2)/df_ridge # From sweep

# Standard errors (submodel)
X_omega = X_[,-omega]; I_ridge_omega = diag(m+1-length(omega)); I_ridge_omega[1,1] = 0
(std.err2 = as.numeric(sigma2_epsilon) * solve(t(X_omega) %*% X_omega + lambda*I_ridge_omega) %*% t(X_omega) %*% X_omega %*% solve(t(X_omega) %*% X_omega + lambda*I_ridge_omega))

# Standard error - classic formula
I2 = diag(m+1-length(omega))[,-1]
sigma2_epsilon * ginv(t(I) %*% t(R_l) %*% R_l %*% I) %*% (diag(m+1-length(omega)) - lambda*I2 %*% t(I2) %*% ginv(t(I) %*% t(R_l) %*% R_l %*% I))

########################### Confidence/prediction intervals - cmb-ridge ###########################

x0 = X_[1,omega_c] # First sample as example
dim(x0) = c(1, length(omega_c))

add = qt(p = 1-0.05/2, df = df_ridge)*sqrt(x0 %*% std.err2 %*% t(x0))
fit = x0 %*% theta_ridge_omega
c(fit = fit, lwr = fit - add, upr = fit + add)
add = qt(p = 1-0.05/2, df = df_ridge)*sqrt(x0 %*% std.err2 %*% t(x0) + sigma2_epsilon)
c(fit = fit, lwr = fit - add, upr = fit + add)


