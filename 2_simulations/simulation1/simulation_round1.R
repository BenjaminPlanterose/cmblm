############################################################################
############################################################################
###########                                                      ###########
###########          Simulation study, round 1 (cmb-lm)          ###########
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
library(matrixcalc)
library(MASS)
library(gplots)
library(mice)
library(gtools)
library(matrixStats)
#library(imputeMissings)
library(ggplot2)
# library(reshape2)
library(reshape)
library(Matrix)
library(fastmatrix)



# Load functions
# Global simulation wrapper
simulation <- function(mu, Cov, theta, p_h, typeNA, subtype, N_iter, conditions, dep)
{
  mus = matrix(0, ncol = 4, nrow = nrow(conditions))
  colnames(mus) = paste("mu", c("full", "mean", "MICE", "cmb_LM"), sep = ":")
  conditions = cbind(conditions, mus)
  
  for(i in 1:nrow(conditions))
  {
    print(paste(i, "out of", nrow(conditions)))
    res = matrix(ncol = 4, nrow = N_iter)
    for(l in 1:N_iter)
    {
      data = simulate_data(conditions[i, 1], mu, Cov, theta, conditions[i, 2], 
                           p_h, conditions[i, 3], typeNA, subtype)
      res[l,] = test_all(data)
    }
    conditions[i, 4:ncol(conditions)] = sapply(1:ncol(res), function(x) paste(round(res[,x], 3), collapse = ", "))
  }
  
  #return(conditions)
  return(prep_results(conditions, typeNA, subtype, dep, N_iter))
}

# Global data generation wrapper
simulate_data <- function(n, mu, Cov, theta, R2, p_h, pNA, typeNA, subtype)
{
  # 1. Data generation
  data = data_generation(n = n, mu = mu, Cov = Cov, theta = theta, R2 = R2)
  
  # 2. Splitting
  data = split(data = data, p_h = p_h, n = n)
  
  # 3. N/A injection
  data = inject.NA(data = data, p = pNA, type = typeNA, subtype = subtype)
  return(data)
}

# 
data_generation <- function(n, mu, Cov, theta, R2)
{
  m = nrow(Cov)
  X = MASS::mvrnorm(n, mu, Cov)
  y_h = as.numeric(cbind(1, X) %*% theta)
  sd_e = sqrt((1-R2)/R2*var(y_h))
  e = rnorm(n); e = (e-mean(e))/sd(e); e = sd_e*e
  y = y_h + e
  return(list(X = X, y = y))
}

# It performs hold-out split X into Xref and Xapp with a proportion p_h
split <- function(data, p_h, n)
{
  indices = sample(1:n, floor(n*p_h))
  ref = list(X = data$X[indices, ], y = data$y[indices])
  app0 = list(X = data$X[!(1:n %in% indices), ], y = data$y[!(1:n %in% indices)])
  return(list(ref = ref, app0 = app0))
}

# It injects missing values in Xapp with a proportion pNA.
inject.NA <- function(data, p, type, subtype = NULL)
{
  if(type == "MCAR")
  {
    X = MCAR(data$app0$X, p)
    if(!is.null(subtype)) (warning("no need to indicate subtype if type = MCAR"))
  }
  else if(type == "MNAR")
  {
    X = MNAR(data$app0$X, p, subtype)
    if(is.null(subtype)) (stop("if type = MNAR or MAR, subtype is required"))
  }
  else if(type == "MAR")
  {
    X = MAR(data$app0$X, p, subtype)
    if(is.null(subtype)) (stop("if type = MNAR or MAR, subtype is required"))
  }
  data$app = list(X = X, y = data$app0$y)
  return(data)
}

# MCAR missing values injecting module
MCAR <- function(X, p)
{
  n = nrow(X); m = ncol(X)
  #message(paste("pNA_max is", 1-1/m))
  k = 0:(m-1)
  probs = choose(m, k)*p^k*(1-p)^(m-k)/(1-p^m)
  n_missing = sample(k, n, replace = T, prob = probs)
  pos_missing = lapply(1:n, function(x) sample(1:m, n_missing[x], replace = F))
  for(i in 1:n)
  {
    X[i, pos_missing[[i]]] = NA
  }
  return(X)
}

# MNAR missing values injecting module
MNAR <- function(X, p, subtype)
{
  x = as.numeric(X); n = nrow(X); m = ncol(X)
  #message(paste("pNA_max is", 1-1/m))
  if(subtype == "R")
  {
    k = uniroot(f = function(k) mean(logistic(x + k)) - p, interval = c(-50, 50))$root
    A = as.numeric(logistic(X + k))
  }
  else if(subtype == "L")
  {
    k = uniroot(f = function(k) mean(logistic(-x + k)) - p, interval = c(-50, 50))$root
    A = as.numeric(logistic(-X + k))
  }
  else if(subtype == "LR")
  {
    k = uniroot(f = function(k) mean(logistic(abs(x) + k)) - p, interval = c(-50, 50))$root
    A = as.numeric(logistic(abs(X) + k))
  }
  else if(subtype == "C")
  {
    k = uniroot(f = function(k) mean(1 - logistic(abs(x) + k)) - p, interval = c(-50, 50))$root
    A = as.numeric(1 - logistic(abs(X) + k))
  }
  #print(k)
  res = matrix(sapply(1:length(A), function(x) sample(c(1,0), size = 1, prob = c(A[x], 1 - A[x]))), nrow(X), ncol(X))
  
  # Rescue pattern - all missing values
  S = which(rowMeans(res) == 1)
  if(length(S) > 0)
  {
    res = rescue_pattern(res, A, S)
  }
  X[as.logical(res)] = NA
  
  return(X)
}

# For MNAR/MCAR, samples with all samples missing need to be substituted by
# an allowed missing value pattern (respecting the imposed probabilities).
rescue_pattern <- function(res, A, S)
{
  A = matrix(A, nrow = nrow(res), ncol = ncol(res))
  pattern = permutations(2, ncol(res), c(1,0), repeats.allowed = T) # NOTE: 1=N/A
  pattern = pattern[!(rowSums(pattern) == ncol(pattern)),] # remove all missing variables pattern
  
  for(i in 1:length(S))
  {
    probs = sapply(1:nrow(pattern), function(x) sum(log(A[S[i],]^pattern[x,]*(1-A[S[i],])^(1 - pattern[x,]))))
    #probs = exp(probs - min(probs))
    probs = exp(probs)
    probs = probs/sum(probs)
    res[S[i],] = pattern[sample(1:nrow(pattern), 1, prob = probs),]
  }
  return(res)
}

# MAR missing values injecting module
MAR <- function(X, p, subtype)
{
  n = nrow(X); m = ncol(X)
  #message(paste("pNA_max is", 1-1/m))
  X_mut = matrix(NA, nrow(X), ncol(X))
  for(j in 1:ncol(X)) # Cyclic permutation
  {
    t_vals = (1-j):(m-j)
    t_vals = t_vals[!(t_vals == 0)]
    t = sample(t_vals, 1)
    X_mut[, j] = X[, j + t]
  }
  D = MNAR(X_mut, p, subtype)
  X[is.na(D)] = NA
  return(X)
}

# Logistic function
logistic <- function(x)
{
  1/(1 + exp(-x))
}

# Test all
test_all <- function(data)
{
  c(test_full(data), 
    tryCatch(test_mean_mode(data), error=function(e){NA}), 
    tryCatch(test_mice(data), error=function(e){NA}),
    tryCatch(test_cmb_lm(data), error=function(e){NA})) 
}

# Test full
test_full <- function(data)
{
  df = as.data.frame(data$ref$X)
  df$y = data$ref$y
  lm_mod = lm(formula = y ~ ., data = df)
  y_pred = predict(lm_mod, as.data.frame(data$app0$X))
  metric(y_pred, data$app0$y)
}

# Test mean-mode
test_mean_mode <- function(data)
{
  df = as.data.frame(data$ref$X)
  df$y = data$ref$y
  lm_mod = lm(formula = y ~ ., data = df)
  
  X_impute = data$app$X
  I = is.na(X_impute)
  NA_count = colSums(I)
  # means = colMeans(data$ref$X) # population means from the reference dataset
  means = colMeans(data$app$X, na.rm = T) # population means from the application dataset
  u = sapply(1:length(means), function(x) rep(means[x], NA_count[x]))
  u = unlist(u)
  X_impute[I] = u
  y_pred = predict(lm_mod, as.data.frame(X_impute))
  metric(y_pred, data$app$y)
}

# Test mice
test_mice <- function(data)
{
  df = as.data.frame(data$ref$X)
  df$y = data$ref$y
  lm_mod = lm(formula = y ~ ., data = df)
  mice_out = mice(data$app$X, m = 1, printFlag = F, method = "pmm", remove.collinear = F) # in Iris maxit=50
  #mice_out = mice(data$app$X, m = 1, printFlag = F, method = "norm", remove.collinear = F) # in Iris maxit=50
  #mice_out = mice(data$app$X, m = 1, printFlag = F, method = "mean", remove.collinear = F) # in Iris maxit=50
  #mice_out = mice(data$app$X, m = 1, printFlag = F, method = "quadratic", remove.collinear = F) # in Iris maxit=50
  data$app$X = complete(mice_out, 1)
  y_pred = predict(lm_mod, as.data.frame(data$app$X))
  metric(y_pred, data$app0$y)
}

test_cmb_lm <- function(data)
{
  QR = qr(cbind(1, data$ref$X))
  R = qr.R(QR, complete = F)
  theta = backsolve(r = R, x = t(qr.Q(QR)) %*% data$ref$y)
  m = length(theta)-1
  yty = t(data$ref$y) %*% data$ref$y
  XtX = rbind(cbind(t(R) %*% R, t(R) %*% R %*% theta), c(t(R) %*% R %*% theta, yty))
  
  X_ = cbind(1, data$app$X)
  Indicator = is.na(X_)
  omegaC_list = lapply(1:nrow(Indicator), function(x) (1:(m+1))[!(1:(m+1) %in% which(Indicator[x,]))])
  y_pred = sapply(1:nrow(X_), function(x) na.omit(X_[x,]) %*% theta_omega(XtX, omegaC_list[[x]]))
  metric(y_pred, data$app0$y)
}

theta_omega <- function(XtX, omega_c)
{
  m = nrow(XtX)-2
  sweep.operator(XtX, k = omega_c)[omega_c,m+2]
}

metric <- function(pred, obs)
{
  rho = cor(pred, obs)
  #atanh(rho)
  rho
}

prep_results <- function(RES, typeNA, subtype, dep, N_iter)
{
  full = Reduce(rbind, strsplit(RES$`mu:full`, ", "))
  colnames(full) = paste("full", 1:ncol(full), sep = "_")
  mean = Reduce(rbind, strsplit(RES$`mu:mean`, ", "))
  colnames(mean) = paste("mean", 1:ncol(mean), sep = "_")
  MICE = Reduce(rbind, strsplit(RES$`mu:MICE`, ", "))
  colnames(MICE) = paste("MICE", 1:ncol(MICE), sep = "_")
  cmbLM = Reduce(rbind, strsplit(RES$`mu:cmb_LM`, ", "))
  colnames(cmbLM) = paste("cmbLM", 1:ncol(cmbLM), sep = "_")
  
  res = Reduce(cbind, list(RES[,1:3], full, mean, MICE, cmbLM))
  res.m = melt(res, id.vars = c("n","R2", "pNA"))
  res.m$value = as.numeric(res.m$value)
  res.m$variable = sapply(strsplit(as.vector(res.m$variable), "_"), function(x) x[1])
  g1 = ggplot(data = res.m, mapping = aes(x = n, y = value^2, col = variable)) +
    facet_wrap(R2 ~ pNA, ncol = 6, labeller = labeller(pNA = label_both, R2 = label_both)) + 
    geom_point(alpha = 0.1, size = 0.1) + geom_smooth(alpha = 0.5, method = "loess", level = 0.95, span = 10) +
    theme(legend.text=element_text(size=8),
          strip.text = element_text(size = 7),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste("Cov = ", dep, ", N_iter = ", N_iter, 
                  ", typeNA = ", typeNA, ", subtypeNA = ", c("NULL", subtype)[1+!is.null(subtype)],
                  sep = "")) +
    ylab("cor^2(y_hat, y)") + ylim(c(0,1))
  res.m$variable = relevel(factor(res.m$variable), ref = "cmbLM")
  #mod1 = lm(atanh(value) ~ n + R2 + pNA + variable, data = res.m)
  mod1 = tryCatch(lm(atanh(value) ~ n + R2 + pNA + variable, data = res.m), error=function(e){NA})
  #print(summary(mod1))
  #print(g1)
  return(list(graphic = g1, mod = mod1, res.m = res.m))
}


# (Discontinued)


# (Discontinued) Test cmb-lm
test_cmb_lm0 <- function(data)
{
  create_Ip <- function(p)
  {
    diag(length(p))[, !p]
  }
  
  X_ = cbind(1, data$ref$X)
  QR = qr(X_)
  R = qr.R(QR, complete = F)
  theta = ginv(X_) %*% data$ref$y
  Qt_y = R %*% theta
  X_t = cbind(1, data$app$X)
  Indicator = is.na(X_t)
  I_list = lapply(1:nrow(Indicator), function(x) create_Ip(Indicator[x,]))
  y_pred = sapply(1:nrow(X_t), function(x) na.omit(X_t[x,]) %*% ginv(R %*% I_list[[x]]) %*% Qt_y)
  metric(y_pred, data$app0$y)
}


########################## Define Parameters ##########################

m = 10
p_h = 0.8
n_min = ceiling((m+1)/p_h)
n = n_min*seq(2, 20, 2)
R2 = seq(0.2, 0.95, 0.25)
pNA = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7)
conditions = expand.grid(n, R2, pNA)
colnames(conditions) = c("n", "R2", "pNA")
mu = c(2, 1, 0.2, -2, 3, 0.1, 1.4, -2, 3, 0.001)
theta = c(0.7, -0.4, 2, 0.1, 0.5, -2, 0.5, 3, 0.8, -0.9, 0.01)
N_iter = 50

# Independent
sigma = c(1.42, 0.82, 1.22, 2.22, 1.67, 2, 0.71, 1.8, 0.73, 1.36)
Cov_ind = diag(sigma)
heatmap.2(Cov_ind/(sqrt(diag(Cov_ind)) %*% t(sqrt(diag(Cov_ind)))), trace = "n", 
          breaks = seq(-1,1, length.out = 30))

# Medium Dependency
Cov_wd = matrix(c(1.42, -0.16, 0.32, 0.49, 0.3, 0.11, 0.24, -0.64, 0.02, -0.8,
                  -0.16, 0.82, 0.28, -0.82, 0.3, -0.39, 0.23, -0.09, -0.24, -0.14,
                  0.32, 0.28, 1.22, -0.36, 0.78, 0.19, 0.46, -0.34, -0.2, 0.28,
                  0.49, -0.82, -0.36, 2.22, 0.27, 0.75, 0.33, 0.03, 0.31, -0.26,
                  0.3, 0.3, 0.78, 0.27, 1.67, 0.84, 0.44, 0.42, -0.17, -0.16,
                  0.11, -0.39, 0.19, 0.75, 0.84, 2, -0.03, 1.14, 0.03, -0.48,
                  0.24, 0.23, 0.46, 0.33, 0.44, -0.03, 0.71, -0.14, -0.06, -0.05,
                  -0.64, -0.09, -0.34, 0.03, 0.42, 1.14, -0.14, 1.8, -0.03, -0.31,
                  0.02, -0.24, -0.2, 0.31, -0.17, 0.03, -0.06, -0.03, 0.73, 0.09,
                  -0.8, -0.14, 0.28, -0.26, -0.16, -0.48, -0.05, -0.31, 0.09, 1.36),
                nrow = 10, byrow = T)
is.positive.definite(Cov_wd)
heatmap.2(Cov_wd/(sqrt(diag(Cov_wd)) %*% t(sqrt(diag(Cov_wd)))), trace = "n", 
          breaks = seq(-1,1, length.out = 30))

# Strong Dependency
Cov_sd = matrix(c(1.42, -0.06, 0.58, 0.36, 0.92, -0.12, 0.51, 0.61, -0.02, -0.02,
                  -0.06, 0.82, -0.15, -0.19, -0.09, 0.58, -0.13, -0.14, 0.01, -0.09,
                  0.58, -0.15, 1.22, 0.84, 1.03, -0.58, 0.79, 1.10, -0.05, 0.11,
                  0.36, -0.19, 0.84, 2.22, 0.99, -0.31, 0.73, 1.02, 0.01, 0.28,
                  0.92, -0.09, 1.03, 0.99, 1.67, -0.21, 0.96, 1.17, -0.15, 0.15,
                  -0.12, 0.58, -0.58, -0.31, -0.21, 2, -0.41, -0.42, -0.03, -0.19,
                  0.51, -0.13, 0.79, 0.73, 0.96, -0.41, 0.71, 0.91, -0.11, 0.13,
                  0.61, -0.14, 1.10, 1.02, 1.17, -0.42, 0.91, 1.8, -0.17, 0.06,
                  -0.02, 0.01, -0.05, 0.01, -0.15, -0.03, -0.11, -0.17, 0.73, -0.13,
                  -0.02, -0.09, 0.11, 0.28, 0.15, -0.19, 0.13, 0.06, -0.13, 1.36),
                nrow = 10, byrow = T)
is.positive.definite(Cov_sd)
heatmap.2(Cov_sd/(sqrt(diag(Cov_sd)) %*% t(sqrt(diag(Cov_sd)))), trace = "n",
          breaks = seq(-1, 1, length.out = 30))


# Quantifying dependence - Normalized eigenvalues of the correlation matrix
l_ind = eigen(Cov_ind/(sqrt(diag(Cov_ind)) %*% t(sqrt(diag(Cov_ind)))))$values/m; sum(l_ind)
l_wd = eigen(Cov_wd/(sqrt(diag(Cov_wd)) %*% t(sqrt(diag(Cov_wd)))))$values/m; sum(l_wd)
l_sd = eigen(Cov_sd/(sqrt(diag(Cov_sd)) %*% t(sqrt(diag(Cov_sd)))))$values/m; sum(l_sd)
plot(1:m, l_ind, type = 'o', ylim = c(0, 0.5), 
     pch = 19, col = alpha('black', 0.7), xlab = c('PC'), ylab = 'Proportion of variance explained (standarized)')
lines(1:m, l_wd, type = 'o', col = alpha('green3', 0.7), pch = 19)
lines(1:m, l_sd, 'o', col = alpha('blue3', 0.7), pch = 19)
legend('right', legend = c('Indep', 'wd', 'sd'), fill = c('black', 'green3', 'blue3'), bty = 'n', cex = 1.3)

# Normalized entropy
-sum(l_ind*log(l_ind)/log(m))
-sum(l_wd*log(l_wd)/log(m))
-sum(l_sd*log(l_sd)/log(m))

# Normalized Gini-Simpson index
l_max = rep(1/m, m)
(1-sum(l_ind^2))/(1-sum(l_max^2))
(1-sum(l_wd^2))/(1-sum(l_max^2))
(1-sum(l_sd^2))/(1-sum(l_max^2))

m/(m-1) - m/(m-1)*sum(l_wd^2)


########################## Simulation 1: Independent+MCAR ##########################
start_time <- Sys.time()
set.seed(123)
out1 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MCAR", subtype = NULL, conditions = conditions, 
                  dep = "independent")
########################## Simulation 2: wd+MCAR ##########################
set.seed(124)
out2 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MCAR", subtype = NULL, conditions = conditions, 
                  dep = "WD")
########################## Simulation 3: sd+MCAR ##########################
set.seed(125)
out3 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MCAR", subtype = NULL, conditions = conditions, 
                  dep = "SD")
########################## Simulation 4: Independent+MAR+L ##########################
set.seed(126)
out4 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "L", conditions = conditions, 
                  dep = "independent")
########################## Simulation 5: wd+MAR+L ##########################
set.seed(127)
out5 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "L", conditions = conditions, 
                  dep = "wd")
########################## Simulation 6: sd+MAR+L ##########################
set.seed(128)
out6 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "L", conditions = conditions, 
                  dep = "sd")
########################## Simulation 7: Independent+MAR+LR ##########################
set.seed(129)
out7 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "LR", conditions = conditions, 
                  dep = "independent")
########################## Simulation 8: wd+MAR+LR ##########################
set.seed(130)
out8 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "LR", conditions = conditions, 
                  dep = "wd")
########################## Simulation 9: sd+MAR+LR ##########################
set.seed(131)
out9 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "LR", conditions = conditions, 
                  dep = "sd")
########################## Simulation 10: Independent+MAR+C ##########################
set.seed(132)
out10 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "C", conditions = conditions, 
                  dep = "independent")
########################## Simulation 11: wd+MAR+C ##########################
set.seed(133)
out11 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "C", conditions = conditions, 
                  dep = "wd")
########################## Simulation 12: sd+MAR+C ##########################
set.seed(134)
out12 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MAR", subtype = "C", conditions = conditions, 
                  dep = "sd")
########################## Simulation 13: Independent+MNAR+L ##########################
set.seed(135)
out13 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "L", conditions = conditions, 
                  dep = "independent")
########################## Simulation 14: wd+MNAR+L ##########################
set.seed(136)
out14 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "L", conditions = conditions, 
                  dep = "wd")
########################## Simulation 15: sd+MNAR+L ##########################
set.seed(137)
out15 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "L", conditions = conditions, 
                  dep = "sd")
########################## Simulation 16: Independent+MNAR+LR ##########################
set.seed(138)
out16 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "LR", conditions = conditions, 
                  dep = "independent")
########################## Simulation 17: wd+MNAR+LR ##########################
set.seed(139)
out17 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "LR", conditions = conditions, 
                  dep = "wd")
########################## Simulation 18: sd+MNAR+LR ##########################
set.seed(140)
out18 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "LR", conditions = conditions, 
                  dep = "sd")
########################## Simulation 19: Independent+MNAR+C ##########################
set.seed(141)
out19 = simulation(mu = mu, Cov = Cov_ind, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "C", conditions = conditions, 
                  dep = "independent")

# Error in diag(length(p)) : 
#   could not find symbol "names" in environment of the generic function
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# create_Ip <- function(p)
# {
#   diag(length(p))[, !p]
# }
# [1] "235 out of 240"


########################## Simulation 20: wd+MNAR+C ##########################
set.seed(142)
out20 = simulation(mu = mu, Cov = Cov_wd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "C", conditions = conditions, 
                  dep = "wd")
########################## Simulation 21: sd+MNAR+C ##########################
set.seed(143)
out21 = simulation(mu = mu, Cov = Cov_sd, theta = theta, p_h = p_h, N_iter = N_iter,
                  typeNA = "MNAR", subtype = "C", conditions = conditions, 
                  dep = "sd")
end_time <- Sys.time()
end_time - start_time # Time difference of 21.84888 hours

################################## Save results ##################################

out = list(out1 = out1,   out2 = out2,   out3 = out3,   out4 = out4,   out5 = out5,
           out6 = out6,   out7 = out7,   out8 = out8,   out9 = out9,   out10 = out10, 
           out11 = out11, out12 = out12, out13 = out13, out14 = out14, out15 = out15, 
           out16 = out16, out17 = out17, out18 = out18, out19 = out19, out20 = out20, 
           out21 = out21)
setwd("/media/ultron/2tb_disk2/0_startallover/CMB_LM/0_simulations/results/round1/")
saveRDS(out, "results_1_0.Rds")

################################## Solve errors in lm ##################################

# round(value, 3)
# if rho = +-1, atanh = +Inf. Linear model spits error:
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# NA/NaN/Inf in 'y'
where = which(sapply(out, function(x) tryCatch(is.na(x$mod$rank), error=function(e){T})))
for(i in 1:length(where))
{
  A = out[[where[i]]]$res.m
  A[which(abs(A$value) == 1),]$value = 1-1e-04
  out[[where[i]]]$mod = lm(atanh(value) ~ n + R2 + pNA + variable, data = A)
}
where = which(sapply(out, function(x) tryCatch(is.na(x$mod$rank), error=function(e){T})))
setwd("/media/ultron/2tb_disk2/0_startallover/CMB_LM/0_simulations/results/round1/")
saveRDS(out, "results_1.Rds")


##
setwd("/media/ultron/2tb_disk2/0_startallover/CMB_LM/0_simulations/results/round4/")
saveRDS(out2, "result_out2_quadraticMICE.Rds")
##

################################## Session Info ##################################

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] fastmatrix_0.3-8196 Matrix_1.3-4        reshape_0.8.8       ggplot2_3.3.5       matrixStats_0.61.0 
# [6] gtools_3.9.2        mice_3.13.0         gplots_3.1.1        MASS_7.3-54         matrixcalc_1.0-5   
# 
# loaded via a namespace (and not attached):
# [1] backports_1.3.0             BiocFileCache_2.0.0         plyr_1.8.6                  splines_4.1.2              
# [5] BiocParallel_1.26.2         GenomeInfoDb_1.28.4         digest_0.6.28               foreach_1.5.1              
# [9] fansi_0.5.0                 magrittr_2.0.1              memoise_2.0.0               tzdb_0.2.0                 
# [13] limma_3.48.3                Biostrings_2.60.2           readr_2.1.0                 annotate_1.70.0            
# [17] askpass_1.1                 siggenes_1.66.0             prettyunits_1.1.1           colorspace_2.0-2           
# [21] blob_1.2.2                  rappdirs_0.3.3              dplyr_1.0.7                 crayon_1.4.2               
# [25] RCurl_1.98-1.5              genefilter_1.74.1           lme4_1.1-27.1               GEOquery_2.60.0            
# [29] survival_3.2-13             iterators_1.0.13            glue_1.5.0                  gtable_0.3.0               
# [33] zlibbioc_1.38.0             XVector_0.32.0              DelayedArray_0.18.0         Rhdf5lib_1.14.2            
# [37] BiocGenerics_0.38.0         HDF5Array_1.20.0            scales_1.1.1                jomo_2.7-2                 
# [41] abind_1.4-5                 mvtnorm_1.1-3               DBI_1.1.1                   rngtools_1.5.2             
# [45] Rcpp_1.0.7                  xtable_1.8-4                progress_1.2.2              bumphunter_1.34.0          
# [49] bit_4.0.4                   mclust_5.4.8                preprocessCore_1.54.0       stats4_4.1.2               
# [53] httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2              farver_2.1.0               
# [57] pkgconfig_2.0.3             XML_3.99-0.8                dbplyr_2.1.1                locfit_1.5-9.4             
# [61] utf8_1.2.2                  tidyselect_1.1.1            rlang_0.4.12                AnnotationDbi_1.54.1       
# [65] munsell_0.5.0               tools_4.1.2                 cachem_1.0.6                generics_0.1.1             
# [69] RSQLite_2.2.8               broom_0.7.10                stringr_1.4.0               mvmeta_1.0.3               
# [73] fastmap_1.1.0               yaml_2.2.1                  bit64_4.0.5                 beanplot_1.2               
# [77] caTools_1.18.2              scrime_1.3.5                purrr_0.3.4                 KEGGREST_1.32.0            
# [81] nlme_3.1-153                doRNG_1.8.2                 sparseMatrixStats_1.4.2     nor1mix_1.3-0              
# [85] xml2_1.3.2                  biomaRt_2.48.3              compiler_4.1.2              rstudioapi_0.13            
# [89] filelock_1.0.2              curl_4.3.2                  png_0.1-7                   tibble_3.1.6               
# [93] stringi_1.7.5               GenomicFeatures_1.44.2      minfi_1.38.0                lattice_0.20-45            
# [97] nloptr_1.2.2.3              multtest_2.48.0             vctrs_0.3.8                 pillar_1.6.4               
# [101] lifecycle_1.0.1             rhdf5filters_1.4.0          data.table_1.14.2           bitops_1.0-7               
# [105] micemd_1.8.0                rtracklayer_1.52.1          mixmeta_1.2.0               GenomicRanges_1.44.0       
# [109] R6_2.5.1                    BiocIO_1.2.0                KernSmooth_2.23-20          IRanges_2.26.0             
# [113] codetools_0.2-18            boot_1.3-28                 assertthat_0.2.1            rhdf5_2.36.0               
# [117] SummarizedExperiment_1.22.0 openssl_1.4.5               rjson_0.2.20                withr_2.4.2                
# [121] GenomicAlignments_1.28.0    Rsamtools_2.8.0             S4Vectors_0.30.2            GenomeInfoDbData_1.2.6     
# [125] parallel_4.1.2              hms_1.1.1                   quadprog_1.5-8              grid_4.1.2                 
# [129] tidyr_1.1.4                 base64_2.0                  minqa_1.2.4                 DelayedMatrixStats_1.14.3  
# [133] illuminaio_0.34.0           MatrixGenerics_1.4.3        Biobase_2.52.0              restfulr_0.0.13  

########################## Read results ##########################


setwd("/media/ultron/2tb_disk2/0_startallover/CMB_LM/0_simulations/results/round1/")
out = readRDS("results_1.Rds")
sapply(1:length(out), function(x) out[[x]]$graphic$layers[[1]]$aes_params$alpha <<- 0.025)
sapply(1:length(out), function(x) out[[x]]$graphic$layers[[2]]$aes_params <<- list(alpha = 0.2, size = 0.3, linetype = 1))

setwd('/media/ultron/2tb_disk2/0_startallover/CMB_LM/0_simulations/results/')
for(i in 1:length(out))
{
  message(i)
  ggsave(plot = out[[i]]$graphic, paste('round1_simulation', i, '.pdf', sep = ''))
  ggsave(plot = out[[i]]$graphic, paste('round1_simulation', i, '.tiff', sep = ''), dpi = 300)
  capture.output(print(paste('############################## Simulation', i, '##############################')), 
                 file = 'linear_models.txt', append = T)
  capture.output(summary(out[[i]]$mod), file = 'linear_models.txt', append = T)
}

