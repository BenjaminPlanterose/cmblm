############################################################################
############################################################################
###########                                                      ###########
###########            Application: Epigenetic ageing            ###########
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
library(data.table)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(wateRmelon)
library(fastmatrix)
library(mice)
library(tidyverse)
library(miceadds)
library(micemd)
library(glmnet)
library(penalized)
library(gplots)

# Load functions
# Perform mean imputation
mean_impute <- function(data)
{
  mu = colMeans(data, na.rm = T)
  for(i in 1:ncol(data))
  {
    data[is.na(data[,i]), i] = mu[i]    
  }
  return(data)
}

# Horvarth transform
transform <- function(age, age_adult = 20)
{
  (age <= age_adult)*log((age+1)/(age_adult + 1)) + (age > age_adult)*(age-age_adult)/(age_adult+1)
}

# Inverse Horvath transform
inverse_transform <- function(age_trans, age_adult = 20)
{
  (age_trans <= 0)*(exp(age_trans)*(age_adult + 1)-1) + (age_trans > 0)*(age_trans*(age_adult+1) + age_adult)
}

# theta_omega, sweep operator version
theta_omega <- function(XtX, omega_c)
{
  m = nrow(XtX)-2
  sweep.operator(XtX, k = omega_c)[omega_c,m+2]
}

# cmb-lm prediction modelling
cmb_lm <- function(R, theta, X_app, yty)
{
  m = length(theta)-1
  XtX = rbind(cbind(t(R) %*% R, t(R) %*% R %*% theta), c(t(R) %*% R %*% theta, yty))
  Indicator = is.na(X_app)
  omegaC_list = lapply(1:nrow(Indicator), function(x) (1:(m+1))[!(1:(m+1) %in% which(Indicator[x,]))])
  y_pred = sapply(1:nrow(X_), function(x) na.omit(X_[x,]) %*% theta_omega(XtX, omegaC_list[[x]]))
  y_pred
}

# cmb-lm prediction modelling (parallelized version)
par_cmb_lm <- function(R, theta, X_app, nThread = NULL, yty)
{
  pred <- function(x, XtX, omega_c)
  {
    na.omit(x) %*% theta_omega(XtX, omega_c)
  }
  if(is.null(nThread))
  {
    nThread <- detectCores(logical = FALSE)
  }
  m = length(theta)-1
  XtX = rbind(cbind(t(R) %*% R, t(R) %*% R %*% theta), c(t(R) %*% R %*% theta, yty))
  Indicator = is.na(X_app)
  omegaC_list = lapply(1:nrow(Indicator), function(x) (1:(m+1))[!(1:(m+1) %in% which(Indicator[x,]))])
  y_pred = parallel::mclapply(1:nrow(X_app), 
                              function(x) pred(X_app[x,], XtX, omegaC_list[[x]]), 
                              mc.cores = nThread)
  y_pred = unlist(y_pred)
  return(y_pred)
}

# Transform cmb-lm model into cmb-ridge model
regularize_L2_model <- function(lambda, model)
{
  I = diag(length(model$theta_ols)); I[1,1] = 0
  theta_ridge = solve(t(model$R) %*% model$R + lambda*I) %*% t(model$R) %*% model$R %*% model$theta_ols
  R_m = chol(t(model$R) %*% model$R + lambda*I)
  return(list(R = R_m, theta_ols = as.numeric(theta_ridge), yTy = model$yTy, lambda = lambda))
}

# Inject missing values (MCAR scheme)
inject.na = function(X, pNA)
{
  p0 = mean(is.na(X))
  if(pNA < p0)
  {
    stop('pNA < p0')
  }
  p = pNA - p0
  DIM = dim(X)
  x = as.numeric(X)
  where = which(is.na(x))
  left = (1:length(x))[-where]
  what = sample(left, floor(p*(length(x))))
  x[what] = NA
  dim(x) = DIM
  return(x)
}

# Testing cmblm and mean-imputation (correlation mode); MICE was too inefficient for this stage
test_cmblm_impute <- function(model, X, y)
{
  # Mice imputation
  # mice_out = mice(X, m = 1, printFlag = F, method = "pmm", remove.collinear = F)
  # X_impute2 = complete(mice_out, 1)
  # y_h = inverse_transform(X_impute2 %*% model$theta_ols)
  # cor_MICE = tryCatch(cor(y_h, y, use = 'complete.obs'), error = function(e) NA)
  
  # mean-imputation
  X_impute = mean_impute(X)
  y_h = inverse_transform(X_impute %*% model$theta_ols)
  cor_impute = tryCatch(cor(y_h, y, use = 'complete.obs'), error = function(e) NA)
  
  # Cmb-lm
  y_h = inverse_transform(par_cmb_lm(R = model$R, theta = model$theta_ols, yty = model$yTy, X_app = X, nThread = 4))
  cor_cmblm = cor(y_h, y)
  
  return(c(cor_impute = cor_impute, cor_cmblm = cor_cmblm))
}

# Testing cmblm and mean-imputation (bias mode); MICE was too inefficient for this stage
test_cmblm_impute2 <- function(model, X, y)
{
  # Mice imputation
  # mice_out = mice(X, m = 1, printFlag = F, method = "pmm", remove.collinear = F)
  # X_impute2 = complete(mice_out, 1)
  # y_h = inverse_transform(X_impute2 %*% model$theta_ols)
  # cor_MICE = tryCatch(cor(y_h, y, use = 'complete.obs'), error = function(e) NA)
  
  # mean-imputation
  X_impute = mean_impute(X)
  y_h = inverse_transform(X_impute %*% model$theta_ols)
  bias_impute = tryCatch(mean(y_h - y), error = function(e) NA)
  
  # Cmb-lm
  y_h = inverse_transform(par_cmb_lm(R = model$R, theta = model$theta_ols, yty = model$yTy, X_app = X, nThread = 4))
  bias_cmblm = mean(y_h - y)
  return(c(bias_impute = bias_impute, bias_cmblm = bias_cmblm))
}

# Wrapper for missing value titration assay (correlation mode)
injection_assay <- function(model, X, y, p_vec, N_iter)
{
  replicate(n = N_iter, expr = {
    sapply(1:length(p_vec), function(x) test_cmblm_impute(model, inject.na(X, p_vec[x]), y)) # correlation mode
  })
}

# Wrapper for missing value titration assay (bias mode)
injection_assay2 <- function(model, X, y, p_vec, N_iter)
{
  replicate(n = N_iter, expr = {
    sapply(1:length(p_vec), function(x) test_cmblm_impute2(model, inject.na(X, p_vec[x]), y)) # Bias mode
  })
}

############################### Section 1: Model building ###############################

# Read phenotype and DNA methylation data from EWAS Data Hub
pheno = fread('age_methylation_v1.txt', nThread = 4, sep = '\t', nrows = 2)
age = as.numeric(pheno[1,-1])
names(age) = colnames(pheno)[-1]
tissue = as.character(pheno[2,-1])
names(tissue) = colnames(pheno)[-1]
data = fread('horvath_age_methylation_v1.txt', nThread = 4, sep = '\t') # Horvath CpGs
CpGs = data$V1
data = as.matrix(data[,-1])
rownames(data) = CpGs
colnames(data) = colnames(pheno)[-1]
data[1:3, 1:3]
dim(data) # 353 8374
sum(is.na(data)) # 123992

# Split into training/testing set, add intercept and transform dependent variable
set.seed(12431)
train_samples = sample(1:ncol(data), floor(0.9*ncol(data)))
test_samples = (1:ncol(data))[!(1:ncol(data) %in% train_samples)]
train_set = cbind(1, t(data[,train_samples]))
train_y = transform(age[colnames(data)[train_samples]])
test_set = cbind(1, t(data[, test_samples]))
test_y = transform(age[colnames(data)[test_samples]])

######## Model 1: NA-robust statistical inference

df = as.data.frame(cbind(y = train_y, train_set[,-1]))
mod = micemd::mice.par(df, m = 5, maxit = 5, method = 'pmm', seed = 500, nnodes = 4) # Very intensive and slow!
Mu = rowMeans(sapply(1:5, function(x) colMeans(complete(mod, x))))[-1]
Mu = c(1, Mu)
Cov = micombine.cov(mod) # This function applies Rubin's rule to combine results from different imputations
Cov = attr(Cov, 'cov_matrix')[-1, -1]
Cov[1:5, 1:5]; cov(complete(mod, 1))[2:6, 2:6]; cov(complete(mod, 2))[2:6, 2:6]
Cov = rbind(0, cbind(0, Cov))
fit = with(data = mod, exp = lm(as.formula(paste('y ~ ', paste(colnames(df)[-1], collapse = ' + '), sep = ''))))
combFit = pool(fit)
theta_ols = combFit$pooled$estimate
n = nrow(train_set)
R = chol((n-1)*Cov + n*Mu %*% t(Mu))
model1 = list(R = R, theta_ols = theta_ols, yTy = t(train_y) %*% train_y)
data(coef); plot(model1$theta_ols, coef[c('(Intercept)', colnames(train_set)[-1])])
saveRDS(model1, 'model1_MICEpooling.Rds')
rm(R, theta_ols)

######## Model 2: Mean impute training set and train cmb-lm model
train_set_imputed = mean_impute(train_set)
sum(is.na(train_set_imputed)) # 0
QR = qr(train_set_imputed)
Q = qr.Q(QR)
R = qr.R(QR)
theta_ols = backsolve(r = R, x = t(Q) %*% train_y)
model2 = list(R = R, theta_ols = theta_ols, yTy = t(train_y) %*% train_y)
plot(model2$theta_ols, coef[c('(Intercept)', colnames(train_set)[-1])])
saveRDS(model2, 'model2_meanImpute.Rds')
rm(R, theta_ols)

######## Regularization

# glmnet can be confusing to interpret coefficients as it standarizes y and transforms lambda.
# But it has a very efficient implementation of CV. So it is worth employing and translating.
# Analytical solution
# lambda = 2
# I = diag(ncol(train_set)); I[1,1] = 0
# reg_mod3 = solve(t(train_set_imputed) %*% train_set_imputed + lambda*I) %*% t(train_set_imputed) %*% train_y
# penalized
# reg_mod2 = penalized(train_y ~ train_set_imputed[,-1], lambda1 = 0, lambda2 = lambda)
# plot(coef(reg_mod2), reg_mod3); abline(0, 1, lty = 2)
# glmnet
# reg_mod1 = glmnet(x = train_set_imputed[,-1], y = train_y, alpha = 0, standardize = F, intercept = T, thresh = 1e-20)
# n = nrow(train_set_imputed); sd_y <- sqrt(var(train_y)*(n-1)/n)
# cofs = as.vector(coef(x = train_set_imputed[,-1], y = train_y, reg_mod1, s = sd_y*lambda/n, exact = TRUE))
# plot(cofs, reg_mod3); abline(0, 1, lty = 2)
# rm(reg_mod1, reg_mod2, reg_mod3, lambda, I, n, cofs)

# Find lambda.min via 10-fold CV
reg_mod1 = cv.glmnet(x = train_set_imputed[,-1], y = train_y, alpha = 0, standardize = F, intercept = T, thresh = 1e-20) # RIDGE
plot(reg_mod1)
lambda.min = reg_mod1$lambda.min
n = nrow(train_set_imputed); sd_y <- sqrt(var(train_y)*(n-1)/n)
cofs = as.vector(coef(x = train_set_imputed[,-1], y = train_y, reg_mod1, s = lambda.min, exact = TRUE))
lambda = lambda.min*n/sd_y
I = diag(ncol(train_set)); I[1,1] = 0
reg_mod3 = solve(t(train_set_imputed) %*% train_set_imputed + lambda*I) %*% t(train_set_imputed) %*% train_y
plot(cofs, reg_mod3); abline(0, 1, lty = 2)
FIT = glmnet(x = train_set_imputed[,-1], y = train_y, alpha = 0, thresh = 1e-20); plot(FIT)

# Building regularized models
model3 = regularize_L2_model(lambda, model1)
model4 = regularize_L2_model(lambda, model2)
saveRDS(model3, 'model3_MICEpooling_L2reg.Rds')
saveRDS(model4, 'model4_meanImpute_L2reg.Rds')

### Testing influence of lambda
train_set_imputed = mean_impute(train_set)
l_vec = c(0, exp(seq(0, 13, 2)))
cor_vec = numeric(length = length(l_vec))
for(i in 1:length(l_vec))
{
  modeli = regularize_L2_model(l_vec[i], model1)
  y_h = inverse_transform(train_set_imputed %*% modeli$theta_ols)
  cor_vec[i] = cor(y_h, inverse_transform(train_y))
  print(i/length(l_vec))
}
plot(l_vec, cor_vec^2, ylim = c(0,1), type = 'l')

#####################################################################################################

######## Read Models
model1 = readRDS('model1_MICEpooling.Rds')
model2 = readRDS('model2_meanImpute.Rds')
model3 = readRDS('model3_MICEpooling_L2reg.Rds')
model4 = readRDS('model4_meanImpute_L2reg.Rds')

######## Model Comparison (QC)
TH = Reduce(cbind, list(model1$theta_ols, model2$theta_ols, model3$theta_ols, model4$theta_ols))
colnames(TH) = paste('mod', 1:4)
heatmap.2(TH, trace = 'n')
plot(model1$theta_ols, model2$theta_ols); abline(0,1, lty = 2)
plot(abs(model1$R), abs(model2$R)); abline(0,1, lty = 2) # Abs because of chol-multiplicity
plot(model1$theta_ols, model3$theta_ols); abline(0,1, lty = 2)
plot(abs(model1$R), abs(model3$R)); abline(0,1, lty = 2) # Abs because of chol-multiplicity
plot(model2$theta_ols, model4$theta_ols); abline(0,1, lty = 2)
plot(abs(model2$R), abs(model4$R)); abline(0,1, lty = 2) # Abs because of chol-multiplicity
plot(model3$theta_ols, model4$theta_ols); abline(0,1, lty = 2)

######## Testing cmb-model on testing set

# Model1
y_hat_test = par_cmb_lm(R = model1$R, theta = model1$theta_ols, yty = model1$yTy, X_app = test_set, nThread = 4)
age_pred_test = inverse_transform(y_hat_test)
plot(age_pred_test, inverse_transform(test_y)); abline(0, 1, lty = 2)
cor(age_pred_test, inverse_transform(test_y)) # 0.9625407

# Model2
y_hat_test = par_cmb_lm(R = model2$R, theta = model2$theta_ols, yty = model2$yTy, X_app = test_set, nThread = 4)
age_pred_test = inverse_transform(y_hat_test)
plot(age_pred_test, inverse_transform(test_y)); abline(0,1,lty = 2)
cor(age_pred_test, inverse_transform(test_y)) # 0.960805

# Model3
y_hat_test = par_cmb_lm(R = model3$R, theta = model3$theta_ols, yty = model3$yTy, X_app = test_set, nThread = 4)
age_pred_test = inverse_transform(y_hat_test)
plot(age_pred_test, inverse_transform(test_y)); abline(0,1,lty = 2)
cor(age_pred_test, inverse_transform(test_y)) # 0.9533582

# Model4
y_hat_test = par_cmb_lm(R = model4$R, theta = model4$theta_ols, yty = model4$yTy, X_app = test_set, nThread = 4)
age_pred_test = inverse_transform(y_hat_test)
plot(age_pred_test, inverse_transform(test_y)); abline(0, 1, lty = 2)
cor(age_pred_test, inverse_transform(test_y)) # 0.9532407

######## Robustness to missing values (correlation mode)

# Model 1
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341)
RES = injection_assay(model1, test_set, inverse_transform(test_y), p_vec, 10)
RES.m = melt(RES)
RES.m$Var2 = p_vec[RES.m$Var2]
RES.m$Var1 = sapply(strsplit(as.vector(RES.m$Var1), 'cor_'), function(x) x[2])
colnames(RES.m) = c('method', 'pNA', 'N_iter', 'rho')
ggplot(data=RES.m, aes(x=pNA, y=rho^2, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95) + xlim(0,1) + ylim(0,1)

# Model 2
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES2 = injection_assay(model2, test_set, inverse_transform(test_y), p_vec, 10)
RES.m2 = melt(RES2)
RES.m2$Var2 = p_vec[RES.m2$Var2]
RES.m2$Var1 = sapply(strsplit(as.vector(RES.m2$Var1), 'cor_'), function(x) x[2])
colnames(RES.m2) = c('method', 'pNA', 'N_iter', 'rho')
ggplot(data=RES.m2, aes(x=pNA, y=rho^2, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95) + xlim(0,1) + ylim(0,1)

# Model 3
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES3 = injection_assay(model3, test_set, inverse_transform(test_y), p_vec, 10)
RES.m3 = melt(RES3)
RES.m3$Var2 = p_vec[RES.m3$Var2]
RES.m3$Var1 = sapply(strsplit(as.vector(RES.m3$Var1), 'cor_'), function(x) x[2])
colnames(RES.m3) = c('method', 'pNA', 'N_iter', 'rho')
ggplot(data=RES.m3, aes(x=pNA, y=rho^2, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95) + xlim(0,1) + ylim(0,1)

# Model 4
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES4 = injection_assay(model4, test_set, inverse_transform(test_y), p_vec, 10)
RES.m4 = melt(RES4)
RES.m4$Var2 = p_vec[RES.m4$Var2]
RES.m4$Var1 = sapply(strsplit(as.vector(RES.m4$Var1), 'cor_'), function(x) x[2])
colnames(RES.m4) = c('method', 'pNA', 'N_iter', 'rho')
ggplot(data=RES.m4, aes(x=pNA, y=rho^2, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95) + xlim(0,1) + ylim(0,1)

RES.m$model = 'thetaMICE_lambda=0'
RES.m2$model = 'thetaMeanImp_lambda=0'
RES.m3$model = 'thetaMICE_lambdaL2=53.19'
RES.m4$model = 'thetaMeanImp_lambdaL2=53.19'
RES_glob = Reduce(rbind, list(RES.m, RES.m2, RES.m3, RES.m4))
lev = c("cmblm thetaMICE_lambda=0", "cmblm thetaMeanImp_lambda=0", "cmblm thetaMICE_lambdaL2=53.19", 
        "cmblm thetaMeanImp_lambdaL2=53.19", "impute thetaMeanImp_lambdaL2=53.19", "impute thetaMICE_lambdaL2=53.19", 
        "impute thetaMeanImp_lambda=0", "impute thetaMICE_lambda=0")
RES_glob$method_model = factor(paste(RES_glob$method, RES_glob$model), levels = lev)
levels(RES_glob$method_model) = c('trainMICE_testCmblm_lambda=0', 'trainMeanImp_testCmblm_lambda=0',
                                  'trainMICE_testCmblm_lambda=53.19', 'trainMeanImp_testCmblm_lambda=53.19',
                                  'trainMeanImp_testMeanImp_lambda=53.19', 'trainMICE_testMeanImp_lambda=53.19',
                                  'trainMeanImp_testMeanImp_lambda=0', 'trainMICE_testMeanImp_lambda=0')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
(g <- ggplot(data=RES_glob, aes(x=pNA, y=rho^2, col=method_model)) + xlim(0,1) + ylim(0,1) +
  geom_point(alpha = 0.2, size = 0.7) + geom_smooth(method = 'loess', level = 0.95, size = 1, span = 0.4) +
  scale_colour_manual(values=cbPalette))

ggsave(plot = g, 'NA_injection_horvath1.pdf')
ggsave(plot = g, 'NA_injection_horvath1.tiff', dpi = 600)
ggsave(plot = g, 'NA_injection_horvath1_2.tiff', dpi = 600, width = 6, height = 4)

ggplot(data=RES_glob, aes(x=pNA, y=rho^2, col=model)) +
  geom_smooth(method = 'loess', level = 0.95, size = 0.3) + xlim(0,1) + ylim(0,1) +
  facet_wrap( ~ method, ncol = 1) + geom_point(alpha = 0.2, size = 0.7) +
  scale_colour_manual(values=cbPalette)

RES_glob$testing_method = RES_glob$method
RES_glob$training_method = sapply(strsplit(RES_glob$model, split = '_'), function(x) x[1])
RES_glob$training_method = sapply(strsplit(RES_glob$training_method, split = 'theta'), function(x) x[2])
RES_glob$regularized = factor(sapply(strsplit(RES_glob$model, split = '_'), function(x) x[2]))
levels(RES_glob$regularized) = c('F', 'T')
mod = lm(atanh(rho) ~ testing_method*regularized*training_method + pNA, data = RES_glob)
summary(mod)
# Residuals:
# Min       1Q   Median       3Q      Max 
# -0.39527 -0.05688  0.02694  0.08332  0.24494 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                                            2.281097   0.010707 213.052  < 2e-16 ***
#   testing_methodimpute                                  -0.424461   0.012782 -33.208  < 2e-16 ***
#   regularizedT                                          -0.155220   0.012782 -12.144  < 2e-16 ***
#   training_methodMICE                                    0.016237   0.012782   1.270  0.20417    
#   pNA                                                   -1.067648   0.011668 -91.499  < 2e-16 ***
#   testing_methodimpute:regularizedT                      0.200099   0.018077  11.070  < 2e-16 ***
#   testing_methodimpute:training_methodMICE              -0.057627   0.018077  -3.188  0.00146 ** 
#   regularizedT:training_methodMICE                      -0.014533   0.018077  -0.804  0.42154    
#   testing_methodimpute:regularizedT:training_methodMICE -0.006159   0.025564  -0.241  0.80966    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1246 on 1511 degrees of freedom
# Multiple R-squared:  0.8867,	Adjusted R-squared:  0.8861 
# F-statistic:  1478 on 8 and 1511 DF,  p-value: < 2.2e-16

saveRDS(RES_glob, 'NA_injection_horvath_models.Rds')
RES_glob = readRDS('NA_injection_horvath_models.Rds')

#########################################################################

######## Robustness to missing values (bias mode)

# Model 1
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES = injection_assay2(model1, test_set, inverse_transform(test_y), p_vec, 10)
RES.m5 = melt(RES)
RES.m5$Var2 = p_vec[RES.m5$Var2]
RES.m5$Var1 = sapply(strsplit(as.vector(RES.m5$Var1), 'bias_'), function(x) x[2])
colnames(RES.m5) = c('method', 'pNA', 'N_iter', 'bias')
ggplot(data=RES.m5, aes(x=pNA, y=bias, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95)

# Model 2
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES = injection_assay2(model2, test_set, inverse_transform(test_y), p_vec, 10)
RES.m6 = melt(RES)
RES.m6$Var2 = p_vec[RES.m6$Var2]
RES.m6$Var1 = sapply(strsplit(as.vector(RES.m6$Var1), 'bias_'), function(x) x[2])
colnames(RES.m6) = c('method', 'pNA', 'N_iter', 'bias')
ggplot(data=RES.m6, aes(x=pNA, y=bias, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95)

# Model 3
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES = injection_assay2(model3, test_set, inverse_transform(test_y), p_vec, 10)
RES.m7 = melt(RES)
RES.m7$Var2 = p_vec[RES.m7$Var2]
RES.m7$Var1 = sapply(strsplit(as.vector(RES.m7$Var1), 'bias_'), function(x) x[2])
colnames(RES.m7) = c('method', 'pNA', 'N_iter', 'bias')
ggplot(data=RES.m7, aes(x=pNA, y=bias, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95)

# Model 4
p0 = mean(is.na(test_set))
p_vec = seq(p0, 0.99, 0.05)
set.seed(1341) # SAME SEED TO COMPARE BETWEEN MODELS
RES = injection_assay2(model4, test_set, inverse_transform(test_y), p_vec, 10)
RES.m8 = melt(RES)
RES.m8$Var2 = p_vec[RES.m8$Var2]
RES.m8$Var1 = sapply(strsplit(as.vector(RES.m8$Var1), 'bias_'), function(x) x[2])
colnames(RES.m8) = c('method', 'pNA', 'N_iter', 'bias')
ggplot(data=RES.m8, aes(x=pNA, y=bias, col=method)) + 
  geom_smooth(method = 'loess', level = 0.95)

RES.m5$model = 'thetaMICE_lambda=0'
RES.m6$model = 'thetaMeanImp_lambda=0'
RES.m7$model = 'thetaMICE_lambdaL2=53.19'
RES.m8$model = 'thetaMeanImp_lambdaL2=53.19'
RES_glob = Reduce(rbind, list(RES.m5, RES.m6, RES.m7, RES.m8))
lev = c("cmblm thetaMICE_lambda=0", "cmblm thetaMeanImp_lambda=0", "cmblm thetaMICE_lambdaL2=53.19", 
        "cmblm thetaMeanImp_lambdaL2=53.19", "impute thetaMeanImp_lambdaL2=53.19", "impute thetaMICE_lambdaL2=53.19", 
        "impute thetaMeanImp_lambda=0", "impute thetaMICE_lambda=0")
RES_glob$method_model = factor(paste(RES_glob$method, RES_glob$model), levels = lev)
levels(RES_glob$method_model) = c('trainMICE_testCmblm_lambda=0', 'trainMeanImp_testCmblm_lambda=0',
                                  'trainMICE_testCmblm_lambda=53.19', 'trainMeanImp_testCmblm_lambda=53.19',
                                  'trainMeanImp_testMeanImp_lambda=53.19', 'trainMICE_testMeanImp_lambda=53.19',
                                  'trainMeanImp_testMeanImp_lambda=0', 'trainMICE_testMeanImp_lambda=0')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
(g <- ggplot(data=RES_glob, aes(x=pNA, y=bias, col=method_model)) + geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.2, size = 0.7) + geom_smooth(method = 'loess', level = 0.95, size = 1, span = 0.4, alpha = 0.3) +
    scale_colour_manual(values=cbPalette))

ggsave(plot = g, 'NA_injection_horvath1_bias.pdf')
ggsave(plot = g, 'NA_injection_horvath1_bias.tiff', dpi = 600)
ggsave(plot = g, 'NA_injection_horvath1_2_bias.tiff', dpi = 600, width = 6, height = 4)
saveRDS(RES_glob, 'NA_injection_horvath_models_bias.Rds')
RES_glob = readRDS('NA_injection_horvath_models_bias.Rds')

RES_glob$testing_method = RES_glob$method
RES_glob$training_method = sapply(strsplit(RES_glob$model, split = '_'), function(x) x[1])
RES_glob$training_method = sapply(strsplit(RES_glob$training_method, split = 'theta'), function(x) x[2])
RES_glob$regularized = factor(sapply(strsplit(RES_glob$model, split = '_'), function(x) x[2]))
levels(RES_glob$regularized) = c('F', 'T')
mod = lm(bias ~ testing_method*regularized*training_method + pNA, data = RES_glob)
summary(mod)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.7819 -0.1676 -0.0018  0.1904  4.4278 
# 
#   Coefficients:
#                                                         Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                                            0.52040    0.04451  11.691  < 2e-16 ***
#   testing_methodimpute                                  -1.00516    0.05314 -18.915  < 2e-16 ***
#   regularizedT                                          -0.28290    0.05314  -5.324 1.17e-07 ***
#   training_methodMICE                                    0.03773    0.05314   0.710    0.478    
#   pNA                                                   -0.75685    0.04851 -15.602  < 2e-16 ***
#   testing_methodimpute:regularizedT                      0.30620    0.07515   4.074 4.85e-05 ***
#   testing_methodimpute:training_methodMICE               1.19486    0.07515  15.900  < 2e-16 ***
#   regularizedT:training_methodMICE                      -0.02785    0.07515  -0.371    0.711    
#   testing_methodimpute:regularizedT:training_methodMICE -1.06301    0.10628 -10.002  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5179 on 1511 degrees of freedom
# Multiple R-squared:  0.4859,	Adjusted R-squared:  0.4832 
# F-statistic: 178.5 on 8 and 1511 DF,  p-value: < 2.2e-16

#########################################################################

sessionInfo()
# version 4.1.2 (2021-11-01)
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
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] fastmatrix_0.3-8196                                gplots_3.1.1                                      
# [3] penalized_0.9-51                                   survival_3.2-13                                   
# [5] glmnet_4.1-3                                       Matrix_1.3-4                                      
# [7] micemd_1.8.0                                       miceadds_3.11-6                                   
# [9] forcats_0.5.1                                      stringr_1.4.0                                     
# [11] dplyr_1.0.7                                        purrr_0.3.4                                       
# [13] readr_2.1.0                                        tidyr_1.1.4                                       
# [15] tibble_3.1.6                                       tidyverse_1.3.1                                   
# [17] mice_3.13.0                                        wateRmelon_1.36.0                                 
# [19] illuminaio_0.34.0                                  IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0
# [21] ROC_1.68.1                                         lumi_2.44.0                                       
# [23] methylumi_2.38.0                                   minfi_1.38.0                                      
# [25] bumphunter_1.34.0                                  locfit_1.5-9.4                                    
# [27] iterators_1.0.13                                   foreach_1.5.1                                     
# [29] Biostrings_2.60.2                                  XVector_0.32.0                                    
# [31] SummarizedExperiment_1.22.0                        MatrixGenerics_1.4.3                              
# [33] FDb.InfiniumMethylation.hg19_2.2.0                 org.Hs.eg.db_3.13.0                               
# [35] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2            GenomicFeatures_1.44.2                            
# [37] AnnotationDbi_1.54.1                               GenomicRanges_1.44.0                              
# [39] GenomeInfoDb_1.28.4                                IRanges_2.26.0                                    
# [41] S4Vectors_0.30.2                                   scales_1.1.1                                      
# [43] matrixStats_0.61.0                                 limma_3.48.3                                      
# [45] Biobase_2.52.0                                     BiocGenerics_0.38.0                               
# [47] reshape2_1.4.4                                     ggplot2_3.3.5                                     
# [49] RColorBrewer_1.1-2                                 MASS_7.3-54                                       
# [51] data.table_1.14.2                                 
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                tidyselect_1.1.1          lme4_1.1-27.1             RSQLite_2.2.8            
# [5] grid_4.1.2                BiocParallel_1.26.2       munsell_0.5.0             codetools_0.2-18         
# [9] preprocessCore_1.54.0     nleqslv_3.3.2             withr_2.4.2               colorspace_2.0-2         
# [13] filelock_1.0.2            knitr_1.36                rstudioapi_0.13           GenomeInfoDbData_1.2.6   
# [17] bit64_4.0.5               rhdf5_2.36.0              vctrs_0.3.8               generics_0.1.1           
# [21] xfun_0.28                 BiocFileCache_2.0.0       R6_2.5.1                  bitops_1.0-7             
# [25] rhdf5filters_1.4.0        cachem_1.0.6              reshape_0.8.8             DelayedArray_0.18.0      
# [29] assertthat_0.2.1          BiocIO_1.2.0              gtable_0.3.0              affy_1.70.0              
# [33] rlang_0.4.12              genefilter_1.74.1         splines_4.1.2             rtracklayer_1.52.1       
# [37] GEOquery_2.60.0           broom_0.7.10              BiocManager_1.30.16       yaml_2.2.1               
# [41] abind_1.4-5               modelr_0.1.8              backports_1.3.0           tools_4.1.2              
# [45] nor1mix_1.3-0             affyio_1.62.0             ellipsis_0.3.2            siggenes_1.66.0          
# [49] Rcpp_1.0.7                plyr_1.8.6                sparseMatrixStats_1.4.2   progress_1.2.2           
# [53] zlibbioc_1.38.0           RCurl_1.98-1.5            prettyunits_1.1.1         openssl_1.4.5            
# [57] haven_2.4.3               fs_1.5.0                  magrittr_2.0.1            reprex_2.0.1             
# [61] mvtnorm_1.1-3             hms_1.1.1                 xtable_1.8-4              XML_3.99-0.8             
# [65] mclust_5.4.8              readxl_1.3.1              shape_1.4.6               compiler_4.1.2           
# [69] biomaRt_2.48.3            KernSmooth_2.23-20        crayon_1.4.2              minqa_1.2.4              
# [73] mgcv_1.8-38               tzdb_0.2.0                lubridate_1.8.0           DBI_1.1.1                
# [77] dbplyr_2.1.1              mvmeta_1.0.3              rappdirs_0.3.3            boot_1.3-28              
# [81] mitools_2.4               cli_3.1.0                 quadprog_1.5-8            pkgconfig_2.0.3          
# [85] GenomicAlignments_1.28.0  xml2_1.3.2                annotate_1.70.0           rngtools_1.5.2           
# [89] multtest_2.48.0           beanplot_1.2              rvest_1.0.2               doRNG_1.8.2              
# [93] scrime_1.3.5              digest_0.6.28             base64_2.0                cellranger_1.1.0         
# [97] DelayedMatrixStats_1.14.3 restfulr_0.0.13           curl_4.3.2                gtools_3.9.2             
# [101] Rsamtools_2.8.0           jomo_2.7-2                rjson_0.2.20              nloptr_1.2.2.3           
# [105] lifecycle_1.0.1           nlme_3.1-153              jsonlite_1.7.2            Rhdf5lib_1.14.2          
# [109] askpass_1.1               fansi_0.5.0               pillar_1.6.4              lattice_0.20-45          
# [113] KEGGREST_1.32.0           fastmap_1.1.0             httr_1.4.2                glue_1.5.0               
# [117] png_0.1-7                 bit_4.0.4                 stringi_1.7.5             HDF5Array_1.20.0         
# [121] blob_1.2.2                caTools_1.18.2            mixmeta_1.2.0             memoise_2.0.0 