## simulation for empirically motivated simulations

## libraries
library(simex)
library(data.table)
library(nonprobsvy)
library(simputation)

## classification matrices from linear and transformer (test data)
## true (rows) / predicted (columns)
lim_mat <- matrix(c(
  153, 50, 17, 3, 6, 0, 0, 1, 2,
  40, 1094, 67, 29, 24, 0, 16, 1, 2,
  12, 103, 292, 23, 34, 0, 30, 6, 3,
  3, 35, 15, 181, 10, 0, 4, 3, 14,
  5, 23, 29, 7, 323, 0, 4, 0, 8,
  0, 1, 0, 0, 0, 3, 3, 0, 0,
  4, 12, 10, 5, 4, 0, 430, 21, 32,
  0, 1, 6, 5, 1, 0, 49, 183, 13,
  1, 1, 11, 9, 12, 0, 27, 5, 154
), nrow = 9, byrow = TRUE)

lin_mat <- lim_mat[-6, -6]

tra_mat <- matrix(c(
  175, 36, 10, 4, 5, 0, 1, 0, 1,
  35, 1108, 70, 25, 26, 0, 8, 1, 0,
  11, 66, 345, 14, 25, 0, 28, 7, 7,
  3, 28, 16, 186, 6, 0, 1, 3, 22,
  4, 15, 25, 7, 342, 0, 1, 0, 5,
  1, 0, 0, 0, 0, 4, 0, 0, 2,
  3, 9, 10, 1, 1, 0, 457, 18, 19,
  0, 1, 3, 5, 1, 0, 45, 197, 6,
  1, 0, 6, 12, 13, 5, 19, 4, 160
), nrow = 9, byrow = TRUE)

tra_mat <- tra_mat[-6, -6]

## genreate matrices for measurement error (pred in rows, true in columns)
lin_mat_errors <- prop.table(t(lin_mat), margin = 2)
lin_mat_errors[lin_mat_errors < 0.01] <- 0.01
lin_mat_errors <- lin_mat_errors/colSums(lin_mat_errors)

lin_mat_correction <- prop.table(lin_mat, margin = 2)
lin_mat_correction[lin_mat_correction < 0.01] <- 0.01
lin_mat_correction <- lin_mat_correction/colSums(lin_mat_correction)


tra_mat_errors <- prop.table(t(tra_mat), margin = 2)
tra_mat_errors[tra_mat_errors < 0.01] <- 0.01
tra_mat_errors <- tra_mat_errors/colSums(tra_mat_errors)

tra_mat_correction <- prop.table(tra_mat, margin = 2)
tra_mat_correction[tra_mat_correction < 0.01] <- 0.01
tra_mat_correction <- tra_mat_correction/colSums(tra_mat_correction)


## population data from the JVS
pop_data <- fread("data-raw/sim-empirical-pop-data.csv")
pop_data[, size := factor(size, c("Small", "Medium", "Large"))]
## substractig more for 1:3 just for the simulation
#pop_data[occup_code %in% 1:3 & prob_y, prob_y:=ifelse(prob_y+0.05 > 1, 1, prob_y+0.05)]
pop_data_expanded <- pop_data[rep(pop_data[, .I], vac_total), .(private, size, occup_code, prob_y, prob_cbop)]
pop_data_expanded[, occup_code:=as.factor(occup_code)]

colnames(tra_mat_errors) <- colnames(tra_mat_correction) <- colnames(lin_mat_errors) <- 
  colnames(lin_mat_correction) <- levels(pop_data_expanded$occup_code)

## add mis-classification
set.seed(2024)
pop_data_expanded[, y:=rbinom(.N, 1, prob_y)]
occup_lin <- simex::misclass(pop_data_expanded[, .(occup_code)], list(occup_code = lin_mat_errors), k = 1)
occup_tra <- simex::misclass(pop_data_expanded[, .(occup_code)], list(occup_code = tra_mat_errors), k = 1)
pop_data_expanded[, occup_code_lin := occup_lin$occup_code]
pop_data_expanded[, occup_code_tra := occup_tra$occup_code]

y_true <- mean(pop_data_expanded$y)

vcd::assocstats(xtabs(~y+ size, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ private, data = pop_data_expanded))

vcd::assocstats(xtabs(~y+ occup_code, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ occup_code_lin, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ occup_code_tra, data = pop_data_expanded))


## simulation
R <- 250
B_simex <- 50
results <- matrix(0, nrow = R, ncol = 19, 
                  dimnames = list(NULL, 
                                  c("naive", 
                                    "ipw", "mi", "dr",
                                    "ipw_lin", "mi_lin", "dr_lin",
                                    "ipw_tra", "mi_tra", "dr_tra",
                                    "ipw_lin_cor", "mi_lin_cor", "dr_lin_cor",
                                    "ipw_tra_cor", "mi_tra_cor", "dr_tra_cor",
                                    "mi_li", "mi_lo", "mi_q")))

results_ci <- matrix(0, nrow = R, ncol = 18, 
                  dimnames = list(NULL, 
                                  c("ipw", "mi", "dr",
                                    "ipw_lin", "mi_lin", "dr_lin",
                                    "ipw_tra", "mi_tra", "dr_tra",
                                    "ipw_lin_cor", "mi_lin_cor", "dr_lin_cor",
                                    "ipw_tra_cor", "mi_tra_cor", "dr_tra_cor",
                                    "mi_li", "mi_lo", "mi_q")))

for (r in 1:R) {
  set.seed(r)
  print(r)
  pop_data_expanded[, cbop:=rbinom(.N, 1, prob_cbop)] 
  cbop_data <- pop_data_expanded[cbop == 1]
  prob_sample <- pop_data_expanded[sample(1:.N, 1000)]
  prob_sample[, weights:=nrow(pop_data_expanded)/1000]
  prob_sample[, ":="(occup_code_lin=NULL, occup_code_tra=NULL)]
  
  ## imputation based mis-classification matrix 
  occup_lin <- simex::misclass(prob_sample[, .(occup_code)], list(occup_code = lin_mat_correction), k = 1)
  occup_tra <- simex::misclass(prob_sample[, .(occup_code)], list(occup_code = lin_mat_correction), k = 1)
  prob_sample_lin <- copy(prob_sample)
  prob_sample_lin[, occup_code := occup_lin$occup_code]
  prob_sample_tra <- copy(prob_sample)
  prob_sample_tra[, occup_code := occup_tra$occup_code]
  
  prob_sample_svy <- svydesign(ids=~1, data=prob_sample, weights= ~ weights)
  prob_sample_lin_svy <- svydesign(ids=~1, data=prob_sample_lin, weights= ~ weights)
  prob_sample_tra_svy <- svydesign(ids=~1, data=prob_sample_tra, weights= ~ weights)
  
  ## original
  res_ipw_cal <- nonprob(data = cbop_data, 
                         selection = ~ private + size + occup_code, 
                         svydesign = prob_sample_svy, 
                         target = ~y, 
                         #, control_selection = controlSel(est_method_sel = "gee")
                         )
  
  res_mi <- nonprob(data = cbop_data, 
                    outcome = y ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial")
  
  res_dr <- nonprob(data = cbop_data, 
                    outcome = y ~ private + size + occup_code, 
                    selection = ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial",
                    #, control_selection = controlSel(est_method_sel = "gee")
                    )
  
  ## with mis-class (linear)
  res_ipw_cal_lin <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_lin)], 
                         selection = ~ private + size + occup_code, 
                         svydesign = prob_sample_svy, 
                         target = ~y
                         #, control_selection = controlSel(est_method_sel = "gee")
                         )
  
  res_mi_lin <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_lin)], 
                    outcome = y ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_lin <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_lin)], 
                    outcome = y ~ private + size + occup_code, 
                    selection = ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial"
                    #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  ## with mis-class (transformer)
  res_ipw_cal_tra <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_tra)], 
                             selection = ~ private + size + occup_code, 
                             svydesign = prob_sample_svy, 
                             target = ~y
                             #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  res_mi_tra <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_tra)], 
                        outcome = y ~ private + size + occup_code, 
                        svydesign = prob_sample_svy,
                        method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_tra <- nonprob(data = cbop_data[, .(y,private,size,occup_code=occup_code_tra)], 
                        outcome = y ~ private + size + occup_code, 
                        selection = ~ private + size + occup_code, 
                        svydesign = prob_sample_svy,
                        method_outcome = "glm", family_outcome = "binomial"
                        #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  
  ## correct for measurement error using known matrix (generate data based on matrix)
  #occup_lin_corr <- simex::misclass(cbop_data[, .(occup_code_lin)], list(occup_code_lin = lin_mat_for_corr), k = 1)
  #occup_tra_corr <- simex::misclass(cbop_data[, .(occup_code_tra)], list(occup_code_tra = tra_mat_for_corr), k = 1)
  #cbop_data[, occup_code_lin_corr := occup_lin_corr$occup_code_lin]
  #cbop_data[, occup_code_tra_corr := occup_tra_corr$occup_code_tra]
  
  ## estimators after correction (linear)
  res_ipw_cal_lin_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_lin_corr)], 
                             selection = ~ private + size + occup_code, 
                             svydesign = prob_sample_lin_svy, 
                             target = ~y
                             #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  res_mi_lin_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_lin_corr)], 
                        outcome = y ~ private + size + occup_code, 
                        svydesign = prob_sample_lin_svy,
                        method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_lin_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_lin_corr)], 
                        outcome = y ~ private + size + occup_code, 
                        selection = ~ private + size + occup_code, 
                        svydesign = prob_sample_lin_svy,
                        method_outcome = "glm", family_outcome = "binomial"
                        #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  ## estimators after correction (transformer)
  res_ipw_cal_tra_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_tra_corr)], 
                             selection = ~ private + size + occup_code, 
                             svydesign = prob_sample_tra_svy, 
                             target = ~y
                             #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  res_mi_tra_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_tra_corr)], 
                        outcome = y ~ private + size + occup_code, 
                        svydesign = prob_sample_tra_svy,
                        method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_tra_cor <- nonprob(data = cbop_data,#[, .(y,private,size,occup_code=occup_code_tra_corr)], 
                        outcome = y ~ private + size + occup_code, 
                        selection = ~ private + size + occup_code, 
                        svydesign = prob_sample_tra_svy,
                        method_outcome = "glm", family_outcome = "binomial"
                        #, control_selection = controlSel(est_method_sel = "gee")
  )
  
  
  # ## model parameters corrected by MC-SIMEX approach

  mi_lin_m1 <- glm(data = cbop_data[, .(y,private,size,occup_code=occup_code_lin)], 
  formula = y ~ private + size + occup_code, family= "binomial")
  
  res_mi_lin_simex_l <- simex::mcsimex(model = mi_lin_m1, 
                                       mc.matrix = list(occup_code = lin_mat_correction), 
                                       SIMEXvariable = c("occup_code"), 
                                       fitting.method	= "linear", 
                                       asymptotic = FALSE, 
                                       #lambda = c(0.005, 0.01, 0.05, 0.1, 0.2, 0.5), 
                                       jackknife.estimation	= FALSE, B = B_simex)

  res_mi_lin_simex_ll <- refit(object = res_mi_lin_simex_l, 
                              fitting.method = "loglinear", 
                              asymptotic = FALSE, 
                              jackknife.estimation	= FALSE)
  
  res_mi_lin_simex_q <- refit(object = res_mi_lin_simex_l, 
                              fitting.method = "quadratic", 
                              asymptotic = FALSE, 
                              jackknife.estimation	= FALSE)
  
  prob_sample_svy <- update(prob_sample_svy,
                            y_li = predict(res_mi_lin_simex_l, prob_sample_svy$variables, type = "response"),
                            y_lo = predict(res_mi_lin_simex_ll, prob_sample_svy$variables, type = "response"),
                            y_qu = predict(res_mi_lin_simex_q, prob_sample_svy$variables, type = "response"))
  
  res <- svymean(~y_li + y_lo + y_qu, design = prob_sample_svy)
  res_ci <- confint(res)
  
  results[r, 1] <- mean(cbop_data$y)
  
  results[r, 2] <- res_ipw_cal$output$mean
  results[r, 3] <- res_mi$output$mean
  results[r, 4] <- res_dr$output$mean
  results[r, 5] <- res_ipw_cal_lin$output$mean
  results[r, 6] <- res_mi_lin$output$mean
  results[r, 7] <- res_dr_lin$output$mean
  results[r, 8] <- res_ipw_cal_tra$output$mean
  results[r, 9] <- res_mi_tra$output$mean
  results[r, 10] <- res_dr_tra$output$mean
  results[r, 11] <- res_ipw_cal_lin_cor$output$mean
  results[r, 12] <- res_mi_lin_cor$output$mean
  results[r, 13] <- res_dr_lin_cor$output$mean
  results[r, 14] <- res_ipw_cal_tra_cor$output$mean
  results[r, 15] <- res_mi_tra_cor$output$mean
  results[r, 16] <- res_dr_tra_cor$output$mean
  results[r, 17] <- res[1]
  results[r, 18] <- res[2]
  results[r, 19] <- res[3]
  
  
  results_ci[r, 1] <- res_ipw_cal$confidence_interval[1] < y_true & res_ipw_cal$confidence_interval[2] > y_true
  results_ci[r, 2] <- res_mi$confidence_interval[1] < y_true & res_mi$confidence_interval[2] > y_true
  results_ci[r, 3] <- res_dr$confidence_interval[1] < y_true & res_dr$confidence_interval[2] > y_true
  results_ci[r, 4] <- res_ipw_cal_lin$confidence_interval[1] < y_true & res_ipw_cal_lin$confidence_interval[2] > y_true
  results_ci[r, 5] <- res_mi_lin$confidence_interval[1] < y_true & res_mi_lin$confidence_interval[2] > y_true
  results_ci[r, 6] <- res_dr_lin$confidence_interval[1] < y_true & res_dr_lin$confidence_interval[2] > y_true
  results_ci[r, 7] <- res_ipw_cal_tra$confidence_interval[1] < y_true & res_ipw_cal_tra$confidence_interval[2] > y_true
  results_ci[r, 8] <- res_mi_tra$confidence_interval[1] < y_true & res_mi_tra$confidence_interval[2] > y_true
  results_ci[r, 9] <- res_dr_tra$confidence_interval[1] < y_true & res_dr_tra$confidence_interval[2] > y_true
  results_ci[r, 10] <- res_ipw_cal_lin_cor$confidence_interval[1] < y_true & res_ipw_cal_lin_cor$confidence_interval[2] > y_true
  results_ci[r, 11] <- res_mi_lin_cor$confidence_interval[1] < y_true & res_mi_lin_cor$confidence_interval[2] > y_true
  results_ci[r, 12] <- res_dr_lin_cor$confidence_interval[1] < y_true & res_dr_lin_cor$confidence_interval[2] > y_true
  results_ci[r, 13] <- res_ipw_cal_tra_cor$confidence_interval[1] < y_true & res_ipw_cal_tra_cor$confidence_interval[2] > y_true
  results_ci[r, 14] <- res_mi_tra_cor$confidence_interval[1] < y_true & res_mi_tra_cor$confidence_interval[2] > y_true
  results_ci[r, 15] <- res_dr_tra_cor$confidence_interval[1] < y_true & res_dr_tra_cor$confidence_interval[2] > y_true
  results_ci[r, 16] <- res_ci[1,1] < y_true & res_ci[1,2] > y_true
  results_ci[r, 17] <- res_ci[2,1] < y_true & res_ci[2,2] > y_true
  results_ci[r, 18] <- res_ci[3,1] < y_true & res_ci[3,2] > y_true
  
}

boxplot(results[1:(r-1),]- y_true)
abline(h=0,col = "red")

colMeans(results_ci[1:(r-1),])


## 
bias <- colMeans(results[1:(r-1),]) - y_true
var <- apply(results[1:(r-1),], 2, var)
rmse <- sqrt(bias^2+var)
data.frame(bias,var, rmse)
round(bias/y_true*100,2)


