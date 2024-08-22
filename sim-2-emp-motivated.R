## simulation for empirically motivated simulations

## libraries
library(simex)
library(data.table)
library(nonprobsvy)
library(simputation)

## classification matrices from linear and transformer (test data)
lin_vals <- c(
  65.9, 21.6, 7.3, 1.3, 2.6, 0.0, 0.0, 0.4, 0.9,
  3.1, 85.9, 5.3, 2.3, 1.9, 0.0, 1.3, 0.1, 0.2,
  2.4, 20.5, 58.1, 4.6, 6.8, 0.0, 6.0, 1.2, 0.6,
  1.1, 13.2, 5.7, 68.3, 3.8, 0.0, 1.5, 1.1, 5.3,
  1.3, 5.8, 7.3, 1.8, 81.0, 0.0, 1.0, 0.0, 2.0,
  0.0, 14.3, 0.0, 0.0, 0.0, 42.9, 42.9, 0.0, 0.0,
  0.8, 2.3, 1.9, 1.0, 0.8, 0.0, 83.0, 4.1, 6.2,
  0.0, 0.4, 2.3, 1.9, 0.4, 0.0, 19.0, 70.9, 5.0,
  0.5, 0.5, 5.0, 4.1, 5.5, 0.0, 12.3, 2.3, 70.0
)

lin_mat <- matrix(lin_vals, nrow = 9, ncol = 9, byrow = TRUE)
lin_mat <- t(lin_mat[-6, -6]/rowSums(lin_mat[-6, -6]))

tra_vals <- values <- c(
  75.4, 15.5, 4.3, 1.7, 2.2, 0.0, 0.4, 0.0, 0.4, 
  2.7, 87.0, 5.5, 2.0, 2.0, 0.0, 0.6, 0.1, 0.0, 
  2.2, 13.1, 68.6, 2.8, 5.0, 0.0, 5.6, 1.4, 1.4, 
  1.1, 10.6, 6.0, 70.2, 2.3, 0.0, 0.4, 1.1, 8.3, 
  1.0, 3.8, 6.3, 1.8, 85.7, 0.0, 0.3, 0.0, 1.3, 
  14.3, 0.0, 0.0, 0.0, 0.0, 57.1, 0.0, 0.0, 28.6, 
  0.6, 1.7, 1.9, 0.2, 0.2, 0.0, 88.2, 3.5, 3.7, 
  0.0, 0.4, 1.2, 1.9, 0.4, 0.0, 17.4, 76.4, 2.3, 
  0.5, 0.0, 2.7, 5.5, 5.9, 2.3, 8.6, 1.8, 72.7 
)

tra_mat <- matrix(tra_vals, nrow = 9, ncol = 9, byrow = TRUE)
tra_mat <- t(tra_mat[-6, -6]/rowSums(tra_mat[-6, -6]))


## population data from the JVS

pop_data <- fread("data-raw/sim-empirical-pop-data.csv")
pop_data[, size := factor(size, c("Small", "Medium", "Large"))]
pop_data_expanded <- pop_data[rep(pop_data[, .I], vac_total), .(private, size, occup_code, prob_y, prob_cbop)]
pop_data_expanded[, occup_code:=as.factor(occup_code)]

colnames(tra_mat) <- colnames(lin_mat) <- levels(pop_data_expanded$occup_code)

## add mis-classification
set.seed(2024)
pop_data_expanded[, y:=rbinom(.N, 1, prob_y)]
occup_lin <- simex::misclass(pop_data_expanded[, .(occup_code)], list(occup_code = lin_mat), k = 1)
occup_tra <- simex::misclass(pop_data_expanded[, .(occup_code)], list(occup_code = tra_mat), k = 1)
pop_data_expanded[, occup_code_lin := occup_lin$occup_code]
pop_data_expanded[, occup_code_tra := occup_tra$occup_code]

## for correction
lin_mat_for_corr <- prop.table(margin = 2, xtabs(~occup_code+occup_code_lin, pop_data_expanded))
lin_mat_for_corr[lin_mat_for_corr < 0.01] <- 0.01
for (i in 1:ncol(lin_mat_for_corr)) lin_mat_for_corr[, i] <- lin_mat_for_corr[,i]/sum(lin_mat_for_corr[,i])
tra_mat_for_corr <- prop.table(margin = 2, xtabs(~occup_code+occup_code_tra, pop_data_expanded))
tra_mat_for_corr[tra_mat_for_corr < 0.01] <- 0.01
for (i in 1:ncol(tra_mat_for_corr)) tra_mat_for_corr[, i] <- tra_mat_for_corr[,i]/sum(tra_mat_for_corr[,i])


y_true <- mean(pop_data_expanded$y)

vcd::assocstats(xtabs(~y+ size, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ private, data = pop_data_expanded))

vcd::assocstats(xtabs(~y+ occup_code, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ occup_code_lin, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ occup_code_tra, data = pop_data_expanded))


## simulation
R <- 50
B_simex <- 50
results <- matrix(0, nrow = R, ncol = 18, 
                  dimnames = list(NULL, 
                                  c("naive", 
                                    "ipw", "mi", "dr",
                                    "ipw_lin", "mi_lin", "dr_lin",
                                    "ipw_tra", "mi_tra", "dr_tra",
                                    "ipw_lin_cor", "mi_lin_cor", "dr_lin_cor",
                                    "ipw_tra_cor", "mi_tra_cor", "dr_tra_cor",
                                    "mi_li", "mi_q")))

results_ci <- matrix(0, nrow = R, ncol = 17, 
                  dimnames = list(NULL, 
                                  c("ipw", "mi", "dr",
                                    "ipw_lin", "mi_lin", "dr_lin",
                                    "ipw_tra", "mi_tra", "dr_tra",
                                    "ipw_lin_cor", "mi_lin_cor", "dr_lin_cor",
                                    "ipw_tra_cor", "mi_tra_cor", "dr_tra_cor",
                                    "mi_li", "mi_q")))

for (r in 1:R) {
  set.seed(r)
  if (r %% 50 == 0) print(r)
  pop_data_expanded[, cbop:=rbinom(.N, 1, prob_cbop)] 
  cbop_data <- pop_data_expanded[cbop == 1]
  prob_sample <- pop_data_expanded[sample(1:.N, 1000)]
  prob_sample[, weights:=nrow(pop_data_expanded)/1000]
  prob_sample[, ":="(occup_code_lin=NULL, occup_code_tra=NULL)]
  
  ## imputation based mis-classification matrix 
  occup_lin <- simex::misclass(prob_sample[, .(occup_code)], list(occup_code = lin_mat), k = 1)
  occup_tra <- simex::misclass(prob_sample[, .(occup_code)], list(occup_code = tra_mat), k = 1)
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
                                       mc.matrix = list(occup_code = lin_mat_for_corr), 
                                       SIMEXvariable = c("occup_code"), fitting.method	= "linear", asymptotic = FALSE, 
                                       lambda = c(0.01, 0.05, 0.1), 
                                       jackknife.estimation	= FALSE, B = B_simex)

  res_mi_lin_simex_q <- refit(object = res_mi_lin_simex_l, 
                              fitting.method = "quadratic", 
                              asymptotic = FALSE, 
                              jackknife.estimation	= FALSE)
  
  prob_sample_svy <- update(prob_sample_svy,
                    y_li = predict(res_mi_lin_simex_l, prob_sample_svy$variables, type = "response"),
                    y_qu = predict(res_mi_lin_simex_q, prob_sample_svy$variables, type = "response"))
  
  res <- svymean(~y_li + y_qu, design = prob_sample_svy)
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
  
}

boxplot(results[1:(r-1),]- y_true)
abline(h=0,col = "red")

colMeans(results_ci)


## 
bias <- colMeans(results) - y_true
var <- apply(results, 2, var)
rmse <- sqrt(bias^2+var)
data.frame(bias,var, rmse)
round(bias/y_true*100,2)


