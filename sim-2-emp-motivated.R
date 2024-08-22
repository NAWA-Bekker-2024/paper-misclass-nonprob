## simulation for empirically motivated simulations

## libraries
library(simex)
library(data.table)
library(nonprobsvy)

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
lin_mat <- lin_mat/rowSums(lin_mat)
lin_mat <- lin_mat[-6, -6]

rownames(lin_mat) <- colnames(lin_mat) <- c("1", "2", "3", "4", "5", "7", "8", "9")

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
tra_mat <- tra_mat/rowSums(tra_mat)
tra_mat <- tra_mat[-6, -6]
rownames(tra_mat) <- colnames(tra_mat) <- c("1", "2", "3", "4", "5", "7", "8", "9")

## population data from the JVS

pop_data <- fread("data-raw/sim-empirical-pop-data.csv")
pop_data$size <- factor(pop_data$size, c("Small", "Medium", "Large"), ordered = T)

## generate target variable, say working in on a single shift
### more public 
### decreases with size
### for specific occupations: it will be significantly lower: 4, 5, 7, 8, 9

pop_data[, pr_y := plogis(13 - 2*(occup_code) - 2*(size=="Large") + 1.5*(private == "Public"))]
pop_data[, pr_peo := vac_peo/vac_total]
pop_data_expanded <- pop_data[rep(pop_data[, .I], vac_total),.(private, size, occup_code, pr_y, pr_peo)]
pop_data_expanded[, occup_code:=as.factor(occup_code)]

set.seed(2024)
pop_data_expanded[, y:=rbinom(.N, 1, pr_y)]

vcd::assocstats(xtabs(~y+ occup_code, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ size, data = pop_data_expanded))
vcd::assocstats(xtabs(~y+ private, data = pop_data_expanded))

m1 <- glm(y~ occup_code+size+private,data=pop_data_expanded, family = binomial())

## simulation
R <- 500
results <- matrix(0, nrow = R, ncol = 10, 
                  dimnames = list(NULL, 
                                  c("naive", 
                                    "ipw", "mi", "dr",
                                    "ipw_lin", "mi_lin", "dr_lin",
                                    "ipw_tra", "mi_tra", "dr_tra")))

for (r in 1:R) {
  set.seed(r)
  pop_data_expanded[, cbop:=rbinom(.N, 1, pr_peo)] 
  cbop_data <- pop_data_expanded[cbop == 1]
  prob_sample <- pop_data_expanded[sample(1:.N, 1000)]
  prob_sample[, weights:=nrow(pop_data_expanded)/1000]
  prob_sample_svy <- svydesign(ids=~1, data=prob_sample, weights= ~ weights)
  
  cbop_copy_lin <- copy(cbop_data)
  cbop_copy_tra <- copy(cbop_data)
  
  occup_lin <- simex::misclass(cbop_copy_lin[, .(occup_code)], list(occup_code = lin_mat), k = 1)
  occup_tra <- simex::misclass(cbop_copy_tra[, .(occup_code)], list(occup_code = tra_mat), k = 1)
  
  cbop_copy_lin[, occup_code := occup_lin$occup_code]
  cbop_copy_tra[, occup_code := occup_tra$occup_code]
  
  ## original
  res_ipw_cal <- nonprob(data = cbop_data, 
                         selection = ~ private + size + occup_code, 
                         svydesign = prob_sample_svy, 
                         target = ~y, 
                         , control_selection = controlSel(est_method_sel = "gee")
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
                    , control_selection = controlSel(est_method_sel = "gee")
                    )
  
  ## with mis-class (linear)
  res_ipw_cal_lin <- nonprob(data = cbop_copy_lin, 
                         selection = ~ private + size + occup_code, 
                         svydesign = prob_sample_svy, 
                         target = ~y
                         , control_selection = controlSel(est_method_sel = "gee")
                         )
  
  res_mi_lin <- nonprob(data = cbop_copy_lin, 
                    outcome = y ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_lin <- nonprob(data = cbop_copy_lin, 
                    outcome = y ~ private + size + occup_code, 
                    selection = ~ private + size + occup_code, 
                    svydesign = prob_sample_svy,
                    method_outcome = "glm", family_outcome = "binomial"
                    , control_selection = controlSel(est_method_sel = "gee")
  )
  
  ## with mis-class (transformer)
  res_ipw_cal_tra <- nonprob(data = cbop_copy_tra, 
                             selection = ~ private + size + occup_code, 
                             svydesign = prob_sample_svy, 
                             target = ~y
                             , control_selection = controlSel(est_method_sel = "gee")
  )
  
  res_mi_tra <- nonprob(data = cbop_copy_tra, 
                        outcome = y ~ private + size + occup_code, 
                        svydesign = prob_sample_svy,
                        method_outcome = "glm", family_outcome = "binomial")
  
  res_dr_tra <- nonprob(data = cbop_copy_tra, 
                        outcome = y ~ private + size + occup_code, 
                        selection = ~ private + size + occup_code, 
                        svydesign = prob_sample_svy,
                        method_outcome = "glm", family_outcome = "binomial"
                        , control_selection = controlSel(est_method_sel = "gee")
  )
  
  
  
  
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
  
}

boxplot(results[1:(r-1),]- mean(pop_data_expanded$y))
abline(h=0,col = "red")





## nonprob




