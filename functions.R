###############################################################################
## Shared helper functions for all analysis scripts.
##
## Source this file at the top of each .Rmd via:
##     source("functions.R")
##
## Several functions rely on the following variables existing in the
## calling (global) environment at the time they are invoked:
##     Tmin, Tmax, step_size, times, n
## Each analysis script sets these before calling the functions.
###############################################################################


# ---- Utility: install & load packages --------------------------------------
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}


# ---- ODE model of antibody dynamics ----------------------------------------
# Same structural model for mice and humans (only the parameter values differ):
#   phase 1 (t <= 35): sigmoid growth toward carrying capacity C
#   phase 2 (t  > 35): exponential decay
ode_ab_pop <- function(pars) {

  k1 <- as.numeric(pars[1]) # Growth rate in phase 1
  k2 <- as.numeric(pars[2]) # Decay rate in phase 2
  A0 <- as.numeric(pars[3]) # Initial condition
  C  <- as.numeric(pars[4]) # Carrying capacity for sigmoid growth phase

  times <- c(seq(Tmin, Tmax, step_size))

  k_time <- data.frame(times = times, k = rep(0, length(times)))
  k_time <- k_time %>% mutate(k = case_when(times <= 35 ~ k1,   # day 35 = peak
                                            times  > 35 ~ -k2))

  k_t <- approxfun(k_time$times, k_time$k, rule = 2)

  derivs <- function(times, y, pars, k_t) {
    with(as.list(c(pars, y)), {
      k <- k_t(times)
      if (times < 35) {
        dA <- k * A * (1 - A / C)   # sigmoid growth phase
      } else {
        dA <- k * A                 # exponential decay phase
      }
      return(list(c(dA)))
    })
  }

  y <- c(A = A0)
  out <- ode(y = y, parms = pars, times = times, func = derivs, k_t = k_t)
  as.data.frame(out)
}


# ---- Population-CI sampling for mice ---------------------------------------
# Draws `num` parameter sets from each population distribution (per vaccine
# covariate), used to propagate uncertainty through the ODE.
sample_ab <- function(pop, num, vac) {

  inds_mean <- which(rownames(pop) %in% c("k1_pop", "k2_pop", "A0_pop", "peak_pop"))
  inds_sd   <- which(rownames(pop) %in% c("omega_k1", "omega_k2", "omega_A0", "omega_peak"))

  pars <- matrix(0, num + 1, length(inds_mean))

  for (i in 1:length(inds_mean)) {
    mean_par <- pop$value[inds_mean[i]]
    sd_par   <- pop$value[inds_sd[i]]

    if (i == 1) {
      beta_k1_0_1ug <- pop$value[inds_mean[i] + 1]
      beta_k1_0_2ug <- pop$value[inds_mean[i] + 2]
      beta_k1_0_5ug <- pop$value[inds_mean[i] + 3]
      beta_k1_1ug   <- pop$value[inds_mean[i] + 4]
      log_k1 = log(mean_par) +
        (vac == '0.1ug Pfizer') * beta_k1_0_1ug +
        (vac == '0.2ug Pfizer') * beta_k1_0_2ug +
        (vac == '0.5ug Pfizer') * beta_k1_0_5ug +
        (vac == '1ug Pfizer')   * beta_k1_1ug
      pars[1:num, i]   <- exp(rnorm(num, mean = log_k1, sd = sd_par))
      pars[num + 1, i] <- exp(log_k1)

    } else if (i == 2) {
      beta_k2_0_1ug <- pop$value[inds_mean[i] + 1]
      beta_k2_0_2ug <- pop$value[inds_mean[i] + 2]
      beta_k2_0_5ug <- pop$value[inds_mean[i] + 3]
      beta_k2_1ug   <- pop$value[inds_mean[i] + 4]
      log_k2 = log(mean_par) +
        (vac == '0.1ug Pfizer') * beta_k2_0_1ug +
        (vac == '0.2ug Pfizer') * beta_k2_0_2ug +
        (vac == '0.5ug Pfizer') * beta_k2_0_5ug +
        (vac == '1ug Pfizer')   * beta_k2_1ug
      pars[1:num, i]   <- exp(rnorm(num, mean = log_k2, sd = sd_par))
      pars[num + 1, i] <- exp(log_k2)

    } else if (i == 3) {
      log_A0 = log(mean_par)
      pars[1:num, i]   <- exp(rnorm(num, mean = log_A0, sd = sd_par))
      pars[num + 1, i] <- exp(log_A0)

    } else if (i == 4) {
      log_C = log(mean_par)
      pars[1:num, i]   <- exp(rnorm(num, mean = log_C, sd = sd_par))
      pars[num + 1, i] <- exp(log_C)
    }
  }
  return(pars)
}


# ---- Run ODE across sampled parameter sets (relies on global n) ------------
run_ODE_control <- function(pars) {
  total_AL <- matrix(NA, nrow = length(seq(Tmin, Tmax, step_size)), ncol = n)
  for (i in 1:n) {
    out <- ode_ab_pop(pars[i, ])
    total_AL[, i] <- out$A
  }
  return(total_AL)
}


# ---- Individual-fit trajectories: MICE -------------------------------------
# Merges with Vaccine.Group from the original observation data so that the
# S-Fig-3 plot can colour individuals by dose group.
# Assumes S = 100 simulation replicates per individual (matches Monolix run).
ind_fit_mice <- function(Est, Simulated, original) {

  Fit <- list()
  for (i in 1:nrow(Est)) {
    pars <- c(k1 = Est$k1_mode[i],
              k2 = Est$k2_mode[i],
              A0 = Est$A0_mode[i],
              C  = Est$peak_mode[i])
    fitted <- ode_ab_pop(pars)
    d1 <- data.frame(Days = (times), y = (fitted$A))
    Code <- Est$id[i]

    S <- 100   # 100 simulated replicates per individual
    P <- matrix(NA, nrow = length(times), ncol = S)

    for (j in 1:S) {
      pars <- c(k1 = Simulated$k1[j + S * (i - 1)],
                k2 = Simulated$k2[j + S * (i - 1)],
                A0 = Simulated$A0[j + S * (i - 1)],
                C  = Simulated$peak[j + S * (i - 1)])
      out <- ode_ab_pop(pars)
      P[, j] <- out$A
    }

    Min <- apply(P, 1, function(x) { quantile(x, 0.005) })
    Max <- apply(P, 1, function(x) { quantile(x, 0.995) })

    fit <- cbind(d1, Min, Max, Code)
    fit$Max[fit$Max > 100] <- 100
    fit$Min[fit$Min > 100] <- 100

    Fit[[i]] <- data.frame(fit)
  }

  ind_fit_unlist <- map_df(Fit, ~as.data.frame(.x))
  unique_Vacc <- original[!duplicated(original$Code), c("Code", "Vaccine.Group")]
  ind_fit <- merge(ind_fit_unlist, unique_Vacc, by = "Code", all.x = TRUE)
  ind_fit$Code <- as.factor(ind_fit$Code)
  ind_fit <- ind_fit[order(ind_fit$Code, ind_fit$Days), ]
  return(ind_fit)
}


# ---- Individual-fit trajectories: HUMAN ------------------------------------
# Single-dose cohort, so no vaccine-group merge.
# Assumes S = 94 simulation replicates per individual (matches Monolix run).
ind_fit_human <- function(Estimated, Simulated) {

  Fit <- list()
  for (i in 1:nrow(Estimated)) {
    pars <- c(k1 = Estimated$k1_mode[i],
              k2 = Estimated$k2_mode[i],
              A0 = Estimated$A0_mode[i],
              C  = Estimated$peak_mode[i])
    fitted <- ode_ab_pop(pars)
    d1 <- data.frame(Days = (times), y = (fitted$A))
    Code <- Estimated$id[i]

    S <- 94    # 94 simulated replicates per individual
    P <- matrix(NA, nrow = length(times), ncol = S)

    for (j in 1:S) {
      pars <- c(k1 = Simulated$k1[j + S * (i - 1)],
                k2 = Simulated$k2[j + S * (i - 1)],
                A0 = Simulated$A0[j + S * (i - 1)],
                C  = Simulated$peak[j + S * (i - 1)])
      out <- ode_ab_pop(pars)
      P[, j] <- out$A
    }

    Min <- apply(P, 1, function(x) { quantile(x, 0.005) })
    Max <- apply(P, 1, function(x) { quantile(x, 0.995) })

    fit <- cbind(d1, Min, Max, Code)
    fit$Max[fit$Max > 100] <- 100
    fit$Min[fit$Min > 100] <- 100

    Fit[[i]] <- data.frame(fit)
  }

  ind_fit_unlist <- map_df(Fit, ~as.data.frame(.x))
  ind_fit_unlist <- ind_fit_unlist[order(ind_fit_unlist$Code, ind_fit_unlist$Days), ]
  return(ind_fit_unlist)
}


# ---- R^2 helpers for nls / lm models ---------------------------------------
r2_k1 <- function(model, data) {
  k1 <- data$k1
  y_pred <- predict(model, newdata = data)
  ss_total    <- sum((k1 - mean(k1))^2)
  ss_residual <- sum((k1 - y_pred)^2)
  return(1 - (ss_residual / ss_total))
}

r2_k2 <- function(model, data) {
  k2 <- data$k2
  y_pred <- predict(model, newdata = data)
  ss_total    <- sum((k2 - mean(k2))^2)
  ss_residual <- sum((k2 - y_pred)^2)
  return(1 - (ss_residual / ss_total))
}


# ---- Allometric scaling: mouse -> human parameters -------------------------
# Given a table of mouse k1,k2 across doses and a reference human parameter
# set at its validated dose, project human k1,k2 across scaled doses by
# applying the mouse-side percentage changes to the human anchor.
allometric_scaling <- function(mice_data, SF, human_data) {

  human_start_dose <- human_data$dose
  mice_start_dose  <- human_start_dose / SF

  mice_start_row <- mice_data %>% filter(abs(dose - mice_start_dose) < 1e-8)
  if (nrow(mice_start_row) == 0) {
    stop("Starting dose for mice not found in mice data.")
  }

  k1_mice_start <- mice_start_row$k1
  k2_mice_start <- mice_start_row$k2

  mice_data <- mice_data %>%
    mutate(
      percentage_change_k1 = (k1 - k1_mice_start) / k1_mice_start * 100,
      percentage_change_k2 = (k2 - k2_mice_start) / k2_mice_start * 100
    )

  human_data <- mice_data %>%
    mutate(
      dose_human = dose * SF,
      k1_human   = human_data$k1     * (1 + percentage_change_k1 / 100),
      k2_human   = human_data$k2     * (1 + percentage_change_k2 / 100),
      k1_lwr     = human_data$k1_lwr * (1 + percentage_change_k1 / 100),
      k1_upr     = human_data$k1_upr * (1 + percentage_change_k1 / 100),
      k2_lwr     = human_data$k2_lwr * (1 + percentage_change_k2 / 100),
      k2_upr     = human_data$k2_upr * (1 + percentage_change_k2 / 100)
    )

  output <- human_data %>%
    dplyr::select(
      dose_mice = dose,
      dose_human,
      k1_mice = k1,
      k2_mice = k2,
      percentage_change_k1,
      percentage_change_k2,
      k1_human,
      k2_human,
      k1_lwr,
      k1_upr,
      k2_lwr,
      k2_upr
    )

  return(output)
}


# ---- Generate a single-dose human trajectory under a given SF --------------
generate_trajectory <- function(dose, SF, mice_data, human_parameter, pop_params) {

  scaling_result   <- allometric_scaling(mice_data, SF, human_parameter)
  closest_dose_idx <- which.min(abs(round(scaling_result$dose_human, 1) - dose))

  k1_val <- scaling_result$k1_human[closest_dose_idx]
  k2_val <- scaling_result$k2_human[closest_dose_idx]

  par <- c(
    k1 = k1_val,
    k2 = k2_val,
    A0 = as.numeric(pop_params["A0_pop",  "value"]),
    C  = as.numeric(pop_params["peak_pop", "value"])
  )

  trajectory <- ode_ab_pop(par)
  data.frame(
    Days = trajectory$time,
    y    = trajectory$A,
    SF   = rep(SF,   nrow(trajectory)),
    dose = rep(dose, nrow(trajectory))
  )
}
