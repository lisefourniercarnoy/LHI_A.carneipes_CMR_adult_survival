###############################################################################

# authors: Lise Fournier-Carnoy, Alex L. Bond, Jennifer  L. Lavers 
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2011-2024 CMR data

# script objective: run RMark analysis

###############################################################################

rm(list = ls())

library(tidyverse) # for data manipulation
library(stringr)
library(RMark) # for mark-recapture analysis
library(R2ucare) # for mark-recapture analysis
library(corrplot) # to check correlations
library(jagsUI) # for the bayesian MCMC
library(kableExtra) # to make tables
library(flextable) # making publishable tables


## Read in data ---------------------------------------------------------------

sbsh  <- readRDS("data/tidy/LHI_FFSH_capture_histories.rds") %>% 
  rename(age_capture = Age) %>% 
  glimpse()

enso    <- readRDS("data/tidy/ENSO_tidy.rds") %>% glimpse()
ao      <- readRDS("data/tidy/AO_tidy.rds") %>% glimpse()
temp    <- readRDS("data/tidy/temp_tidy.rds") %>% dplyr::filter(time>=2013) %>% glimpse()
#pdo     <- readRDS("data/tidy/pdo_tidy.rds") |> glimpse()

clim_list = list(enso, 
                 ao, 
                 temp)
clim <- clim_list %>% 
  reduce(~ inner_join(.x, .y, by = "time"))

clim <- clim %>% dplyr::filter(time >= 2013) # remove the first few years because very few birds are recaptured and estimating survival in these years eats up a lot of parameters.
glimpse(clim)

# make two time groups for recapture
recap_group <- data.frame(time = 2013:2024) %>% 
  mutate(recap_group = ifelse(time > 2019, "A", "B")) # recapture effort was much higher 2020-onwards, so test whether two groups are good to include as covariates.


## Check for collinearity -----------------------------------------------------

# following method by Grosbois et al. 2008, check the collinearity of variables before testing them.
m <- cor(clim)
test <-  cor.mtest(clim, conf.level = 0.95)
corrplot(m, p.mat = test$p, sig.level = 0.1, order = 'hclust')
# any circle not crossed is a pair of collinear variables.

# several pairs of covariates are collinear. 
# if they are tested together in a model, we won't know whether the predicted effect is because of one, or the other.

# here we don't have the power to make a model with multiple climate covariates so we don't need to deal with it
# see XX_archive.R for example code. 

# we can now move on !


## Process the data -----------------------------------------------------------

# we are analysing an open-population system, and the model we fit has assumptions that need to be checked.

# check CJS assumptions by running a goodness-of-fit test
ch <- sbsh %>% dplyr::select(`2013`:`2024`) %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
n <- rep(1, nrow(sbsh))

# I'm following https://jamesepaterson.github.io/jamespatersonblog/2020-05-20_gof_for_CJS

# test 1: 'is there evidence of animals having equal capture probabilities and equal survival'.
overall_CJS(ch, n) # test is not significant, first assumption checked. We can go ahead with the CJS model (Gimenez et al. 2017)
# if this was significant, our data would not be suited for a CJS model, and results would be very wrong.

# if test 1 is significant, run test 2-4 to understand whether recapture depends on when animal was marked or whether marking affects survival.
# test 2: 'does recapture deend on when an animal was first marked?' 
test2ct(ch, n) # here p-val is not significant. in $details, some tests cannot be conducted because of low sample sizes.
test2cl(ch, n) # again, p-val is not significant. in $details, some tests cannot be conducted because of low sample sizes.

# test 3: 'does marking affect survival?'
test3sr(ch, n) # p-val not significant, same thing for $details

# test 4: 'for animals seen again, does when they are recaptured depend on whether they were marked on or before time t?'
test3sm(ch, n)

# We're all good to go. Gimenez et al. 2017 confirms this.


## Checking Variance Inflation Factor (c-hat) ---------------------------------

# c-hat can help adjust model output based on extra-binomial noise (more variance in capture history data than expected) or overdispersion (larger-than-normal discrepancies between observed and predicted values for the binomial model).
# here, an estimation of c-hat is Chi2 / df, and is valid for time-dependent models. (which i initially tested here)
# which is essentially the Chi2 and df from the overall_CJS call just above.

test <- overall_CJS(ch, n)
test$chi2/test$degree_of_freedom # simple calculation from https://jamesepaterson.github.io/jamespatersonblog/2020-06-30_chat_for_CJS.html

# another estimation from gentle intro to mark, chapter 5, section 5.5.1
t2_chi2 <- test2ct(ch, n)$test2ct[1]
t3_chi2 <- test3sr(ch, n)$test3sr[1]
t2_df <- test2ct(ch, n)$test2ct[2]
t3_df <- test3sr(ch, n)$test3sr[2]

(t2_chi2 + t3_chi2)/(t2_df + t3_df) # chi2 of (test2 + test3) divided by degrees of freedom of (test2 + test3).

# a value ~1 indicates good variance given the model, <1 is underdispersed, and >1 overdispersed
# here we have c-hat << 1, which is likely because of small sample size.
# underdispersion is rare, and there doesn't seem to be either 1. a biological reason for it or 2. a fix for it.
# for this reason, we will set c-hat to 1 and just move on.
# see Lebreton 2012, Lampo 2017 and Robertson 2016 for explanations/examples.


## Test model sets ------------------------------------------------------------

# again, following Grosbois et al. 2008, we are first assessing the general shape of survival. is it constant, it is logarithmic etc.
# in addition, we will add recapture group (pre-2019 and post-2019) in the test formulas. if the formula is not good, we can ignore recapture group for the rest of the analysis.
# with the best shape, we can then add the environmental variables to test.
sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2013,
                             groups = "age_capture"
                             )
sbsh_ddl <- make.design.data(sbsh_process)

sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim)
sbsh_ddl$p   <- merge_design.covariates(sbsh_ddl$p, recap_group)

sbsh_ddl$Phi$time <- as.numeric(as.character(sbsh_ddl$Phi$time))
sbsh_ddl$p$time   <- as.numeric(as.character(sbsh_ddl$p$time))

# fix p to 0 in 2019 - no capture effort that year, so this should not be a
# freely estimated parameter (reviewer point re: parameter counts / Table 3
# artefacts)
sbsh_ddl$p$fix <- ifelse(sbsh_ddl$p$time == 2019, 0, NA)

# core candidate set: just 4 models, as recommended
phi.formulas <- list(
  const = list(formula = ~1),
  time  = list(formula = ~factor(time))
)
p.formulas <- list(
  const = list(formula = ~1),
  time  = list(formula = ~factor(time)),
  recap_group = list(formula = ~factor(recap_group))
)

model.combos <- list(
  list(name = "Phi_dot_p_dot",   file = "phi_dot_p_dot",   phi = phi.formulas$const, p = p.formulas$const),
  list(name = "Phi_dot_p_time",   file = "phi_dot_p_time",  phi = phi.formulas$const, p = p.formulas$time),
  list(name = "Phi_time_p_dot",   file = "phi_time_p_dot",  phi = phi.formulas$time,  p = p.formulas$const),
  list(name = "phi_time_p_time",   file = "phi_time_p_time", phi = phi.formulas$time,  p = p.formulas$time),
  list(name = "phi_time_p_recap",   file = "phi_time_p_recap", phi = phi.formulas$time,  p = p.formulas$recap_group),
  list(name = "phi_dot_p_recap",   file = "phi_dot_p_recap", phi = phi.formulas$const,  p = p.formulas$recap_group)

)

models <- list()
dir_name <- paste0("data/rmark_outputs/reference_models/", today())
dir.create(path = dir_name, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(model.combos)) {
  models[[model.combos[[i]]$name]] <- mark(
    data = sbsh_process,
    ddl = sbsh_ddl,
    model.parameters = list(
      Phi = model.combos[[i]]$phi,
      p   = model.combos[[i]]$p
    ),
    filename = paste0(dir_name, "/", model.combos[[i]]$file),
    model.name = model.combos[[i]]$name,
    delete = FALSE,
    chat = 1
  )
}

reference_model_set <- model.table(models, use.AIC = FALSE)
reference_model_set
saveRDS(reference_model_set, paste0(dir_name, "/reference_models_table.rds"))

# Table 1, formatted properly (Lebreton notation) rather than raw RMark output
table_data <- reference_model_set %>%
  dplyr::select(Phi, p, npar, AICc, DeltaAICc, weight, Deviance) %>%
  mutate(
    AICc      = round(AICc, 2),
    DeltaAICc = round(DeltaAICc, 2),
    weight    = ifelse(weight < 0.001, "<0.001", round(weight, 3)),
    Deviance  = round(Deviance, 2)
  ) |> 
  mutate(model = paste0("\u03c6 (", 
                        case_when(Phi == "~1" ~ "1", 
                                  Phi == "~factor(time)" ~ "t",
                                  TRUE ~ NA
                                  ),
                        ")", 
                        
                        " p (", 
                        case_when(p == "~1" ~ "1", 
                                  p == "~factor(time)" ~ "t",
                                  p == "~factor(recap_group)" ~ "r_group",
                                  TRUE ~ NA
                                  ),
                        ")")) |> 
  dplyr::select(-Phi, -p) |> 
  dplyr::select(model, AICc, DeltaAICc, weight, npar, Deviance)
  
ft <- flextable(table_data) %>%
  set_header_labels(
    model     = "Model",
    AICc      = "AICc",
    DeltaAICc = "\u0394AICc",
    weight    = "AICc weight",
    npar      = "Parameters",
    Deviance  = "Deviance"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = "model", align = "left", part = "all") %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  autofit()

ft

save_as_docx(ft, path = "plots/table_1_model_selection.docx")


## Hypothesis testing ---------------------------------------------------------

# now we need to test whether the top model explains a significant amount of variation in survival
# Grosbois suggests 3 tests based on likelihood-ratio.

# set up what we need first: 
# deviance values
Dev_cst <- reference_model_set[reference_model_set$Phi == "~1" & reference_model_set$p == "~1",]$Deviance   # Fcst
Dev_co  <- reference_model_set[1,]$Deviance        # Fco (row 1 = lowest AICc) - or QAICc/QDeviance if c-hat has been adjusted
Dev_t   <- reference_model_set[reference_model_set$Phi == "~factor(time)" & reference_model_set$p == "~factor(time)",]$Deviance # Ft

# number of parameters
J <- reference_model_set[1,]$npar # number of parameters in top model
n <- reference_model_set[reference_model_set$Phi == "~factor(time)" & reference_model_set$p == "~factor(time)",]$npar # number of survival estimates in Ft (i.e number of parameters because it's the fully-time-dependent model)


# 1. we'll test the following null hypothesis: " the top model fits the data as well as the time-dependent model. "
lrt <- Dev_co - Dev_t # Fco - Ft in Grosbois.
lrt
df <- n - J
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value 
# p>0.05, we can't reject the null hypothesis. this means we can move to the second test


# 2. we'll now test a second null hypothesis: " the top model fits the data as well as the constant model. "
lrt <- Dev_cst - Dev_co # Fcst - Fco in Grosbois
lrt
df <- J - 1
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value
# p<0.05, we reject this null hypothesis, the p-value is highly significant, so the top model is better than the constant. this means we can do the third test.


# 3. this is an alternate null hypothesis to 2., in case the 1. was significant.
c_hat <- (Dev_co - Dev_t) / (n - J)
F_stat <- ((Dev_cst - Dev_co) / (J - 1)) / c_hat
df1 <- J - 1
df2 <- n - J
p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)

cat("F-statistic =", F_stat, "\n")
cat("Degrees of freedom =", df1, "and", df2, "\n")
cat("p-value =", p_value, "\n")
# the output is corrected for multiple testing... if p>0.05, we reject this null hypothesis, and the top model is better than the constant.


# Grosbois 2008 also does an Information Theory selection, which I don't understand and didn't do.


## obtain survival from the top model -----------------------------------------

surv <- mark(
  data = sbsh_process,
  ddl = sbsh_ddl,
  model.parameters = list(
    Phi = phi.formulas$const,
    p   = p.formulas$time
  ),
  filename = paste0(dir_name, "/top_model"),
  model.name = "top_m",
  delete = FALSE,
  chat = 1
)

paste0("Estimated survival: ", round(surv$results$real$estimate[1], 3), 
       ", SE = ", round(surv$results$real$se[1], 3),
       ", CI ", round(surv$results$real$lcl[1], 3), " - ", round(surv$results$real$ucl[1], 3)
)


# -- okay. at this point we have an overall survival estimate. good!
# -- what the top model says is: 'this phi~1 model is a good enough estimation of survival"
# -- but this doesn't mean we can't fit other models.
# -- however we would also like to understand how survival varies year-to-year, and for that, we need a fully-time-dependent model.
# -- problem being: 

test <- mark(
  data = sbsh_process,
  ddl = sbsh_ddl,
  model.parameters = list(
    Phi = phi.formulas$time,
    p   = p.formulas$time
  ),
  filename = paste0(dir_name, "/phi_time_p_time_test"),
  model.name = "full_t",
  delete = FALSE,
  chat = 1
)

test$results$real

# -- a lot of parameters (both Phi and p) are at a boundary (0 or 1).
# -- this essentially means that RMark's estimation of survival and recapture probability is not confident at all, due to low sample size.
# -- we will therefore use a Bayesian model to understand the relationship between climate indices and survival.


## BAYESIAN ANALYSIS ----------------------------------------------------------

# okay so. i use the crashcourse on Bayes' theorem to understand this: https://www.youtube.com/watch?v=9TDjifpGj-k&t=67s
# Bayesian statistics is: updating your beliefs about something when new info comes along. PRIOR beliefs update to your POSTERIOR belief
# we update what we think survival is, based on the new climate information we have.

# to test whether a climate covariate is actually useful, the Markov Chain Monte Carlo is going to test many models with different parameters of covariates.
# for example it'll test A. survival ~ 0.1*temp, and B. survival ~ 0.2*temp, and will tell you how likely A and B are, given your data.
# the MCMC 'samples from the POSTERIOR distribution' of parameter values. If the parameter value has a higher likelihood with our data, those values get higher POSTERIOR probability.
# so we test whether with a small effect of temp (in model A) is likely given the recapture data we have. same with a slighly high effect of temp (like in B) etc.
# this gives us a curve: x being the temp parameter, y being the likelihood.
# if the parameter value of 0.1 is given high posterior probability (like in A.), it means that the data highly supports this value.
# say it tests values -1 to 1, it'll give us a 95% confidence interval (we're 95% sure the real parameter value is between _ and _).
# if the confidence interval contains 0, then we consider the covariate doesn't have statistical support, because there's a good chance there is no effect of the covariate on survival.

glimpse(sbsh)
glimpse(sbsh_ddl)

# obtain the first year of capture for each individual
get_first <- function(ch) min(which(ch == 1))
first_capture <- apply(ch, 1, get_first)
ch_matrix <- as.matrix(ch)

# okay here we go.

years <- 2013:2024
effort <- as.numeric(years != 2019)

# -- the model works across two matrices:
# -- y is the observed recaptures: did we see this bird this year? yes/no
# -- z is the unobserved (latent) states: was this bird alive this year? yes/no

# -- the y matrix is fixed, and the model tries to fill in the z matrix.
# -- so, for a given bird i, and a y matrix of: 110, the z matrix will try to assess how likely this recapture history is with different Phi and p values.
# -- the model is going to test lots of values of Phi to fill the last year.
# -- the probabilities of recapture history 110 would be:

# -- for a 3 year dataset (where you tag only the first year), there are 4 possible recapture histories you can get:
# -- 100, 101, 110, 111, 

# -- given that you recaptured the bird in year 2, what is the likely survival from year 2 to year 3?
# -- let's say we fix recapture probability at 60%, and we test survival 0.25 and 0.5.
# -- 2 scenarios are possible: 

# -- survive - survive (recaptured) - survive (not recaptured) :
# -- survive y1 AND be captured in y2 AND survive y2 AND survive but not be captured in y3 
# --   0.25     *        0.6          *   0.25       *       (1-0.6)
# --   0.5      *        0.6          *   0.5        *       (1-0.6)

# -- OR survive - survive (recaptured) - die

# -- survive y1 AND be captured in y2 AND survive y2 AND die
# --   0.25      *        0.6          *   0.25       *       (1-0.25)
# --   0.5       *        0.6          *   0.5        *       (1-0.5)

# -- since we don't know which one actually happened, we add them together:
# -- Phi_0.25: (0.25 * 0.6 * 0.25 * (1-0.6)) + (0.25 * 0.6 * 0.25 * (1-0.25)) = 0.043125
# -- Phi_0.5:  (0.5 * 0.6 * 0.5 * (1-0.6)) + (0.5 * 0.6 * 0.5 * (1-0.5)) = 0.135

# -- from the two Phi parameter values tested, the one with a higher likelihood is Phi = 0.5. 
# -- the function will test many values of Phi, and many values of p, to maximise likelihood.


# -- the z matrix contains '11?' (doesn't know of survival from year 2 and 3), and the JAGS functions tests the likelihood of 110 and 111.
# -- with longer capture histories it gets a bit complicated to explain but from my understanding it tries to find the likelihood that suits the entire capture history...


writeLines("model {

  # -- we'll start with defining priors - what we think we know about survival and recapture probability
  # -- we don't know much, and also the sample size is small, so we will use flat priors (dunif())

  for (t in 1:(T-1)) {
    phi[t] ~ dunif(0, 1) # fully flat prior, no pull toward any value
  }
  for (t in 1:T) {
    p0[t] ~ dunif(0, 1)
    p[t] <- p0[t] * effort[t] # forces p to exactly 0 in zero-effort years
  }
  
  # -- then, for every bird, calculate the probability of the recapture history

  for (i in 1:N) {
  
    z[i, first[i]] <- 1
    
    for (t in (first[i] + 1):T) {
    
      z[i, t] ~ dbern(z[i, t - 1] * phi[t - 1])

    }
    for (t in first[i]:T) {
      y[i, t] ~ dbern(z[i, t] * p[t])
    }
  }
}", "data/JAGS_models/cjs_model.txt")

set.seed(123)

# setup JAGS to run many simulations.
jags_data <- list(
  N      = nrow(ch_matrix), # number of birds banded
  T      = ncol(ch_matrix), # number of years under study
  y      = ch_matrix,
  first  = first_capture,
  effort = effort
)

inits <- function() {
  z_init <- matrix(NA, nrow = jags_data$N, ncol = jags_data$T)
  for (j in 1:jags_data$N) {
    first_j <- jags_data$first[j]
    if (first_j < jags_data$T) {
      for (t in (first_j + 1):jags_data$T) {
        z_init[j, t] <- 1
      }
    }
  }
  list(phi = runif(jags_data$T - 1, 0.3, 0.9),
       p0  = runif(jags_data$T, 0.3, 0.9),
       z   = z_init)
}

params <- c("phi", "p", "deviance")

# quick sanity run first - confirm p[7] (2019) comes back as exactly 0 before committing to the full 100k-iteration run
fit_test <- jags(data = jags_data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "data/JAGS_models/cjs_model.txt",
                 n.chains = 3,
                 n.iter = 5000,
                 n.burnin = 1000,
                 n.thin = 1)
fit_test$summary["p[7]", ]  # mean, sd, everything should be 0

# full run
fit <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "data/JAGS_models/cjs_model.txt",
            n.chains = 3,
            n.iter = 100000,
            n.burnin = 20000,
            n.thin = 1)

fit  # check phi[1]...phi[11] and p[1]...p[12] are all present, Rhat < 1.1

plot(fit)  # visually check chain mixing (should overlap well) per parameter


### Test the correlation between covariates and survival ----------------------

str(fit$sims.list)

phi_samples <- fit$sims.list$phi # survival estimates for each year, for every simulation in JAGS
covariate_years <- years[-length(years)]

get_covariate_vector <- function(var_name) {
  clim %>%
    dplyr::filter(time %in% covariate_years) %>%
    dplyr::arrange(match(time, covariate_years)) %>%
    dplyr::pull({{ var_name }})
}

ao_jja_vec   <- get_covariate_vector(AO_jun_jul_aug)
ao_ond_vec   <- get_covariate_vector(AO_oct_nov_dec)
ao_fma_vec   <- get_covariate_vector(AO_feb_mar_apr)

enso_jja_vec   <- get_covariate_vector(ENSO_jun_jul_aug)
enso_ond_vec   <- get_covariate_vector(ENSO_oct_nov_dec)
enso_fma_vec   <- get_covariate_vector(ENSO_feb_mar_apr)

sst_jja_vec   <- get_covariate_vector(sst_winter)
sst_ond_vec   <- get_covariate_vector(sst_oct_nov_dec)
sst_fma_vec   <- get_covariate_vector(sst_feb_mar_apr)


# -- here we select some JAGS simulation's survival, and calculate the correlation with the climate index.
compute_cor_draws <- function(phi_samples, covariate_vec) {
  apply(phi_samples, 1, function(phi_draw) cor(phi_draw, covariate_vec))
}

cor_ao_jja   <- compute_cor_draws(phi_samples, ao_jja_vec)
cor_ao_ond   <- compute_cor_draws(phi_samples, ao_ond_vec)
cor_ao_fma   <- compute_cor_draws(phi_samples, ao_fma_vec)

cor_sst_jja  <- compute_cor_draws(phi_samples, sst_jja_vec)
cor_sst_ond  <- compute_cor_draws(phi_samples, sst_ond_vec)
cor_sst_fma  <- compute_cor_draws(phi_samples, sst_fma_vec)

cor_enso_jja <- compute_cor_draws(phi_samples, enso_jja_vec)
cor_enso_ond <- compute_cor_draws(phi_samples, enso_ond_vec)
cor_enso_fma <- compute_cor_draws(phi_samples, enso_fma_vec)



summarise_cor <- function(cor_draws, label) {
  data.frame(
    covariate  = label,
    mean_r     = mean(cor_draws),
    lower95    = quantile(cor_draws, 0.025),
    upper95    = quantile(cor_draws, 0.975),
    prop_positive = mean(cor_draws > 0),
    prop_negative = mean(cor_draws < 0)
  )
}

cor_summary <- dplyr::bind_rows(
  summarise_cor(cor_ao_jja,   "AO_jja"),
  summarise_cor(cor_ao_ond,   "AO_ond"),
  summarise_cor(cor_ao_fma,   "AO_fma"),
  
  summarise_cor(cor_sst_jja,  "SST_jja"),
  summarise_cor(cor_sst_ond,  "SST_ond"),
  summarise_cor(cor_sst_fma,  "SST_fma"),
  
  summarise_cor(cor_enso_jja, "ENSO_jja"),
  summarise_cor(cor_enso_ond, "ENSO_ond"),
  summarise_cor(cor_enso_fma, "ENSO_fma")
  
)
cor_summary


cor_draws_df <- dplyr::bind_rows(
  data.frame(covariate = "AO_jja", r = cor_ao_jja),
  data.frame(covariate = "AO_ond", r = cor_ao_ond),
  data.frame(covariate = "AO_fma", r = cor_ao_fma),
  
  data.frame(covariate = "SST_jja", r = cor_sst_jja),
  data.frame(covariate = "SST_ond", r = cor_sst_ond),
  data.frame(covariate = "SST_fma", r = cor_sst_fma),
  
  data.frame(covariate = "ENSO_jja", r = cor_enso_jja),
  data.frame(covariate = "ENSO_ond", r = cor_enso_ond),
  data.frame(covariate = "ENSO_fma", r = cor_enso_fma)
  
)


### Figure S3. Correlation between survival and covariates --------------------

cor <- ggplot(cor_draws_df, aes(x = r, fill = covariate)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  facet_wrap(~ covariate, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = rep('steelblue', length(unique(cor_draws_df$covariate)))) +
  labs(
    x = "Posterior correlation (φ vs covariate)",
    y = "Posterior density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(cor, filename = "plots/fig_S3_survival_correlation.jpg", height = 5, width = 8, dpi = 500)


# -- check the % of draws that cross zero
n_draws <- length(cor_draws_df$r[cor_draws_df$covariate == "SST_jja"])

n_dr_sst <- length(cor_draws_df$r[cor_draws_df$covariate == "SST_jja" & cor_draws_df$r < 0])

perc_neg_sst <- n_dr_sst/n_draws
perc_neg_sst # 100% of draws are negative.


## check recapture prob

p_samples <- fit$sims.list$p
colMeans(p_samples)

par(mfrow = c(1,1))
plot(x = years, y = colMeans(p_samples))

ggplot(p_samples) +
  geom_boxplot()


df <- as.data.frame(p_samples)
colnames(df) <- paste0("V", seq_len(ncol(df)))   # or your own meaningful names

df_long <- pivot_longer(df, cols = everything(),
                        names_to = "parameter", values_to = "value")

# keep columns in numeric order rather than alphabetical
df_long$parameter <- factor(df_long$parameter, levels = colnames(df))

ggplot(df_long, aes(x = parameter, y = value)) +
  geom_boxplot()
summary(p_samples)


## Run permutation tests ------------------------------------------------------

# -- correlation does not mean causation.
# -- we need to statistically test whether the relationship we found just before is just a lucky hit or not.
# -- for this we'll do a permutation test. 

# -- the data we observed has 1 climate value linked to 1 estimated annual survival (e.g. 2020 SST is linked to 2020 survival)
# -- what the permutation test does is reshuffle the climate values (e.g. 2013 SST is paired with 2019 survival, and 2014 SST with 2024 survival etc.),
# -- and test whether the correlation we observed can be found by recombinating survival and climate.
# -- if less than 5% of permutations find a correlation like we observed (p_perm <0.05), then we consider the correlation to be statistically supported.
# -- if not, the evidence is less strong.

set.seed(123)

# remove the last year of the data from being tested (identifiability problem)
# as well as the years around 2019, where recapture effort was zero (identifiability problem again)
keep <- setdiff(1:11, c(6, 7, 11))

phi_use <- phi_samples[, keep] # this is from the JAGS iterations

# pick fewer of the iterations, so it runs a bit quicker (the original has >200,000 rows)
idx <- round(seq(1, nrow(phi_use), length.out = 2000))
phi_use <- phi_use[idx, ]

# calculate the test statistic: mean over posterior draws of cor(phi_draw, covariate)
stat_fun <- function(covariate) {
  r <- cor(t(phi_use), covariate)
  mean(r, na.rm = TRUE)
}

B <- 5000
covariates <- list(ao_jja_vec, ao_ond_vec, ao_fma_vec,
                   enso_jja_vec, enso_ond_vec, enso_fma_vec, 
                   sst_jja_vec, sst_ond_vec, sst_fma_vec
                   )
names(covariates) <- c("AO_jja", "AO_ond", "AO_fma", 
                       "ENSO_jja", "ENSO_ond", "ENSO_fma", 
                       "SST_jja", "SST_ond", "SST_fma")


### Figure S4. Permutation tests ----------------------------------------------

library(patchwork)

permutation_tests <- list()
plots <- list()

for (i in 1:length(covariates)) {
  use <- covariates[[i]][keep]
  obs <- stat_fun(use) # the correlation between the current climate covariate and all iterations
  null <- replicate(B, stat_fun(sample(use))) # for every recombinated set of climate values, calculate the correlation between climate and survival
  p_perm <- (1 + sum(abs(null) >= abs(obs))) / (B + 1) # the % of recombinated sets that have correlation higher than the observed.
  
  permutation_tests[[i]] <- data.frame(obs = obs, 
                                       p_perm = p_perm)
  
  plots[[i]] <- ggplot(data.frame(null = null), aes(x = null)) +
    geom_histogram(bins = 50, fill = "steelblue") +
    geom_vline(xintercept = obs, colour = "red", linewidth = 1) +
    xlab("Mean posterior r") +
    ylab("Frequency") +
    ggtitle(paste0(names(covariates)[[i]], ", p = ", round(p_perm, 3))) +
    theme_minimal()
}

combined <- wrap_plots(plots, ncol = 3)

ggsave("plots/fig_S4_permutation_tests.jpg", combined,
       width = 7, height = 7, dpi = 500) 

# -- the plots show a histogram (with the permutations and what they find to be the correlation between climate and survival)
# -- and a red line, which is what the observed correlation is.
# -- essentially the closer to the bulk of the bars, the less likely the correlation is meaningful
# -- here, SST is the closest to being meaningful, but p_perm is still over 0.05.


## Figures: unconstrained phi(t) vs climate -----------------------------------

phi_summary <- as.data.frame(fit$summary) %>%
  tibble::rownames_to_column("param") %>%
  dplyr::filter(grepl("^phi\\[", param)) %>%
  dplyr::mutate(
    idx  = as.numeric(gsub("phi\\[|\\]", "", param)),
    year = years[idx]
  ) %>%
  dplyr::select(year, mean, `2.5%`, `97.5%`)

phi_summary


### Figure S1. unconstrained phi(t) v. climate --------------------------------

plots <- list()

for (i in 1:length(names(clim[,-1]))) {

  current_cov <- names(clim[,-1])[i]
  x_axis_name <- paste0(ifelse(
    
    #period
    grepl("feb_mar_apr", current_cov), "Late breeding", 
    ifelse(grepl("oct_nov_dec", current_cov), "Early breeding", "Wintering")),
    
    #covariate
    " ",
    ifelse(grepl("sst", current_cov,), "SST*", ifelse(
      grepl("ENSO", current_cov), "ENSO", "AO"
    ))
  )
  
  df <- clim[clim$time %in% phi_summary$year, c("time", current_cov)]
  min <- min(df[,2])
  max <- max(df[,2])
  
  phi <- phi_summary %>% left_join(df, by = c("year" = "time"))
  
  surv <- ggplot(phi, aes(x = .data[[current_cov]], y = mean)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), linewidth = 1, width = 0, color = "steelblue") +
    geom_point(size = 3, shape = 18, colour = "steelblue") +
    geom_text(aes(label = year), vjust = -1, size = 4) +
    labs(x = x_axis_name, y = "Annual survival probability \n (posterior mean)") +
    theme_classic()
  
  plots[[i]] <- surv
  
}
combined <- wrap_plots(plots, ncol = 3)
ggsave(combined, filename = paste0("plots/fig_S2_survival_climate.jpg"), width = 12, height = 8, dpi = 500)


### Figure 2. unconstrained phi(t) v. SST -------------------------------------

sst_df <- clim[clim$time %in% phi_summary$year, c("time", "sst_winter")]
sst_min <- min(sst_df$sst_winter)
sst_max <- max(sst_df$sst_winter)

phi_sst <- phi_summary %>% left_join(sst_df, by = c("year" = "time"))

sst_surv <- ggplot(phi_sst, aes(x = sst_winter, y = mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 1, width = 0, color = "steelblue") +
  geom_point(size = 3, shape = 18, colour = "steelblue") +
  geom_text(aes(label = year), vjust = -1, size = 4) +
  labs(x = "Boreal summer SST anomaly", y = "Annual survival probability \n (posterior mean)") +
  theme_classic()
sst_surv
ggsave(sst_surv, filename = "plots/fig_2_SST_survival.jpg", width = 5, height = 3)

### END ###