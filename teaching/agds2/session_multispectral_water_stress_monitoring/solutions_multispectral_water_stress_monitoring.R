# Preparations -----------------------------------------------------------------
## Load libraries --------------------------------------------------------------
library(caret)
library(ranger)
library(rlang)  # ← required for sym()
library(purrr)
library(dplyr)
library(tidyr)
library(here)
library(ggplot2)
library(readr)
library(tictoc)
library(recipes)
library(vip)
library(themis)
# remotes::install_github("geco-bern/rgeco")
library(rgeco)

set.seed(1982)

## Load data -------------------------------------------------------------------
url_test  <- "https://raw.githubusercontent.com/geco-bern/drought_predictors_competition/02f5fcd9881ad0e407964b636222da3f7b2e73bd/data/competition2025_testing_data.csv"
url_train <- "https://raw.githubusercontent.com/geco-bern/drought_predictors_competition/9e6d3f500bebe59dff203ad2d40d2383a2a183e2/data/competition2025_training_data.csv"

df_train <- readr::read_csv(url_train)
df_test  <- readr::read_csv(url_test)

# Exercise 3.1 Model prediction              -----------------------------------
## Explore data                  -----------------------------------------------


# We have a first look at our target variable fLUE. We consider values below 1 to
# be drought related and color them red. Note how the different sites in our data 
# set show different ranges and regularities of low fLUE values. This is likely
# related to different climates and different vegetation types at these sites.
ggplot(drop_na(df_train, is_flue_drought), 
       aes(x=date, y=flue, 
           color = is_flue_drought, 
           group = site)) +
  geom_line() +
  scale_colour_manual("Drought", 
                      values = c("TRUE"="red","FALSE"="black")) +
  facet_wrap(~site) +
  theme_classic() +
  scale_x_date("date", date_labels = "%Y", 
               breaks = as.Date(c("2005-01-01", "2015-01-01")))

# We check for missing values in our data set:
fix_vis_miss_axis <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
visdat::vis_miss(slice_sample(df_train, n=10000)) + fix_vis_miss_axis

## Fill in missing values        -----------------------------------------------

# Note that there are some missing values in the target and in the predictor 
# values as illustrated above.
df_train2 <- df_train |>
  # not needed
  select(-cluster) |>
  # impute FALSE to missing is_flue_drought (this is only used for over_sampling)
  mutate(is_flue_drought = if_else(is.na(is_flue_drought), FALSE, is_flue_drought)) |>
  # correct
  mutate(is_flue_drought = as.factor(is_flue_drought),
         vegtype         = as.factor(vegtype))

df_test2 <- df_test |>
  # not needed
  select(-cluster) |>
  # correct
  mutate(is_flue_drought = as.factor(is_flue_drought),
         vegtype         = as.factor(vegtype))


# Note the imbalance in the data set with respect to (drought) stressed conditions.
df_train2 |>
  group_by(is_flue_drought) |> summarise(Number_of_observations = n())

# With the data.frame `df_train2` containing our training data we can thus start the model
# training.


## KNN imputation ------------------------------------------------------------------
# perform the KNN imputation here. This should be permissible because
# - KNN imputation as part of the resampling is too slow.
# - we're not imputing missingness caused by the resampling
# - no structured missingness across sites

# This approach contradicts: https://stackoverflow.com/a/60737155 
# But if we truly want global pre‑imputation:
if (FALSE) {
  # For 10k points, This takes about 12 seconds in serial on UBELIX
  # For 20k points, This takes about 51 seconds seconds in serial on UBELIX
  # For 40k points, This takes about 3.3 minutes on UBELIX
  # For 10k points, This takes about 27 seconds in parallel (79 cores) on UBELIX
  # For 20k points, This takes about 1.15 minutes in parallel (79 cores) on UBELIX
  # For 10k points, This takes about 6 seconds on a M4 MacBook Air
  # For 40k points, This takes about 2 minutes on a M4 MacBook Air
  # For 109k points (full training data), this takes about 12.1 minutes on a M4 MacBook Air

  df_train3 <- df_train2 # |> slice_sample(n=40000)
  imp_rec <- recipe(flue ~ ., data = df_train3) %>%
    update_role(site, date, is_flue_drought, new_role = "ID") |>
    step_impute_knn(all_predictors(), neighbors = 5,
                    options = list(
                      nthread = parallelly::availableCores(omit = 1))
    )
  
  start <- Sys.time()
  imp_prep     <- prep(imp_rec, verbose = TRUE); end <- Sys.time(); print(end-start)
  
  # saveRDS(imp_prep,
  #         file = here::here("teaching/agds2",
  #                           "session_multispectral_water_stress_monitoring",
  #                           "imputation_model_kNN_2025.rds"),
  #         compress = "xz")
} else {
  imp_prep <- readRDS(here::here("teaching/agds2/session_multispectral_water_stress_monitoring/imputation_model_kNN_2025.rds"))
}
df_train_imp <- bake(imp_prep, new_data = NULL)      # fully imputed train dataset
df_test_imp  <- bake(imp_prep, new_data = df_test2)  # fully imputed test dataset

visdat::vis_miss(slice_sample(df_train2,    n=10000)) + fix_vis_miss_axis
visdat::vis_miss(slice_sample(df_train_imp, n=10000)) + fix_vis_miss_axis
visdat::vis_miss(slice_sample(df_test2,     n=10000)) + fix_vis_miss_axis
visdat::vis_miss(slice_sample(df_test_imp,  n=10000)) + fix_vis_miss_axis


## Specify target and predictors -----------------------------------------------
rec <- recipe(
  flue ~ .,
  data = df_train_imp
) |>
  # Role handling
  update_role(site, date, is_flue_drought, new_role = "ID") |>
  update_role_requirements("ID", bake = FALSE)

# As a first step we want to impute the missing values (not done in recipe, but done previously)
rec <- rec # |>
# Missing value imputation
# step_impute_knn(all_predictors(), neighbors = 5)   # options = list(nthread = 1, eps = 1e-08),

# As next step, we want to create additional features by considering 
# normalised difference indices between each band combinations,
# e.g. when applied to near-infrared and red this results in the 
# widely used vegetation stress index 
# [NDVI](https://en.wikipedia.org/wiki/Normalized_difference_vegetation_index).
# We define these features by generating for each band combination an expression
# `nd_exprs` and add these features in our recipe with the step `step_mutate`.

# In our recipe we include other model preprocessing steps such as normalization 
# and numeric encoding of vegetation types. Moreover, in our training data set 
# only about a third of all values show stressed conditions (fLUE < 1). Since our 
# focus is on predicting these drops, we artificially increase their presence to 
# make the model evaluation criteria more sensitive to correctly predicting these 
# stressed conditions.

# construct normalised difference indices
band_names <- df_train_imp |>
  select(starts_with("NR_B")) |>
  names()
band_pairs <- combn(band_names, 2, simplify = FALSE)

# Build named list of expressions using set_names()
nd_names <- band_pairs |>
  map_chr(~ paste0("nd_", .x[1], "_", .x[2]))
nd_exprs <- band_pairs |>
  map(~ expr((!!sym(.x[1]) - !!sym(.x[2])) / (!!sym(.x[1]) + !!sym(.x[2]) + 1e-8)))
nd_exprs <- set_names(nd_exprs, nd_names)


rec <- rec |>
  # Feature engineering
  step_mutate(!!!nd_exprs)

## Specify other preprocessing   -----------------------------------------------
rec <- rec |>
  # Preprocessing (only model-predictors!)
  step_normalize(all_numeric_predictors()) |>
  step_novel() |>
  step_dummy(vegtype) # dummy-encode categorical the predictor vegtype

rec <- rec |>
  # upsample cases when flue is < 1. flue-droughts are now overemphasised and
  # make up (over_ratio) times the cases of cases for which is_flue_drought is
  # FALSE. NA-values are imputed with FALSE and thus are not upsampled.
  # step_mutate(is_flue_drought = as.factor(if_else(is.na(is_flue_drought), FALSE, is_flue_drought))) |>
  step_upsample(is_flue_drought, over_ratio = 1) |>
  step_rm(is_flue_drought)

# check roles
summary(rec)

## Specify cross-validation      -----------------------------------------------

# Cross-validation by site (number of folds corresponds to number of sites) or group of sites
set.seed(0)
folds <- caret::groupKFold(
  df_train_imp$site,
  k = length(unique(df_train_imp$site))
)
traincntrlParams <- caret::trainControl(
  index = folds,
  method = "cv",
  savePredictions = "final"   # predictions on each validation resample are then available as modl$pred$Resample
)
tune_grid <- expand.grid(
  .mtry = 7, # c(3, 5, 7),
  .min.node.size = 15, # c(5, 15, 25),
  .splitrule = "variance"
)
# With the training recipe and other settings 
# (`fold`, `traincntrlParams`, `tune_grid`) set up, we can run the model training with the above specified 41 folds.

# Check the folds
# df_train_imp[traincntrlParams$index$Fold01, ]
# df_train_imp[traincntrlParams$index$Fold02, ]
# df_train_imp[traincntrlParams$index$Fold03, ]
# 
# df_train_imp[!(seq_len(nrow(df_train_imp)) %in% traincntrlParams$index$Fold01), ]
# df_train_imp[!(seq_len(nrow(df_train_imp)) %in% traincntrlParams$index$Fold02), ]
# df_train_imp[!(seq_len(nrow(df_train_imp)) %in% traincntrlParams$index$Fold03), ]


## Train the model -------------------------------------------------------------

# rec_prep <- prep(rec, training = df_train_imp)
# rec_juice <- juice(rec_prep) # check the pre-processing steps
# df_train_imp_baked <- bake(rec_prep, df_train_imp)
# df_test_baked <- bake(rec_prep, df_test |> mutate(is_flue_drought = as.factor(is_flue_drought), vegtype = as.factor(vegtype)))

start <- Sys.time()
model <- train(
  rec, # do not use a prepared recipe here
  data            = df_train_imp,
  metric          = "RMSE",
  method          = "ranger",
  tuneGrid        = tune_grid,
  trControl       = traincntrlParams,
  replace         = TRUE,
  sample.fraction = 0.5,
  num.trees       = 500,         # to be boosted to 2000 for the final model
  importance      = "impurity",  # for variable importance analysis, alternative: "permutation"
  num.threads     = 72            # this is passed down to ranger, do not use in combination with carets own parallelism
)
end <- Sys.time()
print(end - start)
# The combination with the optimal resampling statistic is chosen as the final 
# model and the entire training set is used to fit a final model.

# With step_impute_knn() the above takes about 60 minutes for a single fold (on UBELIX).
# With pre-imputed kNN data (df_train_imp), the above takes about 3 minutes for 3 folds (on M4 MacBook Air)
# With pre-imputed kNN data (df_train_imp), the above takes about 1.2 minutes for 3 folds (on UBELIX with 24 RF threads)
# With pre-imputed kNN data (df_train_imp), the above takes about 44 seconds for 3 folds (on UBELIX with 48 RF threads)
# With pre-imputed kNN data (df_train_imp), the above takes about 31 seconds for 3 folds (on UBELIX with 72 RF threads)
# With pre-imputed kNN data (df_train_imp), the above takes about 56 seconds for 6 folds (on UBELIX with 72 RF threads)
# With pre-imputed knn data (df_train_imp), the above takes about 5 minutes for 35 folds (on UBELIX with 72 RF threads)
# With pre-imputed knn data (df_train_imp), the above takes about 22 minutes for 35 folds (on M4 MacBook Air with 72 RF threads speciefied on 10 core machine)

# saveRDS(model,
#         file = here::here("teaching/agds2",
#                           "session_multispectral_water_stress_monitoring",
#                           "model_rf_2025_v3.rds"),
#         compress = "xz")
# model <- readRDS("~/GitHub/fabern/drought_predictors_competition/data/model_rf_2025_v3.rds")

## Evaluate the model cross-validation  -----------------------------------------

# We can now evaluate the goodness of the trained model on the out-of-sample data 
# from our cross-validation folds applied to the training data. 
# Furthermore, the variable importance plots indicates the relative importance of 
# different predictors in the final model. Note how this includes e.g. 
# `nd_NR_B1_NR_B4` which is a feature that was derived (engineered) from the 
# predictors present in the input data set.

preds <- model$pred
preds$site <- df_train_imp$site[preds$rowIndex]

# evaluate overall skill on pooled data from all sites
out <- rgeco::analyse_modobs2(
  preds,
  mod = "pred",
  obs = "obs",
  type = "hex",
  pal = "magma",
  shortsubtitle = TRUE
)

# variable importance plot
vip(model)




# Exercise 3.2 Take part in the competition  -----------------------------------

# Now let's also evaluate our model on the test data. 
# First, we load the test data and `predict()` with our trained model.
# We can then use again the function `rgeco::analyse_modobs2` to compute model 
# skills (RMSE) by comparing our prediction with the true fLUE values on our test 
# data set.
# 
# For the exercise you will not have access to the true fLUE values. You can 
# evaluate your model by committing your predictions to the repository of the 
# course webpage of AGDS II as described under the tab 'For Exercises'.

# The comparison of the fitted model with the test data below indicates a good 
# model fit. The R squared turns out to be better than on the training data 
# (0.84 > 0.63) which indicates that the generalisability of the model appears to 
# be very good.

#       # # NOTE: for the competition, column 'flue' was removed from 
#       # #       the test data. Thus you can NOT compute the scores 
#       # #       and generate a pred vs obs plot yourself.
#       # #       
#       # #       Instead upload your predictions from your best 
#       # #       model as a pull request to
#       # #       https://github.com/geco-bern/agds2_course.
#       # #       
#       # #       To do so: 
#       # #       a) save your results as .*csv: 
#       # readr::write_csv(
#       #   dplyr::select(df_test, site, date, flue_pred),
#       #   file = "~/Desktop/username_results.csv")
#       # #       a) fork the repository 'agds2_course', 
#       # #       b) add and commit your results as 
#       # #          'data/leaderboard/fLUE_fall_2025/[username]_results.csv', 
#       # #          where you replace [username] with your GitHub username
#       # #       c) and open a pull request (PR) from your forked 
#       # #          repository to the main repository at 
#       # #          'geco-bern/agds2_course'.
#       # #       
#       # #       Your predictions will then be added to the the leader 
#       # #       board on the course website at 
#       # #       https://geco-bern.github.io/agds2_course/leaderboard_fLUE_fall_2025.html


# Predict:
test_results1 <- predict(
  model, 
  newdata = df_test_imp |>
    mutate(is_flue_drought = as.factor(is_flue_drought),
           vegtype = as.factor(vegtype))
)

df_test_filled1 <- df_test_imp |>
  mutate(flue_pred = test_results1)

# Simulate alternative output (without gapfilling)
test_results2 <- predict(
  model, 
  newdata = df_test2 |> drop_na() |>  # This removes (instead of imputes) NA # TODO: keep this
    mutate(is_flue_drought = as.factor(is_flue_drought),
           vegtype = as.factor(vegtype))
)
df_test_filled2 <- 
  df_test2 |> drop_na() |>  # This removes (instead of imputes) NA # TODO: keep this
  mutate(flue_pred = test_results2)

# Simulate alternative output (generate random output)
df_test_filled3 <- df_test |> 
  mutate(flue_pred = runif(n= n(), min = 0, max = 1.5))

# Output
df_test_filled1 |>
  dplyr::select(site, date, flue_pred) |>
  readr::write_csv(file = "~/Desktop/example1_fLUE_results.csv")

df_test_filled2 |>
  dplyr::select(site, date, flue_pred) |>
  readr::write_csv(file = "~/Desktop/example2_fLUE_results.csv")

df_test_filled3 |>
  dplyr::select(site, date, flue_pred) |>
  readr::write_csv(file = "~/Desktop/random_fLUE_results.csv")

















# Compute
#' Analyse modelled values versus observed data.
#'
#' Calculates a set of performance statistics and optionally creates plots of modelled
#' versus observed values.
#'
#' @param df A data frame containing columns with names corresponding to arguments
#' \code{mod} and \code{obs}
#' @param mod A character string specifying the variable name (column) of the
#' modelled (simulated) values in data frame \code{df}.
#' @param obs A character string specifying the variable name (column) of the
#' observed values in data frame \code{df}.
#'
analyse_modobs2_metrics <- function(
    df,
    mod,
    obs,
    ...) {
  
  ## rename to 'mod' and 'obs' and remove rows with NA in mod or obs
  df <- df %>%
    as_tibble() %>%
    ungroup() %>%
    dplyr::select(all_of(c(mod = mod, obs = obs))) %>%
    tidyr::drop_na(mod, obs)
  
  ## get linear regression (coefficients)
  linmod <- lm(obs ~ mod, data = df)
  
  ## construct metrics table using the 'yardstick' library
  df_metrics <- df %>%
    yardstick::metrics(obs, mod) %>%
    dplyr::bind_rows(tibble(.metric = "n", .estimator = "standard", .estimate = summarise(df, numb = n()) %>% unlist())) %>%
    dplyr::bind_rows(tibble(.metric = "slope", .estimator = "standard", .estimate = coef(linmod)[2])) %>%
    # dplyr::bind_rows( tibble( .metric = "nse",      .estimator = "standard", .estimate = hydroGOF::NSE( obs, mod, na.rm=TRUE ) ) ) %>%
    dplyr::bind_rows(tibble(.metric = "mean_obs", .estimator = "standard", .estimate = summarise(df, mean = mean(obs, na.rm = TRUE)) %>% unlist())) %>%
    dplyr::bind_rows(tibble(
      .metric = "prmse", .estimator = "standard",
      .estimate = dplyr::filter(., .metric == "rmse") %>% dplyr::select(.estimate) %>% unlist() /
        dplyr::filter(., .metric == "mean_obs") %>%
        dplyr::select(.estimate) %>%
        unlist()
    )) %>%
    dplyr::bind_rows(tibble(
      .metric = "pmae", .estimator = "standard",
      .estimate = dplyr::filter(., .metric == "mae") %>% dplyr::select(.estimate) %>% unlist() /
        dplyr::filter(., .metric == "mean_obs") %>%
        dplyr::select(.estimate) %>%
        unlist()
    )) %>%
    dplyr::bind_rows(tibble(.metric = "bias", .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod - obs), na.rm = TRUE)) %>% unlist())) %>%
    dplyr::bind_rows(tibble(.metric = "pbias", .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod - obs) / obs, na.rm = TRUE)) %>% unlist()))
  
  return(list(df_metrics = df_metrics))
}



# FOR THE GITHUB DO SOMETHING ALONG:
# check if on github
ON_GIT <- Sys.getenv("GITHUB_ACTION") != ""

# if on git read the test labels from encrypted file
if(!ON_GIT) {
  # # Use stored results as reference (not part of the repository):
  ref_df <- read_csv("../../geco-bern/drought_predictors_competition/data/competition2025_testing_data_full.csv") |>
    select(site, date, flue_ref = flue)
  
  # Use encrypted, stored results as reference (part of the repository)
  # use pw to de-zip the encripted zip file
  # pw <- Sys.getenv("fLUE")
  # ref_df <- read_csv(
  #   archive::archive_read(
  #     # archive = here::here("book/images/competition2025_testing_data_full.csv"),
  #     archive = here::here("../../geco-bern/drought_predictors_competition/data/competition2025_testing_data_full.csv.zip"),
  #     file = "competition2025_testing_data_full.csv",
  #     password = pw)
  #   ) |> tibble() |>
  #   select(site, date, flue_ref)
  
  # # Make up random results as reference
  # ref_df  <- df_test |> select(site,date) |> 
  #   mutate(flue_ref = runif(n= n(), min = 0, max = 1.5))
} else {
  
  # Use encrypted, stored results as reference (part of the repository)
  # use pw to de-zip the encripted zip file
  # pw <- Sys.getenv("fLUE")
  # ref_df <- read.csv(
  #   archive::archive_read(
  #     archive = here::here("book/images/reference_results_fLUE.csv.zip"),
  #     file = "reference_results_fLUE.csv",
  #     password = pw)
  #   ) |> tibble()
  
  # Make up random results as reference (
  ref_df  <- df_test |> select(site,date) |> 
    mutate(flue_ref = runif(n= n(), min = 1.2, max = 1.5))
}



# Compare skills of multiple predictions
predictions <- list(
  "gap_filled" = df_test_filled1,
  "with_gaps" = df_test_filled2,
  "random" = df_test_filled3)

# loop over all files (/predictions)
lb <- lapply(names(predictions), function(username) {
  modelled = predictions[[username]]
  
  # calculate SKILL data.frame()
  skills <- try(
    {
      df_pred_obs <- dplyr::inner_join(
        ref_df,
        modelled,
        by = join_by(site, date)) |>
        dplyr::select(site, date, flue_ref, flue_pred)
      # stopifnot(nrow(df_pred_obs) == 20827)
      
      out <- analyse_modobs2_metrics(
        df_pred_obs,
        mod = "flue_pred",
        obs = "flue_ref")
      out$df_metrics |> mutate(User = username)
    },
  )
})

lb <- dplyr::bind_rows(lb)



lb |>
  select(User, .metric, .estimate) |>
  pivot_wider(names_from = .metric, values_from = .estimate) |>
  arrange(rmse) |>
  select(User, n, rmse, rsq, slope)
# 
# # A tibble: 3 × 11
#   User        rmse       rsq   mae     n    slope mean_obs prmse  pmae    bias  pbias
#   <chr>      <dbl>     <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl> <dbl>   <dbl>  <dbl>
# 1 gap_filled 0.149 0.612     0.106 17777  1.12       0.880 0.169 0.121 -0.0170 0.0614
# 2 with_gaps  0.150 0.606     0.107 16280  1.09       0.882 0.170 0.121 -0.0195 0.0540
# 3 random     0.512 0.0000845 0.426 17777 -0.00500    0.880 0.581 0.483 -0.133  0.0234

