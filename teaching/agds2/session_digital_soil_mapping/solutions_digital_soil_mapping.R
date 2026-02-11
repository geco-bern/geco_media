# Preparations -----------------------------------------------------------------
## Load libraries --------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ranger)
library(here)
library(knitr)
library(parallel)
library(yardstick)

set.seed(1982)

## Load data -------------------------------------------------------------------
df_full <- readRDS(here("data/df_full.rds"))

# Exercise 5.1 Simple model                  -----------------------------------

## Specify target and predictors -----------------------------------------------
# Specify target: The pH in the top 10cm
target <- "waterlog.100"

# Make sure the target variable is encoded as a factor
df_full <- df_full |> 
  mutate(across(starts_with("waterlog"), as.factor))

# Specify predictors_all: Remove soil sampling and observational data
predictors_all <- names(df_full)[14:ncol(df_full)]

## Split data ------------------------------------------------------------------
# Split dataset into training and testing sets
df_train <- df_full |> 
  filter(dataset == "calibration")

df_test <- df_full |> 
  filter(dataset == "validation")

# Filter out any NA to avoid error when running a Random Forest
df_train <- df_train |> 
  select(all_of(c(target, predictors_all))) |> 
  drop_na()

df_test <- df_test |> 
  select(all_of(c(target, predictors_all))) |> 
  drop_na()

## Train model -----------------------------------------------------------------
# ranger() crashes when using tibbles, so we are using the
# base R notation to enter the data
rf_basic <- ranger(
  y = df_train[, target],     # target variable
  x = df_train[, predictors_all], # Predictor variables
  seed = 42, # Specify the seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1
) # Use all but one CPU core for quick model training

# Print a summary of fitted model
print(rf_basic)

# Interpretation: 
# - ranger automatically considered it as a classification problem
# - reports OOB prediction error = 21.32%. According to documentation: "For classification this is accuracy (proportion of misclassified observations)"

## Test model ------------------------------------------------------------------
# Make predictions for validation sites
prediction <- predict(
  rf_basic, # RF model
  data = df_test
)

# Save predictions to validation df
df_test$pred <- prediction$predictions

# Evaluate test results - standard metrics
out_metrics <- yardstick::metrics(df_test, truth = "waterlog.100", estimate = "pred")
out_metrics

# full confusion matrix
out_cm <- yardstick::conf_mat(df_test, truth = "waterlog.100", estimate = "pred")
out_cm

# Exercise 5.2 Variable selection  ---------------------------------------------

# Select predictors using Boruta
# run the algorithm
bor <- Boruta::Boruta(
  y = df_train[, target],
  x = df_train[, predictors_all],
  maxRuns = 50, # Number of iterations. Set to 30 or lower if it takes too long
  num.threads = parallel::detectCores() - 2
)

# obtain results: a data frame with all variables, ordered by their importance
df_bor <- Boruta::attStats(bor) |>
  tibble::rownames_to_column() |>
  arrange(desc(meanImp))

# plot the importance result
ggplot(
  aes(
    x = reorder(rowname, meanImp),
    y = meanImp,
    fill = decision
  ),
  data = df_bor
) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = c("grey30", "tomato", "grey70")) +
  labs(
    y = "Variable importance",
    x = "",
    title = "Variable importance based on Boruta"
  ) +
  theme_classic() +
  coord_flip()

# get retained important variables
predictors_selected <- df_bor |>
  filter(decision == "Confirmed") |>
  pull(rowname)

# re-train Random Forest model
rf_bor <- ranger(
  y = df_train[, target], # target variable
  x = df_train[, predictors_selected], # Predictor variables
  seed = 42, # Specify the seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1
) # Use all but one CPU core for quick model training

# Make predictions for validation sites
prediction <- predict(
  rf_bor, # RF model
  data = df_test, # Predictor data
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation df
df_test$pred_bor <- prediction$predictions

# Evaluate test results - standard metrics
out_merics_bor <- yardstick::metrics(df_test, truth = "waterlog.100", estimate = "pred_bor")
out_merics_bor

# full confusion matrix
out_cm_bor <- yardstick::conf_mat(df_test, truth = "waterlog.100", estimate = "pred_bor")
out_cm_bor

# Interpretation: generalises better!

# Looking at OOB:
rf_basic
rf_bor

# => Considering OOB, the model with reduced predictors would also be identified
# as better. OOB prediction error of rf_bor is 20.66%
# while the one for rf_basic is larger with 21.95%. 
# When selecting the model with lower OOB prediction error we would thus also 
# select rf_bor (in agreement with standard metrics).

# Exercise 5.3 Model optimization-------------------------------------------------
# Reduce dataset for clear formula specification
df_train <- df_train |> 
  select(all_of(c(target, predictors_selected)))

df_test <- df_test |> 
  select(all_of(c(target, predictors_selected)))

# Model and pre-processing formulation, use all variables but LW_IN_F
pp <- recipes::recipe(
  waterlog.100 ~ .,
  data = df_train 
  )

# Fit random forest model
mtry_default <- floor(sqrt(length(predictors_selected)))
rf_caret <- caret::train(
  pp,
  data = df_train,
  method = "ranger",
  trControl = caret::trainControl(
    method = "cv", 
    number = 5, 
    savePredictions = "final"
    ),
  tuneGrid = expand.grid(
    .mtry = c(mtry_default, mtry_default + 4,  mtry_default + 6,  mtry_default + 8),
    .min.node.size = c(7, 8, 9 , 10, 11, 12, 13),
    .splitrule = "gini"
  ),
  seed = 1982 # for reproducibility
)

# Make predictions for validation sites
prediction <- predict(
  rf_caret, # RF model
  newdata = df_test
)

# Save predictions to validation df
df_test$pred_caret <- prediction

# Evaluate test results - standard metrics
out_metrics_opt <- yardstick::metrics(df_test, truth = "waterlog.100", estimate = "pred_caret")
out_metrics_opt

# full confusion matrix
out_cm_opt <- yardstick::conf_mat(df_test, truth = "waterlog.100", estimate = "pred_caret")
out_cm_opt

# Interpretation: generalises better!

# Exercise 5.4 Probabilistic predictions----------------------------------------
# Fit random forest model
rf_prob <- ranger(
  y = df_train[, target], # target variable
  x = df_train[, predictors_selected], # Predictor variables
  mtry = rf_caret$bestTune$mtry,
  min.node.size = rf_caret$bestTune$min.node.size,
  splitrule = rf_caret$bestTune$splitrule,
  seed = 1982,
  probability = TRUE
  )

# Make predictions for validation sites
prediction <- predict(
  rf_prob, # RF model
  data = df_test
)

# Save predictions to validation df
df_test <- bind_cols(
  df_test,
  prediction$predictions |> 
    as_tibble() |> 
    setNames(c("is_false_waterlog.100", "is_true_waterlog.100"))
)

# ROC curve
yardstick::roc_curve(
  data = df_test,
  truth = waterlog.100, # this has two levels: 0 ('first' level) and 1 ('second' level)
  # variant 1 (wrong: results in sensitivity for non-water-logging, since event_level == "first"):
  # event_level = "first", # this is the default
  # is_false_waterlog.100 # probabilities for class 0 (specified by event_level, default='first', i.e. 0)
  
  # variant 2 (right: gives sensitivity for water-logging, since event_level = "second"):
  event_level = "second",
  is_true_waterlog.100 # probabilities for class 1 (specified by event_level, here='second', i.e. 1)
) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_text(data = function(df){df |> slice_sample(n=5)},
            aes(label = sprintf("%.2f",.threshold)), color = "red", hjust=0, vjust=1) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw()
