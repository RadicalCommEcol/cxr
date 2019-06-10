library(tidyverse)

# read lambda/alpha estimates
# TODO: update the reading of the results
lambda.results <- readr::read_delim("./results/lambda_estimates.csv",delim = ";")
lambda.results$model <- as.factor(lambda.results$model)
lambda.results$focal.sp <- as.factor(lambda.results$focal.sp)
lambda.results$optim.method <- as.factor(lambda.results$optim.method)

load("./results/param_estimates")

#####################
# diagnostic plots

# 1 - lambda for different models and methods
lambda.models <- ggplot(lambda.results, aes(x = model, y = lambda, group = optim.method)) + 
  geom_point(aes(fill = optim.method,shape = optim.method),size = 2) +
  # geom_errorbar(aes(ymin = lambda.lower.error, ymax = lambda.upper.error, color = optim.method), width = 0.2)+
  facet_wrap(focal.sp~., scales = "free_y",ncol = 8)+
  scale_shape_manual(values = c(21,22)) +
  scale_fill_manual(values = c("darkgreen","darkorange"))+
  scale_color_manual(values = c("darkgreen","darkorange"))+
  # ylim()
  NULL
lambda.models
