library(tidyverse)

# read lambda/alpha estimates
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

# predicted.values <- test1$`lambdas$co`
# 
# # this value cannot be compared to a real one, as it merges all effects into one alpha
# 
# alpha2.plot <- ggplot(predicted.values, aes(x = splist, y = alpha_est2)) + 
#   geom_errorbar(aes(ymin=alpha.lower.error.2, ymax = alpha.upper.error.2), width = 0.05, size = 0.3) +
#   geom_point(size = 2) + 
#   ggtitle("alpha from model 2") + 
#   NULL
# alpha2.plot

##########
# models 3,4,5
# convert to tidy format
alpha.values <- NULL
for(i.model in c(3,4,5)){
  
  if(i.model == 3){
    my.predicted.values <- test2$alpha_matrix3
    my.upper.errors <- test2$alpha_upper_error3
    my.lower.errors <- test2$alpha_lower_error3
  }else if(i.model == 4){
    my.predicted.values <- test2$alpha_matrix4
    my.upper.errors <- test2$alpha_upper_error4
    my.lower.errors <- test2$alpha_lower_error4
  }else if(i.model == 5){
    my.predicted.values <- test2$alpha_matrix5
    my.upper.errors <- test2$alpha_upper_error5
    my.lower.errors <- test2$alpha_lower_error5
  }
  
  predicted.values <- as.data.frame(my.predicted.values)
  predicted.values <- gather(predicted.values,key = "sp",value = "alpha")
  
  upper.errors <- as.data.frame(my.upper.errors)
  upper.errors <- gather(upper.errors,key = "sp",value = "alpha.upper.error")
  
  lower.errors <- as.data.frame(my.lower.errors)
  lower.errors <- gather(lower.errors,key = "sp",value = "alpha.lower.error")
  
  temp.alpha.values <- left_join(predicted.values,upper.errors)
  temp.alpha.values <- left_join(temp.alpha.values,lower.errors)
  
  temp.alpha.values$model <- i.model
  
  alpha.values <- rbind(alpha.values, temp.alpha.values)
}

alpha.values$focal.sp <- 1
alpha.values$model <- as.factor(alpha.values$model)

# append true values
true.values <- data.frame(sp = c("sp1","sp2","sp3"), alpha = c(a_sp1,a_sp2,a_sp3), alpha.upper.error = 0, alpha.lower.error = 0,
                          model = "truth", focal.sp = 1)

alpha.values <- rbind(alpha.values, true.values)

dodge <- position_dodge(.3)
alpha.plot <- ggplot(alpha.values, aes(x = sp, y = alpha, group = model)) + 
  geom_errorbar(aes(ymin=alpha.lower.error, ymax = alpha.upper.error), width = 0.05, size = 0.3, position = dodge) +
  geom_point(aes(color = model), size = 2, position = dodge) + 
  facet_grid(focal.sp~.) + 
  NULL
alpha.plot


