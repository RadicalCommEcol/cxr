
library(ggplot2)

# change
lambda.values <- read.csv2("./data/lambda_values.csv",stringsAsFactors = FALSE)
head(lambda.values)
str(lambda.values)

lambda.values$lambda <- as.numeric(lambda.values$lambda)

# I use ggplot2, which has a syntax in which you add "layers" to a plot sequentially
lambda.plot <- ggplot(lambda.values) + 
  geom_point(aes(x = species, y = lambda)) + 
  facet_grid(~model)+
  theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust=1, size = 6))+
  NULL
lambda.plot

pdf(file = "./data/lambda_values.pdf",width = 7,height = 6)
lambda.plot
dev.off()