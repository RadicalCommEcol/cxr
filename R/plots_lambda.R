
library(ggplot2)

# # change
# lambda.values <- read.csv2("./data/lambda_values_NL1.csv",stringsAsFactors = FALSE)
# head(lambda.values)
# str(lambda.values)
# 
# lambda.values$lambda <- as.numeric(lambda.values$lambda)
# 
# # I use ggplot2, which has a syntax in which you add "layers" to a plot sequentially
# lambda.plot <- ggplot(lambda.values) + 
#   geom_point(aes(x = species, y = lambda)) + 
#   facet_grid(~model)+
#   theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust=1, size = 6))+
#   NULL
# print(lambda.plot)
# 
# # pdf(file = "./data/lambda_values.pdf",width = 7,height = 6)
# # lambda.plot
# # dev.off()
lambda_cov.values <- read.csv2("./data/lambda_cov_values.csv",stringsAsFactors = FALSE)
lambda_cov_NL2.values <- read.csv2("./data/lambda_cov_values_NL2num.csv",stringsAsFactors = FALSE)
lambda_cov_NL3.values <- read.csv2("./data/lambda_cov_values_NL3num.csv",stringsAsFactors = FALSE)
lambda_cov_NL_NL2.values <- read.csv2("./data/lambda_cov_NL_values_NL2num.csv",stringsAsFactors = FALSE)
lambda_cov_NL_NL3.values <- read.csv2("./data/lambda_cov_NL_values_NL3num.csv",stringsAsFactors = FALSE)
for (i in 22:40){
  print(lambda_cov_NL3.values$lambda.cov[i])
curve(expr=as.numeric(lambda_cov_NL2.values$lambda.cov[i])*(1-exp(-as.numeric(lambda_cov_NL_NL2.values$lambda.cov_NL[i])*x)),col="red", from=0, to= 2,main=lambda.cov.values$species[i])
curve(expr=as.numeric(lambda_cov.values$lambda.cov[i])*x,col="black", from=0, to= 2,add=TRUE)
curve(expr=as.numeric(lambda_cov_NL3.values$lambda.cov[i])*x^2/(as.numeric(lambda_cov_NL_NL3.values$lambda.cov_NL[i])+x^2),col="blue", from=0, to= 2,add= TRUE)
}