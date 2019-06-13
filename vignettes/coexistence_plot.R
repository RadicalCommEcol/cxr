# plot some pairwise relationships
source("R/AvgFitnessRatio.R")
source("R/NicheOverlap.R")

library(tidyverse)
library(scales)

#################
# load data

# germination and seed survival rates
load("data/species_rates.RData")
# lambda and alpha coefficients
load("../temp/results/temp_param_estimates.Rdata")

sp.names <- sort(unique(names(param.matrices)))
alpha.matrix <- matrix(0,nrow = length(sp.names),ncol = length(sp.names))
lambda.values <- numeric(length(sp.names))
names(lambda.values) <- sp.names
rownames(alpha.matrix) <- colnames(alpha.matrix) <- sp.names
# from which model and optimization method are we gathering parameters
estimates.model <- "BH_4"
estimates.method <- "optim_NM"
for(i.sp in 1:length(sp.names)){
  lambda.values[i.sp] <- param.matrices[[i.sp]][[estimates.model]][[estimates.method]]$lambda
  alpha.matrix[i.sp,] <- param.matrices[[i.sp]][[estimates.model]][[estimates.method]]$alpha[sp.names]
}

#################
# prepare data

# generate a dataframe with pairwise data
# function from https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i){
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

pairwise.ratios <- as.data.frame(expand.grid.unique(sp.names,sp.names,include.equals = FALSE),stringsAsFactors = FALSE)
names(pairwise.ratios) <- c("sp1","sp2")
pairwise.ratios$fitness.ratio <- 0
pairwise.ratios$niche.overlap <- 0

# gather fitness ratio/niche overlap between each pair
for(i.pair in 1:nrow(pairwise.ratios)){
  my.sp <- c(pairwise.ratios$sp1[i.pair],pairwise.ratios$sp2[i.pair])
  ########### TEST ONLY!!! remove ABS
  my.matrix <- abs(alpha.matrix[my.sp,my.sp])
  ###########
  pairwise.ratios$fitness.ratio[i.pair] <- AvgFitnessRatio(lambda = lambda.values[my.sp],
                                                           germ.rate = species_rates$germination[species_rates$code %in% my.sp],
                                                           survival.rate = species_rates$`seed survival`[species_rates$code %in% my.sp],
                                                           pair.matrix = my.matrix)[[3]]
  pairwise.ratios$niche.overlap[i.pair] <- NicheOverlap(pair.matrix = my.matrix)
}

#################
# prepare plot
x.1 = seq(1, 0, by=-0.001)
x.2 = seq(0, -1, by=-0.001)
y1.1 = 1/((1-x.1))
y2.1 = 1-x.1
y1.2 = 1/((1-x.2))
y2.2 = 1-x.2
my.lim.x = 0.5
my.lim.y = 1.7

x.complete <- seq(1,-1,by=-0.001)
y1 <- 1/(1-x.complete)
y2 <- 1-x.complete

lines <- data.frame(x = c(x.complete,x.complete), 
                    y = c(y1,y2), 
                    fun = c(rep("y1", times = length(seq(1, -1, by=-0.001))), 
                            rep("y2", times = length(seq(1, -1, by=-0.001)))))

coexistence.simple <- ggplot() + 
  geom_line(data = lines, aes(x = x, y = y,linetype = fun,group = fun), col = "black") +
  # grid lines
  geom_abline(intercept = 0, slope = 0, lty = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, col = "darkgrey") +
  # points
  geom_point(data = pairwise.ratios, aes(x = (1 - niche.overlap), y = log(fitness.ratio))) +
  # log-scale
  scale_y_log10() +
  # # trim to desired coords
  coord_cartesian(expand = c(0, 0),
                  xlim = c(-my.lim.x,my.lim.x),
                  ylim = c(1/my.lim.y, my.lim.y)) +
  # axes, legend, etc
  xlab(expression(paste("(1 - ", rho, " )", sep = ""))) + 
  ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
  theme(legend.position = "none", 
        axis.title.y = element_text(angle = 90), 
        axis.title = element_text(size = 20)) +
  NULL
coexistence.simple

# this is the complete version, with inner fill

# auxiliary dataframes for the boundaries of coexistence regions
lines.pos = data.frame(x = c(x.1, x.1), 
                     y = c(y1.1, y2.1), 
                     fun = c(rep("y1", times = length(seq(1, 0, by=-0.001))), 
                                rep("y2", times = length(seq(1, 0, by=-0.001)))))
ribbon.pos = data.frame(min.dim = y1.1, max.dim = y2.1, x.dim = x.1)

lines.neg = data.frame(x = c(x.2, x.2), 
                     y = c(y1.2, y2.2), 
                     fun = c(rep("y1", times = length(seq(0, -1, by=-0.001))), 
                                rep("y2", times = length(seq(0, -1, by=-0.001)))))
ribbon.neg = data.frame(min.dim = y1.2, max.dim = y2.2, x.dim = x.2)

#################
# plot
coexistence.plot <- ggplot() + 
  # positive side
  geom_ribbon(data = ribbon.pos, aes(x = x.dim, ymin = min.dim, ymax = max.dim), fill = "darkgrey", alpha = 0.6) +
  geom_line(data = lines.pos, aes(x = x, y = y,linetype = fun,group = fun), col = "black") +
  # negative side
  geom_ribbon(data = ribbon.neg, aes(x = x.dim, ymin = min.dim, ymax = max.dim), fill = "darkgrey", alpha = 0.6) +
  geom_line(data = lines.neg, aes(x = x, y = y,linetype = fun,group = fun), col = "black") +
  # grid lines
  geom_abline(intercept = 0, slope = 0, lty = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, col = "darkgrey") +
  # points
  geom_point(data = pairwise.ratios, aes(x = (1 - niche.overlap), y = log(fitness.ratio))) +
  # log-scale
  scale_y_log10() +
  # # trim to desired coords
  coord_cartesian(expand = c(0, 0),
                  xlim = c(-my.lim.x,my.lim.x),
                  ylim = c(1/my.lim.y, my.lim.y)) +
  # axes, legend, etc
  xlab(expression(paste("(1 - ", rho, " )", sep = ""))) + 
  ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
  theme(legend.position = "none", 
        axis.title.y = element_text(angle = 90), 
        axis.title = element_text(size = 20)) +
  NULL
 coexistence.plot
