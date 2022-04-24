source("fit_multi.R")


# read in data
load("../data/weather.RData")
y_mat <- weather[,c(4,5)]
locs0 <- weather[,c(1,2,3)]

# make zero mean
y_mat[,1] <- y_mat[,1]  - mean(y_mat[,1])
y_mat[,2] <- y_mat[,2]  - mean(y_mat[,2])  

# fit model
fit_bivariate_matern(y_mat, locs0)


# read in data
data <- read.csv("../data/3a_1_train.csv")
y_mat <- as.matrix(data[1:10000,3:4])
locs0 <- as.matrix(data[1:10000,1:2])

# make zero mean
y_mat[,1] <- y_mat[,1]  - mean(y_mat[,1])
y_mat[,2] <- y_mat[,2]  - mean(y_mat[,2])  

# fit model
fit_bivariate_matern(y_mat, locs0)

