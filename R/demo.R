source("fit_multi.R")


# weather example
load("../data/weather.RData")
y_mat <- weather[,c(4,5)]
locs0 <- weather[,c(1,2,3)]
n <- nrow(y_mat)

y1 <- y_mat[,1]
X1 <- as.matrix(rep(1, n))

y2 <- y_mat[,2]
X2 <- as.matrix(rep(1, n))

NNarray <- find_ordered_nn(locs0, 20)


# fit model
fit_bivariate_matern(y1, y2, X1, X2, locs0)


# competition example
data <- read.csv("../data/3a_1_train.csv")
y_mat <- as.matrix(data[1:10000,3:4])
locs0 <- as.matrix(data[1:10000,1:2])
n <- nrow(y_mat) 


y1 <- y_mat[,1]
X1 <- as.matrix(rep(1, n))

y2 <- y_mat[,2]
X2 <- as.matrix(rep(1, n))

NNarray <- find_ordered_nn(locs0, 20)


# fit model
fit_bivariate_matern(y1, y2, X1, X2, locs0)
