library(rcompanion)
library(rrcov)
data <- get(data(OsloTransect))
data <- data[complete.cases(data), ]

locs <- as.matrix(data[c(5,6)])

y1 <- as.vector(blom(data$Pb))
y2 <- as.vector(blom(data$Mn))
y <- c(y1,y2)

locs <- rbind(locs,locs)
locs <- cbind(locs, 1)
locs[351:700,3] <- 2

data <- matrix(0, 700,4)
data[,1] <- y
data[,2:4] <- locs


save(data,file =  "../data/new.RData")
