library(quantmod)
library(data.table)
library(stochvol)
#提取数据
ssec <- getSymbols("^SSEC", from = "2017-01-01", to = "2017-01-31", auto.assign = F)
sp500 <- getSymbols("^GSPC", from = "2017-01-01", to = "2017-01-31", auto.assign = F)
ft100 <- getSymbols("^FTSE", from = "2017-01-01", to = "2017-01-31", auto.assign = F)
n225 <- getSymbols("^N225", from = "2017-01-01", to = "2017-01-31", auto.assign = F)
bvsp <- getSymbols("^BVSP", from = "2017-01-01", to = "2017-01-31", auto.assign = F)

ssec <- as.data.table(ssec[, 6])
sp500 <- as.data.table(sp500[, 6])
ft100 <- as.data.table(ft100[, 6])
n225 <- as.data.table(n225[, 6])
bvsp <- as.data.table(bvsp[, 6])

data_handle <- function(x, y) {
  lx <- nrow(x)
  ly <- nrow(y)
  if (lx > ly) {
    x_y <- x[index %in% y$index]
    temp <- y[index %in% x_y$index]
    x_y <- cbind(x_y, temp[, -1])
  }
  else {
    x_y <- y[index %in% x$index]
    temp <- x[index %in% x_y$index]
    x_y <- cbind(x_y, temp[, -1])
  }
}
a <- data_handle(sp500, ssec)
b <- data_handle(ft100, a)
c <- data_handle(n225, b)
new_data <- data_handle(bvsp, c)
write.csv(new_data, "Realtest.csv")
new_data <- read.csv("Realtest.csv")

bvsp_ret <- logret(new_data$BVSP.Adjusted)
n225_ret <- logret(new_data$N225.Adjusted)
ft100_ret <- logret(new_data$FTSE.Adjusted)
sp500_ret <- logret(new_data$GSPC.Adjusted)
ssec_ret <- logret(new_data$SSEC.Adjusted)
new_ret <- data.frame(BVSP_ret = bvsp_ret, SP500_ret = sp500_ret, FT100_ret = ft100_ret,
                      N225_ret = n225_ret, SSEC_ret = ssec_ret)

income_GAS <- sum(as.matrix(new_ret[1:5,]) %*% port_GAS[[1]]) # 5%
income_GAS_l <- sum(as.matrix(new_ret[1:5,]) %*% port_GAS_l[[1]]) # 1%
income_GAS_u <- sum(as.matrix(new_ret[1:5,]) %*% port_GAS_u[[1]]) # 10%

income_sVine <- sum(as.matrix(new_ret[1:5,]) %*% port_sVine[[1]]) # 5%
income_sVine_l <- sum(as.matrix(new_ret[1:5,]) %*% port_sVine_l[[1]])
income_sVine_u <- sum(as.matrix(new_ret[1:5,]) %*% port_sVine_u[[1]])

income_t <- sum(as.matrix(new_ret[1:5,]) %*% port_t[[1]]) # 5%
income_t_l <- sum(as.matrix(new_ret[1:5,]) %*% port_t_l[[1]])
income_t_u <- sum(as.matrix(new_ret[1:5,]) %*% port_t_u[[1]])
