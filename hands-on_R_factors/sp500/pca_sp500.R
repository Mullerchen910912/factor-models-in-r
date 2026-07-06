
#用 S&P500 裡８９家公司的日報酬資料作為例子。
#資料時間為 2013/01/02 - 2022/06/22 (九年半),　
#總共 ２１萬筆資料。
#前３個 Principal Components 
#可以解釋資料 ７３％ 的變動。
#接著建構 統計因子模型（Statistical Factor Model)，
#common factors 數目設定為３個。
#Factor1 影響「醫療產業」
#Factor2 影響「科技產業」
#Factor3 影響「能源產業」


rm(list=ls())


sp500_prices <- read.csv('/Users/chenshengyuan/Desktop/NTU_ECON/Course/Advanced econometric topic/factor_models/hands-on_R_factors/sp500/sp500_prices.csv')

library(tidyverse)

# compute the simple returns (percentage change in prices) and drop NA
rtn <- sp500_prices %>%
  arrange(date) %>%
  mutate(ret = (adjusted - lag(adjusted)) / lag(adjusted)) %>%
  select(date, symbol,ret) %>% 
  drop_na(ret) 

# log returns (continuously compounded returns)
#returns_apple <- log(returns_apple + 1)

# 2013/01/02 - 2022/06/22, NT = 2385*89 = 212,265
rtn_wide <- rtn %>% spread(symbol, ret) %>% drop_na
rtn_wide <- rtn_wide[,-1]


########################################################################
# PCA
pca_rtn <- prcomp(rtn_wide, center = TRUE, scale. = TRUE)

# PC1 (the egienvector corresponding to largest value of eigenvalues)
pca_rtn$rotation
summary(pca_rtn)
# scree plot for eigenvalues
# type = c("barplot", "lines")
screeplot(pca_rtn, type="lines")

# R-square-type interpretation
summary(pca_rtn)

stat_fac <- factanal(rtn_wide, factors=3, 
                     method="mle", rotation="varimax", lower = 0.01)
stat_fac


par(mfrow = c(3, 1))
#adjust plot margins
par(mar = c(3, 3, 3, 3))

loadings_factor_1 <- sort(stat_fac$loadings[, 1], decreasing = TRUE)
barplot(loadings_factor_1[1:10], main="Factor1,  loadings")

loadings_factor_2 <- sort(stat_fac$loadings[, 2], decreasing = TRUE)
barplot(loadings_factor_2[1:10], main="Factor2,  loadings")

loadings_factor_3 <- sort(stat_fac$loadings[, 3], decreasing = TRUE)
barplot(loadings_factor_3[1:10], main="Factor3,  loadings")




