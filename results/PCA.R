install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
data(decathlon2)
decathlon2.active <- decathlon2[1:23, 1:10]
head(decathlon2.active[, 1:6], 4)
test <- read.csv('~/git_repos/upec_mouse_model/test2.csv', row.names = 1)

library("FactoMineR")
res.pca <- PCA(test, graph = FALSE)

eig.val <- get_eigenvalue(res.pca)
eig.val

var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
     

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
     