library(data.table)
library(dplyr)
library(tidyverse)
library(devtools)
library(ggbiplot)
options(scipen = 100)


library(ggpubr)
library(factoextra)

pca_data <- read.table("C:/Users/tjarva/Desktop/orthofinder_selection/pca_data.txt", header = TRUE, sep = "\t")


###############################################################################################

dt <- data.table(#gene  = character(534L), 
  #aveexpress = numeric(534L),
  brain = numeric(543L), liver = numeric(543L), kidney = numeric(543L), ovary = numeric(543L), skin= numeric(543L),
  gravy= numeric(543L), omega = numeric(543L), percent= numeric(543L)) 

#get gene names to set row names
namevec <- pca_data$gene
for(i in 1:length(namevec)){
  rownames(dt)[i] <- namevec[i]
}


#dt$gene <- data$gene
dt$brain <- pca_data$brain
dt$liver <- pca_data$liver
dt$kidney <- pca_data$kidney
dt$ovary <- pca_data$ovary
dt$skin  <- pca_data$skin
dt$gravy <- as.numeric(pca_data$gravy)
dt$omega <- as.numeric(pca_data$omega)
dt$percent <- as.numeric(pca_data$percent)


plot(dt[,1:8])

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(dt, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

set.seed(123)

fviz_nbclust(dt, kmeans, method = "wss")

k2 <- kmeans(dt, centers = 2, nstart = 25)
k3 <- kmeans(dt, centers = 3, nstart = 25)
k4 <- kmeans(dt, centers = 4, nstart = 25)
k5 <- kmeans(dt, centers = 5, nstart = 25)
p1 <- fviz_cluster(k2, geom = "point", data = dt) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = dt) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = dt) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = dt) + ggtitle("k = 5")


library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

#Get principal component vectors using prcomp instead of princomp
pc <- prcomp(dt)
# First for principal components
comp <- data.frame(pc$x[,1:8])
k <- kmeans(comp, 3, nstart=25, iter.max=1000)


# Compute k-means with k = 3
set.seed(123)
res.km <- kmeans(scale(dt), 3, nstart = 25)
# K-means clusters showing the group of each individuals
res.km$cluster


fviz_cluster(res.km, data = dt,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#C83D95"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


# Dimension reduction using PCA
res.pca <- prcomp(dt,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$Gene <- data$gene
# Data inspection
head(ind.coord)


write.table(ind.coord, file = "pca_genes_clusters.txt", sep = '\t', row.names = FALSE, quote = FALSE, append = FALSE)

eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

library(ggpubr)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2",
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cluster", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

####################################################################################################
# other R tutorial
library(factoextra)
pca <- prcomp(dt, scale = TRUE, center = TRUE)
plot(pca)
summary(pca)

fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
g <- ggbiplot(pca, labels = rownames(dt),
              obs.scale = 1,
              var.scale = 1,
              groups = dt$percent,
              circle = TRUE)

g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)

#compuate the average silhouette width
library(cluster)

avg_sil <- function(k) {
  km.res <- kmeans(dt, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(dt))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")
fviz_nbclust(dt, kmeans, method = "silhouette")

set.seed(123)
gap_stat <- clusGap(dt, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

fviz_gap_stat(gap_stat)
