library(dplyr)

setwd("C:/Users/tjarva/Desktop")

files <- c(list.files(pattern = "stats"))

library(data.table)

# Load data sets

a1_kidney <- read.table(files[1])
a1_medulla <- read.table(files[2])
a1_ovary <- read.table(files[3])
a3_kidney <- read.table(files[4])
a4_diencephalon <- read.table(files[5])
a4_liver <- read.table(files[6])
a1_liver <- read.table(files[7])
a1_metencephalon <- read.table(files[8])
a1_telencephalon <- read.table(files[9])
a3_skin <- read.table(files[10])
a4_kidney <- read.table(files[11])


#join expression values from all tissue samples (this is for ease of selecting a different sample if need be)
join1 <- full_join(a1_kidney, a1_medulla, by = "V1", keep = FALSE) 
join1 <- join1 %>% select(-c(V4.x, V4.y, V2.y))
join1 <- join1 %>% rename("length" = V2.x, "a1.kidney" = V3.x,"a1.medulla" = V3.y)
join2 <- full_join(join1, a1_ovary, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A1.ovary" = V3)
join3 <- full_join(join2, a1_liver, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A1.liver" = V3)
join4 <- full_join(join3, a1_metencephalon, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A1.metencephalon" = V3)
join5 <- full_join(join4, a1_telencephalon, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A1.telencephalon" = V3)

join6 <- full_join(join5, a3_kidney, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A3.kidney" = V3)
join7 <- full_join(join6, a3_skin, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A3.skin" = V3)

join8 <- full_join(join7, a4_diencephalon, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A4.diencephalon" = V3)
join9 <- full_join(join8, a4_kidney, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A4.kidney" = V3)

join10 <- full_join(join9, a4_liver, by = "V1", keep = FALSE) %>% select(-c(V2, V4)) %>% rename("A4.liver" = V3)

rm(join1, join2, join3, join4, join5, join6, join7, join8, join9)
rm(a1_telencephalon, a1_kidney, a1_liver, a1_medulla, a1_metencephalon, a1_ovary, a3_kidney, a3_skin, a4_diencephalon, a4_kidney, a4_liver)


gene_file <- read.table("C:/Users/tjarva/Desktop/get2", fill = TRUE, header = FALSE, col.names= paste0("V", seq_len(7)),  sep = "\t")

keycol <- "transcript"
valuecol <- "gene"
gathercols <- c("V2", "V3", "V4", "V5", "V6", "V7")

newdf <- gather(gene_file, keycol, valuecol, gathercols)
newdf <- newdf %>% mutate_all(na_if,"") %>% na.omit() %>% select(-c(keycol))
express <- rename(express, "valuecol" = V1)

all_join <- full_join(newdf, express, by = "valuecol", keep = FALSE)
write.table(all_join, file = "all_catshark_expression", quote = FALSE, row.names = FALSE, col.names = TRUE)

all_join <- read.table("C:/Users/tjarva/Desktop/all_catshark_expression", header=TRUE, sep = " ")


#### normalize catshark expression data
#a1.kidney, a1.medulla, a1.ovary, a1.liver, a1.metencephalon,
#a1. telencephalon, a3.kidney,
#a3.skin, a4.diencephalon, 
#a4.kidney, a4.liver
libsize <- c(59932, 60497, 71234, 27745, 62688, 
             68090, 59687,
             58960, 66717,
             62076, 41583)
#express <- data.frame(expression_file)
scipen = 100
#calculate RPKM of genes of interest, store in data frame to write out
rpkm_vec <- data.frame(name = character(811L), gene = character(811L), length = numeric(811L), a1.kidney = numeric(811L), a1.medulla  = numeric(811L), a1.ovary = numeric(811L),
                       a1.liver  = numeric(811L), a1.metencephalon  = numeric(811L), a1.telencephalon  = numeric(811L), 
                       a3.kidney  = numeric(811L), a3.skin  = numeric(811L), a4.diencephalon = numeric(811L), 
                       a4.kidney  = numeric(811L), a4.liver  = numeric(811L))
for(i in 1:nrow(all_join)){
   gene_length <- as.integer(all_join[i,3])
   name <- as.character(all_join[i,2])
   gene <- as.character(all_join[i,1])
   rpkm <- vector(mode="numeric", length = 11L)
   for(j in 4:14){
      index <- j - 3
      numReads <- as.integer(all_join[i, j])
      rpkm[index] <- numReads/(gene_length/1000*libsize[index]/100000)
      
         }
   rpkm_vec[i,] <- list(name, gene, gene_length, rpkm[1], rpkm[2], rpkm[3], rpkm[4], rpkm[5], rpkm[6], rpkm[7], rpkm[8], rpkm[9], rpkm[10], rpkm[11])
}
final_rpkm <- rpkm_vec

#write.table(final_rpkm, file = "catshark_rpkm_express.txt", quote = FALSE, append = FALSE, row.names = FALSE, col.names = TRUE)

####### choose the longest length for each catshark gene
catGeneList <- final_rpkm$gene %>% unique() %>% na.omit()

results <- data.frame()
for (g in catGeneList){
   current <- final_rpkm %>%
      filter(gene == g) %>%
      arrange(desc(length))
   results <- rbind(results, current[1,])
}

results <- results %>% select(-c(length, a1.metencephalon, a1.telencephalon, a3.kidney, a4.diencephalon, a4.kidney, a4.liver)) %>% 
   rename("kidney"= a1.kidney, "brain" = a1.medulla, "ovary"= a1.ovary, "liver" = a1.liver, "skin" = a3.skin)

# later manually edited to change transcript names to match sawfish
saw_rpkm <- read.table("C:/Users/tjarva/Desktop/542Expression", header = TRUE, sep = "\t")



#####################################################################################################
library(data.table)

data <- read.table("C:/Users/tjarva/Documents/catsaw_expression", header = TRUE, sep = "\t")

#### find the differences in expression for each tissue
data <- data %>% mutate(diff.brain = brain.x - brain.y, diff.liver = liver.x - liver.y, diff.kidney= kidney.x - kidney.y, 
                        diff.ovary= ovary.x - ovary.y, diff.skin = skin.x - skin.y) %>%
   select(-c(brain.x, brain.y, liver.x, liver.y, kidney.x, kidney.y, ovary.x, ovary.y, skin.x, skin.y))

## get gravy, omega, percent sites data from pca_data
pca_data <- read.table("C:/Users/tjarva/Desktop/pca_data.txt", header = TRUE, sep = "\t")

pca_data <- select(pca_data, -c(brain, liver, kidney, ovary, skin))
## join expression change data with pca_data

new_join <- full_join(pca_data, data, by = "gene")

new_join <- new_join %>% select(-c(name.x, name.y))


dt <- data.table(diff.brain = numeric(543L), diff.liver = numeric(543L), diff.kidney = numeric(543L), 
   diff.ovary = numeric(543L), diff.skin= numeric(543L),
   gravy= numeric(543L), omega = numeric(543L), percent= numeric(543L)) 

#get gene names to set row names
namevec <- new_join$gene
for(i in 1:length(namevec)){
   rownames(dt)[i] <- namevec[i]
}


#get gene names to set row names

#dt$gene <- data$gene
dt$diff.brain <- new_join$diff.brain
dt$diff.liver <- new_join$diff.liver
dt$diff.kidney <- new_join$diff.kidney
dt$diff.ovary <- new_join$diff.ovary
dt$diff.skin  <- new_join$diff.skin
dt$gravy <- new_join$gravy
dt$omega <- new_join$omega
dt$percent <- new_join$percent

plot(dt[,1:8])


library(ggpubr)
library(factoextra)
library(data.table)
library(dplyr)
library(tidyverse)
library(devtools)
library(ggbiplot)
options(scipen = 100)

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

############################################################

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

#########################################################################
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


write.table(ind.coord, file = "pca_catshark_k3.txt", sep = '\t', row.names = FALSE, quote = FALSE, append = FALSE)

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


