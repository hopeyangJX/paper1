## Model fit using random forest
library(randomForest)
library(randomForestExplainer)
setwd("D:/Writing/Extreme drought/HST.AMF.2017/Data analyses/results/AMF/Model.fit")
data <- read.csv(file = "env.csv", row.names = 1)
head(data)

forest <- randomForest(AMF ~ Plant.richness + Moisture + ANPP + Bbiomass + C + N + AP + CN + pH, 
                       data = data, localImp = TRUE, importance = TRUE)
forest

min_depth_frame <- min_depth_distribution(forest)
head(min_depth_frame)

plot_min_depth_distribution(min_depth_frame)

importance_frame <- measure_importance(forest)
importance_frame

plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes")

plot_multi_way_importance(importance_frame, x_measure = "mse_increase", 
                          y_measure = "node_purity_increase", 
                          size_measure = "p_value", no_of_labels = 10)

importance(forest, type = 1)