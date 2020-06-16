library(caret)

data(ChemicalManufacturingProcess, package = "AppliedPredictiveModeling")

set.seed(100)
training.samples <- createDataPartition(ChemicalManufacturingProcess, p = 0.8, list = FALSE)
train.data  <- ChemicalManufacturingProcess[training.samples,]
test.data <- ChemicalManufacturingProcess[-training.samples,]

y = train.data$Yield

pp.missing.data.model <- preProcess(train.data, method="medianImpute")
train.data <- predict(pp.missing.data.model, train.data)

pp.range.model <- preProcess(train.data, method=c("scale", "center", "YeoJohnson", "nzv"))
train.data <- predict(pp.range.model, train.data)

m = dim(train.data)[1]
n = dim(train.data)[2]

train.data$Yield = y

library(TDAmapper)
library(igraph)

X <- train.data[,2:n]
D <- dist(X)
mds <- cmdscale(D, k = 2)
f.values <- list(as.double(mds[,1]), as.double(mds[,2]))

mapper <- mapper2D(
  distance_matrix = D,
  filter_values = f.values,
  num_intervals = c(14,14),
  percent_overlap = 80,
  num_bins_when_clustering = 10)

g <- graph.adjacency(mapper$adjacency, mode="undirected")
plot(g , layout = layout.auto(g))

library(networkD3)

nodes <- mapperVertices(mapper, 1:dim(X)[1])
links <- mapperEdges(mapper)

node.yields <- c()
for (i in 1:length(nodes$Nodename)) {
  node <- nodes$Nodename[i]
  x <- strsplit(toString(node),":")[[1]][2]
  indicies <- as.numeric(strsplit(x,split=", ",fixed=TRUE)[[1]])
  node.yields <- c(node.yields, mean(train.data[indicies,]$Yield))
}

range <- max(node.yields) - min(node.yields)
min = as.numeric(min(node.yields))

node.yields[node.yields <= min + (1/3)*range] = "Low"
node.yields[min + (1/3)*range < node.yields & node.yields <= min + (2/3)*range] = "Medium"
node.yields[node.yields != "Low" & node.yields != "Medium"] = "High"

nodes$Nodegroup <- node.yields

network <- forceNetwork(Nodes = nodes, Links = links,
                        Source = "Linksource", Target="Linktarget",
                        Value = "Linkvalue", NodeID="Nodename",
                        Group = "Nodegroup", opacity=0.8,
                        zoom = TRUE, Nodesize="Nodesize",
                        legend = TRUE)

library(magrittr)

saveNetwork(network, file = 'Net1.html')

# From graph extract 2 low yield sub-graphs and 2 high yield sub-graphs

get.data.for.group <- function(group) {
  group.data = as.data.frame(matrix(nrow=0, ncol=dim(train.data)[2]))
  colnames(group.data) <- colnames(train.data)
  for (i in group) {
    node <- nodes$Nodename[i]
    x <- strsplit(toString(node),":")[[1]][2]
    indicies <- as.numeric(strsplit(x,split=", ",fixed=TRUE)[[1]])
    group.data <- rbind(group.data, train.data[indicies,])
  }
  return(group.data)
}

H1.data <- get.data.for.group(c(158, 159, 185, 157, 182, 200))
H2.data <- get.data.for.group(c(59, 58, 60, 61, 20))
L1.data <- get.data.for.group(c(290, 273, 275, 271, 269, 220))
L2.data <- get.data.for.group(c(71, 115, 144, 100, 140))

# Now we need to do KS tests between H1/L1, H1/L2, H2/L1, H2/L2.


results <- as.data.frame(matrix(nrow=n, ncol=2))
colnames(results) <- c("Feature", "KS")
for (i in 1:n) {
  test1 <- ks.test(H1.data[,i], L1.data[,i])
  test2 <- ks.test(H1.data[,i], L2.data[,i])
  test3 <- ks.test(H2.data[,i], L1.data[,i])
  test4 <- ks.test(H2.data[,i], L2.data[,i])
  
  results[i,1] = colnames(train.data)[i]
  results[i,2] = max(test1$statistic, test2$statistic, test3$statistic, test4$statistic)
}
results <- results[2:n,]


# print features with greater than 0.90 KS score

print(results[results$KS >= 0.9,])



