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
  num_bins_when_clustering = 5)

g <- graph.adjacency(mapper$adjacency, mode="undirected")
plot(g , layout = layout.auto(g))

library(networkD3)

nodes <- mapperVertices(mapper, 1:dim(X)[1])
links <- mapperEdges(mapper)

node.yields <- c()
for (i in 1:length(nodes$Nodename)) {
  n <- nodes$Nodename[i]
  x <- strsplit(toString(n),":")[[1]][2]
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
                        zoom = TRUE, Nodesize="Nodesize")

library(magrittr)

saveNetwork(network, file = 'Net1.html')
