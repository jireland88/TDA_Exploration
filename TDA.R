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

forceNetwork(Nodes = nodes, Links = links,
             Source = "Linksource", Target="Linktarget",
             Value = "Linkvalue", NodeID="Nodename",
             Group = "Nodegroup", opacity=0.8,
             zoom = TRUE, Nodesize="Nodesize")


