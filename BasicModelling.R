library(caret)

data(ChemicalManufacturingProcess, package = "AppliedPredictiveModeling")

set.seed(100)
training.samples <- createDataPartition(ChemicalManufacturingProcess, p = 0.8, list = FALSE)
train.data  <- ChemicalManufacturingProcess[training.samples,]
test.data <- ChemicalManufacturingProcess[-training.samples,]

y = train.data$Yield
ytest = test.data$Yield

library(naniar)
for (i in names(ChemicalManufacturingProcess)) {
  Q1 <- fivenum(train.data[[i]])[2]
  Q2 <- fivenum(train.data[[i]])[3]
  Q3 <- fivenum(train.data[[i]])[4]
  lower <- Q1 - 1.5*(Q3-Q1)
  upper <- Q3 + 1.5*(Q3-Q1)
  print(lower)
  train.data = replace_with_na(train.data, replace = list(i <= lower))
  train.data = replace_with_na(train.data, replace = list(i >= upper))
  
  test.data = replace_with_na(test.data, replace = list(i <= lower))
  test.data = replace_with_na(test.data, replace = list(i >= upper))
}

library(skimr)
skimmed <- skim(train.data)

pp.missing.data.model <- preProcess(train.data, method="medianImpute")
train.data <- predict(pp.missing.data.model, train.data)

pp.range.model <- preProcess(train.data, method=c("YeoJohnson", "nzv"))
train.data <- predict(pp.range.model, train.data)

high.corr.columns <- findCorrelation(cor(train.data))
train.data <- train.data[,-high.corr.columns]

train.data$Yield = y

set.seed(100)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x=train.data[,2:47], y=train.data$Yield, rfeControl = ctrl)

lmProfile

# From this we see top 5 variables are: ManufacturingProcess32, ManufacturingProcess13, BiologicalMaterial03, ManufacturingProcess17, ManufacturingProcess09

# Reduced model
library(dplyr)
train.data.reduced <- select(train.data, Yield, ManufacturingProcess32, ManufacturingProcess13, BiologicalMaterial03, ManufacturingProcess17, ManufacturingProcess09)

lm_1 = train(Yield ~., data=train.data, method="lm")
lm_2 = train(Yield ~., data=train.data.reduced, method="lm")

# put test data through pipeline
test.data <- predict(pp.missing.data.model, test.data)
test.data <- predict(pp.range.model, test.data)
test.data <- test.data[,-high.corr.columns]
test.data$Yield <- ytest
test.data.reduced <- select(test.data, Yield, ManufacturingProcess32, ManufacturingProcess13, BiologicalMaterial03, ManufacturingProcess17, ManufacturingProcess09)

p1 <- predict(lm_1, test.data)
p2 <- predict(lm_2, test.data.reduced)

RMSE(p1, test.data$Yield)
RMSE(p2, test.data.reduced$Yield)

R2(p1, test.data$Yield)
R2(p2, test.data.reduced$Yield)

# So our linear model actually works better using fewer variables !!
