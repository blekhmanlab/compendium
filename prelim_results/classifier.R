library(ggfortify)
setwd('~/code/shithouse/prelim_results')

x <- read.csv('studies_consolidated.tsv', sep='\t')
reads <- as.data.frame(x)
# get sample names
rownames(reads) <- reads$X
reads$X <- NULL
reads <- reads[rowSums(reads) > 1000,] # get rid of empty samples
# taxa must have > 70 reads across ALL samples
reads <- reads[,colSums(reads) > 70]
# taxa must appear in a minimum number of samples
#reads = reads[,!sapply(df, function(x) mean(x==0))>0.1]


library(dplyr)    # for data manipulation
library(tidyr)

# CLR TRANSFORM VERSION
# fix zeroes
library(zCompositions)

load(file="reads.nozero.Rda")
#reads.nozero <- cmultRepl(reads, label=0, method='SQ', output="p-counts") # default method won't work
#save(reads.nozero, file="reads.nozero.Rda")
reads.nozero <- unique(reads.nozero) # GET RID OF DUPLICATES

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
clr <- apply(reads.nozero, 1, function(a) log(a/gm_mean(a)))
clr <- as.data.frame(t(clr))

# trim the taxa to show only genus
clr.trimmedname <- clr
colnames(clr.trimmedname) <- gsub('.+\\.+(\\w+)$', '\\1', colnames(clr.trimmedname))
# trim the samples to show only sample name (no project)
rownames(clr.trimmedname) <- gsub('^[^_]+_','', rownames(clr.trimmedname))

# bring in country names
allcountries <- read.csv('sample_countries.csv')

samplenames <- as.data.frame(rownames(clr.trimmedname))
colnames(samplenames) <- c('srr')
annotated <- samplenames %>% left_join(allcountries, by = "srr")



library(ggplot2)  # for data visualization
library(h2o)  # for fitting autoencoders

h2o.init(max_mem_size = "10g")
variance <- sapply(clr.trimmedname, sd)
hist(variance)
sum(variance > 0.35)
# only include the ~550 OTUs with the highest variance
clr.high <- clr.trimmedname[,variance > 0.35]

# set aside test samples
set.seed(123)
index <- sample(1:nrow(clr.high), floor(nrow(clr.high) * 0.2))
sampled_digits <- clr.high[index, ]
training <- clr.high[-index,]
features <- as.h2o(training)

# test all hyperparams together:
hyper_grid <- list(
  sparsity_beta = c(0.01, 0.1, 0.2, 0.5),
  hidden = list(
    c(25),
    c(80),
    c(100),
    c(50,15,50),
    c(100, 50, 100)
  ),
  average_activation = c(-0.5, -0.1, 0, 0.1, 0.5)
)

# Execute grid search
big_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = 'big_grid',
  autoencoder = TRUE,
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  ignore_const_cols = FALSE,
  seed = 123
)
big_grid <- h2o.getGrid('big_grid', sort_by = 'mse', decreasing = FALSE)
best_model_id <- big_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# get the LVs!!!
new_codings <- h2o.deepfeatures(best_model, features, layer = 2)
new_codings <- as.data.frame(new_codings)
rownames(new_codings) <- rownames(training)

library(Rtsne)
set.seed(1234)
tsned <- Rtsne(as.matrix(new_codings), perplexity=45,
               pca=FALSE, max_iter=3000, num_threads=6,
               theta=0.5)

# make a list of all countries for the samples
countries_toplot <- as.data.frame(rownames(new_codings))
colnames(countries_toplot) <- c('srr')
countries_toplot <- countries_toplot %>% left_join(allcountries, by = "srr")

toplot <- as.data.frame(tsned$Y)
ggplot(toplot, aes(x=V1, y=V2, color=countries_toplot$standard)) +
  geom_point() +
  theme_bw()

# add countries to the actual data
rftrain <- new_codings
rftrain$srr <- rownames(new_codings)
rftrain <- rftrain %>% inner_join(allcountries, by = "srr")

# add continent annotation
continents <- read.csv('continents.csv')
rftrain <- rftrain %>% inner_join(continents, by = "standard")
rftrain <- rftrain[rftrain$continent != '',]
rftrain$continent <- as.factor(rftrain$continent)

predictors <- colnames(rftrain)[1:50]
response <- "continent"

rftrain <- as.h2o(rftrain)

rfsplit <- h2o.splitFrame(data = rftrain, ratios = 0.8, seed = 1234)
train <- rfsplit[[1]]
valid <- rfsplit[[2]]

# Build and train the model:
cars_drf <- h2o.randomForest(x = predictors,
                             y = response,
                             ntrees = 50,
                             max_depth = 20,
                             min_rows = 10,
                             balance_classes = TRUE,
                             training_frame = train,
                             validation_frame = valid)

# Eval performance:
perf <- h2o.performance(cars_drf)
perf

################################
# Pull in the REAL validation dataset:
topredict <- as.h2o(sampled_digits)
lockbox <- h2o.deepfeatures(best_model, topredict, layer = 2)
lockbox <- as.data.frame(lockbox)
rownames(lockbox) <- rownames(sampled_digits)

lockbox$srr <- rownames(lockbox)
lockbox <- lockbox %>% inner_join(allcountries, by = "srr")

continents <- read.csv('continents.csv')
lockbox <- lockbox %>% inner_join(continents, by = "standard")
lockbox <- lockbox[lockbox$continent != '',]
lockbox$continent <- as.factor(lockbox$continent)
lockbox <- as.h2o(lockbox)

predict <- h2o.predict(cars_drf, newdata = lockbox)
z <- as.data.frame(predict)
sum(predict$predict == lockbox$continent)
