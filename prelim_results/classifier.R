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

h2o.init(max_mem_size = "10g")#, nthreads=3) # we should specify threads
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
    c(50),
    c(80),
    c(100)
    #c(100, 50, 100)
  ),
  #average_activation = c(-0.5, -0.1, 0, 0.1, 0.5),
  epochs = c(10,25,50)
)

# Execute grid search
big_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = 'big_grid_epochs',
  autoencoder = TRUE,
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  seed = 123,
  export_weights_and_biases = TRUE
)
big_grid <- h2o.getGrid('big_grid_epochs', sort_by = 'mse', decreasing = FALSE)
best_model_id <- big_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

final_weights <- as.data.frame(h2o.weights(best_model))
# most influential taxa
z <- colMeans(abs(final_weights))


# get the LVs!!!
new_codings <- h2o.deepfeatures(best_model, features, layer = 1)
new_codings <- as.data.frame(new_codings)
rownames(new_codings) <- rownames(training)

# clustering:
d <- dist(heatplot.high, method = "euclidean")
clustering <- hclust(d, method='complete')
plot(clustering)

# t-SNE!
library(Rtsne)
set.seed(1234)
tsned <- Rtsne(as.matrix(new_codings), perplexity=45,
               pca=FALSE, max_iter=1000, num_threads=3,
               theta=0.5)

# make a list of all countries for the samples
countries_toplot <- as.data.frame(rownames(new_codings))
colnames(countries_toplot) <- c('srr')
countries_toplot <- countries_toplot %>% left_join(allcountries, by = "srr")

toplot <- as.data.frame(tsned$Y)
ggplot(toplot, aes(x=V1, y=V2, color=countries_toplot$standard)) +
  geom_point() +
  theme_bw()

# add continent annotation
continents <- read.csv('continents.csv')
countries_toplot <- countries_toplot %>% left_join(continents, by = "standard")
ggplot(toplot, aes(x=V1, y=V2, color=countries_toplot$continent)) +
  geom_point() +
  theme_bw()


# heatmap
heatplot <- new_codings
lv_variance <- sapply(heatplot, sd)
sum(lv_variance > 0.1)
heatplot.high <- as.matrix(heatplot[,lv_variance > 0.1])

heatmap(heatplot.high)

library(gplots)
heatmap.2(heatplot.high,dendrogram="none",
          trace="none",Rowv=FALSE, RowSideColors=countries_toplot$continent)

############################3
### TRAIN RANDOM FOREST MODEL ON LVs
####################3
# add countries to the actual data
rftraindata <- new_codings
rftraindata$srr <- rownames(new_codings)
rftraindata <- rftraindata %>% inner_join(allcountries, by = "srr")

# add continent annotation
rftraindata <- rftraindata %>% inner_join(continents, by = "standard")
rftraindata <- rftraindata[rftraindata$continent != '',]
rftraindata$continent <- as.factor(rftraindata$continent)

predictors <- colnames(rftraindata)[1:(ncol(rftraindata)-4)]
response <- "continent"

rftrain <- as.h2o(rftraindata)

rfsplit <- h2o.splitFrame(data = rftrain, ratios = 0.8, seed = 1234)
train <- rfsplit[[1]]
valid <- rfsplit[[2]]

# Build and train the model:
lv_model <- h2o.randomForest(x = predictors,
                             y = response,
                             ntrees = 250,
                             max_depth = 40,
                             min_rows = 10,
                             balance_classes = TRUE,
                             training_frame = train,
                             validation_frame = valid)

# Eval performance:
perf <- h2o.performance(lv_model)
perf

################################
# Pull in the REAL validation dataset:
topredict <- as.h2o(sampled_digits)
lockbox <- h2o.deepfeatures(best_model, topredict, layer = 1)
lockbox <- as.data.frame(lockbox)
rownames(lockbox) <- rownames(sampled_digits)

lockbox$srr <- rownames(lockbox)
lockbox <- lockbox %>% inner_join(allcountries, by = "srr")

continents <- read.csv('continents.csv')
lockbox <- lockbox %>% inner_join(continents, by = "standard")
lockbox <- lockbox[lockbox$continent != '',]
lockbox$continent <- as.factor(lockbox$continent)
lockbox <- as.h2o(lockbox)

lvlockbox_perf <- h2o.performance(lv_model, newdata = lockbox)
lvlockbox_perf



##########################
### TEST USING ASV TABLE INSTEAD OF LVs
##########################

# add countries to the actual data
asvtrain <- training
asvtrain$srr <- rownames(training)
asvtrain <- asvtrain %>% inner_join(allcountries, by = "srr")

# add continent annotation
continents <- read.csv('continents.csv')
asvtrain <- asvtrain %>% inner_join(continents, by = "standard")
asvtrain <- asvtrain[asvtrain$continent != '',]
asvtrain$continent <- as.factor(asvtrain$continent)

predictors <- colnames(asvtrain)[1:(ncol(asvtrain)-4)]
response <- "continent"

asvtrain <- as.h2o(asvtrain)

asvsplit <- h2o.splitFrame(data = asvtrain, ratios = 0.8, seed = 1234)
asvtrain <- asvsplit[[1]]
asvvalid <- asvsplit[[2]]

# Build and train the model:
asv_model <- h2o.randomForest(x = predictors,
                              y = response,
                              ntrees = 250,
                              max_depth = 40,
                              min_rows = 10,
                              balance_classes = TRUE,
                              training_frame = asvtrain,
                              validation_frame = asvvalid)

# Eval performance:
asvperf <- h2o.performance(asv_model)
asvperf

################################
# Pull in the REAL validation dataset:
asvlockbox <- sampled_digits
asvlockbox$srr <- rownames(asvlockbox)
asvlockbox <- asvlockbox %>% inner_join(allcountries, by = "srr")

continents <- read.csv('continents.csv')
asvlockbox <- asvlockbox %>% inner_join(continents, by = "standard")
asvlockbox <- asvlockbox[asvlockbox$continent != '',]
asvlockbox$continent <- as.factor(asvlockbox$continent)
asvlockbox <- as.h2o(asvlockbox)

asvlockbox_perf <- h2o.performance(asv_model, newdata = asvlockbox)
asvlockbox_perf


########################
## DIFFERENTIAL ACTIVATION
#########################

asia <- rftraindata[rftraindata$continent=='Asia',]
na <- rftraindata[rftraindata$continent=='NAmerica',]
x <- wilcox.test(asia$DF.L1.C1, na$DF.L1.C1)

heatmap.2(as.matrix(rbind(asia,na)),dendrogram="none",trace="none",
          Rowv=FALSE, RowSideColors=as.character(as.numeric(dat$GO)))

results = data.frame(matrix(ncol = 3, nrow = 0))
colnames(results) <- c('LV', 'W', 'p')
for(col in colnames(heatplot[,lv_variance > 0.1])) {
  test <- wilcox.test(asia[[col]], na[[col]])
  results <- rbind(results, data.frame(LV=col, W=test$statistic, p=test$p.value))
}
results
lvweights.abs <- t(abs(final_weights[32,]))
lvweights <- t(final_weights[100,])
median(na$DF.L1.C100)
median(asia$DF.L1.C100)

combined <- rbind(asia, na)
library(ggExtra)
#plain <- ggplot(combined, aes(x=DF.L1.C32, y=DF.L1.C41, color=continent)) +
plain <- ggplot(combined, aes(x=DF.L1.C32, y=DF.L1.C41, color=continent)) +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  )
ggExtra::ggMarginal(plain,
  type='density',
  groupColour = TRUE, groupFill=TRUE)

