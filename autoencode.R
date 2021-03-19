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


# make everything relative abundances
rel <- apply(reads, 1, function(a) a/sum(a))
rel <- as.data.frame(t(rel))

library(dplyr)    # for data manipulation
library(tidyr)

# trim the names to be only the genus
colnames(rel) <- gsub('.+\\.+(\\w+)$', '\\1', colnames(rel))

pca <- prcomp(rel)
rel$study <- gsub('_.*','', rownames(rel))
#rel_long <- rel %>% pivot_longer(-study, names_to='taxon', values_to='true')
#length(rel_long[is.infinite(rel_long$true),]$true)

autoplot(pca, data = rel, colour = 'study',
         loadings=TRUE, loadings.label = TRUE,
         loadings.label.colour='black') +
  theme_bw() +
  theme(
    legend.position = 'none'
  )

# CLR TRANSFORM VERSION
# fix zeroes
library(zCompositions)

load(file="reads.nozero.Rda")
#reads.nozero <- cmultRepl(reads, label=0, method='SQ', output="p-counts") # default method won't work
#save(reads.nozero, file="reads.nozero.Rda")


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
clr <- apply(reads.nozero, 1, function(a) log(a/gm_mean(a)))
clr <- as.data.frame(t(clr))

clr.trimmedname <- clr
colnames(clr.trimmedname) <- gsub('.+\\.+(\\w+)$', '\\1', colnames(clr.trimmedname))

pca_clr <- prcomp(clr.trimmedname)
# add study IDs to data matrix
clr.trimmedname$study = gsub('_.*','', rownames(clr.trimmedname))
#clr_long <- clr %>% pivot_longer(-study, names_to='taxon', values_to='true')
autoplot(pca_clr, data = clr.trimmedname, colour = 'study',
     loadings=TRUE, loadings.label = TRUE,
     loadings.label.colour='black') +
  theme_bw() +
  theme(
    legend.position = 'none',
  )

# https://bradleyboehmke.github.io/HOML/autoencoders.html#undercomplete-autoencoders
library(ggplot2)  # for data visualization
library(h2o)  # for fitting autoencoders


h2o.init(max_mem_size = "8g")
clr$study <- NULL

variance <- sapply(clr, sd)
hist(variance)
sum(variance > 0.35)
# only include the ~400 OTUs with the highest variance
clr.high <- clr[,variance > 0.35]

# set aside test samples
index <- sample(1:nrow(clr.high), 2174)
sampled_digits <- clr.high[index, ]
training <- clr.high[-index,]
# we have to make sure there aren't duplicates in the training
# set so T-SNE will do its thing
training <- unique(training)

features <- as.h2o(training)

# test all hyperparams together:
hyper_grid <- list(
  sparsity_beta = c(0.01, 0.1, 0.2, 0.5),
  hidden = list(
    #c(5),
    #c(25,10,25)
    c(25),
    c(50),
    c(80),
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
h2o.getGrid('big_grid', sort_by = 'mse', decreasing = FALSE)
best_model_id <- big_grid@model_ids[[1]] # !!!! CHANGE BACK TO 1
best_model <- h2o.getModel(best_model_id)

new_codings <- h2o.deepfeatures(best_model, features, layer = 1)
new_codings <- as.data.frame(new_codings)


# Plot the LVs using T-SNE

# annotate with countries!!
countries <- read.csv('sample_countries.csv')

# throw out the ones that aren't US/China
newcountries <- countries
newcountries[!newcountries$standard %in% c('BOL','KOR','ECU'),]$standard <- NA
#newcountries[!newcountries$standard %in% c('JPN','CAN'),]$standard <- NA

#samples <- as.data.frame(gsub('^([^_]+)_.+', '\\1', rownames(training)))
samples <- as.data.frame(gsub('^[^_]+_','', rownames(training)))
colnames(samples) <- c('srr')

annotated <- samples %>% left_join(newcountries, by = "srr")

library(Rtsne)
set.seed(1234)
tsned <- Rtsne(as.matrix(new_codings), perplexity=75,
    pca=FALSE, max_iter=3000, num_threads=6,
    theta=0.5)

toplot <- as.data.frame(tsned$Y)
ggplot(toplot, aes(x=V1, y=V2, color=annotated$standard)) +
  geom_point() +
  theme_bw()

# leaving out the ones without a color assigned:
toplot <- as.data.frame(tsned$Y)
toplot <- subset(toplot, !is.na(annotated$standard))
ggplot(toplot, aes(x=V1, y=V2, color=annotated[!is.na(annotated$standard),]$standard)) +
  geom_point() +
  theme_bw()

# plot all the samples using the two biggest codes
codings <- as.data.frame(new_codings)
ggplot(codings, aes(x=DF.L1.C13, y=DF.L1.C12, color=annotated$standard)) +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = 'none'
  )

# all samples, but leave out the ones without colors:
toplot <- subset(codings, !is.na(annotated$standard))
ggplot(toplot, aes(x=DF.L1.C1, y=DF.L1.C4, color=annotated[!is.na(annotated$standard),]$standard)) +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = 'none'
  )