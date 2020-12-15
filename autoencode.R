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
reads.nozero <- cmultRepl(reads, label=0, method='SQ', output="p-counts") # default method won't work

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
# only include the ~400 OTUs with the highest variance
clr.high <- clr[,variance > 0.5]

# set aside test samples
index <- sample(1:nrow(clr.high), 2174)
sampled_digits <- clr.high[index, ]
training <- clr.high[-index,]

features <- as.h2o(training)



# Find the BEST model:

# DEPTH hyperparameter search grid
hyper_grid <- list(hidden = list(
  c(5),
  c(10),
  c(20),
  c(50),
  c(50,10,50),
  c(100, 50, 100)
))
# Execute grid search
ae_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = 'clrautoencoder_grid_new',
  autoencoder = TRUE,
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  ignore_const_cols = FALSE,
  seed = 123
)
# Print grid details
h2o.getGrid('clrautoencoder_grid_new', sort_by = 'mse', decreasing = FALSE)

# Measure how well it reconstructs the data
best_model_id <- ae_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Plot reconstructed pixel values
reconstructed_digits <- predict(best_model, as.h2o(sampled_digits))
reconstructed <- as.data.frame(reconstructed_digits)
colnames(reconstructed) <- colnames(sampled_digits)
reconstructed$sample <- rownames(sampled_digits)
reconstructed_long <- reconstructed %>% pivot_longer(-sample, names_to='taxon', values_to='decoded')
sampled_digits$sample <- rownames(sampled_digits)
sampled_long <- sampled_digits %>% pivot_longer(-sample, names_to='taxon', values_to='true')
newtable <- inner_join(reconstructed_long,sampled_long, by=c('sample','taxon'))
stats::cor(newtable$true, newtable$decoded, method='pearson')




# SPARSITY Hyperparameter search grid
hyper_grid <- list(sparsity_beta = c(0.01, 0.05, 0.1, 0.2, 0.5))
# Execute grid search
ae_sparsity_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = '1clrsparsity_grid',
  autoencoder = TRUE,
  hidden = 25, # include best model from depth grid search HERE
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  average_activation = -0.1,
  ignore_const_cols = FALSE,
  seed = 123
)
h2o.getGrid('1clrsparsity_grid', sort_by = 'mse', decreasing = FALSE)


best_model_id <- ae_sparsity_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)



# test all hyperparams together:
hyper_grid <- list(
  sparsity_beta = c(0.01, 0.05, 0.1, 0.2, 0.5),
  hidden = list(
    c(5),
    c(10),
    c(50),
    c(50,10,50),
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
best_model_id <- big_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

results <- data.frame(
  model=character(),
  clr_pearson=double(),
  rel_pearson=double()
)

library(compositions) # for clrInv()

for(i in seq(5,length(big_grid@model_ids))) {
  model_id <- big_grid@model_ids[[i]]
  model <- h2o.getModel(model_id)
  # make predictions
  reconstructed_digits <- predict(model, as.h2o(sampled_digits))
  # evaluate CLR results
  reconstructed <- as.data.frame(reconstructed_digits)
  colnames(reconstructed) <- colnames(training)
  reconstructed$sample <- rownames(sampled_digits)
  reconstructed_long <- reconstructed %>% pivot_longer(-sample, names_to='taxon', values_to='decoded')
  sampled_digits$sample <- rownames(sampled_digits)
  sampled_long <- sampled_digits %>% pivot_longer(-sample, names_to='taxon', values_to='true')
  clrtable <- inner_join(reconstructed_long,sampled_long, by=c('sample','taxon'))

  # evaluate relative abundance results
  reconstructed.rel <- as.data.frame(reconstructed_digits)
  colnames(reconstructed.rel) <- colnames(training)
  reconstructed.rel <- as.data.frame(clrInv(reconstructed.rel))
  reconstructed.rel$sample <- rownames(sampled_digits)
  reconstructed.rel_long <- reconstructed.rel %>% pivot_longer(-sample, names_to='taxon', values_to='decoded')

  sampled_digits.rel <- sampled_digits
  sampled_digits.rel$sample <- NULL
  sampled_digits.rel <- as.data.frame(clrInv(sampled_digits.rel))

  sampled_digits.rel$sample <- rownames(sampled_digits)
  sampled.rel_long <- sampled_digits.rel %>% pivot_longer(-sample, names_to='taxon', values_to='true')
  reltable <- inner_join(reconstructed.rel_long,sampled.rel_long, by=c('sample','taxon'))

  # append results
  results <- rbind(results, data.frame(
    model = model_id,
    clr_pearson = stats::cor(clrtable$true, clrtable$decoded),
    rel_pearson = stats::cor(reltable$true, reltable$decoded)
  ))
}


##########################################3

best_model <- h2o.getModel('big_grid_model_1') # from results table

# Plot reconstructed pixel values using CLR-transformed numbers
reconstructed_digits <- predict(best_model, as.h2o(sampled_digits))
reconstructed <- as.data.frame(reconstructed_digits)
names <- gsub('_.*','', rownames(sampled_digits))
colnames(reconstructed) <- colnames(training)
reconstructed$study <- names
reconstructed_long <- reconstructed %>% pivot_longer(-study, names_to='taxon', values_to='decoded')
sampled_digits$study <- namesj
sampled_digits$sample <- NULL
sampled_long <- sampled_digits %>% pivot_longer(-study, names_to='taxon', values_to='true')
newtable <- inner_join(reconstructed_long,sampled_long, by=c('study','taxon'))
ggplot(newtable, aes(x=true, y=decoded, color=sample)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, color='red') +
  theme_bw() +
  labs(x='Observed',y='Predicted') +
  theme(
    legend.position='none'
  )

# Plot same stuff using RELATIVE ABUNDANCE
reconstructed.rel <- as.data.frame(reconstructed_digits)
colnames(reconstructed.rel) <- colnames(training)

library(compositions)
reconstructed.rel <- as.data.frame(clrInv(reconstructed.rel))
reconstructed.rel$sample <- rownames(sampled_digits)
reconstructed.rel_long <- reconstructed.rel %>% pivot_longer(-sample, names_to='taxon', values_to='decoded')

sampled_digits.rel <- sampled_digits
sampled_digits.rel$sample <- NULL
sampled_digits.rel <- as.data.frame(clrInv(sampled_digits.rel))

sampled_digits.rel$sample <- rownames(sampled_digits)
sampled.rel_long <- sampled_digits.rel %>% pivot_longer(-sample, names_to='taxon', values_to='true')
newtable.rel <- inner_join(reconstructed.rel_long,sampled.rel_long, by=c('sample','taxon'))
library(scales)
ggplot(newtable.rel, aes(x=true, y=decoded, color=sample)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, color='red') +
  theme_bw() +
  labs(x='Observed',y='Predicted') +
  theme(
    legend.position='none'
  ) +
  #scale_x_log10(labels = percent) +
  #scale_y_log10(labels = percent) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent)

stats::cor(newtable.rel$true, newtable.rel$decoded)



new_codings <- h2o.deepfeatures(best_model, features, layer = 1)
# plot all the samples using the two biggest codes
ggplot(as.data.frame(new_codings), aes(x=DF.L1.C1, y=DF.L1.C2)) +
  geom_point(aes(color=gsub('_.*','', rownames(training)))) +
  theme_bw() +
  theme(
    legend.position = 'none'
  )

reconstruction_errors <- h2o.anomaly(best_model, features)
reconstruction_errors <- as.data.frame(reconstruction_errors)
ggplot(reconstruction_errors, aes(Reconstruction.MSE)) +
  geom_histogram(bins=100)
