# ----------- RAREFACTION

REPLICATES <- 30
todo <- c(10, 50, seq(from=100, to=2000, by=100), seq(from=3000,to=10000,by=1000))

rarefaction_sampling <- function(taxtable) {
  richness <- taxtable
  richness <- richness %>% mutate_if(is.numeric, ~1 * (. > 0))
  found <- data.frame(colSums(richness)) %>% mutate_if(is.numeric, ~1 * (. > 0))
  found$taxon <- rownames(found)
  colnames(found) <- c('observed','taxon')
  observed <- sum(found$observed)
  rarefaction<-data.frame(nrow(richness), observed)
  names(rarefaction)<-c("scount","observed")
  
  for(scount in todo) {
    print(paste('SAMPLE SIZE:', scount))
    for(iter in seq(0,REPLICATES)) {
      sub <- dplyr::sample_n(richness, scount)
      found <- data.frame(colSums(sub)) %>% mutate_if(is.numeric, ~1 * (. > 0))
      found$taxon <- rownames(found)
      colnames(found) <- c('observed','taxon')
      observed <- sum(found$observed)
      result <- data.frame(scount, observed)
      names(result)<-c("scount","observed")
      rarefaction <- rbind(rarefaction, result)
    }
  }

  print(paste('SAMPLE SIZE: ALL'))
  sub <- richness
  found <- data.frame(colSums(sub)) %>% mutate_if(is.numeric, ~1 * (. > 0))
  found$taxon <- rownames(found)
  colnames(found) <- c('observed','taxon')
  observed <- sum(found$observed)
  result <- data.frame(nrow(richness), observed)
  names(result)<-c("scount","observed")
  rarefaction <- rbind(rarefaction, result)

  return(rarefaction)
}

set.seed(42)
rarefaction.phylum <- rarefaction_sampling(taxphylum)
rarefaction.phylum$level <- 'phylum'
rarefaction.class <- rarefaction_sampling(taxclass)
rarefaction.class$level <- 'class'
rarefaction.order <- rarefaction_sampling(taxorder)
rarefaction.order$level <- 'order'
rarefaction.family <- rarefaction_sampling(taxfamily)
rarefaction.family$level <- 'family'
rarefaction.genus <- rarefaction_sampling(taxgenus)
rarefaction.genus$level <- 'genus'

rarefaction <- rbind(rarefaction.phylum, rarefaction.class,
                     rarefaction.order, rarefaction.family,
                     rarefaction.genus)
rarefaction$level <- as.factor(rarefaction$level)
#rarefaction$scount <- as.factor(rarefaction$scount)
#ggplot(rarefaction, aes(x=scount, y=observed,
#                        fill=level)) +
#  geom_boxplot(outlier.shape = NA, width=1) +
#  #scale_y_continuous(limits=c(0, 1400)) +
#  theme_bw() +
#  labs(x='Samples', y='Taxa observed')

summary <- rarefaction %>%
  group_by(level, scount) %>%
  summarise(
    sd = sd(observed, na.rm = TRUE),
    observed = mean(observed)
  )
summary

ggplot(summary, aes(x=scount, y=observed,
                    ymin = observed-sd, ymax = observed+sd,
                    color=level)) +
  geom_errorbar(width = 50) +
  geom_line() +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(x='Samples', y='Unique taxa') +
  scale_x_continuous(
    limits=c(0,2000),
    labels = comma_format()
    #trans="sqrt"
  )
