# SAMBA
CRISPR Screen analysis with moderated Bayesian statistics and adaptive gene aggregation scoring

# Install
install.packages(devtools)

devtools::install_github('Prenauer/SAMBA')

# Create dummy screen data
samples.screen <- sapply(1:4, function(x){ rnbinom(80000, mu = 1000, size = 0.2) })

samples.ctrl <- sapply(1:4, function(x){ rnbinom(80000, mu = 1000, size = 1) })

library.sgrna <- sapply(1:80000, function(x) { paste0('sgRNA_',x) })

library.gene <- c(lapply(1:(79000/4), function(x) { rep(paste0('Gene_',x),4) }) %>% unlist(),
                  rep('NTC',1000))
                  
counts <- data.frame(sgRNA = library.sgrna,
                     Gene = library.gene,
                     samples.ctrl,
                     samples.screen)
                     
colnames(counts)[3:10] <- c(paste0('Ctrl_',1:4), paste0('Screen_',1:4))

# Create design matrix
Screen <- c(0,0,0,0,1,1,1,1)

design <- model.matrix(~ Screen)


# Run analysis
results <- Samba(data = counts, design = design, coefficient = 'Screen', score.method = 'GeneScore')

View(results$GeneResults)


