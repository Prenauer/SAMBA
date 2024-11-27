## SAMBA
CRISPR Screen analysis with moderated Bayesian statistics and aggregated gene scoring (SAMBA). 
This analysis algorithm employs negative binomial generalized linear models (GLMs) to fit CRISPR 
screen data, and it analyzes guide-level enrichment/depletion using moderated Bayesian statistics 
and quasi-likelihood F tests. Gene-level scores are then calculated as a modified weighted sum 
and compared with randomized null distribution to calculate enrichment and depletion statistics 
for the screen.


## Install
```{r}
## Install packages
install.packages('devtools')  
devtools::install_github('Prenauer/SAMBA')

## Load SAMBA
library(SAMBA)
```

## Setting up the data
#### All you need to run SAMBA is a (1) dataframe of count data and (2) sample names. 
The example below demonstrates a simple screen analysis. However, more complex analyses
are possible by supplying a design matrix, and even a contrast matrix, to take full 
advantage of the GLM capabilities of SAMBA. 

##### (1) Get dataframe of sgRNA counts
This needs the following columns: sgRNA, Gene, counts. See the example "counts" 
dataset with four control and four screen samples below.

```{r}
## Create 4 screen and 4 control samples, each with random counts of 80,000 sgRNA.
samples.screen <- sapply(1:4, function(x) rnbinom(80000, mu = 1000, size = 0.2))  
samples.ctrl <- sapply(1:4, function(x) rnbinom(80000, mu = 1000, size = 1))  

## Generate names for the sgRNAs and genes included in the screen library.
library.sgrna <- sapply(1:80000, function(x) paste0('sgRNA_',x))  
library.gene <- c(unlist(lapply(1:(79000/4), function(x) rep(paste0('Gene_',x),4))), 
                  rep('NTC',1000))

## Create a dataframe with the CRISPR library information, followed by the sgRNA counts
##    for each sample.
counts <- data.frame(sgRNA = library.sgrna,
                     Gene = library.gene,
                     samples.ctrl,
                     samples.screen)
## Name the samples (optional).
colnames(counts)[3:10] <- c(paste0('Ctrl_',1:4), paste0('Screen_',1:4))
head(counts)
```

##### (2) Get sample names
Create character vectors of the screen samples and control samples.
```{r}
screen.names <- paste0('Screen_',1:4)
ctrl.names <- paste0('Ctrl_',1:4)
```

## Run the analysis
The simplest way to run SAMBA is to use the all-in-one "Samba" function, which is 
demonstrated below.
```{r}
## Run SAMBA using the count data, the sample 
##    Note that Samba works better when a non-targeting control is used. In the
##    example dataset, these are named "NTC". Also note that the output is a list 
##    of sgRNA-level results and gene-level results.
results <- Samba(data = counts,
    screen.names=screen.names,
    ctrl.names=ctrl.names,
    control.gene='NTC',
    coefficient = 'screen')

## View Gene-level results.
head(results$GeneResults)
```

