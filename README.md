## SAMBA
CRISPR Screen analysis with moderated Bayesian statistics and aggregated gene scoring (SAMBA). 
Note that this is a working version, but a more finalized version is coming soon.


## Install
```{r}
## Install packages
install.packages('devtools')  
devtools::install_github('Prenauer/SAMBA')

## Load SAMBA
library(SAMBA)
```

## Setting up the data
#### All you need to run SAMBA is a (1) dataframe of count data and (2) a design matrix. Each are described below.
##### (1) Dataframe of sgRNA counts
This needs the following columns: sgRNA, Gene, counts. See the example "counts" dataset with four control and four screen samples below.
```{r}
## Create 4 screen and 4 control samples, each with random counts of 80,000 sgRNA.
samples.screen <- sapply(1:4, function(x) rnbinom(80000, mu = 1000, size = 0.2))  
samples.ctrl <- sapply(1:4, function(x) rnbinom(80000, mu = 1000, size = 1))  

## Generate names for the sgRNAs and genes included in the screen library.
library.sgrna <- sapply(1:80000, function(x) paste0('sgRNA_',x))  
library.gene <- c(lapply(1:(79000/4), function(x) rep(paste0('Gene_',x),4)) %>% unlist(), 
                  rep('NTC',1000))

## Create a dataframe with the CRISPR library information, followed by the sgRNA counts for each sample.
counts <- data.frame(sgRNA = library.sgrna,
                     Gene = library.gene,
                     samples.ctrl,
                     samples.screen)
## Name the samples (optional).
colnames(counts)[3:10] <- c(paste0('Ctrl_',1:4), paste0('Screen_',1:4))
head(counts)
```

##### (2) Design matrix
This matrix simply states which samples are the control or screen samples. You can borrow the code below to designate your samples with a 0 = control or 1 = screen. Note that the order of the samples in the design matrix must match the order of the samples in the "counts" dataset. You can also use more complicated design matrices, similar to edgeR analyses.
```{r}
screen <- c(0,0,0,0,1,1,1,1)  
design <- model.matrix(~ screen)
```

## Run the analysis
The simplest way to run SAMBA is to use the all-in-one "Samba" function, which is demonstrated below.
```{r}
## Run SAMBA using the count data, design matrix, and the screen samples as the coefficient.
##    Note that the coefficient is a character vector that indicates a column name in the design matrix.
##    Also note that the output is a list of sgRNA-level results and gene-level results.
results <- Samba(data = counts, design = design, control.gene='NTC', coefficient = 'screen')

## View Gene-level results.
View(results$GeneResults)
```

