#' Samba
#'
#' This is an all-in-one function to preprocesses and analyze read count data from a CRISPR
#' screen. The structure of the input data is assumed to have columns labeled "sgRNA" and
#' "Gene", followed by read counts from each sample in appropriately labeled columns.
#'
#' @param data data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param null.dist.gene logical, indicating the gene ID for control guides. This is used to generate a null distribution for gene-level analyses. Defaults to "NTC", but this can parameter can also be set to NULL to use all guides for the null distribution.
#' @param score.method character, indicating whether to use the 'KS' (Kolmogorov-Smirnov), 'MetaAnalysis', or 'GeneScore' method to perform gene-level analysis.
#' @param test.method character, indicating whether to use 'QLF' or 'LRT' tests.
#' @param GuideMap data.frame that maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.
#' @param file.prefix Prefix for writing output files.
#' @return DGEList object of the edgeR package.
#' @export
Samba <- function(data, design, coefficient = NULL, contrast = NULL, null.dist.gene = 'NTC', normalization.method = 'TMMwsp', uq2.f=100,
                  score.method = 'KS', test.method = 'QLF', GuideMap = NULL, file.prefix = NULL){
    
    dge <- Preprocess_Samba(data = data, design = design, normalization.method = normalization.method, uq2.f=uq2.f)
    sgRes <- Analyze_Samba_Guides(dge = dge, coefficient = coefficient, method = test.method, contrast = contrast, file.prefix = file.prefix)
    geneRes <- Analyze_Samba_Genes(sgRes = sgRes, null.dist.gene = null.dist.gene, score.method = score.method, file.prefix = file.prefix)
    
    return(list(GuideResults = sgRes, GeneResults = geneRes))
}


#' Preprocess_Samba
#'
#' This function preprocesses read count data from CRISPR screen data. The structure of the
#' input data is assumed to have columns labeled "sgRNA" and "Gene", followed by read counts
#' from each sample in appropriately labeled columns.
#'
#' @param data data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @return DGEList object of the edgeR package.
#' @export
Preprocess_Samba <- function(data, design, min.guides = 1, group = NULL, normalization.method = 'TMMwsp', dispersion.estimate.method='GLMTagwise', uq2.f=100){
    cat('Starting preprocessing.\n')
    # Default settings
    hdr.gene = 'Gene'
    hdr.guide = 'sgRNA'
    
    # prep input
    sgInfo <- data.frame(data[,1:2])
    data <- data.frame(data[,-c(1:2)])
    rownames(data) <- sgInfo[, hdr.guide]
    rownames(sgInfo) <- sgInfo[, hdr.guide]
    design <- as.matrix(design)
    
    # check input
    if (sum(colnames(sgInfo) %in% c(hdr.gene, hdr.guide)) != 2) 
        warning("Count data is missing the columns: \"Guide\" and \"Gene\"!")
    if (nrow(design) != ncol(data)) 
        warning("Design matrix has incorrect number of rows!")
    if ((mode(design) != "numeric") & (mode(design) != "integer")) 
        warning("Design matrix does not have integer/numeric format!")
    
    # get sample groups
    if(is.null(group)) group <- as.character(design[,ncol(design)]) 
    
    # filter data
    if(is.integer(min.guides) | is.numeric(min.guides)){
        if (min.guides > 0) {
            rm.guide <- lapply(unique(group), function(g) which(rowSums(data[,which(group == g), drop=F]) < min.guides) %>% names()) %>% unlist() %>% unique()
            keep.guide <- setdiff(rownames(data), rm.guide)
            cat(paste0("   Removed ", length(rm.guide), " from ", nrow(data), " guides.\n"))
        } else {
            keep.guide <- rownames(data)
        }
    }
    data <- data[keep.guide,]
    sgInfo <- sgInfo[keep.guide,]
    
    # normalize data
    cat("Performing normalization: ",normalization.method,"\n")
    data <- switch(normalization.method,
                   uq2=UpperQuantNorm2Step(data.frame(data), uq2.f),
                   total=apply(data, 2, function(x) 1e6*x/sum(x)),
                   zscore=scale(data, center=F),
                   data
    )
    if(normalization.method %in% c('uq2','total','zscore'))
        normalization.method <- 'none'

    #  Make DGE object and filter counts 
    dge <- edgeR::DGEList(data, remove.zeros = T, group = group)
    
    #  Apply median normalization, if applicable
    if(normalization.method == 'median') dge$samples$norm.factors <- EBSeq::MedianNorm(dge$counts)
    #  Apply default TMMwsp or no normalization
    if(normalization.method != 'median') dge <- edgeR::calcNormFactors(dge, method = normalization.method)
    
    # calculate dispersion
    cat('Calculating dispersions ...\n')
    if(dispersion.estimate.method == 'GLMTagwise'){
        dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
        cat('   Calculated common dispersion\n')
        dge <- edgeR::estimateGLMTrendedDisp(dge, design = design, method = 'bin.loess', span = 1/3)
        cat('   Calculated trended dispersion\n')
        dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
        cat('   Calculated sgRNA-wise dispersion\n')
    }

    # return DGE object
    dge$Design <- design
    dge$GuideMap <- sgInfo
    cat('Preprocessing complete!\n')
    return(dge)
}

#' Analyze_Samba_Guides
#'
#' This function performs guide-level screen enrichment analysis of SAMBA, using a pre-processed DGEList object as an input.
#'
#' @param dge DGEList object with normalization and dispersion already calculated. This should be the output from the "Preprocess_Samba" function.
#' @param method character, indicating which statistical model and test to use for guide-level analysis. Currently, only "QLF" and "LRT" are supported.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param GuideMap data.frame that maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.
#' @param file.prefix (optional) file.path and/or character prefix to use for exporting the result tables
#' @return data.frame of the results from the guide-level screen enrichment analysis.
#' @export
Analyze_Samba_Guides <- function(dge, method = 'QLF', design = NULL, coefficient = NULL, contrast = NULL, GuideMap = NULL, file.prefix = NULL){
    # set defaults
    hdr.gene = 'Gene'
    hdr.guide = 'sgRNA'
    fdr.threshold = 0.1
    
    # check inputs
    if(!xor(is.null(coefficient), is.null(contrast))) warning('Please input a valid coefficient or contrast!')
    if(!(method %in% c('QLF','LRT'))) warning('Please choose either "LRT" or "QLF" as an analysis method!')
    if(is.null(GuideMap)) GuideMap <- dge$GuideMap
    if(is.null(GuideMap)) warning('Please provide a dataframe that maps guides to genes!')
    
    # get design info
    if(is.null(design)) design <- dge$Design
    
    # sgRNA-level analysis
    cat('Performing sgRNA-level analysis ...\n')
    if(method == 'QLF'){
        fit <- edgeR::glmQLFit(dge, design = design, robust = T)
        sgRes <- edgeR::glmQLFTest(fit, coef = coefficient, contrast = contrast)
    }
    if(method == 'LRT'){
        fit <- edgeR::glmFit(dge, design = design, robust = T)
        sgRes <- edgeR::glmLRT(fit, coef = coefficient, contrast = contrast)
    }
    cat('Done!\n')
    
    # Format sg-level results
    sgRes <- edgeR::topTags(sgRes, n = Inf, sort.by = 'PValue')$table
    colnames(sgRes)[3] <- 'Statistic'
    sgRes$sgRNA <- rownames(sgRes)
    sgRes <- merge(GuideMap, sgRes, by = 'sgRNA', all.x = F, all.y = T)
    sgRes$zscore <- scale(sgRes$logFC, center = F, scale = T)[,1]
    sgRes <- sgRes[with(sgRes, order(Statistic,decreasing = T)),]
    
    # save data
    if(!is.null(file.prefix)){
        saveRDS(fit, file = paste0(file.prefix, '_FittedData.rds'))
        write.table(sgRes, file = paste0(file.prefix, '_GuideLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }
    
    # return data
    return(sgRes)
}

#' Analyze_Samba_Genes
#'
#' This function performs the gene-level screen enrichment analysis of SAMBA, using guide-level SAMBA results as an input.
#'
#' @param sgRes data.frame of the results from the guide-level screen enrichment analysis. This is the output from the "Analyze_Samba_Guides" function.
#' @param null.dist.gene logical, indicating the gene ID for control guides. This is used to generate a null distribution for gene-level analyses. Defaults to "NTC", but this can parameter can also be set to NULL to use all guides for the null distribution.
#' @param score.method character, indicating whether to use the 'MetaAnalysis' or 'GeneScore' method to perform gene-level analysis.
#' @param file.prefix (optional) file.path and/or character prefix to use for exporting the result tables
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
Analyze_Samba_Genes <- function(sgRes, null.dist.gene = 'NTC', score.method = 'KS', min.guide.per.gene = 2, file.prefix = NULL){
    # Set defaults
    control.gene <- null.dist.gene
    hdr.gene <- 'Gene'
    fdr.threshold <- 0.1
    
    # Check inputs
    if(!(score.method %in% c('GeneScore','KS'))) warning('Please choose either "GeneScore" or "MetaAnalysis" as a score.method!')
    
    ## Filter guide-level data (>= 2sgRNA/gene and not NA)
    sgPerGene <- plyr::ddply(sgRes, 'Gene', nrow)
    keep.genes <- sgPerGene[which(sgPerGene[,2] >= min.guide.per.gene), 'Gene']
    sgRes <- sgRes[which(sgRes$Gene %in% keep.genes),]
    sgRes <- sgRes[!is.na(sgRes$Gene ),]
    
    ## set data for null dist
    cat('Generating null distribution.\n')
    if(sum(sgRes$Gene == control.gene) == 0) {
        nullData <- sgRes
        nGeneGuides = sgPerGene[,'V1']
    }
    if(sum(sgRes$Gene == control.gene) > 0) {
        nullData <- sgRes[which(sgRes$Gene  == control.gene),]
        sgRes <- sgRes[which(sgRes$Gene != control.gene),]
        nGeneGuides = sgPerGene[which(sgPerGene$Gene != control.gene),'V1'] 
    }
    
    rIndices <- RandomIndexGenerator(sgPerGene = nGeneGuides, nTotalGuides = nrow(nullData))
    nullData <- nullData[unlist(rIndices),]
    i <- lapply(rIndices, length)
    nullData[,hdr.gene] <- paste0('NTC', lapply(1:length(i), function(x){ rep(x,i[[x]]) }) %>% unlist())
    
    # perform enrichment and depletion analyses
    if(score.method == 'GeneScore'){
        cat('Calculating enrichment.\n')
        enriched <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1)
        cat('Calculating depletion.\n')
        depleted <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = -1)
        cat('Done!\n')
        output <- merge(enriched,depleted,by = 'Gene', all = T)
    }
    if(score.method == 'KS'){
        cat('Calculating KS-stats\n')
        enriched <- KStest_Samba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1, filter_beta = T)
        depleted <- KStest_Samba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = -1, filter_beta = T)
        output <- merge(enriched,depleted,by = 'Gene', all = T)
        output <- output[order(output$pval_pos, -output$score_pos, decreasing = F),]
        cat('Done!\n')
    }
    output <- output[order(output$pval_pos, decreasing = F),]
    
    # save data
    if(!is.null(file.prefix)){
        write.table(output, file = paste0(file.prefix, '_GeneLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }
    
    return(output)
}


#' GeneScoreSamba
#'
#' This function performs SAMBA gene-level aggregation, using the gene score method. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of guide-level screen enrichment data for the null distribution.
#' @param fdr.threshold numeric between 0 and 1, indicating the FDR to use as a threshold.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
GeneScoreSamba <- function(data, nullData, fdr.threshold, direction = 1){
    data$logFC <- direction * data$logFC
    nullData$logFC <- direction * nullData$logFC
    sg.threshold <- quantile(nullData$logFC, 1-fdr.threshold)
    if(direction == 1) cat(paste0('FDR threshold for guides: ', sg.threshold,'\n'))
    
    # Get scores
    nullScores <- WtSumScore(nullData, sg.threshold, direction)
    targetScores <- WtSumScore(data, sg.threshold, direction)
    
    # Calculate z-score from null distribution
    targetScores$score <- scale(targetScores$score, 
                                median(nullScores$score), 
                                sd(nullScores$score))[,1]
    
    # z-test
    if(direction==  1) pval <- pnorm(q=targetScores$score, lower.tail=F)
    if(direction== -1) pval <- pnorm(q=targetScores$score, lower.tail=T)
    qval <- p.adjust(pval, 'fdr')
    
    # Generate output
    output <- data.frame(targetScores, pval = pval, qval = qval)
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output)[-1] <- paste0(colnames(output)[-1],'_',direction)
    return(output)
}


#' KStest_Samba
#'
#' This function performs SAMBA gene-level aggregation, using the Kolmogorov-Smirnov test. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of guide-level screen enrichment data for the null distribution.
#' @param fdr.threshold numeric between 0 and 1, indicating the FDR to use as a threshold.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @param filter_beta logical, describing whether to use only half of the guides for each gene.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
KStest_Samba <- function(data, nullData, fdr.threshold, direction = 1, filter_beta = F){
    data$logFC <- direction * data$logFC
    data <- data[order(data$logFC, decreasing = T),]
    nullData$logFC <- direction * nullData$logFC

    # Filter coefficients to include only the top half
    if(filter_beta) {
        data <- dplyr::reframe(data, .by= 'Gene', n_guides = max(length(Gene),3), tophalf = quantile(logFC,seq(.5,1,.125))) 
        data$logFC <- data$tophalf
        nullData <- dplyr::reframe(data, .by= 'Gene', n_guides = max(length(Gene),3), tophalf = quantile(logFC,seq(.5,1,.125))) 
        nullData$logFC <- nullData$tophalf
    }
    
    # Calculate z-score
    data$scale <- scale(data$logFC, median(nullData$logFC), sd(nullData$logFC))[,1]
    
    # Calculate p & q values
    ifelse(direction==1, alt <- 'less', alt <- 'greater')
    runKS <- function(x) return(ks.test(x = x, 'pnorm', alternative = alt, simulate.p.value = F, B = 10000)$p.value)
    output <- dplyr::reframe(data, score = quantile(scale, 0.75), .by = 'Gene')
    output$pval <- dplyr::reframe(data, pval = runKS(scale), .by = 'Gene')[,2]
    output$qval <- p.adjust(output$pval, 'fdr')
    
    # Generate output
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output)[-1] <- paste0(colnames(output)[-1],'_',direction)
    return(output)
}


#' RandomIndexGenerator
#'
#' This function creates a list of reproducible, random integers from which to create a null distribution.
#'
#' @param sgPerGene integer vector, representing the number of guides for each gene.
#' @param nTotalGuides integer, indicating the total number of guides in the filtered library.
#' @return list of integer vectors, each of random numbers from which to construct a null data distribution.
#' @export
RandomIndexGenerator <- function(sgPerGene, nTotalGuides){
    # make sure there is between 10,000 and 100,000 iterations
    nIterations <- min(max(length(sgPerGene), 10000), 100000) 
    # make sure there is a min of 2 grna/gene
    sgPerGene[which(sgPerGene < 2)] <- 2
    # if n-genes < n-iterations, bootstrap up to n-iterations
    set.seed(42)
    if(length(sgPerGene) < nIterations) sgPerGene <- dqrng::dqsample(sgPerGene, size = nIterations, replace = T)
    
    set.seed(42)
    randomSeeds <- dqrng::dqsample(1:nIterations, replace = F)
    randData <- lapply(1:nIterations, function(i) {
        set.seed(randomSeeds[i])
        dqrng::dqsample(1:nTotalGuides, size = min(nTotalGuides,sgPerGene[i]), replace = F)
    })
    return(randData)
}


#' WtSumScore
#'
#' This function calculates a single gene score for a given set of guide log-FC values.
#'
#' @param scores numeric vector of the log-FC values for all guides of a single gene.
#' @param sg.threshold numeric, giving the logFC value of the null data that represents the FDR threshold
#' @return numeric score of a gene.
#' @export
WtSumScore <- function(data, fdr.thresh, direction){
    wt.increment <- 1/4
    # sort grna-level results by logfc
    data <- data[with(data, order(Gene, -logFC)),]
    # get the mean of the top half of grna log-FCs, using quantiles to account for difference in the #s of gRNAs/Gene
    score_tophalf <- dplyr::reframe(data, .by= 'Gene', n.grna = max(length(Gene),3), 
                                    score_tophalf = sum(quantile(logFC,seq(0.5,1,.25)))/n.grna) 

    # calculate weights for filtered grna
    if(nrow(data[which(data$logFC > fdr.thresh),]) > 0){
        score_fdr <- data[which(data$logFC > fdr.thresh),] %>% 
            dplyr::group_by(Gene) %>% 
            dplyr::mutate(wt = 1 + (wt.increment * (1:length(logFC)))) 
        score_fdr[which(score_fdr$wt > 2),'wt'] <- 2 # cap wts
    } else {
        score_fdr <- data %>% 
            dplyr::group_by(Gene) %>% 
            dplyr::mutate(wt = 0) 
    }

    # summarize weighted logFCs
    score_fdr <- score_fdr %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(n_fdr_guides = length(logFC), score_fdr = sum(logFC * wt)/length(logFC))
    # merge the score_tophalf, score_fdr, and count_fdr
    scores <- merge(score_fdr, score_tophalf, by = 'Gene', all = T)
    scores[is.na(scores)] <- 0
    
    scores <- data.frame(Gene = scores$Gene, 
                         score = (scores$score_fdr + scores$score_tophalf), 
                         n_fdr_guides = scores$n_fdr_guides)
    
    return(scores)
}


#' UpperQuantNorm
#'
#' This function normalizes a count matrix by the 75th percentile of each sample.
#'
#' @param X Count matrix with genes as rows and samples as columns.
#' @return Normalized matrix.
#' @export
UpperQuantNorm <- function(X){
    X<-X+0.1
    upperQ<-apply(X,2,function(y) quantile(y, 0.75))
    f<-upperQ/mean(upperQ) # calculate normalization factor
    res<-scale(X,center=FALSE,scale=f) 
    return(res)
}

#' UpperQuantNorm2Step
#'
#' This function performs a 2-step normalization, first with quantile normalization
#' of samples, then by median normalization of gRNA values.
#'
#' @param X Count matrix with genes as rows and samples as columns.
#' @param f Multiplication factor for nornalized data.
#' @return Normalized matrix.
#' @export
UpperQuantNorm2Step <-function(X, f = 50){
    # per gene normalization by Median: pgQ2
    # X: a matrix of data with the multiplication of factor (f) as:
    # f=100 as a default
    
    uq.res <- UpperQuantNorm(X) #  perform UQ normalization per sample
    
    ## perform per gene nomralization
    m<-apply(uq.res,1, median)
    
    idx<-which(m==0) # avoid 0 median 
    m[idx]=0.1
    
    si<-m/f # calculate normalization factor
    X1<-scale(t(uq.res),center=FALSE,scale=si)
    res<-t(X1)
    return(res)  
}




