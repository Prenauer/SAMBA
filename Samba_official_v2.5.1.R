#' @import {plyr}
#' @import {dplyr}
#' @import {stringr}
#' @import {edgeR}
#' @import {limma}
#' @import {dqrng}
#' Samba
#'
#' This is an all-in-one function to preprocesses and analyze read count data from a CRISPR
#' screen. The structure of the input data is assumed to have columns labeled "sgRNA" and
#' "Gene", followed by read counts from each sample in appropriately labeled columns.
#'
#' @param data Path to the input file.
#' @param design Design matrix for fitting objects.
#' @param coefficient Coefficient to use for comparing data models.
#' @param contrast Contrast matrix to use for comparing data models.
#' @param ntc.as.null.dist Use NTCs or all guides to generate a null distribution for gene-level analyses.
#' @param score.method Use 'MetaAnalysis' or 'GeneScore' method to perform gene-level analysis.
#' @param test.method Use 'QLF' or 'LRT' tests.
#' @param GuideMap (optional) data.frame with sgRNA and gene info.
#' @param file.prefix Prefix for writing output files.
#' @return DGEList object of the edgeR package.
#' @export
Samba <- function(data, design, coefficient = NULL, contrast = NULL, control.gene = NULL,
                  score.method = 'GeneScore', test.method = 'QLF', GuideMap = NULL, file.prefix = NULL){
    
    dge <- Preprocess_Samba(data = data, design = design)
    sgRes <- Analyze_Samba_Guides(dge = dge, coefficient = coefficient, method = test.method, contrast = contrast, file.prefix = file.prefix)
    geneRes <- Analyze_Samba_Genes(sgRes = sgRes, score.method = score.method, control.gene=control.gene, file.prefix = file.prefix)
    
    return(list(GuideResults = sgRes, GeneResults = geneRes))
}


#' Preprocess_Samba
#'
#' This function preprocesses read count data from CRISPR screen data. The structure of the
#' input data is assumed to have columns labeled "sgRNA" and "Gene", followed by read counts
#' from each sample in appropriately labeled columns.
#'
#' @param data Path to the input file.
#' @param design Design matrix for fitting objects.
#' @return DGEList object of the edgeR package.
#' @export
Preprocess_Samba <- function(data, design){
    cat('Starting preprocessing.\n')
    # Default settings
    hdr.gene = 'Gene'
    hdr.guide = 'sgRNA'
    min.guides = 0
    pseudocount = 4
    dispersion.method = 'GLM'
    
    # check input
    data <- data.frame(data)
    design <- as.matrix(design)
    if(sum(colnames(data) %in% c(hdr.gene, hdr.guide)) != 2) warning('Count data is missing the columns: "Guide" and "Gene"!')
    if((nrow(design) + 2) != ncol(data)) warning('Design matrix has incorrect number of rows!')
    if((mode(design) != 'numeric') & (mode(design) != 'integer')) warning('Design matrix does not have integer/numeric format!')
    
    # get sample groups
    if(ncol(design) == 1){
        design <- cbind(baseline=1, design)
    } 
    if(ncol(design) == 2){
        group <- as.integer(design[,2])
    }else {
        group <- as.integer(rowSums(design[,-1]))
    }
    if(length(unique(group)) >= (nrow(design))) {
        group = NULL
        cat('No sample grouping will be used!\n')
    }
    # make/filter DGE object
    rownames(data) <- data[,hdr.guide]
    if(min.guides > 0){
        dge <- edgeR::DGEList(data[,-c(1:2)], remove.zeros = F, group = group)
        keep.exprs <- FilterCountData(dge = dge, group = group, min.guides = min.guides)
        dge <- edgeR::DGEList(data[keep.exprs,-c(1:2)] + pseudocount, remove.zeros = T, group = group)
    } else {
        dge <- edgeR::DGEList(data[,-c(1:2)] + pseudocount, remove.zeros = T, group = group)
    }
    # generate guide weights based on detection across samples
    dge <- WeighGuides(dge = dge)
    
    # calculate size factors
    cat('Calculating TMM-wsp size factors\n')
    dge <- edgeR::calcNormFactors(dge, method = 'TMMwsp')
    cat('   Size factor summary:\n')
    print(summary(dge$samples$norm.factors))
    
    # calculate size factors and dispersion, if more than 2 samples
    if(nrow(dge$samples) > 2){
        cat('Calculating dispersions ...\n')
        dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
        cat('   Calculated common dispersion\n')
        dge <- edgeR::estimateGLMTrendedDisp(dge, design = design, 
                                             method = 'bin.loess', span = 1/3)
        cat('   Calculated trended dispersion\n')
        dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
        cat('   Calculated sgRNA-wise dispersion\n')
    }
    
    # return DGE object
    dge$Design <- design
    dge$GuideMap <- data[,c(hdr.guide, hdr.gene)]
    cat('Preprocessing complete!\n')
    return(dge)
}

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
    
    # sgRNA-level analysis, if more than 2 samples
    if(nrow(dge$samples) > 2){
        cat('Performing sgRNA-level analysis ...\n')
        if(method == 'QLF'){
            fit <- edgeR::glmQLFit(dge, design = design, robust = T)
            sgRes <- edgeR::glmQLFTest(fit, coef = coefficient, 
                                       contrast = contrast)
        }
        if(method == 'LRT'){
            fit <- edgeR::glmFit(dge, design = design, robust = T)
            sgRes <- edgeR::glmLRT(fit, coef = coefficient, contrast = contrast)
        }
        cat('Done!\n')
        sgRes <- edgeR::topTags(sgRes, n = Inf, sort.by = 'PValue')$table
    }
    # create result table for dge with 2samples
    # it is assumed that the screen data has coefficient==1, and 
    #    control coefficient==0
    if(nrow(dge$samples) == 2){
        samp.screen <- which(dge$Design[,coefficient]==1)
        samp.ctrl <- which(dge$Design[,coefficient]==0)
        normdata <- edgeR::cpm(dge, log=T)
        rownames(normdata) <- rownames(dge$counts)
        sgRes <- data.frame(logFC=normdata[,samp.screen]-normdata[,samp.ctrl],
                            logCPM=aveLogCPM(dge), Statistic=NA, PValue=NA,
                            FDR=NA)
    }
    # Format sg-level results
    colnames(sgRes)[3] <- 'Statistic'
    sgRes$sgRNA <- rownames(sgRes)
    sgRes <- merge(GuideMap, sgRes, by = 'sgRNA', all.x = F, all.y = T)
    sgRes$zscore <- scale(sgRes$logFC, center = F, scale = T)[,1]
    sgRes <- sgRes[with(sgRes, order(Statistic, logFC, decreasing = T)),]
    
    # save data
    if(!is.null(file.prefix)){
        saveRDS(fit, file = paste0(file.prefix, '_FittedData.rds'))
        write.table(sgRes, file = paste0(file.prefix, '_GuideLevelResults.txt'), 
                    sep = '\t', quote = F, row.names = F)
    }
    
    # return data
    return(sgRes)
}

#' @export
Analyze_Samba_Genes <- function(sgRes, control.gene = NULL, 
                                score.method = 'GeneScore', file.prefix = NULL){
    # Get FDR threshold for genes based on control guides
    if(!is.null(control.gene)){
        if(sum(sgRes$Gene  == control.gene) >= 50){
            ntc <- sgRes$logFC[which(sgRes$Gene  == control.gene)]
            sg.threshold.up <- quantile(ntc, 0.9)
            sg.threshold.dn <- quantile(ntc, 0.1)
        }
        if(sum(sgRes$Gene  == control.gene) < 50){
            control.gene <- NULL
        }
        
    }
    if(is.null(control.gene)){
        control.gene <- ''
        sg.threshold.up <- quantile(sgRes$logFC, 0.9)
        sg.threshold.dn <- quantile(sgRes$logFC, 0.1)
    }
    

    # Check inputs
    if(!(score.method %in% c('GeneScore','KS','Riger'))) 
        warning('Please choose either "GeneScore", "KS", or "Riger" or as a score.method!')
    
    ## Filter guide-level data (>= 2sgRNA/gene and not NA)
    sgPerGene <- plyr::ddply(sgRes, 'Gene', nrow)
    sgRes <- sgRes[!is.na(sgRes$Gene ),]
    sgRes <- sgRes[which(sgRes$Gene != control.gene),]
    
    
    ## set data for null dist
    cat('Generating null distribution.\n')
    nullData <- sgRes
    nGeneGuides = sgPerGene[which(sgPerGene$Gene != control.gene),'V1']
    rIndices <- RandomIndexGenerator(sgPerGene = nGeneGuides, nTotalGuides = nrow(sgRes))
    nullData <- sgRes[unlist(rIndices),]
    i <- lapply(rIndices, length)
    nullData[,'Gene'] <- paste0('Null', lapply(1:100000, function(x){ 
        rep(x,i[[x]]) }) %>% unlist())
    
    # perform enrichment and depletion analyses
    if(score.method == 'GeneScore'){
        cat('Calculating enrichment.\n')
        enriched <- GeneScoreSamba(data = sgRes, nullData = nullData, 
                                   sg.threshold = sg.threshold.up, 
                                   direction = 1)
        cat('Calculating depletion.\n')
        depleted <- GeneScoreSamba(data = sgRes, nullData = nullData, 
                                   sg.threshold = sg.threshold.dn, 
                                   direction = -1)
        cat('Done!\n')
        output <- merge(enriched,depleted,by = 'Gene')
    }
    if(score.method == 'KS'){
        cat('Calculating KS-stats\n')
        enriched <- KStest_Samba(data = sgRes, nullData = nullData, 
                                 fdr.threshold = sg.threshold.up, 
                                 direction = 1, filter_beta = T)
        depleted <- KStest_Samba(data = sgRes, nullData = nullData, 
                                 fdr.threshold = sg.threshold.dn, 
                                 direction = -1, filter_beta = T)
        output <- merge(enriched,depleted,by = 'Gene', all = T)
        output <- output[order(output$pval_pos, -output$score_pos, 
                               decreasing = F),]
        cat('Done!\n')
    }
    if(score.method == 'Riger'){
        cat('Calculating Riger-stats\n')
        output <- RigerTest_Samba(data = sgRes, nullData = nullData)
        output <- output[order(output$pval_pos, -output$adjScore, decreasing = F),]
        cat('Done!\n')
    }
    # sort results
    output <- output[order(output$pval_pos),]
    
    # save data
    if(!is.null(file.prefix)){
        write.table(output, file = paste0(file.prefix, '_GeneLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }
    
    return(output)
}


#' @export
GeneScoreSamba <- function(data, nullData, sg.threshold, direction = 1){
    data$logFC <- direction * data$logFC
    nullData$logFC <- direction * nullData$logFC

    # Get Null scores
    nullData <- nullData[with(nullData, order(Gene, -logFC)),]
    nullScores <- WtSumScore(nullData, sg.threshold)
    nullScores <- nullScores$score
    
    # Get gene scores
    data <- data[with(data, order(Gene, -logFC)),]
    targetScores <- WtSumScore(data, sg.threshold)
    
    ## Adjust scores (normalize to mean of the random positive or negative scores)
    positiveScoreAdjustmentFactor = mean(nullScores[which(nullScores >= 0)])
    negativeScoreAdjustmentFactor = -1 * mean(nullScores[which(nullScores < 0)])
    targetScores$score[which(targetScores$score > 0)] <- 
        targetScores$score[which(targetScores$score > 0)]/positiveScoreAdjustmentFactor
    targetScores$score[which(targetScores$score < 0)] <- 
        targetScores$score[which(targetScores$score < 0)]/negativeScoreAdjustmentFactor
    
    # Calculate p and q values
    targetScores$pval <- pnorm(q = as.numeric(targetScores$score),  lower.tail = F)#, 
    targetScores$pval <- targetScores$pval/length(nullScores)
    targetScores$fdr <- p.adjust(targetScores$pval, method='fdr')

    # Generate output
    output <- targetScores[,c('Gene','score','pval','fdr','n_fdr_guides')]
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output) <- c('Gene',paste0('score_',direction),
                          paste0('pval_',direction), paste0('fdr_',direction), 
                          paste0('fdrGuides_',direction))
    return(output)
}


## FilterCountData function generates list of guides whose combined counts 
##   across screen samples is greater than the desired threshold (min.guides).
#' @export
FilterCountData <- function(# input the DGE object, group info and min.guides
        dge, group = NULL, min.guides, verbose=T){
    # print message, if verbose setting is true
    if(verbose) {
        cat('Filtering guides with low representation across screen samples.\n')
    }
    # if group info provided, choose samples with coefficient > 0
    if(!is.null(group)){
        samples.screen <- which(group != 0)
    } else { # if no group info, choose all samples
        samples.screen <- colnames(dge$counts)
    }
    # get sum of guide counts across samples
    keep.exprs = dge$counts[,samples.screen] %>% rowSums(.)
    # get names of guides whose total counts are greater than min.guides
    keep.exprs <- keep.exprs[keep.exprs > min.guides] %>% names()
    # if verbose, output number of guides excluded
    if(verbose) {
        cat(paste0('   Removed ',nrow(dge) - length(keep.exprs),' of ',nrow(dge) ,' guides.\n'))
    }
    # return filtered guide names
    return(keep.exprs)
}

## WeighGuides function adds a weight matrix to the DGE object, based on the 
##   sigmoid function of the guide-detection rate across samples.
#' @export
WeighGuides <- function(dge){ # input DGE object
    # extract counts
    df <- edgeR::cpm(dge, log = F)
    # set min guide detection rate to 5 percentile of all counts or at least 1
    min_guide_for_det <- max(quantile(rowMeans(df),.05),1)
    # for each guide, get proportion of samples detected > min. detection rate
    wt <- apply(df,1, function(x) {length(x[x > min_guide_for_det])/ncol(df)})
    # add pseudocount of 0.1 
    wt <- ((wt) + 0.1)
    # calculate wt from sigmoid function
    wt <- 1/(1+exp(-wt))
    # generate matrix of wts with same dimnames as DGE object
    wt <- matrix(rep(wt, ncol(dge)), nrow = nrow(dge), byrow = F)
    dimnames(wt) <- dimnames(dge)
    # add wt matrix to DGE object
    dge$weights = wt
    # return updated DGE object
    return(dge)
}



#' @export
RandomIndexGenerator <- function(sgPerGene, nTotalGuides){
    # choose number of null guide genes as the number of genes in library,
    #   or at least 100,000
    nIterations <- max(length(sgPerGene), 100000)
    # randomly select # of guides per null genes, using a seed
    set.seed(42)
    if(length(sgPerGene) < nIterations) sgPerGene <- dqrng::dqsample(sgPerGene, size = nIterations, replace = T)
    # make sure there are at least 2 guides per null gene
    sgPerGene[which(sgPerGene < 2)] <- 2
    set.seed(42)
    
    randomSeeds <- dqrng::dqsample(1:nIterations, replace = F)
    randData <- lapply(1:nIterations, function(i) {
        set.seed(randomSeeds[i])
        dqrng::dqsample(1:nTotalGuides, size = min(nTotalGuides,sgPerGene[i]), replace = F)
    })
    return(randData)
}

#' KStest_Samba
#'
#' This function performs SAMBA gene-level aggregation, using the 
#'   Kolmogorov-Smirnov test. Note that this is an internal function.
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
    sg.threshold <- quantile(nullData$logFC, 1-fdr.threshold)
    #if(direction == 1) cat(paste0('FDR threshold for guides: ', sg.threshold,'\n'))
    nullScores <- as.numeric(nullData$logFC)
    # Calculate p & q values
    ifelse(direction==1, alt <- 'less', alt <- 'greater')
    runKS <- function(x) return(ks.test(x = x, 'pnorm', alternative = alt, 
                                        simulate.p.value = F, 
                                        B = 10000)$p.value)
    output <- dplyr::reframe(data, .by = 'Gene', 
                             score = quantile(scale, 0.75), 
                             pval = runKS(scale))
    qval <- function(x) qvalue::lfdr(x, pi0=limma::propTrueNull(x))
    output$qval <- qval(output$pval)
    # Generate output
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output)[-1] <- paste0(colnames(output)[-1],'_',direction)
    return(output)
}

#' RigerTest_Samba
#'
#' This function performs riger gene-level aggregation, using a modified 
#'   Kolmogorov-Smirnov test. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of guide-level screen enrichment data for the null distribution.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
RigerTest_Samba <- function(data, nullData){
    ScoreGeneSet_Riger <- function(scores, n_total_guides, score.rank, wts=NULL, alpha=1) {
        n_guides <- length(scores)
        if(is.null(wts)) wts <- rep(1, length(scores))
        nonTargetSetScoreWeight = -1 / (n_total_guides - length(scores))
        weightedTargetSetScores = abs(scores * wts)
        sumOfWeightedTargetSetScores = sum(weightedTargetSetScores^alpha)
        # proportion of wtScore
        targetHairpinIndexToScoreMap <- weightedTargetSetScores / sumOfWeightedTargetSetScores        
        orderedscore.rank = sort(score.rank)
        cumulativeScore = 0
        maxCumulativeScore = 0
        minCumulativeScore = 0
        lastTargetHairpinIndex = -1
        for (i in 1:length(orderedscore.rank)) {
            targetHairpinIndex = orderedscore.rank[i]
            numSkippedNonTargetSetScores = targetHairpinIndex - lastTargetHairpinIndex - 1
            # because nonTargetSetScoreWeight is negative, no need to update maxCumulativeScore here
            cumulativeScore = cumulativeScore + (numSkippedNonTargetSetScores * nonTargetSetScoreWeight)
            minCumulativeScore = min(minCumulativeScore, cumulativeScore)
            # because the target hairpin score is positive, no need to update minCumulativeScore here
            cumulativeScore = cumulativeScore + targetHairpinIndexToScoreMap[i]
            maxCumulativeScore = max(maxCumulativeScore, cumulativeScore)
            lastTargetHairpinIndex = targetHairpinIndex
        }
        numSkippedNonTargetSetScores = n_total_guides - lastTargetHairpinIndex - 1
        # because nonTargetSetScoreWeight is negative, no need to update maxCumulativeScore here
        cumulativeScore = cumulativeScore + (numSkippedNonTargetSetScores * nonTargetSetScoreWeight)
        minCumulativeScore = min(minCumulativeScore, cumulativeScore)
        if (abs(maxCumulativeScore) > abs(minCumulativeScore)){
            return(signif(maxCumulativeScore, 5))
        } else{
            return(signif(minCumulativeScore, 5))
        }
    }
    GeneScoreAdjuster_Riger <- function(x, nullScores) {
        positiveScoreAdjustmentFactor = mean(nullScores[which(nullScores >= 0)])
        negativeScoreAdjustmentFactor = -1 * mean(nullScores[which(nullScores < 0)])
        x[which(x >= 0)] <- x[which(x >= 0)]/positiveScoreAdjustmentFactor
        x[which(x < 0)] <-  x[which(x < 0)]/negativeScoreAdjustmentFactor
        return(x)
    }
    # for each gene
    # calculate score
    n_total_guides <- nrow(data)
    data$score.rank <- rank(data$logFC)
    gs <- dplyr::reframe(data, .by='Gene', geneScore=
                             ScoreGeneSet_Riger(logFC, n_total_guides, score.rank))
    # calculate scores for nullGenes
    n_total_guides <- nrow(nullData)
    nullData$score.rank <- rank(nullData$logFC)
    ns <- dplyr::reframe(nullData, .by='Gene', geneScore=
                             ScoreGeneSet_Riger(logFC, n_total_guides, score.rank))
    # adjust score
    gs$adjScore <- GeneScoreAdjuster_Riger(gs$geneScore, ns$geneScore)
    # calculate proportion of random scores that are better
    ns <- sort(ns$geneScore)
    gs$pval_pos <- (lapply(gs$geneScore, function(s) sum(s > ns)) %>% unlist())/length(ns)
    gs$pval_neg <- (lapply(gs$geneScore, function(s) sum(s < ns)) %>% unlist())/length(ns)
    # calculate q values
    qval <- function(x) qvalue::lfdr(x, pi0=limma::propTrueNull(x))
    gs$qval_pos <- qval(gs$pval_pos)
    gs$qval_neg <- qval(gs$pval_neg)
    return(gs)
}


#' WtSumScore
#'
#' This function calculates a single gene score for a given set of guide log-FC values.
#'
#' @param scores numeric vector of the log-FC values for all guides of a single gene.
#' @param sg.threshold numeric, giving the logFC value of the null data that represents the FDR threshold
#' @return numeric score of a gene.
#' @export
WtSumScore <- function(data, fdr.thresh){
    wt.increment <- 1/4
    # sort grna-level results by logfc
    data <- data[with(data, order(Gene, -logFC, decreasing = FALSE)),]
    # get the sum of the top half of grna log-FCs, 
    #   using quantiles to account for difference in the #s of gRNAs/Gene
    data <- mutate(data, .by='Gene', n_guides=length(logFC))
    median_n_guide <- median(slice_head(data, by='Gene', n=1)$n_guides)
    data$sf <- median_n_guide/data$n_guides
    score_tophalf <- dplyr::reframe(data, .by= 'Gene', n_guides=n_guides[1],
                                    score_tophalf = sum(quantile(logFC,c(0.5,0.75,1)), na.rm=T))
    
    # calculate weighted sum of guides that pass FDR threshold
    if(nrow(data[which(data$logFC > fdr.thresh),]) > 0){
        score_fdr <- data[which(data$logFC > fdr.thresh),] %>% 
            dplyr::mutate(.by='Gene', n_fdr_guides = length(logFC))
        score_fdr <-  dplyr::mutate(score_fdr, .by='Gene', wt = 
                                        (wt.increment * sf[1] * (1:length(logFC)))) 
    } else {
        score_fdr <- data %>% 
            dplyr::group_by(Gene) %>% 
            dplyr::mutate(wt = 0) 
    }
    # create size factor to normalize the weighted sum score, 
    #  relative to the number of guides for that gene
    score_fdr <- dplyr::reframe(score_fdr, .by='Gene', 
                                score_fdr=sum(logFC * wt),
                                n_fdr_guides=n_fdr_guides[1])
    
    # merge the score_tophalf, score_fdr, and count_fdr
    scores <- merge(score_fdr, score_tophalf, by = 'Gene', all = T)
    scores[is.na(scores)] <- 0
    
    scores <- data.frame(Gene = scores$Gene, 
                         score_fdr = scores$score_fdr, 
                         score_tophalf = scores$score_tophalf,
                         score = (scores$score_fdr + scores$score_tophalf), 
                         n_fdr_guides = scores$n_fdr_guides,
                         n_guides=scores$n_guides)
    
    return(scores)
}

