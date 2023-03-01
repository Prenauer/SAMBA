#' @import {plyr}
#' @import {dplyr}
#' @import {stringr}
#' @import {edgeR}
#' @import {limma}


#' SAMBA
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
SAMBA <- function(data, design, coefficient = NULL, contrast = NULL, ntc.as.null.dist = T,
                  score.method = 'MetaAnalysis', test.method = 'QLF', GuideMap = NULL, file.prefix = NULL){

    dge <- Preprocess_Samba(data = data, design = design)
    sgRes <- Analyze_Samba_Guides(dge = dge, coefficient = coefficient, method = test.method, contrast = contrast, file.prefix = file.prefix)
    geneRes <- Analyze_Samba_Genes(sgRes = sgRes, score.method = score.method, file.prefix = file.prefix)

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
    min.guides = 1
    pseudocount = 4
    dispersion.method = 'GLM'

    # check input
    data <- data.frame(data)
    design <- as.matrix(design)
    if(sum(colnames(data) %in% c(hdr.gene, hdr.guide)) != 2) warning('Count data is missing the columns: "Guide" and "Gene"!')
    if((nrow(design) + 2) != ncol(data)) warning('Design matrix has incorrect number of rows!')
    if((mode(design) != 'numeric') & (mode(design) != 'integer')) warning('Design matrix does not have integer/numeric format!')

    # get sample groups
    if(ncol(design) == 2){
        group <- as.integer(design[,-1])
    } else {
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

    # calculate size factors and dispersion
    cat('Calculating dispersions ...\n')
    dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
    cat('   Calculated common dispersion\n')
    dge <- edgeR::estimateGLMTrendedDisp(dge, design = design, method = 'bin.loess', span = 1/3)
    cat('   Calculated trended dispersion\n')
    dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
    cat('   Calculated sgRNA-wise dispersion\n')

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

#' @export
Analyze_Samba_Genes <- function(sgRes, ntc.as.null.dist = T, score.method = 'MetaAnalysis', file.prefix = NULL){
    # Set defaults
    control.gene <- 'NTC'
    hdr.gene <- 'Gene'
    fdr.threshold <- 0.1

    # Check inputs
    if(!(score.method %in% c('GeneScore','MetaAnalysis'))) warning('Please choose either "GeneScore" or "MetaAnalysis" as a score.method!')

    ## Filter guide-level data (>= 2sgRNA/gene and not NA)
    sgPerGene <- plyr::ddply(sgRes, 'Gene', nrow)
    keep.genes <- sgPerGene[which(sgPerGene[,2] > 1), 'Gene']
    sgRes <- sgRes[which(sgRes$Gene %in% keep.genes),]
    sgRes <- sgRes[!is.na(sgRes$Gene ),]

    ## set data for null dist
    cat('Generating null distribution.\n')
    ifelse(ntc.as.null.dist, nullData <- sgRes[which(sgRes$Gene  == control.gene),], nullData <- sgRes)
    sgRes <- sgRes[which(sgRes$Gene != control.gene),]
    nGeneGuides = sgPerGene[which(sgPerGene$Gene != control.gene),'V1']
    rIndices <- RandomIndexGenerator(sgPerGene = nGeneGuides, nTotalGuides = nrow(nullData))
    nullData <- nullData[unlist(rIndices),]
    i <- lapply(rIndices, length)
    nullData[,hdr.gene] <- paste0('NTC', lapply(1:100000, function(x){ rep(x,i[[x]]) }) %>% unlist())

    # perform enrichment and depletion analyses
    if(score.method == 'GeneScore'){
        cat('Calculating enrichment.\n')
        enriched <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1)
        cat('Calculating depletion.\n')
        depleted <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = -1)
        cat('Done!\n')
        output <- merge(enriched,depleted,by = 'Gene')
    }
    if(score.method == 'MetaAnalysis'){
        cat('Calculating enrichment.\n')
        output <- MetaAnalysisSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1)
        cat('Done!\n')
    }

    # save data
    if(!is.null(file.prefix)){
        write.table(sgRes, file = paste0(file.prefix, '_GeneLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }

    return(output)
}


MetaAnalysisSamba <- function(data, nullData, fdr.threshold, direction = 1){
    data$zscore <- direction * data$zscore
    sg.threshold <- quantile(nullData$zscore, 1-fdr.threshold)

    data <- data[with(data, order(zscore, decreasing = T)),]
    bottom <- data[duplicated(data$Gene),]
    bottom <- bottom[duplicated(bottom$Gene),]
    top2 <- setdiff(data$sgRNA, bottom$sgRNA)
    topfdr <- data[which(data$zscore > sg.threshold),'sgRNA']
    sg.keep <- union(topfdr, top2)
    data <- data[which(data$sgRNA %in% sg.keep),]
    fdr.sg <- plyr::ddply(data[which(data$sgRNA %in% topfdr),],'Gene',nrow)
    colnames(fdr.sg) <- c('Gene','FDRGuides')

    ## Gene-level aggregation
    geneRes <- plyr::ddply(data, 'Gene', summarise, score = mean(zscore), fdrGuides = length(Gene), pval = suppressWarnings(metap::sumlog(PValue)$p))
    geneRes$FDR <- p.adjust(geneRes$pval, method = 'fdr')
    geneRes <- merge(geneRes, fdr.sg, by = 'Gene', all.x = T, all.y = F)
    geneRes[which(is.na(geneRes$FDRGuides)),'FDRGuides'] <- 0
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(geneRes) <- c('Gene',paste0('score_',direction),paste0('FilteredGuides_',direction), paste0('pval_',direction),paste0('qval_',direction), paste0('fdrGuides_',direction))

    return(geneRes)
}

GeneScoreSamba <- function(data, nullData, fdr.threshold, direction = 1){
    data$logFC <- direction * data$logFC
    sg.threshold <- quantile(nullData$logFC, 1-fdr.threshold)

    # Get Null scores
    split.data <- split(nullData, factor(nullData$Gene))
    nullScores <- lapply(split.data, function(x) { WtSumScore(scores = x[,'logFC'], sg.threshold = sg.threshold) }) %>% unlist()
    fdr.sg <- lapply(split.data, function(x){ sum(x[,'logFC'] > sg.threshold) })

    # Get gene scores
    split.data <- split(data, factor(data$Gene))
    targetScores <- lapply(split.data, function(x) { WtSumScore(scores = x[,'logFC'], sg.threshold = sg.threshold) }) %>% unlist()
    fdr.sg <- lapply(split.data, function(x){ sum(x[,'logFC'] > sg.threshold) })

    ## Adjust scores (normalize to mean of the random positive or negative scores)
    positiveScoreAdjustmentFactor = sum(nullScores[which(nullScores >= 0)])/length(nullScores[which(nullScores >= 0)])
    negativeScoreAdjustmentFactor = -1 * sum(nullScores[which(nullScores < 0)])/length(nullScores[which(nullScores < 0)])
    adjTargetScores <- sapply(targetScores, function(unadjustedScore) {
        ifelse (unadjustedScore >= 0,
                return(unadjustedScore / positiveScoreAdjustmentFactor),
                return(unadjustedScore / negativeScoreAdjustmentFactor)
        )})

    # Calculate p/q values
    pval <- pnorm(q = as.numeric(targetScores), lower.tail = F, mean = mean(as.numeric(nullScores)), sd = sd(as.numeric(nullScores)))
    qval <- p.adjust(pval, 'fdr')


    # Generate output
    output <- data.frame(Gene = names(split.data), score = adjTargetScores, pval = pval, qval = qval, fdrGuides = unlist(fdr.sg))
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output) <- c('Gene',paste0('score_',direction),paste0('pval_',direction),paste0('qval_',direction), paste0('fdrGuides_',direction))
    return(output)
}

FilterCountData <- function(dge, group = NULL, min.guides){
    cat('Filtering guides with low representation across screen samples.\n')
    if(!is.null(group)){
        samples.screen <- which(group != 0)
    } else {
        samples.screen <- colnames(dge$counts)
    }
    keep.exprs = dge$counts[,samples.screen] %>% rowSums(.)
    keep.exprs <- keep.exprs[keep.exprs > min.guides] %>% names()
    cat(paste0('   Removed ',nrow(dge) - length(keep.exprs),' of ',nrow(dge) ,' guides.\n'))
    return(keep.exprs)
}

WeighGuides <- function(dge){
    df <- edgeR::cpm(dge, log = F)
    upper_rate_limit <- 1
    min_guide_for_det <- max(quantile(rowMeans(df),.05),1)
    wt <- apply(df,1, function(x) {min(length(x[x > min_guide_for_det])/ncol(df),upper_rate_limit)})
    wt <- ((wt) + .1)
    wt <- 1/(1+exp(-wt))
    wt <- matrix(rep(wt, ncol(dge)), nrow = nrow(dge), byrow = F)
    dimnames(wt) <- dimnames(dge)
    dge$weights = wt
    return(dge)
}

RandomIndexGenerator <- function(sgPerGene, nTotalGuides){
    nIterations <- max(length(sgPerGene), 100000) # at least 10000 genes
    set.seed(42)
    if(length(sgPerGene) < nIterations) sgPerGene <- sample(sgPerGene, size = nIterations, replace = T)
    sgPerGene[which(sgPerGene < 2)] <- 2
    set.seed(42)
    randomSeeds <- sample(1:nIterations, replace = F)
    randData <- lapply(1:nIterations, function(i) {
        set.seed(randomSeeds[i])
        sample(1:nTotalGuides, size = min(nTotalGuides,sgPerGene[i]), replace = F)
    })
    return(randData)
}


WtSumScore <- function(scores,sg.threshold){
    fdr.sg <- which(scores > sg.threshold)
    fdr.sg <- suppressWarnings(scores[fdr.sg])

    ifelse(length(scores) >= 4, min.sgPerGene <- 4, min.sgPerGene <- 2)
    n.sg <- max(length(scores[scores >= quantile(scores,.5)]),min.sgPerGene)

    if(length(fdr.sg) > 0){
        fdr.sg <- sum(sort(fdr.sg, decreasing = T) * seq(0.25,(length(fdr.sg)*0.25),0.25))
    } else {
        fdr.sg <- 0
    }
    top.sg <- sort(scores) %>% tail(.,n.sg)
    top.sg <- sum(top.sg)
    return(top.sg + fdr.sg)
}




