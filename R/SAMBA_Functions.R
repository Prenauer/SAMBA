#' Samba
#'
#' This is an all-in-one function to preprocesses and analyze read count data from a CRISPR
#' screen. The structure of the input data is assumed to have columns labeled "sgRNA" and
#' "Gene", followed by read counts from each sample in appropriately labeled columns.
#'
#' @param data data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @param screen.names character vector of the screen sample names.
#' @param ctrl.names character vector of the control sample names.
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param control.gene character, gene label of the internal control guides (e.g. "NTC").
#' @param score.method character, indicating whether to use the 'KS' (Kolmogorov-Smirnov), 'MetaAnalysis', or 'GeneScore' method to perform gene-level analysis.
#' @param test.method character, indicating whether to use 'QLF' or 'LRT' tests.
#' @param GuideMap data.frame, maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.
#' @param file.prefix character, prefix for writing output files.
#' @param verbose logical, indicating whether to display verbose output information.
#' @return list, contains the guide-level result data.frame and the gene-level result data.frame.
#' @export
Samba <- function(data, design=NULL, screen.names=NULL, ctrl.names=NULL,
                  coefficient = NULL, contrast = NULL,
                  control.gene = NULL,score.method = 'GeneScore',
                  test.method = 'QLF', GuideMap = NULL, file.prefix = NULL,
                  verbose=T){

    # if both screen.names and ctrl.names are given...
    if(!is.null(screen.names) & !is.null(ctrl.names)){
        # set the coefficient to 'screen'
        coefficient <- 'screen'
    }

    # preprocess data to generate DGE object
    dge <- Preprocess_Samba(data = data, design = design,
                            screen.names=screen.names, ctrl.names=ctrl.names)
    # perform guide-level analysis
    sgRes <- Analyze_Samba_Guides(dge = dge, coefficient = coefficient,
                                  method = test.method, contrast = contrast,
                                  file.prefix = file.prefix)
    # perform gene-level analysis
    geneRes <- Analyze_Samba_Genes(sgRes = sgRes, score.method = score.method,
                                   control.gene=control.gene,
                                   file.prefix = file.prefix)
    # return guide-level and gene-level results in a list
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
#' @param screen.names character vector of the screen sample names.
#' @param ctrl.names character vector of the control sample names.
#' @param hdr.gene character, name of column with gene names. Defaults to "Gene".
#' @param hdr.guide character, name of column with guide names. Defaults to "sgRNA".
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param verbose logical, indicating whether to display verbose output information.
#' @return DGEList object of the edgeR package.
#' @export
Preprocess_Samba <- function(data, design=NULL, screen.names=NULL,
                             ctrl.names=NULL, hdr.gene = 'Gene',
                             hdr.guide = 'sgRNA', coefficient=NULL,
                             verbose=T){
    ## Check that there is either design matrix or sample names
    # if there is no design matrix...
    if(is.null(design)){
        # if either screen.names or ctrl.names is not provided...
        if(is.null(screen.names) | is.null(ctrl.names)){
            # print warning
            warning('Please provide either design matrix or sample names.')
        }
        # if both screen.names and ctrl.names are given...
        if(!is.null(screen.names) & !is.null(ctrl.names)){
            # create a design matrix from sample names.
            design <- CreateDesignMatrix(colnames(data)[3:ncol(data)],
                                         screen.names, ctrl.names)
        }
    }

    # if verbose, print message
    if(verbose) cat('Starting preprocessing.\n')

    # set minimum guides threshold to zero, removing no guides in preprocessing.
    min.guides = 0
    # set pseudocount to 4. This will be added to all counts to create noise
    #   where small guide-counts are similar to the 0 counts. This method
    #   performs better than using a hard threshold
    pseudocount = 4

    # convert input data to data.frame
    data <- data.frame(data)
    # convert input design to matrix format
    design <- as.matrix(design)

    # check that data has gene and guide columns
    if(sum(colnames(data) %in% c(hdr.gene, hdr.guide)) != 2)
        warning('Count data is missing the columns: "Guide" and "Gene"!')
    # check that design matrix has the correct # of rows
    if((nrow(design) + 2) != ncol(data))
        warning('Design matrix has incorrect number of rows!')
    # check if design matrix has integer or numeric format
    if((mode(design) != 'numeric') & (mode(design) != 'integer'))
        warning('Design matrix does not have integer/numeric format!')

    # if design matrix has one column, add a baseline column
    if(ncol(design) == 1){
        design <- cbind(baseline=1, design)
    }

    ## determine groups of samples for performing guide-filtering
    # if coefficient name is given, match it to column-name of design.
    if(!is.null(coefficient)){
        if(sum(colnames(design) %in% coefficient) > 0){
            # set groups based on coefficient column
            group <-
                as.integer(design[,which(colnames(design)==coefficient)])
        }
    }
    # if there are two columns and the first is baseline...
    if((ncol(design) == 2) & length(unique(design[,1]))==1){
        # set groups based on the non-baseline column
        group <- as.integer(design[,2])
    } else { # otherwise do not use grouping for filtering.
        group = NULL
        if(verbose) cat('No sample grouping will be used!\n')
    }

    ## Make/filter DGE object
    # name data rows with guide-names
    rownames(data) <- data[,hdr.guide]
    # if min.guides is set (hard threshold), filter guides with low counts.
    if(min.guides > 0){
        # create DGE object with count data, and set groups
        dge <- edgeR::DGEList(data[,-c(1:2)], remove.zeros = F, group = group)
        # run FilterCountData function to get names of filtered guides
        keep.exprs <- FilterCountData(dge = dge, group=group,
                                      min.guides=min.guides)
        # filter guides, add pseudocount, and remove zero counts.
        dge <- edgeR::DGEList(data[keep.exprs,-c(1:2)] + pseudocount,
                              remove.zeros = T, group = group)
    # if min.guides isn't set...
    } else {
        # create DGE object with count data, add pseudocount,
        #   remove zero counts, and set groups
        dge <- edgeR::DGEList(data[,-c(1:2)] + pseudocount,
                              remove.zeros = T, group = group)
    }

    # generate guide weights using WeighGuides function, and add to DGE object
    dge <- WeighGuides(dge = dge)

    # calculate size factors
    if(verbose) cat('Calculating TMM-wsp size factors\n')
    dge <- edgeR::calcNormFactors(dge, method = 'TMMwsp')

    # estimate dispersion, if more than 2 samples
    if(nrow(dge$samples) > 2){
        if(verbose) cat('Calculating dispersions ...\n')
        # estimate common dispersion
        dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
        if(verbose) cat('   Calculated common dispersion\n')
        # estimate trended dispersion
        dge <- edgeR::estimateGLMTrendedDisp(dge, design = design,
                                             method = 'bin.loess', span = 1/3)
        if(verbose) cat('   Calculated trended dispersion\n')
        # estimate guide-wise dispersion
        dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
        if(verbose) cat('   Calculated guide-wise dispersion\n')
    }

    # add design matrix to DGE
    dge$Design <- design
    # add guide mapping info to DGE
    dge$GuideMap <- data[,c(hdr.guide, hdr.gene)]
    if(verbose) cat('Preprocessing complete!\n')

    # return DGE object
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
#' @param verbose logical, indicating whether to display verbose output information.
#' @return data.frame of the results from the guide-level screen enrichment analysis.
#' @export
Analyze_Samba_Guides <- function(dge, method = 'QLF', design = NULL,
                                 coefficient = NULL, contrast = NULL,
                                 GuideMap = NULL, file.prefix = NULL,
                                 verbose = T){
    ## Get design matrix and guide info
    # if no guidemap is provided, get it from the DGE
    if(is.null(GuideMap)) GuideMap <- dge$GuideMap
    # if no guidemap provided or in the DGE, print warning.
    if(is.null(GuideMap))
        warning('Please provide a dataframe that maps guides to genes!')
    # get design info from DGE
    if(is.null(design)) design <- dge$Design

    ## Check for valid inputs
    # check that there is either a coefficient or a contrast
    if(!xor(is.null(coefficient), is.null(contrast))){
        # try to designate a coefficient from design matrix
        coefficient <- colnames(design)[ncol(design)]
        warning(paste0('No coefficient provided. Using "',coefficient,
                       '" as the coefficient.'))
    }
    # check that QLF or LRT is indicated
    if(!(method %in% c('QLF','LRT')))
        warning('Please choose either "LRT" or "QLF" as an analysis method!')


    ## sgRNA-level analysis
    # if more than 2 samples...
    if(nrow(dge$samples) > 2){
        if(verbose) cat('Performing sgRNA-level analysis ...\n')
        # if quasi-likelihood F test selected (default)
        if(method == 'QLF'){
            # fit NB GLM with quasi-dispersion estimates to data
            fit <- edgeR::glmQLFit(dge, design = design, robust = T)
            # run QL F-test on fitted data
            sgRes <- edgeR::glmQLFTest(fit, coef = coefficient,
                                       contrast = contrast)
        }
        # if likelihood ratio test chosen...
        if(method == 'LRT'){
            # fit NB GLM to data
            fit <- edgeR::glmFit(dge, design = design, robust = T)
            # run LRT
            sgRes <- edgeR::glmLRT(fit, coef = coefficient, contrast = contrast)
        }
        if(verbose) cat('Done!\n')
        # extract results as data.frame
        sgRes <- edgeR::topTags(sgRes, n = Inf, sort.by = 'PValue')$table
    }

    ## create a result table when only 2 samples given
    # if only 2 samples provided...
    if(nrow(dge$samples) == 2){
        # identify screen sample (coefficient == 1)
        samp.screen <- which(dge$Design[,coefficient]==1)
        # identify control sample (coefficient == 0)
        samp.ctrl <- which(dge$Design[,coefficient]==0)
        # get log-normalized count-data
        normdata <- edgeR::cpm(dge, log=T)
        # add guide names to normdata row names
        rownames(normdata) <- rownames(dge$counts)
        # create data.frame of results
        sgRes <- data.frame(logFC=normdata[,samp.screen]-normdata[,samp.ctrl],
                            logCPM=aveLogCPM(dge), Statistic=NA, PValue=NA,
                            FDR=NA)
    }

    ## Format guide-level results
    # change 3rd column name from "F" (might cause problems) to "Statistic"
    colnames(sgRes)[3] <- 'Statistic'
    # add guide names to row names
    sgRes$sgRNA <- rownames(sgRes)
    # add guidemap to sgRes
    sgRes <- merge(GuideMap, sgRes, by = 'sgRNA', all.x = F, all.y = T)
    # create z-score of results. Not actually used for analysis
    sgRes$zscore <- scale(sgRes$logFC, center = F, scale = T)[,1]
    # sort rows by the test statistic, then logFC
    sgRes <- sgRes[with(sgRes, order(Statistic, logFC, decreasing = T)),]

    # if file prefix given...
    if(!is.null(file.prefix)){
        # save fitted data to RDS file
        saveRDS(fit, file = paste0(file.prefix, '_FittedData.rds'))
        # write results to file
        write.table(sgRes, file = paste0(file.prefix, '_GuideLevelResults.txt'),
                    sep = '\t', quote = F, row.names = F)
    }

    # return sgRes data.frame
    return(sgRes)
}

#' Analyze_Samba_Genes
#'
#' This function performs the gene-level screen enrichment analysis of SAMBA, using guide-level SAMBA results as an input.
#'
#' @param sgRes data.frame of the results from the guide-level screen enrichment analysis. This is the output from the "Analyze_Samba_Guides" function.
#' @param control.gene character, gene label of the internal control guides (e.g. "NTC").
#' @param score.method character, indicating whether to use the 'KS' (Kolmogorov-Smirnov) or 'GeneScore' method to perform gene-level analysis.
#' @param fdr.threshold numeric, indicating the proportion of the control guides, or total guides if control.gene = NULL, will be used as a threshold for the weighted sum score.
#' @param file.prefix (optional) file.path and/or character prefix to use for exporting the result tables
#' @param verbose logical, indicating whether to display verbose output information.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
Analyze_Samba_Genes <- function(sgRes, control.gene = NULL,
                                score.method = 'GeneScore', fdr.threshold=0.1,
                                file.prefix = NULL,
                                verbose = T){
    ## Check inputs
    # if a supported score.method is not set, print a warning
    if(!(score.method %in% c('GeneScore','KS','Riger')))
        warning('Please choose either "GeneScore" or "KS" as a score.method!')

    # if an internal-control gene is specified...
    if(!is.null(control.gene)){
        # if there are at least 50 internal-control null guides...
        if(sum(sgRes$Gene  == control.gene) >= 50){
            # get logFC's of internal-control null guides
            ntc <- sgRes$logFC[which(sgRes$Gene  == control.gene)]
            # set enrichment threshold. Default is 90th percentile (fdr=10%)
            sg.threshold.up <- quantile(ntc, (1-fdr.threshold))
            # set depletion threshold. Default is 10th percentile (fdr=10%)
            sg.threshold.dn <- quantile(ntc, fdr.threshold)
        }
        # if there are less than 50 internal-control null guides...
        if(sum(sgRes$Gene  == control.gene) < 50){
            # set control.gene to NULL, so that it isn't used.
            control.gene <- NULL
        }
    }
    # if control.gene is NULL...
    if(is.null(control.gene)){
        # set control.gene to ''
        control.gene <- ''
        # set enrichment threshold, based on targeting-guide data
        sg.threshold.up <- quantile(sgRes$logFC, (1-fdr.threshold))
        # set depletion threshold, based on targeting-guide data
        sg.threshold.dn <- quantile(sgRes$logFC, fdr.threshold)
    }


    ## Filter sgRes guides for Gene-level analysis
    # filter sgRes for NA genes
    sgRes <- sgRes[!is.na(sgRes$Gene ),]
    # remove internal-control guides, if present
    sgRes <- sgRes[which(sgRes$Gene != control.gene),]


    ## Create random data to use for null distribution calculations
    if(verbose) cat('Generating null distribution.\n')
    # get # of guides/gene
    sgPerGene <- dplyr::reframe(sgRes, .by='Gene', freq=length(Gene))
     # generate randomized, reproducible guide-gene data for a null distribution
    rIndices <- RandomIndexGenerator(sgPerGene = sgPerGene$freq,
                                     nTotalGuides = nrow(sgRes))
    # rearrange sgRes with random indices
    nullData <- sgRes[unlist(rIndices),]
    # calculate length of null genes
    i <- lapply(rIndices, length)
    # create labels for null genes
    nullData[,'Gene'] <- paste0('Null', lapply(1:100000, function(x){
        rep(x,i[[x]]) }) %>% unlist())


    ## perform enrichment and depletion analyses
    # if 'GeneScore' is used...
    if(score.method == 'GeneScore'){
        if(verbose) cat('Calculating enrichment.\n')
        # analyze enrichment results
        enriched <- GeneScoreSamba(data = sgRes, nullData = nullData,
                                   sg.threshold = sg.threshold.up,
                                   direction = 1)
        if(verbose) cat('Calculating depletion.\n')
        # analyze depletion results
        depleted <- GeneScoreSamba(data = sgRes, nullData = nullData,
                                   sg.threshold = sg.threshold.dn,
                                   direction = -1)
        if(verbose) cat('Done!\n')
        # merge enrichment and depletion result data.frames
        output <- merge(enriched,depleted,by = 'Gene')
        # sort output based on enrichment results
        output <- output[order(output$pval_pos),]
    }
    # if 'KS' is used...
    if(score.method == 'KS'){
        if(verbose) cat('Calculating KS-stats\n')
        # analyze enrichment results
        enriched <- KStest_Samba(data = sgRes, direction = 1)
        # analyze depletion results
        depleted <- KStest_Samba(data = sgRes, direction = -1)
        # merge enrichment and depletion result data.frames
        output <- merge(enriched,depleted,by = 'Gene', all = T)
        # sort output based on enrichment results
        output <- output[order(output$pval_pos, -output$score_pos,
                               decreasing = F),]
        if(verbose) cat('Done!\n')
    }


    ## write and/or return analysis results
    # if a prefix is provided...
    if(!is.null(file.prefix)){
        # write results to file
        write.table(output, file = paste0(file.prefix, '_GeneLevelResults.txt'),
                    sep = '\t', quote = F, row.names = F)
    }
    # return data.frame of results
    return(output)
}


#' GeneScoreSamba
#'
#' This function performs SAMBA gene-level aggregation, using the gene score method. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of randomly arranged null-genes. This is used to create the null score distribution.
#' @param sg.threshold numeric, indicating the internal-control threshold to use, based on FDR.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
GeneScoreSamba <- function(data, nullData, sg.threshold, direction = 1){
    ## Adjust data and nullData for depletion, if needed
    # if depletion is indicated (direction = -1), flip sign of data logFC
    data$logFC <- direction * data$logFC
    # if depletion is indicated (direction = -1), flip sign of nullData logFC
    nullData$logFC <- direction * nullData$logFC
    # if depletion is indicated (direction = -1), flip sign of sg.threshold
    sg.threshold <- direction * sg.threshold

    ## Get Null scores
    # sort nullData by 'gene', then 'logFC'
    nullData <- nullData[with(nullData, order(Gene, -logFC)),]
    # get weighted-sum score
    nullScores <- WtSumScore(nullData, sg.threshold)
    # keep only the scores for the nullData
    nullScores <- nullScores$score

    ## Get gene scores
    # sort data by 'gene', then 'logFC'
    data <- data[with(data, order(Gene, -logFC)),]
    # get weighted-sum scores
    targetScores <- WtSumScore(data, sg.threshold)

    ## normalize results without altering sign
    # pos scale factor = mean positive null scores
    pos.sf = mean(nullScores[which(nullScores > 0)])
    # neg scale factor = absolute mean negative null scores
    neg.sf = abs(mean(nullScores[which(nullScores < 0)]))
    # positive scores are normalized to mean of the null positive scores
    targetScores$score[which(targetScores$score > 0)] <-
        targetScores$score[which(targetScores$score > 0)]/pos.sf
    # negative scores are normalized to mean of the null negative scores
    targetScores$score[which(targetScores$score < 0)] <-
        targetScores$score[which(targetScores$score < 0)]/neg.sf


    ## Calculate p and q values
    # calculate pval from normal probability distribution function
    targetScores$pval <- pnorm(q=as.numeric(targetScores$score), lower.tail=F)
    # calculate fdr-adjusted p values
    targetScores$fdr <- p.adjust(targetScores$pval, method='fdr')

    ## Generate output
    # rearrange data.frame to get output format
    output <- targetScores[,c('Gene','score','pval','fdr','n_fdr_guides')]
    # get direction-based suffix
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    # add direction-suffix to appropriate column names
    colnames(output) <- c('Gene',paste0('score_',direction),
                          paste0('pval_',direction), paste0('fdr_',direction),
                          paste0('fdrGuides_',direction))
    # return output
    return(output)
}


#' FilterCountData
#'
#' This function filters DGEList of raw count data, based on a minimum number of total counts.
#'
#' @param dge DGEList object with raw count data.
#' @param group factor, indicating the grouping of all samples. Default is NULL.
#' @param min.guides numeric, indicating the minimum number of counts for each guide across all samples. Default is 0.
#' @return character vector of the guides that passed filtering.
#' @export
FilterCountData <- function(
        dge, group = NULL, min.guides, verbose=T){
    # print message, if verbose setting is true
    if(verbose) {
        cat('Filtering guides with low representation across screen samples.\n')
    }

    # if group info provided, choose samples with coefficient != 0
    if(!is.null(group)){
        samples.screen <- which(group != 0)
    } else { # if no group info, choose all samples instead of group-filtering
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


#' WeighGuides
#'
#' This function calculates a weight for each guide, based on the number of screen samples that have detected counts.
#'
#' @param dge DGEList object with raw count data.
#' @return DGEList object that contains guide weights as a numeric vector.
#' @export
WeighGuides <- function(dge){
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


#' RandomIndexGenerator
#'
#' This function creates a list of reproducible, random integers from which to create a null distribution.
#'
#' @param sgPerGene integer vector, representing the number of guides for each gene.
#' @param nTotalGuides integer, indicating the total number of guides in the filtered library.
#' @return list of integer vectors, each of random numbers from which to construct a null data distribution.
#' @export
RandomIndexGenerator <- function(sgPerGene, nTotalGuides){
    # choose number of null guide genes as the number of genes in library,
    #   or at least 100,000
    nIterations <- max(length(sgPerGene), 100000)
    # randomly select # of guides per null genes, using a seed
    set.seed(42)
    if(length(sgPerGene) < nIterations)
        sgPerGene <- dqrng::dqsample(sgPerGene, size = nIterations, replace = T)
    # make sure there are at least 2 guides per null gene
    sgPerGene[which(sgPerGene < 2)] <- 2
    # randomly select seeds for each null-gene, using a starting seed
    set.seed(42)
    randomSeeds <- dqrng::dqsample(1:nIterations, replace = F)
    # iterate through null-genes
    randData <- lapply(1:nIterations, function(i) {
        # set null-gene seed, then select random sgPerGene
        set.seed(randomSeeds[i])
        dqrng::dqsample(1:nTotalGuides,
                        size = min(nTotalGuides,sgPerGene[i]),
                        replace = F)
    })
    # return randomized indices
    return(randData)
}


#' KStest_Samba
#'
#' This function performs SAMBA gene-level aggregation, using the
#'   Kolmogorov-Smirnov test. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of randomly arranged null-genes. This is used to create the null score distribution.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
KStest_Samba <- function(data, nullData, direction=1){
    # sort data by logFC
    data <- data[order(data$logFC, decreasing = T),]

    # get null distribution as nullData logFC
    nullDist <- nullData$logFC

    ## create function for KS analysis
    # set alternative hypothesis, based on enrichment or depletion analysis
    ifelse(direction==1, alt <- 'less', alt <- 'greater')
    # runKS function will perform two-sample KS test, comparing the logFC of
    #   a gene to the null distribution. One-tailed p-value will be given for
    #   enrichment or depletion.
    runKS <- function(x) return(ks.test(x = x, y=nullDist,
                                        alternative=alt,
                                        simulate.p.value = F,
                                        B = 10000)$p.value)
    # use runKS to get pvals for each gene set.
    output <- dplyr::reframe(data, .by = 'Gene',
                             pval = runKS(logFC))

    ## Arrange output data.frame
    # get rank of results
    output$rank <- rank(output$pval)
    # calculate fdr-adjusted p value.
    output$fdr <- p.adjust(output$pval)
    # rearrange order of data.frame columns
    output <- output[,c('Gene','score','pval','fdr')]
    # rename score, pval, and fdr columns with either pos for enrichment or
    #   neg for depletion.
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output)[-1] <- paste0(colnames(output)[-1],'_',direction)

    # return results
    return(output)
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
    # set the increment of the weighted sum
    wt.increment <- 1/4
    # sort grna-level results by logfc
    data <- data[with(data, order(Gene, -logFC, decreasing = FALSE)),]

    # add the # of guides/gene to the results
    data <- mutate(data, .by='Gene', n_guides=length(logFC))
    # calculate the median # of guides/gene
    median_n_guide <- median(slice_head(data, by='Gene', n=1)$n_guides)
    # calculate a size-factor for the # of guides/gene
    data$sf <- median_n_guide/data$n_guides


    # get the sum of the top half of grna log-FCs,
    #   using quantiles to account for difference in the #s of gRNAs/Gene
    score_tophalf <- dplyr::reframe(data, .by= 'Gene', n_guides=n_guides[1],
                                    score_tophalf =
                                        sum(quantile(logFC,c(0.5,0.75,1)),
                                            na.rm=T))

    ## calculate weighted sum of guides that pass the FDR threshold
    # if any guides pass the threshold...
    if(nrow(data[which(data$logFC > fdr.thresh),]) > 0){
        # filter guides that pass FDR threshold
        score_fdr <- data[which(data$logFC > fdr.thresh),] %>%
            # add # of guides that pass FDR threshold
            dplyr::mutate(.by='Gene', n_fdr_guides = length(logFC))
        # calculate wt as step-wise increment multiplied by size-factor
        score_fdr <-  dplyr::mutate(score_fdr, .by='Gene', wt =
                                        (wt.increment*sf[1]*(1:length(logFC))))
    } else { # if no guides pass the threshold...
        # set all wts to 0.
        score_fdr <- data %>%
            dplyr::group_by(Gene) %>%
            dplyr::mutate(wt = 0)
    }
    # Calculate weighted fdr score by multiplying weights and logFC
    score_fdr <- dplyr::reframe(score_fdr, .by='Gene',
                                score_fdr=sum(logFC * wt),
                                n_fdr_guides=n_fdr_guides[1])


    # merge the score_tophalf and score_fdr
    scores <- merge(score_fdr, score_tophalf, by = 'Gene', all = T)
    # pad the NA values of the score_fdr with 0's
    scores[is.na(scores)] <- 0

    # create summary data.frame. Add the score_fdr to the score_tophalf
    scores <- data.frame(Gene = scores$Gene,
                         score_fdr = scores$score_fdr,
                         score_tophalf = scores$score_tophalf,
                         score = (scores$score_fdr + scores$score_tophalf),
                         n_fdr_guides = scores$n_fdr_guides,
                         n_guides=scores$n_guides)

    # return results
    return(scores)
}


#' CreateDesignMatrix
#'
#' This function generates a simple design matrix when the user inputs a count file and a
#' list of names for the screen and control samples.
#'
#' @param hdr.names character vector of the count data column names.
#' @param screen.names character vector of the screen sample names.
#' @param ctrl.names character vector of the control sample names.
#' @return data.frame of the design matrix.
#' @export
CreateDesignMatrix <- function(hdr.names, screen.names, ctrl.names){
    # check if the input names are in the column-names for the count data
    if(sum(hdr.names %in% c(screen.names, ctrl.names)) != length(hdr.names)){
        warning('The provided sample names do not match the count data header.')
    }

    # create design matrix with intercept column, and a "screen" column, where
    #   screen samples have a 1.
    design <- data.frame(intercept=1,
                         screen=as.integer(hdr.names %in% screen.names))

    # return design matrix
    return(design)
}

