library("Seurat")
library("pROC")
library("parallel")


F1_score <- function(x, y) {

    if (sum(dim(x) == dim(y)) < 2) {
        stop("x and y have different dimensions\n")
    }
    if (sum(rownames(x) %in% rownames(y)) != nrow(x)) {
        stop("x and y have different rownames\n")
    }
    if (sum(colnames(x) %in% colnames(y)) != ncol(x)) {
        stop("x and y have different colnames\n")
    }
    y <- y[rownames(x),colnames(x)]
    F1 <- sapply(1:ncol(x), function(j) (x[,j]*y[,j])/(x[,j]+y[,j]))
    return(F1)
}


# compute Jaccard index between all gene set pairs
Jaccard <- function(gene.list.1, gene.list.2) {
    
    jaccard <- matrix(NA, nrow = length(gene.list.1), ncol = length(gene.list.2))
    rownames(jaccard) <- names(gene.list.1)
    colnames(jaccard) <- names(gene.list.2)
    i <- 1
    for (a in gene.list.1) {
        j <- 1
        for (b in gene.list.2) {
            jaccard[i,j] <- length(a[a %in% b]) / length(unique(c(a,b)))
            # jaccard[i,j] <- length(intersect(a,b)) / length(union(a,b)) # more elegant like this
            j <- j + 1
        }
        i <- i + 1
    }
    
    return(jaccard)
}


mean_abs_logFC_nhood_groups_from_data_frame <- function(df1, df2, param1, param2) {

    nhg1_logFC <- df1$logfc[df1$param == param1]
    nhg2_logFC <- df2$logfc[df2$param == param2]
    mean_logFC <- matrix(NA, nrow = length(nhg1_logFC), ncol = length(nhg2_logFC))
    groups1 <- df1$group[df1$param == param1]
    groups2 <- df2$group[df2$param == param2]
    rownames(mean_logFC) <- names(nhg1_logFC) <- groups1[!is.na(groups1)]
    colnames(mean_logFC) <- names(nhg2_logFC) <- groups2[!is.na(groups2)]
    for (i in names(nhg1_logFC)) { 
        mean_logFC[i,] <- sapply(names(nhg2_logFC), function(j) abs(mean(c(nhg1_logFC[i],nhg2_logFC[j]))))
    }
    
    return(mean_logFC)
}


delta_logFC_nhood_groups_from_data_frame <- function(df1, df2, param1, param2) {

    nhg1_logFC <- df1$logfc[df1$param == param1]
    nhg2_logFC <- df2$logfc[df2$param == param2]
    delta_logFC <- matrix(NA, nrow = length(nhg1_logFC), ncol = length(nhg2_logFC))
    groups1 <- df1$group[df1$param == param1]
    groups2 <- df2$group[df2$param == param2]
    rownames(delta_logFC) <- names(nhg1_logFC) <- groups1[!is.na(groups1)]
    colnames(delta_logFC) <- names(nhg2_logFC) <- groups2[!is.na(groups2)]
    for (i in names(nhg1_logFC)) {
        delta_logFC[i,] <- sapply(names(nhg2_logFC), function(j) abs(nhg1_logFC[i]-nhg2_logFC[j]))
    }
    
    return(delta_logFC)
}


mean_abs_logFC_nhood_groups <- function(nhg1, nhg2, param1, param2) {

    nhg1_logFC <- logFC_nhood_groups(nhg1, param1)
    nhg2_logFC <- logFC_nhood_groups(nhg2, param2)
    mean_logFC <- matrix(NA, nrow = length(nhg1_logFC), ncol = length(nhg2_logFC))
    groups1 <- unique(nhg1[[param1]])
    groups2 <- unique(nhg2[[param2]])
    rownames(mean_logFC) <- groups1[!is.na(groups1)]
    colnames(mean_logFC) <- groups2[!is.na(groups2)]
    for (i in names(nhg1_logFC)) { 
        mean_logFC[i,] <- sapply(names(nhg2_logFC), function(j) abs(mean(c(nhg1_logFC[i],nhg2_logFC[j]))))
    }
    
    return(mean_logFC)
}


logFC_similarity_nhood_groups <- function(nhg1, nhg2, param1, param2) {

    nhg1_logFC <- logFC_nhood_groups(nhg1, param1)
    nhg2_logFC <- logFC_nhood_groups(nhg2, param2)
    delta_logFC <- matrix(NA, nrow = length(nhg1_logFC), ncol = length(nhg2_logFC))
    groups1 <- unique(nhg1[[param1]])
    groups2 <- unique(nhg2[[param2]])
    rownames(delta_logFC) <- groups1[!is.na(groups1)]
    colnames(delta_logFC) <- groups2[!is.na(groups2)]
    for (i in names(nhg1_logFC)) {
        delta_logFC[i,] <- sapply(names(nhg2_logFC), function(j) abs(nhg1_logFC[i]-nhg2_logFC[j]))
    }
    
    return(apply(delta_logFC, 2, function(x) 1-x/max(delta_logFC)))
}
    

logFC_nhood_groups <- function(nhg, param) {

    groups <- unique(nhg[[param]])
    groups <- groups[!is.na(groups)]
    nhg_logFC <- c()
    for (i in groups) {
        val <- mean(nhg$logFC[nhg[[param]] == i], na.rm = TRUE)
        nhg_logFC <- c(nhg_logFC, val)
    }
    # OR: nhg_logFC <- unlist(lapply(groups, function(i) mean(nhg$logFC[nhg[[param]] == i], na.rm = TRUE)))
    names(nhg_logFC) <- groups
    
    return(nhg_logFC)
}


# compute the AUC separately for every group, using all markers for that group
mean_markers_auc <- function(obj1, obj2, param1, param2, gene.list.1, gene.list.2) {

    auc1 <- markers_auc_multi(object = obj1, param = param1, gene.list = gene.list.1)
    auc2 <- markers_auc_multi(object = obj2, param = param2, gene.list = gene.list.2)
    mean_auc <- matrix(NA, nrow = length(gene.list.1), ncol = length(gene.list.2))
    rownames(mean_auc) <- rownames(auc1)
    colnames(mean_auc) <- colnames(auc1)
    for (i in 1:nrow(mean_auc)) { 
        mean_auc[i,] <- sapply(auc2, function(x) mean(c(auc1[i],x)))
    }
    
    return(mean_auc)
}


# compute the AUC on markers shared between two groups
mean_shared_markers_auc <- function(obj1, obj2, param1, param2, gene.list.1, gene.list.2) {

    mean_auc <- matrix(NA, nrow = length(gene.list.1), ncol = length(gene.list.2))
    rownames(mean_auc) <- names(gene.list.1)
    colnames(mean_auc) <- names(gene.list.2)
    for (i in names(gene.list.1)) { 
        for (j in names(gene.list.2)) {
            shared_markers <- intersect(gene.list.1[[i]], gene.list.2[[j]])
            if (length(shared_markers) > 0) {
                auc1 <- markers_auc(object = obj1, meta.field = param1, group = i, gene.list = shared_markers)
                auc2 <- markers_auc(object = obj2, meta.field = param2, group = j, gene.list = shared_markers)
                mean_auc[i,j] <- mean(c(auc1, auc2))
            }
        }
    }
    
    return(mean_auc)
}


mean_shared_markers_auc_parallel <- function(obj1, obj2, param1, param2, gene.list.1, gene.list.2) {

    mean_auc <- mclapply(names(gene.list.1), function(i) { 
        unlist(lapply(names(gene.list.2), function(j) {
            shared_markers <- intersect(gene.list.1[[i]], gene.list.2[[j]])
            if (length(shared_markers) > 0) {
                # this calculation is the bottleneck; use multiple cores whenever possible
                auc1 <- markers_auc(object = obj1, meta.field = param1, group = i, gene.list = shared_markers)
                auc2 <- markers_auc(object = obj2, meta.field = param2, group = j, gene.list = shared_markers)
                mean(c(auc1, auc2))
            } else {
                NA
            }
        }))
    })
    mean_auc <- t(matrix(unlist(mean_auc), ncol = length(gene.list.1), nrow = length(gene.list.2)))
    rownames(mean_auc) <- names(gene.list.1)
    colnames(mean_auc) <- names(gene.list.2)
    
    return(mean_auc)
}


shared_marker_matrix <- function(gene.list.1, gene.list.2) {

    shared_markers <- unlist(lapply(names(gene.list.2), function(j) { 
        unlist(lapply(names(gene.list.1), function(i) {
            paste(intersect(gene.list.1[[i]], gene.list.2[[j]]), collapse = ",")
        }))
    }))
    shared_markers <- matrix(shared_markers, nrow = length(gene.list.1), ncol = length(gene.list.2))
    rownames(shared_markers) <- names(gene.list.1)
    colnames(shared_markers) <- names(gene.list.2)

    return(shared_markers)
}


# param is the identifier of the parameters used to build the nhood group
# gene.list is a named list of vectors, each vector corresponding to the marker for a single group
# the name of each vector is the same as the group ID in the meta.data 
# TODO: test this!!!
markers_auc_multi <- function(object, meta.field, gene.list) {

    # the object should be annotated with the nhood groups
    if (!(meta.field %in% colnames(object@meta.data))) {
        stop(paste(meta.field, "field was not found in object@meta.data\n"))
    }
    groups <- names(gene.list)
    if (sum(groups %in% unique(object@meta.data$param)) < length(groups)) {
        stop("Groups in gene.list and object@meta.data$param are not the same\n")
    }
    object <- AddModuleScore(object = object, features = gene.list, assay = "RNA")
    # compute the AUC for each group
    auc_markers <- rep(NA, length(groups))
    names(auc_markers) <- groups
    for (i in 1:length(groups)) {
        meta.name <- paste0("Cluster",i)
        pred <- object@meta.data[[meta.name]]
        res <- object@meta.data$param == groups[i]
        res[is.na(res)] <- rep(FALSE, sum(is.na(res)))
        if (length(unique(res)) == 2) {
            auc_markers[i] <- suppressMessages(auc(response = res, predictor = pred, direction = "<"))
        }
    }
    
    return(auc_markers)
}


# param is the identifier of the parameters used to build the nhood group
# gene.list is a vector
markers_auc <- function(object, meta.field, group, gene.list) {

    # the object should be annotated with the nhood groups
    if (!(meta.field %in% colnames(object@meta.data))) {
        stop(paste(meta.field, "field was not found in object@meta.data\n"))
    }
    object <- AddModuleScore(object = object, features = list(gene.list), assay = "RNA")
    # compute the AUC
    meta.name <- "Cluster1"
    pred <- object@meta.data[[meta.name]]
    res <- object@meta.data[[meta.field]] == group
    res[is.na(res)] <- rep(FALSE, sum(is.na(res)))
    if (length(unique(res)) == 2) {
        auc_markers <- suppressMessages(auc(response = res, predictor = pred, direction = "<"))
    } else {
        auc_markers <- NA
    }
    
    return(auc_markers)
}



# return a named list of vectors
# each vector contains the names of the markers
# names are nhood groups
# TODO: test this!!!
extract_nhood_group_markers <- function(df, param, LOGFC_MARKER = 1) {

    df_sub <- df[df$NhoodGroupParams == param,]
    groups <- as.character(unique(df_sub$NhoodGroup))
    l <- list()
    for (i in groups) {
        v <- df_sub$gene[df_sub$NhoodGroup == i & df_sub$avg_log2FC >= LOGFC_MARKER]
        l[[i]] <- v
    }
    
    return(l)
}



