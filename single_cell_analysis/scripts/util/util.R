# Qiwen Hu 2018
# core functions

kmeans_eval <- function(feature, celltype, iter = 50, seed = 1234){
 # This function is used to performame iterative k-means clustering based on projected 
 # features for single cell data and then evaluate the performance according to 
 # Nomalized mutual information (NMI) and adjusted rand index (ARI)
 #
 # Args:
 #  feature: a data.frame contains projected features (n dimensional space, n = 2, 3 ...), the columns are
 #    projected features in n dimensional space for a sample (cell), rows are samples
 #  celltype: vector contains cell type informantion
 #  iter: number of interation for k-means clustering
 #  seed: seed for k-means clustering
 # Returns:
 #   mean and stard deviation of NMI and ARI
  
  set.seed(seed)
  
  # remove the first column that contain sample ID
  if("id" %in% names(feature)){
    feature <- feature[, -1]
  } 
  
  # convert celltype into numbers
  celltype_num <- as.numeric(celltype)
  sample_id <- seq(1:nrow(feature))
  
  # iterative k-means
  nmi_score_all <- c()
  ari_score_all <- c()
  all_cluster <- list()
  
  for(i in 1:iter){
    # set k equal to the number of celltypes in the dataset
    k <- length(unique(celltype))
    # perform k means clustering
    km <- kmeans(feature, k)
    # get predicted cluster
    feature$cluster <- km$cluster
    
    # true clusters
    orignal_data <- data.frame(sample_id, celltype_num)
    
    # predicted clusters
    cl_data <- data.frame(sample_id, feature$cluster)
    # cacluate NMI and ARI score
    nmi_score <- NMI::NMI(orignal_data, cl_data)$value
    ari_score <- mclust::adjustedRandIndex(feature$cluster, celltype_num)
    
    nmi_score_all <- c(nmi_score_all, nmi_score)
    ari_score_all <- c(ari_score_all, ari_score)
  }
  
  return(data.frame(nmi = round(mean(nmi_score_all), 2), ari = round(mean(ari_score_all), 2), 
                    nmi_var = round(sd(nmi_score_all), 2), ari_var = round(sd(ari_score_all), 2)))
}

cal_performance <- function(pred_class, true_class, num_class_labels){
  # Performance evaluation for classification
  # Args:
  #  pred_class: predicted class label
  #  true_class: true class label
  #  num_class_labels: number of class labels in the training set
  # Returns:
  #  evaluation metrics (accruacy, precision, recall and f1)
  
  confusion <- as.matrix(table(true_class, pred_class, deparse.level = 0))
  
  if(num_class_labels > 2){
    n <- sum(confusion) # number of instances
    nc <- nrow(confusion) # number of classes
    diag <- diag(confusion) # number of correctly classified instances per class
    accuracy <- sum(diag)/sum(confusion)
    rowsums <- apply(confusion, 1, sum) # number of instances per class
    colsums <- apply(confusion, 2, sum) # number of predictions per class
    p <- rowsums / n # distribution of instances over the actual classes
    q <- colsums / n # distribution of instances over the predicted classes
    precision <- diag / colsums
    recall <- diag / rowsums
    f1 <- 2 * precision * recall / (precision + recall)
    macroPrecision <- mean(precision, na.rm = T)
    macroRecall <- mean(recall, na.rm = T)
    macroF1 <- mean(f1, na.rm = T)
    return(data.frame(accuracy = accuracy, macroPrecision = macroPrecision,
                      macrorecall = macroRecall, macrof1 = macroF1))
  } else {
    accuracy <- (confusion[1,1] + confusion[2,2]) / sum(confusion)
    precision <- confusion[1,1] / (confusion[1,1] + confusion[1,2])
    recall <- confusion[1,1] / (confusion[1,1] + confusion[2,1])
    f1 <- 2 * precision * recall / (precision + recall)
    return(data.frame(accuracy = accuracy, Precision = precision,
                      recall = recall, f1 = f1))
  }
}

knn_eval <- function(feature, celltype, k = 5, seed = 1234){
  # This function performs knn based evaluation 
  # Args:
  #  feature: a data.frame contains projected features (n dimensional space, n = 2, 3 ...), the columns are
  #    projected features in n dimensional space for a sample (cell), rows are samples
  #  celltype: vector contains cell type informantion
  #  k: cross validation fold
  # Returns:
  #  performance metrics and confusion matrix
  
  cv <- cvTools::cvFolds(nrow(feature), K = k, R = 1)
  feature$celltype <- celltype
  
  # remove id in feature
  if("id" %in% names(feature)){
    feature <- feature[, -1]
  } 
  
  perf.eval <- list()
  confusion.matrix <- 0
    
  for(i in 1:k){
    train <- feature[cv$subset[-which(cv$which == i)], ]
    test <- feature[cv$subset[which(cv$which == i)], ]
    knn_fit <- caret::train(celltype ~., data = train, method = "knn",
                     trControl = trainControl(method = "cv", number = 3),
                     preProcess = c("center", "scale"),
                     tuneLength = 10)
    knn_pred <- predict(knn_fit, newdata = subset(test, select = -c(celltype)))
    perf.eval[[i]] <- round(cal_performance(knn_pred, test$celltype, option = 3), 2)
    matrix <- as.matrix(table(test$celltype, knn_pred, deparse.level = 0))
    confusion.matrix <- matrix + confusion.matrix
  }
  
  # get mean performance of cross validation
  perf.eval <- dplyr::bind_rows(perf.eval)
  perf.mean <- as.data.frame(t(colMeans(perf.eval)))
  acc.sd <- sd(perf.eval$accuracy)
  perf.mean$acc.sd <- acc.sd
  
  return(list(performance = perf.mean, confusion = confusion.matrix))
}

knn_eval_cell <- function(feature, celltype, k = 5, seed = 1234){
  # This function performs knn based evaluation - validation is based on holdout cells
  # Args:
  #  feature: a data.frame contains projected features (n dimensional space, n = 2, 3 ...), the columns are
  #    projected features in n dimensional space for a sample (cell), rows are samples
  #  celltype: vector contains cell type informantion
  #  k: cross validation fold
  # Returns:
  #  performance metrics and confusion matrix
  
  # get cell information
  cell.id <- feature$id
  cell.id <- gsub("[.].*", "", cell.id)

  feature$celltype <- celltype
  
  #assign cell group information to feature table
  feature$group <- cell.id
  cell.id <- unique(cell.id)
  
  # cross validation - holdout for cells
  cv <- cvTools::cvFolds(length(cell.id), K = k, R = 1)
  
  perf.eval <- list()
  confusion.matrix <- 0
  
  for(i in 1:k){
    # get cell ids for training set
    cell.id.training <- cell.id[cv$subset[-which(cv$which == i)]]
    # get cell ids for testing set
    cell.id.testing <- cell.id[cv$subset[which(cv$which == i)]]
    
    train <- feature[feature$group %in% cell.id.training, c(2:4)]
    test <- feature[feature$group %in% cell.id.testing, c(2:4)]
    
    knn_fit <- caret::train(celltype ~., data = train, method = "knn",
                            trControl = trainControl(method = "cv", number = 3),
                            preProcess = c("center", "scale"),
                            tuneLength = 10)
    knn_pred <- predict(knn_fit, newdata = subset(test, select = -c(celltype)))
    perf.eval[[i]] <- round(cal_performance(knn_pred, test$celltype, option = 3), 2)
    matrix <- as.matrix(table(test$celltype, knn_pred, deparse.level = 0))
    confusion.matrix <- matrix + confusion.matrix
  }
  
  # get mean performance of cross validation
  perf.eval <- dplyr::bind_rows(perf.eval)
  perf.mean <- as.data.frame(t(colMeans(perf.eval)))
  acc.sd <- sd(perf.eval$accuracy)
  perf.mean$acc.sd <- acc.sd
  
  return(list(performance = perf.mean, confusion = confusion.matrix))
}


dim_reduc <- function(data, option = "SIMLR", seed = 12345){
  # Perform SIMLR or PCA dimension reduction 
  # Args:
  #   data: matrix used to do dimension reduction
  #   option: SIMLR/PCA: SIMLR is an statistical approach to find subpopulation of cells 
  #           based on similarity. For detail description of the algorithm, see
  #           https://www.bioconductor.org/packages/release/bioc/html/SIMLR.html
  #
  # Returns:
  #   features projected into two dimensional spaces
  #
  set.seed(seed)
  
  if(option == "SIMLR"){
    simlr_result <- SIMLR::SIMLR(X = data, c = 20, cores.ratio = 0)
    return(simlr_result$ydata)
  } 
  
  else if(option == "PCA"){
    gene_exp_pca <- as.data.frame(t(data))
    gene.pca <- prcomp(gene_exp_pca, center = TRUE, scale. = FALSE) 
    pca_result <- data.frame(gene.pca$x[, 1], gene.pca$x[, 2])
    return(pca_result)
  }
}

ave_sil <- function(feature, celltype){
  # This function computes the average silhouette score for projected features 
  # Args:
  #  feature: a data frame contains projected features (generated by an algorithm, e.g. PCA/ZIFA/VAE),
  #           and their correpodent cluster (columns), rows are samples.
  #  celltype: a vector contains cell type information for each sample (cell) 
  # Returns:
  #  average silhouette score
  
  # convert celltype label into numeric variable
  cluster <- as.numeric(celltype)
  
  # remove the first column that contain sample ID
  if("id" %in% names(feature)){
    feature <- feature[, -1]
  } 
  
  # computer silhouette value for each point in a specific cluster
  ss <- cluster::silhouette(cluster, dist(feature))
  
  # average silhouette value
  ss.mean <- mean(ss[, 3])
  return(round(ss.mean, 2))
}

coranking_metric <- function(count.matrix, dm.matrix){
  # Compute coranking metric to evaluate the performance of dimension reduction.
  # In general, coranking matrix is based on the distance matrix of high- and low-
  #  dimensional feature matrix. For each vetor i in the feature matrix, the rank of vector 
  #  j with respect to i is computed based on distance. It compares the rank of each element
  #  in the high dimensional space with the rank of the same element in the low dimensional space.
  # for detail definition of coranking matrix, see 
  # Lee, J.A., Lee, J.A., Verleysen, M., 2009. Quality assessment of dimensionality reduction:
  #  Rank-based criteria. Neurocomputing 72
  #
  # Args:
  #   count.matrix: orignal count.matrix (normalized) before dimension reduction
  #   dm.matrix: low dimension matrix after dimension reduction
  # Returns:
  #   q.local and q.global
  
  coranking.matrix <- coRanking::coranking(t(count.matrix), dm.matrix)
  lcmc <- numeric(nrow(coranking.matrix))
  lcmc <- coRanking::LCMC(coranking.matrix)
  
  # get the position of max lcmc value
  Kmax <- which.max(lcmc)
  
  #q.local is the mean over the values left of the maximum
  q.local <- round(mean(lcmc[1:(Kmax-1)]), 2)
 
  #q.global is the mean over the values right of the maximum
  q.global <- round(mean(lcmc[(Kmax+1):length(lcmc)]), 2)
  coranking.eval <- c(q.local, q.global)
  
  names(coranking.eval) <- c("Q.local", "Q.global")
  return(coranking.eval)
}

get_param <- function(filename, option = "simulated"){
  # This function is used to extract parameters/dataset information from simulated/real single cell datasets
  # 
  # Args:
  #   filename: file name
  #   option: simulated/real - simulated or real datasets
  # Returns
  #   parameters/dataset information
  
  if(option == "simulated"){
    parameter <- unlist(strsplit(filename, split = "[.]"))
    return(paste(parameter[2], parameter[3], sep = "."))
  } else if(option == "real"){
    parameter <- unlist(strsplit(tybalt_d1_files[i], split = "[.]"))[1]
  }
  return(parameter)
}
