## (c) Vlado Uzunangelov 2017 ## uzunangelov@gmail.com

#' AKLIMATE : Algorithm for Kernel Learning with Approximating Tree Ensembles
#' @title aklimate
#' @param dat  samples x features data frame where columns might be of different type
#' @param dat_grp  a list of vectors, each consisting of suffixes for data types that match the ones used in dat. Each vector corresponds to a particular combination of data types that will be tested for each component RF. Only the data type combination with the
#' best performance for a given feature set is retained.
#' The data type suffixes should be distinct from one another so that none is a proper substring of another - i.e. c('cnv','cnv_gistic') is not OK, but c('MUTA:HOT','MUTA:NONSENSE') is. This argument is considered experimental - we recommend supplying a list of length 1, with the list entry a vector of all possible suffixes.
#' @param fsets  list of prior knowledge feature sets
#' @param lbls  vector of training data labels
#' @param always_add  vector of dat column names that are to be included with each fset
#' @param rf_pars list of parameters for RF base kernels run
#' \describe{
#' \item{ntree}{Number of trees for RF kernel construction. Default is 1000.}
#' \item{min_node_prop}{Minimal size of leaf nodes (unit is proportion of training set size). Default is 0.01.}
#' \item{min_nfeat}{Minimal size of feature set (across all data modalities) for an RF to be constructed. Default is 15.}
#' \item{mtry_prop}{Proportion of features to be considered for each splitting decision. Default is 0.25}
#' \item{regression_q}{ For regression predictions only. Quantile of the per-sample empirical distribution of absolute differences between RF sample predictions and sample label. Used for binarization of sample predictions during best RF selection. Default 0.05.}
#' \item{replace}{TRUE/FALSE. Is subsampling to be done with replacement? Default is FALSE.}
#' \item{sample_frac}{Fraction of training data points to subsample for each tree. Default is 0.5 for sampling without replacement and 1 for bootstrapping.}
#' \item{ttype}{Type of learning task - choices are "binary","multiclass", and "regression". Default is "binary".}
#' \item{split_rule}{Type of splitting criteria- choices are "gini","hellinger","variance",and "beta". See ranger documentation for more details. Default is "gini".}
#' \item{importance}{Rule for calculating feature and feature set importance - choices are "impurity_corrected","permutation",and "impurity". Default is "impurity_corrected".}
#' \item{metric}{Metric for ranking RF base learner performance used in the selection of best RFs. Choices are "roc","pr","acc","bacc","mar","rmse","rsq","mae","pearson", and "spearman". Default is "roc".}
#' \item{unordered_factors}{How to treat unordered factors. Choices are "order","ignore", and "partition". See ranger for more details. Default is "order".}
#' \item{oob_cv}{A data frame of parameters to tune during trainings of all RF base learners, with OOB metric performance (from choices above) used to select best combination. Each row of the data frame includes a different combination of RF hyperparameters. The data frame has to contain at least two columns, with one column being "ntree". Having too many hyperparameter combinations can lead to significant slowdown in computation time. Default is a data frame of 1 row using the "min_node_prop","mtry_prop", and "ntree"/2 values of the rf_pars list. This argument is experimental - we recommend using the default setting.}
#' }
#' @param akl_pars list of parameters for RF best kernel selection and MKL meta-learner
#' \describe{
#' \item{topn}{number of RF kernels (ranked by metric specified in rf_pars) that correctly predict a given sample to be included in best RF list. Default is 5.}
#' \item{cvlen}{Number of random MKL hyperparameter combinations to be tested during MKL CV step. Default is 100.}
#' \item{nfold}{Number of folds to be used in MKL CV. Default is 5.}
#' \item{lamb}{Interval bounds from which random MKL hyperparameter combinations are drawn (log2 units). Default is (-20,0).}
#' \item{subsetCV}{ TRUE/FALSE. When TRUE, the MKL CV step also randomly varies the number of RF kernels in addition to the MKL regularization hyperparameters. It does so by training on a subset of kernels of size K, randomly selected on the (0,number best RF kernels) interval. Once K is selected, the top K kernels (ranked by metric specified in rf_pars) are included in current CV run. Default is TRUE.}
#' \item{type}{Type of predictions - possible choices are "response" and "probability". Default is "response".}
#' \item{celnet}{Hyperparameters for MKL elastic net run. Should be a vector of length 2. Default is NULL - hyperparameters are tuned via internal cross-validation.}
#' }
#' @param store_kernels TRUE/FALSE. Should the model store the training RF kernels. Default is FALSE.
#' @param verbose TRUE/FALSE. Should the model print verbose progress statements. Default is FALSE.
#' @return a model of class AKLIMATE with the following fields:
#' \describe{
#' \item{rf_stats}{List of metrics and predictions from training run on all RF base learners.}
#' \item{kernels}{ RF kernels used in MKL training step. NULL id store_kernels is set to FALSE. }
#' \item{kern_cv}{if akl_pars$celnet is NULL, heyperparameter vectors examined during MKL cross-validation, along with matching metric scores.}
#' \item{rf_models}{Set of RF base learners used to produce RF kernels for stacked MKL.}
#' \item{akl_model}{Trained spicer MKL model, with either user-supplied elastic net hyperparameters, or the hyperparameters selected via CV tuning.}
#' \item{rf_pars_global}{rf_pars argument}
#' \item{rf_pars_local}{optimal RF parameters for each RF base learner. Those will be the same (with the exception of ntree) as the rf_pars_global parameters except if rf_pars$oob_cv was specified by the user.}
#' \item{akl_pars}{akl_pars argument}
#' \item{dat_grp}{dat_grp argument}
#' \item{idx_train}{Vector of training data instances.}
#' \item{preds_train}{AKLIMATE predictions on training set.}
#' }
#' @references
#'   \itemize{
#'   \item  V. Uzunangelov, C. K. Wong, and J. Stuart. Highly Accurate Cancer Phenotype Prediction with AKLIMATE, a Stacked Kernel Learner Integrating Multimodal Genomic Data and Pathway Knowledge. bioRxiv, July 2020.
#'   }
#' @export
aklimate <- function(dat, dat_grp, lbls, fsets, always_add = NULL, rf_pars = list(), akl_pars = list(),
    store_kernels = FALSE, verbose = FALSE) {

    rf_pars <- rf_pars_default(rf_pars)

    ## for sorage purposes only
    rf_pars$always_add <- always_add

    akl_pars <- aklimate_pars_default(akl_pars)


    ##########################################################

    switch(rf_pars$ttype, multiclass = , binary = {
        lbls <- data.frame(labels = factor(lbls))
    }, regression = {
        lbls <- data.frame(labels = lbls)

    }, error("Trying a method that is not implemented yet!"))


    idx_train <- rownames(lbls)

    #######################################

    rf_out <- train_forest_stats(dat, dat_grp, fsets, lbls, rf_pars, NULL, "_", verbose)

    ## overall
    switch(rf_pars$ttype,
           regression=,
           binary={

             idx <- rownames(rf_out$predictions_match)[sort(unique(unlist(lapply(1:ncol(rf_out$predictions_match),function(x) head(which(rf_out$predictions_match[, x]), n = akl_pars$topn)))))]


           },
           multiclass={

             idx <- rownames(rf_out$predictions_match)[sort(unique(unlist(lapply(1:ncol(rf_out$predictions_match),function(x) head(which(rf_out$predictions_match[, x]), n = akl_pars$topn)))))]

             ####################
             lvls <- levels(lbls[, 1])
             probs<-rf_out$probabilities
             probs[!rf_out$predictions_match]<- 0


              #######################################
             mult <- foreach(i = 1:length(lvls)) %do% {
               clvl <- which(lbls[, 1] == lvls[i])
               po<-rowSums(probs[,clvl],na.rm=TRUE) -
                 sapply(1:nrow(probs),function(x) {
                   sum(probs[x,lbls[,1]!=lvls[i] &
                               rf_out$predictions[x,]==lvls[i] &
                               !rf_out$predictions_match[x,]])
                   })
               oo<-order(po,decreasing=TRUE)
               oopick <- unique(c(rownames(probs)[oo[1:ceiling(akl_pars$topn/2)]],
                     idx))

             }
             names(mult) <- lvls

             #################

             mult <- c(list(combined = idx), mult)

             um <- unique(unlist(mult))
             krel <- lapply(um, function(x) names(mult)[sapply(mult, function(y) x %in% y)])
             names(krel) <- um

             idx <- rownames(rf_out$predictions)[rownames(rf_out$predictions)%in%um]

           })



    ############
    idx <- rf_out$pars_local[idx, , drop = FALSE]

    ###########################################################

    rf_models <- train_forest_kernels(dat, dat_grp, fsets, lbls, idx, rf_pars, always_add, "_", verbose)
    if (rf_pars$ttype == "multiclass") {
        for (k in 1:length(krel)) {
            rf_models[[names(krel)[k]]]$multic_rel <- krel[[k]]
        }
    }



    ## guarding against some forests that have predominantly empty trees
    rf_models <- rf_models[sapply(rf_models,function(x) is_valid_forest(x,minl=3))]

    if (is.null(akl_pars$celnet)) {
        ## you have to subset dat and idx_train, otherwise it will construct the test kernels!!
        kout <- forest_to_kernel_oob(rf_models, dat[idx_train, , drop = FALSE], dat_grp, fsets, always_add,
            idx_train, "_", verbose)



        mkl_pars <- cv_grid(nkern = dim(kout)[3], len = akl_pars$cvlen, lam_b = akl_pars$lamb)
        if (!akl_pars$subsetCV)
            mkl_pars$nkern <- dim(kout)[3]

        kcv <- kernel_cv(kout, lbls[idx_train, 1], mkl_pars, akl_pars$nfold, rf_pars$ttype, rf_pars$metric)

        if (akl_pars$subsetCV) {
            if (rf_pars$ttype == "multiclass") {
                ll <- sapply(rf_models, function(x) length(x$multic_rel))
                sel <- 1:which(cumsum(ll) >= kcv$pars[kcv$best_id, "nkern"])[1]
            } else {
                sel <- 1:kcv$pars[kcv$best_id, "nkern"]
            }

        } else {
            sel <- 1:length(rf_models)
        }


    } else {

        sel <- 1:length(rf_models)

    }


    kout <- forest_to_kernel(rf_models[sel], dat[idx_train, , drop = FALSE], dat_grp, fsets, always_add,
        idx_train, "_", verbose)


    if (is.null(akl_pars$celnet)) {
        akl_model <- SPICER::spicer(kout[idx_train, idx_train, 1:kcv$pars[kcv$best_id, "nkern"], drop = FALSE],
            lbls[idx_train, 1], C = c(kcv$pars[kcv$best_id, "lam1"], kcv$pars[kcv$best_id, "lam2"]),
            opt = list(regname = "elasticnet", display = 1))

    } else {
        akl_model <- SPICER::spicer(kout[idx_train, idx_train, , drop = FALSE], lbls[idx_train, 1], C = akl_pars$celnet,
            opt = list(regname = "elasticnet", display = 1))
    }


    preds <- predict(akl_model, kout[idx_train, idx_train, , drop = FALSE], type = akl_pars$type)

    res <- list(rf_stats = rf_out, kernels = if (store_kernels) kout else NULL, kern_cv = if (is.null(akl_pars$celnet)) kcv else NULL,
        rf_models = rf_models, akl_model = akl_model, rf_pars_global = rf_pars,
        always_add = always_add, rf_pars_local = idx, akl_pars = akl_pars, dat_grp = dat_grp, idx_train = idx_train,
        preds_train = preds)

    class(res) <- c("aklimate", class(res))

    return(res)

}  ## end of aklimate

############################################################

#' Compute predictions from an aklimate model
#' @title predict.aklimate
#' @param akl_obj  an AKLIMATE model
#' @param dat  samples x features data frame where columns might be of different type
#' @param fsets  list of prior knowledge feature sets
#' @param kernels  Pre-computed kernels for the MKL prediction part. By default the value is NULL, i.e. test kernels are computed on the fly from the RF models stored in the aklimate object.
#' @param store_kernels Boolean variable indicating whether test_kernels should be stored in output. Default is FALSE.
#' @return List with the following elements:
#' \describe{
#' \item{preds}{Predictions. \cr
#' If akl_obj$akl_pars$type is "response", this is a vector class predictions or regression estimates. \cr
#' If akl_obj$akl_pars$type is "probability" (classification tasks only), it is a matrix with predicted probabilities of membership for each class.}
#' \item{kernels}{If store_kernels is TRUE, the test kernels. Default is FALSE.}
#' }
#' @export
predict.aklimate <- function(akl_obj, dat, fsets, kernels = NULL, store_kernels = FALSE) {

    idx_train <- akl_obj$idx_train
    ## check all idx_train indices are in the data provided
    stopifnot(length(idx_train) == length(intersect(rownames(dat), idx_train)))

    idx_test <- setdiff(rownames(dat), idx_train)
    ## if no idx_test , then stop
    if (length(idx_test) == 0)
        stop("No test samples to make predictions on!")



    if (is.null(kernels)) {
        active <- if (akl_obj$rf_pars_global$ttype == "multiclass")
            unique(unlist(lapply(akl_obj$akl_model, function(x) names(x$sorted_kern_weight)))) else names(akl_obj$akl_model$sorted_kern_weight)

        active <- names(akl_obj$rf_models)[unname(sapply(gsub("\\+", "\\\\+", names(akl_obj$rf_models)),
            function(x) sum(grepl(x, active, fixed = TRUE)) > 0))]

        kernels <- forest_to_kernel(akl_obj$rf_models[active], dat, akl_obj$dat_grp, fsets, akl_obj$always_add,
            idx_train, "_", TRUE)

    }

    out <- predict(akl_obj$akl_model, kernels[idx_train, idx_test, , drop = FALSE], type = akl_obj$akl_pars$type)

    return(list(preds = out, kernels = if (store_kernels) kernels else NULL))
}

