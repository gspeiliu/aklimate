## (c) Vlado Uzunangelov 2017 ## uzunangelov@gmail.com

#' AKLIMATE : Algorithm for Kernel Learning with Approximating Tree Ensembles
#' @param dat - samples x features data frame where columns might be of different type   #' @param dat_grp - a list of vectors, each consisting of suffixes for data types that match the ones used in dat. Each vector corresponds to a particular combination of data types that will be tested for each component RF. Only the data type combination with the
#' best performance for a given feature set is retained.
#' The data type suffixes should be distinct from one another so that none is a proper substring of another - i.e. c('cnv','cnv_gistic') is not OK, but c('MUTA:HOT','MUTA:NONSENSE') is
#' @param fsets - list of prior knowledge feature sets
#' @param lbls - vector of training data labels
#' @param always_add - vector of dat column names that are to be included with each fset
#' @return a model of class "aklimate" with the following fields:
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


    idx <- rownames(rf_out$predictions_match)[sort(unique(unlist(lapply(1:ncol(rf_out$predictions_match),function(x) head(which(rf_out$predictions_match[, x]), n = akl_pars$topn)))))]
    ## overall
    switch(rf_pars$ttype,
           binary={
             probs<-rf_out$probabilities
             probs[rf_out$predictions_match]<-NA
             idx_1 <- unique(unlist(lapply(1:ncol(probs),function(x) rownames(probs)[order(probs[, x],decreasing=TRUE)[1:akl_pars$topn]])))
             idx<-c(idx,idx_1)
             idx<-rownames(probs)[rownames(probs)%in%idx]

           },
           multiclass={

             ####################
             lvls <- levels(lbls[, 1])
             probs<-rf_out$probabilities
             probs[!rf_out$predictions_match]<-NA
             # lpm <- foreach(j = 1:nrow(rf_out$predictions), .combine = rbind) %docomb% {
             #   confM <- caret::confusionMatrix(factor(rf_out$predictions[j, ], levels = levels(lbls[, 1])),
             #                                   lbls[, 1])
             #
             #   unname(confM$byClass[, "Balanced Accuracy"])
             #
             # }
             # rownames(lpm) <- rownames(rf_out$predictions)

              #######################################
             mult <- foreach(i = 1:length(lvls)) %do% {
               clvl <- which(lbls[, 1] == lvls[i])
               po<-rowSums(probs[,clvl],na.rm=TRUE)
               oo<-order(po,decreasing=TRUE)
               oopick <- unique(c(rownames(probs)[oo[1:akl_pars$topn]],
                     idx[which(po[idx] > quantile(po, 0.95,na.rm=TRUE))]))

             }
             names(mult) <- lvls

             #################

             mult <- c(list(combined = idx), mult)

             um <- unique(unlist(mult))
             krel <- lapply(um, function(x) names(mult)[sapply(mult, function(y) x %in% y)])
             names(krel) <- um

             idx <- rownames(rf_out$predictions)[rownames(rf_out$predictions)%in%um]

           },
           regression={})



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
#' @param akl_obj - an aklimate model
#' @param dat - samples x features data frame where columns might be of different type
#' @param fsets - list of prior knowledge feature sets
#' @param kernels - (Optional) pre-computed kernels for the MKL prediction part. By default the test kernels are computed on the fly from the RF models stored in the aklimate object.
#' @return A vector of predictions.
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

