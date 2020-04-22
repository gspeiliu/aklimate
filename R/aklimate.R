## (c) Vlado Uzunangelov 2017 ## uzunangelov@gmail.com

## dat is samples by features data frame where columns might be of different type (convert to data
## table later on?)  dat.grp - a list of the suffixes used to differentiate dat columns - can be more
## than one together - they should be distinct from one another so that making them regexes does not
## cause unwanted behavior i.e. c('cnv','cnv_gistic') is not OK, but c('MUTA:HOT','MUTA:NONSENSE') is

## always add - vector of dat column names that are to be included with each fset
#' @export
aklimate <- function(dat, dat.grp, lbls, fsets, always.add = NULL, rf.pars = list(), akl_pars = list(),
    store.kernels = FALSE, verbose = FALSE) {

    rf.pars <- rf.pars.default(rf.pars)

    ## for sorage purposes only
    rf.pars$always.add.vars <- always.add

    akl_pars <- aklimate_pars_default(akl_pars)


    ##########################################################

    switch(rf.pars$ttype, multiclass = , binary = {
        lbls <- data.frame(labels = factor(lbls))
    }, regression = {
        lbls <- data.frame(labels = lbls)

    }, error("Trying a method that is not implemented yet!"))


    idx.train <- rownames(lbls)

    #######################################

    rf.out <- train.forest.stats(dat, dat.grp, fsets, lbls, rf.pars, NULL, "_", verbose)



    ## overall
    idx <- rownames(rf.out$predictions.match)[sort(unique(unlist(lapply(1:ncol(rf.out$predictions.match),
        function(x) head(which(rf.out$predictions.match[, x]), n = akl_pars$topn)))))]

    ## ##multiclass extensions
    if (rf.pars$ttype == "multiclass") {
        lvls <- levels(lbls[, 1])

        lpm <- foreach(j = 1:nrow(rf.out$predictions), .combine = rbind) %docomb% {
            confM <- caret::confusionMatrix(factor(rf.out$predictions[j, ], levels = levels(lbls[, 1])),
                lbls[, 1])

            unname(confM$byClass[, "Balanced Accuracy"])

        }
        rownames(lpm) <- rownames(rf.out$predictions)


        mult <- foreach(i = 1:length(lvls)) %do% {
            oo <- order(lpm[, i], rf.out$msums, decreasing = TRUE)

            clvl <- which(lbls[, 1] == lvls[i])
            oopick <- sort(unique(unlist(lapply(1:length(clvl), function(x) head(which(rf.out$predictions.match[oo,
                clvl[x]]), n = akl_pars$topn/2)))))

            unique(c(rownames(lpm)[oo[oopick]], idx[which(lpm[idx, i] > quantile(lpm[, i], 0.95))]))

        }
        names(mult) <- lvls
        mult <- c(list(combined = idx), mult)

        um <- unique(unlist(mult))
        krel <- lapply(um, function(x) names(mult)[sapply(mult, function(y) x %in% y)])
        names(krel) <- um

        idx <- sort(match(um, rownames(rf.out$predictions)))
    }

    ############
    idx <- rf.out$pars.local[idx, , drop = FALSE]

    ###########################################################

    rf.models <- train.forest.kernels(dat, dat.grp, fsets, lbls, idx, rf.pars, always.add, "_", verbose)
    if (rf.pars$ttype == "multiclass") {
        for (k in 1:length(krel)) {
            rf.models[[names(krel)[k]]]$multic.rel <- krel[[k]]
        }
    }



    ## guarding against some forests that have predominantly empty trees
    medianl <- sapply(1:length(rf.models), function(x) median(sapply(rf.models[[x]]$forest$child.nodeIDs,
        function(x) length(x[[1]]))))
    rf.models <- rf.models[medianl >= 3]

    if (is.null(akl_pars$c.elnet)) {
        ## you have to subset dat and idx.train, otherwise it will construct the test kernels!!
        k.out <- forest.to.kernel.oob(rf.models, dat[idx.train, , drop = FALSE], dat.grp, fsets, always.add,
            idx.train, "_", verbose)



        mkl.pars <- cv_grid(nkern = dim(k.out)[3], len = akl_pars$cvlen, lam.b = akl_pars$lamb)
        if (!akl_pars$subsetCV)
            mkl.pars$nkern <- dim(k.out)[3]

        kcv <- kernel.cv(k.out, lbls[idx.train, 1], mkl.pars, akl_pars$nfold, rf.pars$ttype, rf.pars$metric)

        if (akl_pars$subsetCV) {
            if (rf.pars$ttype == "multiclass") {
                ll <- sapply(rf.models, function(x) length(x$multic.rel))
                sel <- 1:which(cumsum(ll) > kcv$pars[kcv$best.id, "nkern"])[1]
            } else {
                sel <- 1:kcv$pars[kcv$best.id, "nkern"]
            }

        } else {
            sel <- 1:length(rf.models)
        }


    } else {

        sel <- 1:length(rf.models)

    }


    k.out <- forest.to.kernel(rf.models[sel], dat[idx.train, , drop = FALSE], dat.grp, fsets, always.add,
        idx.train, "_", verbose)


    if (is.null(akl_pars$c.elnet)) {
        akl_model <- SPICER::spicer(k.out[idx.train, idx.train, 1:kcv$pars[kcv$best.id, "nkern"], drop = FALSE],
            lbls[idx.train, 1], C = c(kcv$pars[kcv$best.id, "lam1"], kcv$pars[kcv$best.id, "lam2"]),
            opt = list(regname = "elasticnet", display = 1))

    } else {
        akl_model <- SPICER::spicer(k.out[idx.train, idx.train, , drop = FALSE], lbls[idx.train, 1], C = akl_pars$c.elnet,
            opt = list(regname = "elasticnet", display = 1))
    }


    preds <- predict(akl_model, k.out[idx.train, idx.train, , drop = FALSE], type = akl_pars$type)

    res <- list(rf.stats = rf.out, kernels = if (store.kernels) k.out else NULL, kernel.cv = if (is.null(akl_pars$c.elnet)) kcv else NULL,
        rf.models = rf.models, akl_model = akl_model, fsets = rownames(idx), rf.pars.global = rf.pars,
        always.add = always.add, rf.pars.local = idx, akl_pars = akl_pars, dat.grp = dat.grp, idx.train = idx.train,
        preds.train = preds)

    class(res) <- c("aklimate", class(res))

    return(res)

}  ## end of aklimate

############################################################

## dat - same input as aklimate
#' @export
predict.aklimate <- function(akl_obj, dat, fsets, kernels = NULL, store.kernels = FALSE) {

    idx.train <- akl_obj$idx.train
    ## check all idx.train indices are in the data provided
    stopifnot(length(idx.train) == length(intersect(rownames(dat), idx.train)))

    idx.test <- setdiff(rownames(dat), idx.train)
    ## if no idx.test , then stop
    if (length(idx.test) == 0)
        stop("No test samples to make predictions on!")



    if (is.null(kernels)) {
        active <- if (akl_obj$rf.pars.global$ttype == "multiclass")
            unique(unlist(lapply(akl_obj$akl_model, function(x) names(x$sorted_kern_weight)))) else names(akl_obj$akl_model$sorted_kern_weight)

        active <- names(akl_obj$rf.models)[unname(sapply(gsub("\\+", "\\\\+", names(akl_obj$rf.models)),
            function(x) sum(grepl(x, active, fixed = TRUE)) > 0))]

        kernels <- forest.to.kernel(akl_obj$rf.models[active], dat, akl_obj$dat.grp, fsets, akl_obj$always.add,
            idx.train, "_", TRUE)

    }

    out <- predict(akl_obj$akl_model, kernels[idx.train, idx.test, , drop = FALSE], type = akl_obj$akl_pars$type)

    return(list(preds = out, kernels = if (store.kernels) kernels else NULL))
}

