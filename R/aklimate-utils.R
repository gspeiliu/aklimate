# (c) Vlado Uzunangelov 2017                                                             ## uzunangelov@gmail.com

#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom iterators iter


##reads a file in a listt format,
##e.g. Name1 member1 member2 ...
##     Name2 member1 member2 ...
##returns a list with each line representing a separate vector object named by first entry on line and having as members subsequent entries
read_set_list <- function(file, delim="\t") {
  l <- readLines(file)
  lflds <- strsplit(l, delim)
  r <- lapply(lflds, function(x) as.vector(x[-1]))
  names(r) <- sapply(lflds, "[[", 1)
  return(r)
}


######################################################
## reverse operation of read_set_list

write_set_list <- function(setList,outfile,delim="\t"){
    #since we are appending to a file,
    #we need to check if it exists and remove it
    if(outfile!=stdout() && file.exists(outfile)){
        file.remove(outfile)
    }
    for (i in 1:length(setList)){
        write(paste(names(setList)[i],paste(setList[[i]],collapse=delim),sep=delim),file=outfile,append=TRUE)
    }

}


##Computes a statistic for each pair of entries in two lists
##list1 - can be a data.frame (list of columns), will represent rows in output matrix
##list2 - can be a data.frame , will represent columns in output matrix
##compfunc - a function whose arguments are source for the (full) list and target for the individual entry (vector) being compared against each element of the list
##output is a matrix of size length(list1)xlength(list2) with each column representing
##the statistics computed by compfunc(list1,list2[[x]],...)
get_cross_table <- function(list1,list2,compfunc,parallel=FALSE,...) {
    if (parallel) {
        cores <- if(foreach::getDoParRegistered()) foreach::getDoParWorkers() else parallel::detectCores()-2
        out <- parallel::mcmapply(compfunc,target=list2,MoreArgs=c(list(source=list1),list(...)),mc.cores=cores)
    } else {
        out <- mapply(compfunc,target=list2,MoreArgs=c(list(source=list1),list(...)))
    }
    return(out)
}

##source is a list with the rows of the output matrix
##target is a vector that will produce one column of the output matrix
sample_sim_vect <- function(source,target) {
    sapply(source,function(x)  {
        inds<-(!is.na(x))&(!is.na(target))
        sum(x[inds]==target[inds])/sum(inds)
    })
}

sample_sim_vect_mask <- function(source,target,mask) {
    mapply(function(x,y,target) sum(x[y]==target[y])/sum(y),source,mask,MoreArgs = list(target=target) )
}


sample_sim_vect_eucl <- function(source,target) {
    sapply(source,function(x) {
        inds<-(!is.na(x))&(!is.na(target))
        sqrt(sum((x[inds]-target[inds])^2))
    })
}

sample_sim_vect_eucl_mask <- function(source,target,mask) {
    mapply(function(x,y,target) sqrt(sum((x[y]-target[y])^2)),source,mask,MoreArgs = list(target=target))
}

#'Relative contribution of data types based on cumulative importance of features
#' @param suffs - vector of suffixes of participating data types
#' @param ranked - vactor of ranked features produced by rank_features()
#' @return A vector of controbutions for each data type, normalized to a sum of 1.
#' @export
rank_importance_type <- function(suffs,ranked,intron="_") {

    ranked_pos <- ranked[ranked>0]
    res <- sapply(suffs,function(x) {

        ## ##This breaks down if you have suffixes that are a substring of another suffix         ## at beginning - i.e. (cnv, cnv_gistic)
        ##use gistic_cnv, cnv instead

        sum(ranked_pos[grepl(paste0(intron,x),names(ranked_pos),perl=TRUE)])
    })

    res <- sort(res,decreasing=TRUE)
    res <- res/sum(ranked_pos)
    return(res)
}


expand_names <- function(names,suffs,intron="_"){
    if(is.null(suffs)) return(names)
    res <- apply(expand.grid(names,suffs),1,function(x) paste(x,collapse=intron))
    return(res)
}

##cleans up names by replacing dots with underscores
clean_names <- function(names) {
    ##make names converts "-" to ".", so end result is "-" -> "_"!! Be very careful!!
    res <- gsub('[.]', '_',make.names(names,unique=TRUE))
    ##res <- gsub('\\.{1}$','',res)
    return(res)
}

#' Combine entries in a vector that match a particular pattern
#' @param patterns - vector of patterns to match
#' @param weights - named numeric vector that needs to be collapsed based on the patterns
#' @return A named vector of scores for each pattern, sorted from largest to smallest.
#' @export
collapse_weights<-function(patterns, weights){

  res<-rep(0,length(patterns))
  names(res)<-patterns
  for(i in patterns){
    res[i]<-sum(weights[grepl(i,names(weights))])
  }

  res<-sort(res[res>0],decreasing = TRUE)

  if (sum(res)!=sum(weights)) warning("There were unmatched or overlapping patterns in your weight aggregation call!")

  return(res)
}

##Produces a combined data frame of properly named features from all data types used
##dat - a list of data.frames sample x feature that will be collapsed into one sample x feature df
## the names of dat will become the suffixes of the features in a given data type ( to differentiate entries with the same name from different data types)
## output is a sample x feature data frame

collapse_data <- function(dat,intron="_"){

    if(length(dat)>1){
        same_names <- sapply(2:length(dat),function(x) identical(rownames(dat[[x]]),
                                                             rownames(dat[[x-1]])))

        stopifnot(sum(same_names)!=(length(same_names)-1))
    }
    dat_adj <- lapply(names(dat),function(x) {
        res <- dat[[x]]
        colnames(res) <- paste(colnames(res),x,sep=intron)
        res
    })

    out <- do.call(cbind,dat_adj)
    out <- data.frame(out)
    colnames(out) <- clean_names(colnames(out))
    return(out)
}


create_grp_suff <- function(dat_grp,intron="_") {
    paste0(intron,sapply(dat_grp,paste,collapse=intron))
    ##paste0(intron,paste0(dat_grp,collapse=intron))
}


split_set_suff <- function(combined,dat_grp,...){
    grp_suff <- create_grp_suff(dat_grp,...)

    fset_ind <- strsplit(combined,paste(grp_suff,collapse="|"))[[1]][1]
    feat_ind <- which(grp_suff==strsplit(combined,fset_ind,fixed=TRUE)[[1]][2])
    feat_suff <- if(length(feat_ind)>0) dat_grp[[feat_ind]] else NULL

    return(list(fset_ind=fset_ind,feat_suff=feat_suff))
}

extract_features<-function(allfeats,infeats,ids){
    allfeats[grepl(paste(infeats,collapse="|"),allfeats,ignore.case = TRUE)&
             grepl(paste(ids,collapse="|"),allfeats,ignore.case = TRUE)]
}


##the rows of metrics nees to be ranked from most to least desirable
which_keep <- function(metrics,comm_names,suffs){
    ind_keep <- c()
    for(i in comm_names){

        idx <- which(rownames(metrics)%in%paste0(i,suffs))
        if (length(idx)>0) ind_keep <- c(ind_keep,idx[1])
    }

    return(rownames(metrics)[sort(ind_keep)])
}

cv_grid <- function(nkern=500,len=250,lam_b=c(-5,5)){
        lam1 <- 2^runif(len, lam_b[1], lam_b[2])
        ##lam1 <- c(runif(len-2),0,1)
        lam2 <- 2^runif(len, lam_b[1], lam_b[2])

        nkern <- sample(1:nkern, size = len, replace = TRUE)
        newGrid <- data.frame(lam1=lam1,lam2=lam2,nkern=nkern)
        newGrid<-newGrid[order(newGrid$lam1/newGrid$lam2),]
        return(newGrid)
}


##pars - a data.frame with columns named lam1, lam2 & nkern giving the parameters to be tried with MKL
kernel_cv <- function(kerns,lbls,pars=cv_grid(nkern=dim(kerns)[3]),nfolds=5,type="binary",measure="bacc") {

    ############################
    kfolds <- caret::createFolds(if(type=="regression") lbls else factor(lbls),
                          k=nfolds,
                          list=TRUE,
                          returnTrain = TRUE)

    .loss <- switch(type,
                    multiclass=,
                    binary="logit",
                    regression="square")

    res <- foreach(i=1:nrow(pars),.options.multicore=list(preschedule=FALSE))%docomb%{

        simn<-10
        if(i<2*simn){
            Sys.sleep(10*floor(i/simn))
        }

        res1 <- foreach(j=1:length(kfolds))%do%{
            mod <- tryCatch({
                SPICER::spicer(kerns[kfolds[[j]],kfolds[[j]],1:pars[i,"nkern"] ,drop=FALSE],
                           lbls[kfolds[[j]]],
                           C=c(pars[i,"lam1"],pars[i,"lam2"]),
                       opt=list(loss=.loss,regname="elasticnet",
                           display=1))
            },error=function(err){
                return(NULL)
            })
            .rtype <- switch(type,
                             multiclass=,
                             binary="probability",
                             regression="response")

            switch(type,
                   binary=,
                   regression={
                       if(is.null(mod) || length(mod$sorted_kern_weight)==0) {
                           return(c(metric=NA,nused=NA))
                       }
                   },
                   multiclass={
                       if( is.null(mod) ||
                          any(sapply(mod,is.null)) ||
                          any(sapply(mod,function(x) length(x$sorted_kern_weight))==0) ) {
                           return(c(metric=NA,nused=NA))
                       }
                   })

            preds <- predict(mod,
                             kerns[kfolds[[j]],-kfolds[[j]],
                                   1:pars[i,"nkern"],drop=FALSE],
                             type=.rtype)

            if(any(is.na(preds))) {
                metric <- NA
            } else {
                ##big metrics are considered good, small metrics are inferior

                metric <- select_metric(preds,lbls[-kfolds[[j]]],type,measure)

            }

            ##for multiclass problems, kernel weights are stored in an attribute
            nused <- max(length(mod$sorted_kern_weight),
                         length(attributes(mod)$sorted_kern_weight))
            rm(mod,preds)
            return(c(metric=metric,
                   nused=nused))
        }
        res1 <-do.call(rbind,res1)
        c(colMeans(res1),sd=sd(res1[,"metric"])/sqrt(nfolds))

    }

    res <- do.call(rbind,res)
    res <- data.frame(res)
    cc <- complete.cases(res)
    res <- res[cc,]
    pars <- pars[cc,]

    #####################
    o <- order(res$metric,decreasing=TRUE)
    res <- res[o,]
    pars <- pars[o,]

    pars_id <- round(0.1*nrow(res))
    return(list(res=res,pars=pars,best_id=pars_id))
}



is_valid_forest <- function(rf,minl=3) {
   medianl <- median(sapply(rf$forest$child.nodeIDs, function(x) length(x[[1]])))
   return(medianl>=minl)

}

correct_regression <- function(preds,lbls) {

      sq_err_f <- function(par) {
                    sum((lbls- par*preds )^2)
                }

       oo <- tryCatch({
                    optim(1,sq_err_f,method="L-BFGS-B",lower=0)
                }, error=function(cond) { return(list(par=0,value=Inf)) } )

      if(is.infinite(oo$value)) {
                    return (preds)
                } else {
                    return(oo$par*preds)
                }

}
## the binary/multiclass predictions need to be of type probability (i.e. dim=2)
select_metric <- function(preds,lbls,ttype,measure="roc") {

    metric <- switch(ttype,
           binary={

               ##auroc
               rocr_pred <- ROCR::prediction(preds[,2],lbls)

               confM <- caret::confusionMatrix(
                   factor(apply(preds,1,function(x) colnames(preds)[which.max(x)]),
                          levels=levels(lbls)),lbls)

               msr <- switch(measure,
                             roc={
                                 ROCR::performance(rocr_pred,"auc")@y.values[[1]]
                                },
                             pr={
                                 perf <- ROCR::performance(rocr_pred,"prec","rec")
                                 perf@y.values[[1]][1] <- 1
                                 pracma::trapz(perf@x.values[[1]],perf@y.values[[1]])
                             },
                             acc={
                                 unname(confM$overall["Accuracy"])

                                  },
                             bacc={
                                 unname(confM$byClass["Balanced Accuracy"])
                             },
                             mar={
                                           unname(mean(
                                               confM$byClass[,"Sensitivity"]))
                                       })

               msr


           },
                     regression={
                                msr <- switch(measure,
                                       rmse={
                                           1/sqrt(sum((lbls-preds)^2)/length(preds))
                                       },
                                       mae={
                                           length(preds)/sum(abs(lbls-preds))
                                       },
                                       pearson={
                                           cor(lbls,preds,method="pearson")
                                       },
                                       spearman={
                                           cor(lbls,preds,method="spearman")
                                       },
                                       rsq={
                                           1- (sum((lbls-preds)^2)/sum((lbls-mean(lbls))^2))
                                       })
                         msr

           },
                     multiclass={
                         confM <- caret::confusionMatrix(
                             factor(apply(preds,1,function(x) colnames(preds)[which.max(x)]),
                                    levels=levels(lbls)),lbls)
                         msr <- switch(measure,
                                       acc={
                                           unname(confM$overall["Accuracy"])
                                       },
                                       bacc={
                                           unname(mean(
                                               confM$byClass[,"Balanced Accuracy"]))
                                       },
                                       mar={
                                           unname(mean(
                                               confM$byClass[,"Sensitivity"]))
                                       })
                         msr
                     },
           error("Trying a method that is not implemented yet!"))

    return(metric)
}


compose_features <- function(dat,suff="composite",deg=3,topn=ncol(dat),idx=rownames(dat)) {

     comp_list <- foreach(i=iter(2:deg))%do%{
         combos <- combn(colnames(dat),i)
         out <- apply(combos,2,function(x) Emcdf::emcdf(dat[idx,x,drop=FALSE],dat[idx,x,drop=FALSE]))
         ##out <- apply(combos,2,function(x) copula::C.n(pobs(dat[idx,x,drop=FALSE]),dat[idx,x,drop=FALSE],ties.method="average"))
         rownames(out) <- idx
         colnames(out) <- paste0(apply(combos,2,function(x) paste0(x,collapse="_")),
                                    "_",suff)
         out
     }

     res <- do.call(cbind,comp_list)

     idx_test <- setdiff(rownames(dat),idx)
     res1 <- NULL

     if(length(idx_test)>0) {

     comp_list <- foreach(i=iter(2:deg))%do%{
         combos <- combn(colnames(dat),i)
         out <- apply(combos,2,function(x) Emcdf::emcdf(dat[idx,x,drop=FALSE],dat[idx_test,x,drop=FALSE]))
         ##out <- apply(combos,2,function(x) copula::C.n(pobs(dat[idx_test,x,drop=FALSE]),dat[idx,x,drop=FALSE],ties.method="average"))
         rownames(out) <- idx_test
         colnames(out) <- paste0(apply(combos,2,function(x) paste0(x,collapse="_")),
                                    "_",suff)
         out
     }
     res1 <- do.call(cbind,comp_list)
     }


     ##cutoff <- quantile(apply(dat[idx,],2,sd),0.5)
     sds <- apply(res,2,sd)
     sds <- sds[order(sds,decreasing=TRUE)]
     res <- rbind(res,res1)
     ##return(res[,names(sds)[sds>cutoff]])
     return(res[,names(sds)[1:min(topn,ncol(res))]])
 }

##nbr - number of bins
##idx - indices of training data
##if idx == rownames(dat) (default) then all data is considered simultaneously for binning
## if idx!=rownames(dat) the training set is idx, which is what you use to construct the bins,
##the test set is rownames(dat)-idx, which uses the bins constructed on idx to bin data

quantize_data <- function(dat,nbr=5,idx=rownames(dat)) {

    idx_test <- setdiff(rownames(dat),idx)
    res <- foreach(i=1:ncol(dat),.combine=cbind)%docomb%{
        breaks <- unique(quantile(dat[idx,i],seq(0,1,1/nbr)))
        breaks[length(breaks)] <- Inf
        breaks[1] <- -Inf
        as.numeric(cut(dat[idx,i],breaks=breaks))
        }

    rownames(res) <- idx

    if(length(idx_test)>0){

        res1 <- foreach(i=1:ncol(dat),.combine=cbind)%docomb%{
            breaks <- unique(quantile(dat[idx,i],seq(0,1,1/nbr)))
            breaks[length(breaks)] <- Inf
            breaks[1] <- -Inf
            as.numeric(cut(dat[idx_test,i],breaks=breaks))
        }
        rownames(res1) <- idx_test
        res <- rbind(res,res1)
    }


    colnames(res) <- colnames(dat)
    res <- res[rownames(dat),]

    return(res)
}


##default parameters for RF selected individually for each feature set

rf_pars_default <- function(rf_pars=list()) {

    ##rf params defaults
    if(is.null(rf_pars$ntree)) rf_pars$ntree <- 1000
    if(is.null(rf_pars$min_node_prop)) rf_pars$min_node_prop <- 0.01
    if(is.null(rf_pars$min_nfeat)) rf_pars$min_nfeat <- 15
    if(is.null(rf_pars$mtry_prop)) rf_pars$mtry_prop <- 0.25
    if(is.null(rf_pars$regression_q)) rf_pars$regression_q <- 0.05

    rf_pars$replace <- if(is.null(rf_pars$replace)) FALSE else as.logical(match.arg(as.character(rf_pars$replace),c("TRUE","FALSE")))

    if(is.null(rf_pars$sample_frac)) rf_pars$sample_frac <- ifelse(rf_pars$replace,1,0.5)



    rf_pars$ttype <- if(is.null(rf_pars$ttype)) "binary" else match.arg(rf_pars$ttype,c("binary","regression","multiclass"))

    rf_pars$split_rule <- if(is.null(rf_pars$split_rule)) "gini" else match.arg(rf_pars$split_rule,c("gini","hellinger","variance","beta"))

    rf_pars$importance <- if(is.null(rf_pars$importance)) "impurity_corrected" else match.arg(rf_pars$importance,c("impurity_corrected","permutation","impurity"))

    rf_pars$metric <- if(is.null(rf_pars$metric)) "roc" else match.arg(rf_pars$metric,c("roc","pr","acc","bacc","mar","rmse","rsq","mae","pearson","spearman"))

    rf_pars$unordered_factors <- if(is.null(rf_pars$unordered_factors)) "order" else match.arg(rf_pars$unordered_factors,c("order","ignore","partition"))

    ##paramaters for oob_cv
    ##the list should include ntree as a minimum, since the CV part will be run on forests with fewer trees
    ## for speed of computation
    ## The list should include at least one more parameter, otherwise RR ranking is done with the parameters
    ## specified above (min_node_prop, min_nfeat,mtry_prop)
    if (is.null(rf_pars$oob_cv)) rf_pars$oob_cv <- data.frame(min_node_prop=rf_pars$min_node_prop,
                                                              mtry_prop=rf_pars$mtry_prop,
                                                              ntree=rf_pars$ntree/2,
                                                              check.names=FALSE,
                                                              stringsAsFactors = FALSE)

    if (!"ntree"%in%colnames(rf_pars$oob_cv)||ncol(rf_pars$oob_cv)<2) stop("Misspecified rf_pars$oob_cv parameters!")

    return(rf_pars)
}

aklimate_pars_default<-function(akl_pars=list()) {

  if(is.null(akl_pars$topn)) akl_pars$topn <- 5
  if(is.null(akl_pars$cvlen)) akl_pars$cvlen <- 100
  if(is.null(akl_pars$nfold)) akl_pars$nfold <- 5
  if(is.null(akl_pars$lamb)) akl_pars$lamb <- c(-20,0)


  akl_pars$subsetCV <- if(is.null(akl_pars$subsetCV)) TRUE else as.logical(match.arg(as.character(akl_pars$subsetCV),c("TRUE","FALSE")))

  if(is.null(akl_pars$type)) akl_pars$type <- "response" else match.arg(akl_pars$type,c("response","probability"))
  return(akl_pars)

}



##Helper function that trains RF for each feature set

##lbls -a data.frame of 1 column
##dat -a data.frame of predictors
##rf_pars - list of parameters for RF
##always_split - vector of variables that will be added to the ones considered for splitting for each tree
## returns a trained RF model
rf_train <- function(dat,lbls,rf_pars,always_split=NULL) {

    colnames(lbls) <- c("labels")

    idx_train <- rownames(lbls)

    ##dat <- mlr::createDummyFeatures(dat)

    dat <- cbind(lbls,dat[idx_train,,drop=FALSE])
    switch(rf_pars$ttype,
           multiclass=,
           binary={


               rf <- ranger::ranger(data=dat[idx_train,,drop=FALSE],
                            dependent.variable.name="labels",
                            splitrule = rf_pars$split_rule,
                            always.split.variables=always_split,
                            classification = TRUE,
                            sample.fraction = rf_pars$sample_frac,
                            num.trees=rf_pars$ntree,
                            mtry=ceiling((ncol(dat)-1)*rf_pars$mtry_prop),
                            min.node.size=ceiling(nrow(dat)*rf_pars$min_node_prop),
                            class.weights = rf_pars$class_weights,
                            num.threads=1,
                            probability=TRUE,
                            respect.unordered.factors = rf_pars$unordered_factors,
                            importance=rf_pars$importance,
                            write.forest=TRUE,
                            keep.inbag=TRUE,
                            replace=rf_pars$replace)

           },
           regression={


               rf <- ranger::ranger(data=dat[idx_train,,drop=FALSE],
                            dependent.variable.name="labels",
                            splitrule = rf_pars$split_rule,
                            always.split.variables=always_split,
                            sample.fraction = rf_pars$sample_frac,
                            num.trees=rf_pars$ntree,
                            mtry=ceiling((ncol(dat)-1)*rf_pars$mtry_prop),
                            min.node.size=ceiling(nrow(dat)*rf_pars$min_node_prop),
                            num.threads=1,
                            importance=rf_pars$importance,
                            write.forest=TRUE,
                            respect.unordered.factors = rf_pars$unordered_factors,
                            keep.inbag=TRUE,
                            replace=rf_pars$replace)
           },
           error("Trying a method that is not implemented yet!"))


    return(rf)
}

##Helper function for selection of most relevant RFs from all RFs trained
##dat - a samples x featres df with all features together
## dat_grp - list of vectors, each vector containing a group of data type suffixes to be considered together
##fsets - a list of feature sets
##lbls -a df of one column - a factor - for classification and  muti-class analysis
## (status column should be 0-1 with 1's corresponding to deaths)
##always_add - vector of dat columns that should be added to each fset before RF trainins

## returns a ranking feature sets based on RF out-of-bag prediction performance on trainin set
train_forest_stats <- function(dat,dat_grp,fsets,lbls,rf_pars_global=rf_pars_default(list()),always_add=NULL,sep="_",verbose=FALSE){

    rf_pars_global <- rf_pars_default(rf_pars_global)

        idx_train <- rownames(lbls)

    grp_suff <- create_grp_suff(dat_grp,intron=sep)
    ##train models


    if(verbose) {
        print(date())
        print("Starting feature set forests")
    }

    zz <- foreach(k=1:nrow(rf_pars_global$oob_cv))%do%{

        rf_pars_local <- rf_pars_global
        rf_pars_local[colnames(rf_pars_global$oob_cv)] <- rf_pars_global$oob_cv[k,]

        vv <- foreach(j=iter(dat_grp)) %do% {
            tt <- foreach(i=iter(fsets)) %docomb% {

                feats<-extract_features(colnames(dat),i,j)

                if(length(feats)<rf_pars_global$min_nfeat) return(list(metric=NA,imps=NA,preds=NA))

                curr_dat <- dat[idx_train,unique(c(feats,always_add)),drop=FALSE]

                rf <- rf_train(curr_dat,
                               lbls,
                               rf_pars_local,
                               NULL)

                imps <- ranger::importance(rf)
                imps <- imps[order(imps,decreasing=TRUE)]
                imps <- imps[imps>0]
                imps <- imps


                if(!is_valid_forest(rf,3)) return(list(metric=NA,imps=NA,preds=NA))
                metric <- select_metric(rf$predictions,lbls[,1],rf_pars_global$ttype,rf_pars_global$metric)


                ############

                if(length(dim(rf$predictions))==2) {
                  ##binary and multiclass
                    preds<-apply(rf$predictions,1,
                                 function(x) levels(lbls$labels)[which.max(x)])
                    probs<-apply(rf$predictions,1,max)
                    names(probs)<-rownames(lbls)
                } else {
                    preds<-rf$predictions
                    probs<-NULL
                }

                if(all(is.na(preds))) return(list(metric=NA,imps=NA,preds=NA))

                names(preds) <- rownames(lbls)
                return(list(metric=metric,imps=imps,preds=preds,probs=probs))

            }

            names(tt) <- paste0(names(fsets),paste0(sep,paste(j,collapse=sep)))
            tt <- tt[!sapply(tt,function(x) any(is.na(x)))]

            mms_in <- sapply(tt,function(x) x$metric)

            imps_in <- lapply(tt,function(x) x$imps)

            preds_in <- lapply(tt,function(x) x$preds)

            probs_in <- lapply(tt,function(x) x$probs)
            list(mms_in=mms_in,imps_in=imps_in,preds_in=preds_in,probs_in=probs_in)
        }

        mms <- unlist(lapply(vv,function(x) x$mms_in))

        imps <- do.call(c,lapply(vv, function(x) x$imps_in))

        preds <- do.call(c,lapply(vv, function(x) x$preds_in))

        probs <- do.call(c,lapply(vv, function(x) x$probs_in))
        list(mms=mms,imps=imps,preds=preds,probs=probs)
    }

    if(verbose) {
        print(date())
        print("Finished feature set forests")
    }

    metrics <- purrr::reduce(lapply(1:length(zz),
                             function(x) {

                                 out <- data.frame(id=names(zz[[x]]$mms),
                                               metric=zz[[x]]$mms,
                                               stringsAsFactors = FALSE)
                                 colnames(out)[2] <- paste0("metric_",x)
                                 out
                             }),merge,by="id")
    rownames(metrics) <- metrics[,1]
    metrics <- as.matrix(metrics[,-1,drop=FALSE])


    ################################

    ranks <- metrics


    #############################
    rankmax <- sapply(1:nrow(ranks),function(x) {
        o <- order(ranks[x,],decreasing=TRUE)
        o[1]
        ##o[max(1,floor(0.05*ncol(ranks)))]
    })
    names(rankmax) <- rownames(ranks)


    pars <- data.frame(t(sapply(rankmax,
                                function(x) {
                                    unlist(rf_pars_global$oob_cv[x,colnames(rf_pars_global$oob_cv),drop=FALSE])
                                })))

    pars <- pars[,setdiff(colnames(pars),"ntree"),drop=FALSE]

    imps <- lapply(1:length(rankmax),function(x) {
        zz[[rankmax[x]]]$imps[[names(rankmax)[x]]]
    })
    names(imps) <- names(rankmax)

    preds <- lapply(1:length(rankmax),function(x) {
        zz[[rankmax[x]]]$preds[[names(rankmax)[x]]]
    })
    names(preds) <- names(rankmax)

    probs <- lapply(1:length(rankmax),function(x) {
      zz[[rankmax[x]]]$probs[[names(rankmax)[x]]]
    })
    names(probs) <- names(rankmax)


    msums <- rowMeans(ranks,na.rm=TRUE)


    o <- order(msums,decreasing = TRUE)

    msums <- msums[o]
    metrics <- metrics[o,,drop=FALSE]
    pars <- pars[o,,drop=FALSE]
    ranks <- ranks[o,,drop=FALSE]
    imps <- imps[o]
    preds <- preds[o]
    probs<-probs[o]


    ##############################
    if(length(dat_grp)>1){
        ##keep only highest scoring RF corresponding to a given feature set
        kept <- which_keep(metrics,names(fsets),grp_suff)
        metrics <- metrics[kept,,drop=FALSE]
        msums <- msums[kept]
        pars <- pars[kept,,drop=FALSE]
        ranks <- ranks[kept,,drop=FALSE]
        imps <- imps[kept]
        preds <- preds[kept]
        probs<-probs[kept]
    }
    ##################################

    preds <- t(sapply(preds,function(x) x))

    ##needs to be extended for regression (refactored into a separate function)
    if(rf_pars_global$ttype=="regression"){
        preds_match <- t(apply(preds,1,function(x) x-lbls[idx_train,1]))
        preds_match <- apply(abs(preds_match),2,function(x) x<quantile(x,rf_pars_global$regression_q))
        probs<-NULL

    } else {
        preds_match <- t(apply(preds,1,function(x) x==levels(lbls[idx_train,1])[lbls[idx_train,1]]))
        probs<-t(sapply(probs,function(x) x))
    }
    return(list(pars_global=rf_pars_global,pars_local=pars,res=zz,stats=metrics,msums=msums,rankmax=rankmax,importance=imps,predictions=preds,probabilities=probs,predictions_match=preds_match))
}



##lbls - a data.frame of 1 or 2 columns
##rf_pars_local - the individual parameters for each forest (determined by oob error)- a df with colnames corresponding to parameters in rf_pars_default - rows correspond to each forest
train_forest_kernels <- function(dat,dat_grp,fsets,lbls,rf_pars_local,rf_pars_global=rf_pars_default(list()),always_add=NULL,sep="_",verbose=FALSE) {

    idx_train <- rownames(lbls)

    if(verbose) {
        print(date())
        print("Started training forest")
    }

    ii <- if(is.null(rf_pars_local)) names(fsets) else rownames(rf_pars_local)

    rf_out <- foreach(i=iter(ii),.options.multicore=list(preschedule=TRUE))%docomb%{

        ##fset_ind is the name of the set
        ##feat_suff is a vector of suffices that denominate data types
        SPICER::expand(split_set_suff(i,dat_grp),nms=c("fset_ind","feat_suff"))

        feats<-extract_features(colnames(dat),fsets[[fset_ind]],feat_suff)

        ##feats <- colnames(dat)[colnames(dat)%in%expand_names(fsets[[fset_ind]],feat_suff,sep)]

        curr_dat <- dat[idx_train,unique(c(feats,always_add)),drop=FALSE]

        if (!is.null(rf_pars_local)) rf_pars_global[colnames(rf_pars_local)] <- rf_pars_local[i,]
        rf <- rf_train(curr_dat,
                       lbls,
                       rf_pars_global,
                       NULL)

        return(rf)

    }
    names(rf_out) <- ii

    ##gc()
    if(verbose) {
        print(date())
        print("Finished training forests")
    }

    return(rf_out)
}

#####
##dat - same input as aklimate; by default idx_train are presumed to be all indices in rows of dat - if some of those are idx_test, specify accordingly
##dat rownames (samples) should be ordered accordingly

forest_to_kernel <- function(rf_models,dat,dat_grp,fsets,always_add=NULL,idx_train=rownames(dat),sep="_",verbose=FALSE){

    ##check all idx_train indices are in the data provided
    stopifnot(length(idx_train)==length(intersect(rownames(dat),idx_train)))


    idx_test <- setdiff(rownames(dat),idx_train)
    ##if no idx_test , then we ar econstructing trianing kernels
    if (length(idx_test)==0) idx_test <- idx_train

    if(verbose) {
        print(date())
        print("Started kernel construction")
    }

    medianl <- sapply(1:length(rf_models),function(x) median(sapply(rf_models[[x]]$forest$child.nodeIDs, function(x) length(x[[1]]))))
    names(medianl) <- names(rf_models)


    kerns <- foreach(i=iter(names(rf_models)),.options.multicore=list(preschedule=FALSE))%docomb%{

        simn<-10
        if(match(i,names(rf_models))<2*simn){
            Sys.sleep(10*floor(match(i,names(rf_models))/simn))
        }


        SPICER::expand(split_set_suff(i,dat_grp),nms=c("fset_ind","feat_suff"))

        feats<-extract_features(colnames(dat),fsets[[fset_ind]],feat_suff)


        curr_dat <- dat[,unique(c(feats,always_add)),drop=FALSE]
        ##curr_dat <- mlr::createDummyFeatures(curr_dat)

        q <- 1
        p <- 2

        if(medianl[i]>1){

            preds <- predict(rf_models[[i]],
                             curr_dat,
                             predict.all = TRUE)$predictions

            rownames(preds) <- rownames(curr_dat)

            multic<-length(dim(preds))>2 && dim(preds)[2]>2

            if(multic){
                extensions<-rf_models[[i]]$multic_rel
                counter<-which(rf_models[[i]]$forest$levels%in%rf_models[[i]]$multic_rel)
            } else {
                extensions<-"combined"
            }


            knames<-c(paste(i,extensions,sep=sep))
            ##just for names, we change "_combined" to ""
            ##this way non-multicalss weight names are not plastered with "_combined at end
            knames<-gsub(paste0(sep,"combined$"),"",knames)


            proximity<-replicate(length(extensions),matrix(1,length(idx_train),length(idx_test)))

            dimnames(proximity)[[1]] <- idx_train
            dimnames(proximity)[[2]] <- idx_test
            dimnames(proximity)[[3]] <- knames

            if(multic && length(counter)>0) {
                for(j in 1:length(counter)){

                    pA <- Similarity::distance(preds[idx_train,counter[j],],
                                               preds[idx_test,counter[j],],
                                               method="euclidian",p=p)

                    if(identical(idx_train,idx_test)){
                        sA <- quantile(pA,probs=q)
                    } else {
                        sA <- quantile(Similarity::distance(preds[idx_train,counter[j],],
                                                            preds[idx_train,counter[j],],
                                                            method="euclidian",p=p),
                                       probs=q)
                    }

                    ##first matrix is reserved for combined matrix
                    if("combined"%in%extensions) {
                        proximity[,,j+1]<-1-(pA/sA)
                    } else {
                        proximity[,,j]<-1-(pA/sA)
                    }

                    rm(pA)
                }
            }

            if("combined"%in%extensions) {

                if(length(dim(preds))>2) {
                    if(dim(preds)[2]==2){
                        ##binary classification task with probability estimates prediction -
                        ##returns a 3D array with second dimension of 2 (number of classes)
                        ##using Euclidean distance on the probabilities of either class
                        ##produces the same distance matrix
                        preds <- preds[,1,]
                    } else {
                        ##multiclass classification -
                        ##output is again 3D array but second dimension corresponds to more
                        ##classes and you can't just pick one to compute distance on
                        ##will have to concatenate them
                        preds <- t(apply(preds,1,c))
                    }
                }
                ##no need for else since preds have right format in regression case

             prox <- Similarity::distance(preds[idx_train,],
                               preds[idx_test,],
                               method="euclidian",p=p)

             if(identical(idx_train,idx_test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma <- quantile(Similarity::distance(preds[idx_train,],
                                            preds[idx_train,],
                                            method="euclidian",p=p),probs=q)
             }

             proximity[,,1]<-1-(prox/sigma)


             ###############################
             preds1 <- predict(rf_models[[i]],
                               curr_dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr_dat)

             prox <- Similarity::distance(preds1[idx_train,],preds1[idx_test,],method="euclidian",p=p)

             if(identical(idx_train,idx_test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma <- quantile(Similarity::distance(preds1[idx_train,],
                                  preds1[idx_train,],
                                  method="euclidian",p=p),probs=q)

             }

             proximity[,,1]<-proximity[,,1]*(1-(prox/sigma))

             #####################################################

             preds1 <- as.data.frame(t(preds1))

             prox  <-  get_cross_table(preds1[,idx_train],
                                       preds1[,idx_test],
                                       sample_sim_vect)

             ##Similarity package ranger interactions seems to be broken now
             ##prox <- Similarity::proximityMatrixRanger(curr_dat[idx_train,],
             ##                             curr_dat[idx_test,],
             ##                             rf_models[[i]])

             proximity[,,1]<-proximity[,,1]*prox
             proximity[,,1]<-proximity[,,1]^(1/3)

             for(p in 1:dim(proximity)[3]){
                 proximity[,,p][is.na(proximity[,,p])]<-0
             }

            }

        ############################################


        } else {
            warning("No kernel was computed because the forest has mostly empty trees!")
            proximity <- NULL
        }


        rm(curr_dat,preds,preds1,prox)
        ##gc()
        proximity
    }

    kerns <- abind::abind(kerns,along=3)

    if(verbose) {
        print(date())
        print("Finished kernel construction")
    }

    rm(rf_models)

    return(kerns)
}

forest_to_kernel_oob <- function(rf_models,dat,dat_grp,fsets,always_add=NULL,idx_train=rownames(dat),sep="_",verbose=FALSE){

    ##check all idx_train indices are in the data provided
    stopifnot(length(idx_train)==length(intersect(rownames(dat),idx_train)))


    idx_test <- setdiff(rownames(dat),idx_train)
    ##if no idx_test , then we ar econstructing trianing kernels
    if (length(idx_test)==0) idx_test <- idx_train

    if(verbose) {
        print(date())
        print("Started kernel construction")
    }

    medianl <- sapply(1:length(rf_models),function(x) median(sapply(rf_models[[x]]$forest$child.nodeIDs, function(x) length(x[[1]]))))
    names(medianl) <- names(rf_models)

    kerns <- foreach(i=iter(names(rf_models)),.options.multicore=list(preschedule=FALSE))%docomb%{

        simn<-10
        if(match(i,names(rf_models))<2*simn){
            Sys.sleep(10*floor(match(i,names(rf_models))/simn))
        }


        SPICER::expand(split_set_suff(i,dat_grp),nms=c("fset_ind","feat_suff"))

        feats<-extract_features(colnames(dat),fsets[[fset_ind]],feat_suff)

        curr_dat <- dat[,unique(c(feats,always_add)),drop=FALSE]
        ##curr_dat <- mlr::createDummyFeatures(curr_dat)

        q <- 1
        p <- 2

        if(medianl[i]>1){

            mask <- rf_models[[i]]$inbag.counts
            mask <- t(matrix(unlist(mask),ncol=length(mask[[1]]),byrow=TRUE))
            mask <- mask!=0


            preds <- predict(rf_models[[i]],
                             curr_dat,
                             predict.all = TRUE)$predictions

            rownames(preds) <- rownames(curr_dat)



            multic<-length(dim(preds))>2 && dim(preds)[2]>2

            if(multic){
                extensions<-rf_models[[i]]$multic_rel
                counter<-which(rf_models[[i]]$forest$levels%in%rf_models[[i]]$multic_rel)
            } else {
                extensions<-"combined"
            }

            knames<-c(paste(i,extensions,sep=sep))
            ##just for names, we change "_combined" to ""
            ##this way non-multicalss weight names are not plastered with "_combined at end
            knames<-gsub(paste0(sep,"combined$"),"",knames)

            proximity<-replicate(length(extensions),matrix(1,length(idx_train),length(idx_test)))

            dimnames(proximity)[[1]] <- idx_train
            dimnames(proximity)[[2]] <- idx_test
            dimnames(proximity)[[3]] <- knames

            if(multic && length(counter)>0) {
                for(j in 1:length(counter)){
                    pl <- preds[,counter[j],]

                    pl[mask]<-NA

                    pl<-as.data.frame(t(pl))

                    pA<-get_cross_table(pl[,idx_train],
                                             pl[,idx_test],
                                             sample_sim_vect_eucl)

                    if(identical(idx_train,idx_test)){
                        sA <- quantile(pA,probs=q)
                    } else {

                        sA<-quantile(get_cross_table(pl[,idx_train],
                                                 pl[,idx_train],
                                                 sample_sim_vect_eucl),
                                     probs=q)
                    }

                    ##first matrix is reserved for combined matrix
                    if("combined"%in%extensions) {
                        proximity[,,j+1]<-1-(pA/sA)
                    } else {
                        proximity[,,j]<-1-(pA/sA)
                    }
                    rm(pA)
                }

            }


            if("combined"%in%extensions) {

                if(length(dim(preds))>2) {
                    if(dim(preds)[2]==2){
                    ##binary classification task with probability estimates prediction -
                    ##returns a 3D array with second dimension of 2 (number of classes)
                    ##using Euclidean distance on the probabilities of either class
                    ##produces the same distance matrix
                        preds <- preds[,1,]

                        preds[mask] <- NA
                    } else {
                        rps <- dim(preds)[2]
                    ##multiclass classification -
                    ##output is again 3D array but second dimension corresponds to more
                    ##classes and you can't just pick one to compute distance on
                    ##will have to concatenate them
                        preds <- t(apply(preds,1,c))
                    ##need to replicate mask matrix to match
                        mask_extend <- matrix(data=apply(mask,2,rep,rps),
                                   ncol=ncol(mask)*rps)

                        preds[mask_extend]<-NA
                    }
                } else {

                    preds[mask]<-NA

                }

                preds<-as.data.frame(t(preds))

                prox<-get_cross_table(preds[,idx_train],
                                         preds[,idx_test],
                                         sample_sim_vect_eucl)

             if(identical(idx_train,idx_test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma<-quantile(get_cross_table(preds[,idx_train],
                                          preds[,idx_train],
                                          sample_sim_vect_eucl),
                                 probs=q)
             }


                proximity[,,1]<-1-(prox/sigma)


             ###############################
             preds1 <- predict(rf_models[[i]],
                               curr_dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr_dat)


                preds1[mask] <- NA

                preds1<-as.data.frame(t(preds1))

                prox<-get_cross_table(preds1[,idx_train],
                                           preds1[,idx_test],
                                           sample_sim_vect_eucl)


             if(identical(idx_train,idx_test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma<-quantile(get_cross_table(preds1[,idx_train],
                                                      preds1[,idx_train],
                                                      sample_sim_vect_eucl),
                                 probs=q)

             }


             proximity[,,1]<-proximity[,,1]*(1-(prox/sigma))

             #####################################################


             preds1 <- predict(rf_models[[i]],
                               curr_dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr_dat)
                preds1[mask] <- NA

              preds1 <- as.data.frame(t(preds1))

                prox  <-  get_cross_table(preds1[,idx_train],
                                               preds1[,idx_test],
                                               sample_sim_vect)


             proximity[,,1]<-proximity[,,1]*prox
             proximity[,,1]<-proximity[,,1]^(1/3)

                for(p in 1:dim(proximity)[3]){
                    proximity[,,p][is.na(proximity[,,p])]<-0
                }
            }


        } else {
            warning("No kernel was computed because the forest has mostly empty trees!")
            proximity <- NULL
        }


        rm(curr_dat,preds,preds1,prox,mask)
        proximity
    }

    kerns <- abind::abind(kerns,along=3)

    if(verbose) {
        print(date())
        print("Finished kernel construction")
    }

    rm(rf_models)

    return(kerns)
}

#' Feature importance calculation from an aklimate model
#' @param akl_obj - an aklimate model
#' @return A ranked vector of all features with non-zero importance.
#' @export
rank_features <- function(akl_obj) {

    if(akl_obj$rf_pars_global$ttype=="multiclass"){
        imps_list <- foreach(k=iter(akl_obj$akl_model))%do%{
                    wghts <- k$sorted_kern_weight
                    imps <- foreach(i=iter(names(wghts)))%do%{
                        ##you need to anchor the patterns at the start of the string
                        ## in case there are pathway names that are exact substrings of
                        ##other pathway names
                        ##this will become a problem if the matching substring is
                        ##at the beginning of the longer name
                        ind<-which(stringr::str_detect(i,paste0("^",names(akl_obj$rf_models))))

                        res <- sort(ranger::importance(akl_obj$rf_models[[ind]]),decreasing=TRUE)
                        res <- res[res>0]
                        res
                    }
                    names(imps) <- names(wghts)

                    scaled_imps <- lapply(names(wghts), function(x) imps[[x]]*wghts[x]/sum(imps[[x]]))
                    names(scaled_imps) <- names(wghts)

                    all_names <- unique(unlist(lapply(scaled_imps,names)))
                    comb_imps <- rep(0,length(all_names))
                    names(comb_imps) <- all_names

                    for(i in names(scaled_imps)) {
                        comb_imps[names(scaled_imps[[i]])] <- comb_imps[names(scaled_imps[[i]])] + scaled_imps[[i]]
                    }
                    comb_imps <- sort(comb_imps,decreasing=TRUE)
        }

        full_names <- unique(names(unlist(imps_list)))
        out <- rep(0,length(full_names))
        names(out) <- full_names

        for(l in 1:length(imps_list)){
            out[names(imps_list[[l]])] <-  out[names(imps_list[[l]])] + imps_list[[l]]
        }
        out <- sort(out,decreasing=TRUE)/sum(out)


    } else {
        wghts <- akl_obj$akl_model$sorted_kern_weight
        imps <- foreach(i=iter(names(wghts)))%do%{
            ind<-which(stringr::str_detect(i,names(akl_obj$rf_models)))
            res <- sort(ranger::importance(akl_obj$rf_models[[ind]]),decreasing=TRUE)
            res <- res[res>0]
            res
        }
        names(imps) <- names(wghts)

        scaled_imps <- lapply(names(wghts), function(x) imps[[x]]*wghts[x]/sum(imps[[x]]))
        names(scaled_imps) <- names(wghts)

        all_names <- unique(unlist(lapply(scaled_imps,names)))
        comb_imps <- rep(0,length(all_names))
        names(comb_imps) <- all_names

        for(i in names(scaled_imps)) {
            comb_imps[names(scaled_imps[[i]])] <- comb_imps[names(scaled_imps[[i]])] + scaled_imps[[i]]
        }
        out <- sort(comb_imps,decreasing=TRUE)


    }

    return(out)
}

##dat - a data frame samples x features with all data types combined (can have factors)
##clean up names of pathways beforehand
##you should pre-filter set_wghts and feat_wghts to the # you want plotted!!
##the first entry in should_scale is the one whose column heatmap is used if cluster_col is TRUE
##the last column in lbls is the one that we use to group columns
plot_heatmap <- function(feat_wghts,dat_grp,dat,lbls,pdf_title,fset_wghts=NULL,fsets=NULL,should_scale=NULL,idx_ordered=1,fsize=6,cluster_rows=TRUE,cluster_cols=TRUE,hsplit=NULL){


    library(ComplexHeatmap)
    library(circlize)

    ##require(colorRamps)
    require(RColorBrewer)

    ##color pallettes
    qual_col_pals  <-  RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    colors <-  unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))


    cont_pals <- lapply(c("RdYlBu","RdBu","RdYlGn","PuOr","RdGy","PRGn","BrBG"),
                        function(x) RColorBrewer::brewer.pal(n = 11,name = x))

    seq_col_pals <- lapply(rownames(RColorBrewer::brewer.pal.info)[RColorBrewer::brewer.pal.info$category == 'seq'],
                        function(x) RColorBrewer::brewer.pal(n = 9,name = x))

    ##prep data
    factor_feats <- names(feat_wghts)[sapply(dat[,names(feat_wghts)],is.factor)]
    dat[,factor_feats] <- apply(dat[,factor_feats,drop=FALSE],2,function(x) as.numeric(addNA(x)))

    include_fsets <- !is.null(fsets) && length(fset_wghts) >1

    stopifnot(length(intersect(rownames(dat),rownames(lbls)))==nrow(lbls))

    dat <- dat[rownames(lbls),]
    ##annotation labels should be in a sample x annotation labels df
    oc <- do.call(order,lbls[,idx_ordered,drop=FALSE])
    lbls <- lbls[oc,,drop=FALSE]
    dat <- dat[oc,,drop=FALSE]

    ##Data Preprocessing
    prsd_feats <- t(sapply(names(feat_wghts),function(x) {
        split_set_suff(x,unique(unlist(dat_grp)))
    }))

    if(!is.null(fsets)){

        nfset <- gsub(paste0(paste0("_",sapply(dat_grp,
                                               function(x) paste0(x,collapse ="_"))),
                             collapse="|"),
                      "",names(fset_wghts))
        names(fset_wghts) <- nfset

        fsets_trimmed <- lapply(fsets[names(fset_wghts)],function(x) {
            rownames(prsd_feats)[prsd_feats[,"fset_ind"]%in%x]
        })


    fset_wghts <- fset_wghts[names(fsets_trimmed)]
    fset_annot <- list.to.matrix.membership(fsets_trimmed)

    ##if just showing one feature set, trim features to match it
    if(length(fsets_trimmed)==1) prsd_feats <- prsd_feats[fsets_trimmed[[1]],,drop=FALSE]

    }




    #######################
    uniq_dtypes <- unlist(unique(prsd_feats[,"feat_suff"]))
    uniq_feats <-   unlist(unique(prsd_feats[,"fset_ind"]))

    if(is.null(should_scale)) should_scale <- setNames(rep(FALSE,length(uniq_dtypes)),
                                                       uniq_dtypes)
    num_feats <- table(as.character(prsd_feats[,2]))
    num_feats <- num_feats[uniq_dtypes]

    ##add pseudocounts for small groups
    num_feats <- num_feats+4
    ##account for annotations
    num_feats[1] <- 2*num_feats[1]

    if(include_fsets) {
        ##this is for the features that did not fall in any of the top pathways
        empty_annot <- matrix(0,nrow=nrow(prsd_feats)-nrow(fset_annot),ncol=ncol(fset_annot))
        rownames(empty_annot) <- setdiff(rownames(prsd_feats),rownames(fset_annot))
        fset_annot <- rbind(fset_annot,empty_annot)
        o <- order(fset_wghts,decreasing=TRUE)
        fset_annot <- fset_annot[,names(fset_wghts)[o],drop=FALSE]

        fset_annot_col <- sapply(1:ncol(fset_annot),function(x) {
            out <- rep(0,nrow(fset_annot))
            out[fset_annot[,x]>0] <- x
            out
        })
        rownames(fset_annot_col) <- rownames(fset_annot)
        colnames(fset_annot_col) <- colnames(fset_annot)

        fset_annot <- fset_annot_col
        rm(fset_annot_col)
    }



    ##############################
    ##Heatmap column annotations

    ##color mapping
    cc <- 0
    cc_seq <- 1
    lbl_col <- list()

    for (x in colnames(lbls)){
        if(!is.factor(lbls[,x])) {
                    br <- unname(quantile(lbls[,x],probs=seq(0,1,0.12),na.rm=TRUE))
                    br[1] <- floor(br[1])
                    br[length(br)] <- ceiling(br[length(br)])

                    lbl_col <- c(lbl_col,list(cols=circlize::colorRamp2(br,seq_col_pals[[cc_seq]])))
                    cc_seq <- cc_seq+1

         } else  {
            lbl_col <-  c(lbl_col,list(cols=setNames(colors[cc+1:length(levels(lbls[,x]))], levels(lbls[,x]))))
             cc <- cc+length(levels(lbls[,x]))
         }
    }
    names(lbl_col) <- colnames(lbls)


    ##unit(rep(3,ncol(lbls)),"mm")
    ##since the annotation height is included in the top heatmap height,
    ##best to keep the units relative since this is how we determine the
    ## relative size of viewports on the layout
    ha_top <- HeatmapAnnotation(df=lbls,
                                   col=lbl_col,
                                   na_col = "whitesmoke",
                                   show_annotation_name = TRUE,
                                   annotation_name_offset = unit(1,"mm"),
                                   annotation_name_gp=gpar(fontsize=fsize*1.5),
                                   annotation_name_side = "left",
                                   height=unit(0.05*ncol(lbls),"npc"),
                                annotation_height=if(ncol(lbls)<5) unit(rep(0.05,ncol(lbls)),"npc") else unit(rep(0.5/ncol(lbls),ncol(lbls)),"npc"),
                                   show_legend=FALSE)
    legend_list_top <- lapply(colnames(lbls),function(x) color_mapping_legend(ha_top@anno_list[[x]]@color_mapping,nrow=NULL,ncol=2,plot = FALSE))

    names(legend_list_top) <- colnames(lbls)

    legend_list_bottom <- list()
    ############################

    ##Plotting and Display

    pdf(pdf_title,paper="a4",
        width=unit(8,"inches"),
        height=unit(12*sum(num_feats)/100,"inches"))
    grid.newpage()
    main_vp <- viewport(name="main_vp",
                        layout = grid.layout(nr = length(uniq_dtypes),
                              nc = 1,
                              height=c(num_feats/sum(num_feats),0.05)))
    pushViewport(main_vp)

    if(include_fsets) {
        col_bplot_annot <- HeatmapAnnotation(fset_weights=column_anno_barplot(fset_wghts[colnames(fset_annot)],
                                                 baseline=floor(0.7*min(fset_wghts)),
                                             bar_width=0.5,
                                                ylim=c(floor(0.7*min(fset_wghts)),
                                                    max(fset_wghts)),
                                                 border=FALSE,axis=FALSE,
                                                 axis_side="left",
                                             gp=gpar(fill="#CCCCCC",col="#CCCCCC")),
                                             which="column",
                                         annotation_height = unit(15,"mm"))
    }


    for(i in 1:length(uniq_dtypes)) {
        ind <- uniq_dtypes[i]
        curr_feat <- rownames(prsd_feats)[prsd_feats[,2]==ind]
        names(curr_feat) <- as.character(prsd_feats[prsd_feats[,2]==ind,1])

        if(should_scale[ind]) {
            dat_curr <- scale(dat[rownames(lbls),curr_feat,drop=FALSE])
        } else {
            dat_curr <- dat[rownames(lbls),curr_feat,drop=FALSE]
        }

            dat_curr <- t(dat_curr)
            rownames(dat_curr) <- names(curr_feat)


        ##annotation_width=unit(2,"cm"),
            row_bplot_annot <-  HeatmapAnnotation(feat_weights=row_anno_barplot(feat_wghts[curr_feat],
                                             baseline=0,
                                             size=unit(min(fsize/4,2),"mm"),
                                                  ylim=c(0,1.1*max(feat_wghts)),
                                             border=FALSE,axis=FALSE,axis_side="top",
                                             gp=gpar(fill="#CCCCCC",col="#CCCCCC")),
                                              which="row",
                                              width=unit(0.067,"npc"))
            row_bplot_annot@anno_list[[1]]@name <- paste0(ind,"_feat_weights")
            ############

        ## if(i==1) {

        ##     chclust <- group_hclust(1-cor(dat_curr,method="spearman"),lbls)
        ## }


        ##viridis::inferno(11)),
                                     ##column_dend_reorder=TRUE,
                                     ##show_column_dend=if(i==1) TRUE else FALSE,

        breaks <- unname(quantile(dat_curr,probs=seq(0,1,0.1)))
        breaks[1] <- floor(breaks[1])
        breaks[length(breaks)] <- ceiling(breaks[length(breaks)])

        hm1 <- Heatmap(dat_curr,
                                     col=colorRamp2(breaks=breaks,rev(cont_pals[[i]])),
                                     name=ind,
                                     clustering_distance_rows = "spearman",
                                     clustering_method_rows = "average",
                                     show_row_dend=cluster_rows,
                                     row_dend_reorder = cluster_rows,
                                     row_dend_side = "left",
                                     row_dend_width=unit(0.1,"npc"),
                                     cluster_rows=cluster_rows,
                                     width= 5,
                                     cluster_columns=FALSE,
                                     show_row_names=TRUE,
                                     row_names_gp=gpar(fontsize=fsize),
                                     row_names_max_width = unit(0.01,"npc"),
                                     show_column_names=FALSE,
                                     show_heatmap_legend=FALSE,
                                     heatmap_legend_param=list(title=ind,
                                         color_bar="continuous",
                                         legend_direction="horizontal",
                                         nrow=1,
                                         gp=list(fontsize=fsize/1.5),
                                         title_position="topleft",
                                         legend_height=unit(0.005,"cm")),
                                     top_annotation= if(i==1) ha_top else NULL)

        legend_list_top <- c(legend_list_top,
                             list(color_mapping_legend(hm1@matrix_color_mapping,
                                                       ncol=2,
                                                   legend_direction="horizontal",
                                                   plot = FALSE)))


        layer <- hm1+row_bplot_annot

        if(include_fsets){


            hm2 <- Heatmap(name="fset_membership",
                            fset_annot[curr_feat,,drop=FALSE],
                           col=if((max(fset_annot[curr_feat,])-min(fset_annot[curr_feat,]))==0) "whitesmoke" else c("whitesmoke",colors[(length(colors)-ncol(fset_annot)+1):length(colors)]),
                            cluster_rows=FALSE,
                            cluster_columns=FALSE,
                            width=1,
                            top_annotation= if(i==1) col_bplot_annot else NULL,
                            show_row_names=FALSE,
                            show_column_names=FALSE,
                            show_heatmap_legend=FALSE,
                            heatmap_legend_param=list(
                                title="pathway",
                                at=1:ncol(fset_annot),
                                color_bar="discrete",
                                labels=colnames(fset_annot)))

            ## hm2@matrix_color_mapping
            ## you need to set the legend colors with the full heatmap set of colors!
            ## if you just extract hm2@matrix_color_mapping from the last dtype,
            ## if not all pathways have members in that dtype, weird things start
            ## happening with the legend colors
            if (i==length(uniq_dtypes) && include_fsets) legend_list_bottom <- c(legend_list_bottom,
                         list(color_mapping_legend(ColorMapping(colors=setNames(
                                                                    colors[(length(colors)-ncol(fset_annot)+1):length(colors)],
                                                               colnames(fset_annot))),
                                                   title="pathway",
                                                   at=colnames(fset_annot),
                                                   color_bar="discrete",

                                                   legend_direction="horizontal",
                                                   labels=sapply(colnames(fset_annot),function(x) substr(x,1,60)),
                                                   plot = FALSE)))

                    layer <- hm1 + hm2 + row_bplot_annot
        }

        ###################################################################

        pushViewport(viewport(layout.pos.row =i, layout.pos_col = 1))

        .spl <- if(nrow(dat_curr)>10) hsplit else NULL
        draw(layer,
             split=.spl,
             gap=unit(c(1.25,0.25),units="cm"),
             newpage=FALSE)


        ##separator lines
        if(i < length(uniq_dtypes))  grid.lines(x=unit(c(0.05,0.95),"npc"),
                                              y=unit(c(0,0),"npc"),
                                              gp=gpar(lty="dashed",lwd=2))



        if(include_fsets) {
            decorate_annotation("fset_weights", { grid.yaxis(at =seq(floor(0.7*min(fset_wghts)),
                                                             round(1.1*(max(fset_wghts)),digits=2),length.out=5),
                                                   label = seq(floor(0.7*min(fset_wghts)),
                                                       round(1.1*(max(fset_wghts)),digits=2),length.out=5),
                                                   main=TRUE,
                                                   gp=gpar(fontsize=fsize,
                                                       lineheight=0.6))
                                                 grid.lines(c(0, 1),
                                                            unit(c(0.04, 0.04), "npc"),
                                                            gp = gpar(col = "black"))
                                                  grid.text("Feature Set\nWeights",
                                                            y=unit(1,"npc"),
                                                            gp=gpar(fontsize=fsize*1.5))
                                              })
        }

                    if (i==1) {
                decorate_annotation(paste0(ind,"_feat_weights_1"), { grid.xaxis(at =seq(0,round(1.33*(max(feat_wghts)),digits=2),length.out=5),
                                                   label = seq(0,round(1.33*(max(feat_wghts)),digits=2),length.out=5),
                                                   main=FALSE,
                                                   gp=gpar(fontsize=fsize,
                                                       lineheight=0.75))
                                                  grid.text("Feature\nWeights",
                                                            y=unit(0.9,"npc"),
                                                            x=unit(0.8,"npc"),
                                                            gp=gpar(fontsize=fsize*1.5))


                                             })
            }


        if(is.null(.spl)){
                decorate_annotation(paste0(ind,"_feat_weights"),{
                                                 grid.lines(c(0.04,0.04),
                                                            unit(c(1, 0), "npc"),
                                                            gp = gpar(col = "black"))
                                    })
        } else {
            for(i in 1:.spl){
                decorate_annotation(paste0(ind,"_feat_weights_",i),{
                                                 grid.lines(c(0.04,0.04),
                                                            unit(c(1, 0), "npc"),
                                                            gp = gpar(col = "black"))
                                    })
            }
        }

        upViewport()

    }

        seekViewport("main_vp")
    ##pushViewport(viewport(layout.pos.row =1, layout.pos_col = 1))
    ##necessary because you cannot plot the list of legends directly
    if(length(legend_list_top)>5)   {

        draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend_list_top[1:5],
                                          heatmap_legend_side="left",
             newpage=TRUE)



        draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend_list_top[6:length(legend_list_top)],
                                          heatmap_legend_side="left",
             newpage=TRUE)


   }  else   {
       draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend_list_top,
                                          heatmap_legend_side="left",
                                          newpage=TRUE)

   }
    ##upViewport()

    if(include_fsets) draw(Heatmap(matrix(nrow=0,ncol=1),
                 show_row_names=FALSE,
                 show_column_names=FALSE,
                 show_heatmap_legend = FALSE),
        heatmap_legend_list=legend_list_bottom,
       heatmap_legend_side="left",
             newpage=TRUE)

    dev.off()
}

################################
##need to fix labels for feature sets so that they appear as numbers in graph and have their separate legend

plot_network <- function(akl_obj,fsets,title,sep="_",topn=c(10,30)) {

    require(RColorBrewer)
    require(scales)
    require(igraph)
    qual_col_pals  <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
    colors <-  unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

    set_wghts <- akl_obj$akl_model$sorted_kern_weight
    feat_wghts <- rank_features(akl_obj)

    set_wghts <- set_wghts[1:min(length(set_wghts),topn[1])]
    feat_wghts <- feat_wghts[1:min(length(feat_wghts),topn[2])]


    adj_df <- foreach(i=iter(names(set_wghts)))%docomb%{
        SPICER::expand(split_set_suff(i,akl_obj$dat_grp),nms=c("fset_ind","feat_suff"))
        feats <- akl_obj$rf_models[[i]]$forest$independent.variable.names[akl_obj$rf_models[[i]]$forest$independent.variable.names%in%names(feat_wghts)]
        if(length(feats)>0){
            return(data.frame(from=toupper(fset_ind),to=toupper(feats),weight=set_wghts[i]*feat_wghts[feats],color="grey",stringsAsFactors = FALSE,check.names=FALSE) )
        } else {
            return (NULL)
        }
    }

    adj_df <- adj_df[!sapply(adj_df,is.null)]
    adj_df <- do.call(rbind,adj_df)


    nodes_set <- foreach(i=iter(names(set_wghts)),.combine=rbind)%docomb%{
        SPICER::expand(split_set_suff(i,akl_obj$dat_grp),nms=c("fset_ind","feat_suff"))
        data.frame(nodeID=toupper(fset_ind),
                   nodeID_clean=gsub("genesigdb_|_up$|_down$|_dn$|^go_","",fset_ind,ignore.case = TRUE),
                   type=TRUE,
                   weight=set_wghts[i],
                   color=colors[1],stringsAsFactors = FALSE,check.names=FALSE)
        ##colors.sets[which(names(set_wghts)%in%i)]

    }

    unique_grps <- unique(unlist(akl_obj$dat_grp))
    colors_feats <- colors[2:(length(unique_grps)+1)]
    nodes_feat <- foreach(i=iter(names(feat_wghts)),.combine=rbind)%docomb%{
        SPICER::expand(split_set_suff(i,unique_grps),nms=c("feat","feat_suff"))

        data.frame(nodeID=toupper(i),
                   nodeID_clean=feat,
                   type=FALSE,
                   weight=feat_wghts[i],
                   color=colors_feats[unique_grps%in%feat_suff],
                   stringsAsFactors = FALSE,check.names=FALSE)

    }



    nodes <- rbind(nodes_set,nodes_feat)
    ig <- graph_from_data_frame(d=adj_df,vertices=nodes,directed=FALSE)

    V(ig)$degree <- degree(ig,mode="all")
    ig <- delete_vertices(ig,V(ig)[V(ig)$degree==0])
    V(ig)$shape <- ifelse(V(ig)$type, "square", "circle")

   ## edge.coff <- mean(E(ig)$weight)
   ## ig <- delete_edges(ig,E(ig)[E(ig)$weight<edge.coff])


    E(ig)$width <- 0.1+0.5*log(E(ig)$weight/min(E(ig)$weight))
    V(ig)$size <- 3+log(V(ig)$weight/min(V(ig)$weight),1.2)

    l <- layout.davidson.harel(ig)
    l <- norm_coords(l, ymin=-2, ymax=2, xmin=-2, xmax=2)

    color_sets <- colors[(length(unique_grps)+2):(length(unique_grps)+length(set_wghts)+1)]
    ##this sets the transparency of the shades to 50%, seems reasonable for now
    color_sets <- scales::alpha(color_sets,alpha=0.5)
        ends <- ends(ig,E(ig))
    grp_mark <- lapply(unique(ends[,1]),function(x) {
        which(names(V(ig))%in%ends[ends[,1]==x,2])
    })
    names(grp_mark) <- unique(ends[,1])

    pdf(title)
    plot(ig,mark.groups=grp_mark,mark.color=color_sets,mark.border=NA,layout=l,vertex.label.font=10,vertex.label=V(ig)$nodeID_clean,vertex.label.cex=0.7,vertex.label.dist=1.1,vertex.label.degree=pi/2)
    dev.off()
}


