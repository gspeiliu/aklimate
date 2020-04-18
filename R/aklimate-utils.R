# (c) Vlado Uzunangelov 2017                                                             ## uzunangelov@gmail.com

##library(foreach)
library(doParallel)
library(ranger)
library(ROCR)
library(caret)
library(Similarity)
library(proxy)
library(purrr)
library(pracma)
library(stringr)
`%docomb%` <- if(getDoParRegistered()) `%dopar%` else `%do%`


##reads a file in a listt format,
##e.g. Name1 member1 member2 ...
##     Name2 member1 member2 ...
##returns a list with each line representing a separate vector object named by first entry on line and having as members subsequent entries
read_set_list <- function(file, delim="\t") {
  l <- readLines(file)
  l.fields <- strsplit(l, delim)
  r <- lapply(l.fields, function(x) as.vector(x[-1]))
  names(r) <- sapply(l.fields, "[[", 1)
  return(r)
}


######################################################
## reverse operation of read_set_list

write_set_list <- function(setList,out.file,delim="\t"){
    #since we are appending to a file,
    #we need to check if it exists and remove it
    if(out.file!=stdout() && file.exists(out.file)){
        file.remove(out.file)
    }
    for (i in 1:length(setList)){
        write(paste(names(setList)[i],paste(setList[[i]],collapse=delim),sep=delim),file=out.file,append=TRUE)
    }

}

##Description: Wrapper around write.table that adds the name of the first column
##Args:
##df - the data frame to write
##row.names.id - name for the first column
##out.file - name of file to write

write_df <- function(df,row.names.id='',out.file){
  output <- cbind(rownames(df),df)
  colnames(output)[1] <-row.names.id
  write.table(output ,file=out.file,quote=FALSE,append=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}

##Computes a statistic for each pair of entries in two lists
##list.1 - can be a data.frame (list of columns), will represent rows in output matrix
##list.2 - can be a data.frame , will represent columns in output matrix
##comp.func - a function whose arguments are source for the (full) list and target for the individual entry (vector) being compared against each element of the list
##output is a matrix of size length(list.1)xlength(list.2) with each column representing
##the statistics computed by comp.func(list.1,list.2[[x]],...)
get_cross_table <- function(list.1,list.2,comp.func,parallel=FALSE,...) {
    if (parallel) {
        cores <- if(getDoParRegistered()) getDoParWorkers() else detectCores()-2
        out <- mcmapply(comp.func,target=list.2,MoreArgs=c(list(source=list.1),list(...)),mc.cores=cores)
    } else {
        out <- mapply(comp.func,target=list.2,MoreArgs=c(list(source=list.1),list(...)))
    }
    return(out)
}

##source is a list with the rows of the output matrix
##target is a vector that will produce one column of the output matrix
sample.sim.vect <- function(source,target) {
    sapply(source,function(x)  {
        inds<-(!is.na(x))&(!is.na(target))
        sum(x[inds]==target[inds])/sum(inds)
    })
}

sample.sim.vect.mask <- function(source,target,mask) {
    mapply(function(x,y,target) sum(x[y]==target[y])/sum(y),source,mask,MoreArgs = list(target=target) )
}


sample.sim.vect.eucl <- function(source,target) {
    sapply(source,function(x) {
        inds<-(!is.na(x))&(!is.na(target))
        sqrt(sum((x[inds]-target[inds])^2))
        ##sqrt(sum((x[inds]-target[inds])^2)/sum(inds))
    })
}

sample.sim.vect.eucl.mask <- function(source,target,mask) {
    mapply(function(x,y,target) sqrt(sum((x[y]-target[y])^2)),source,mask,MoreArgs = list(target=target))
}



##Computes relative contribution of each data type to the prediction accuracy of a junkle model
##suffs - vector of suffixes of participating feature types
##ranked - vactor of ranked features produced by rank.features.jklm()$all.weights
rank.importance.type <- function(suffs,ranked,intron="_") {

    ranked.pos <- ranked[ranked>0]
    res <- sapply(suffs,function(x) {

        ## ##need to add $ at end to protect against suffs that are included in other suffs
        ## ##(cnv, cnv_gistic)
        ## ## then need to add first part of regex to pick up one hot encoded variables

        ## ##sum(ranked.pos[grepl(paste0(intron,x,"\\.","|",intron,x,"$"),names(ranked.pos),perl=TRUE)])

        ##correction to previous comments - not necessary if you use proper suffixes
        sum(ranked.pos[grepl(paste0(intron,x),names(ranked.pos),perl=TRUE)])
    })

    res <- sort(res,decreasing=TRUE)
    res <- res/sum(ranked.pos)
    return(res)
}

expand.names <- function(names,suffs,intron="_"){
    if(is.null(suffs)) return(names)
    res <- apply(expand.grid(names,suffs),1,function(x) paste(x,collapse=intron))
    return(res)
}

##cleans up names by replacing dots with underscores
clean.names <- function(names) {
    ##make names converts "-" to ".", so end result is "-" -> "_"!! Be very careful!!
    res <- gsub('[.]', '_',make.names(names,unique=TRUE))
    ##res <- gsub('\\.{1}$','',res)
    return(res)
}

##Produces a combined data frame of properly named features from all data types used
##dat - a list of data.frames sample x feature that will be collapsed into one sample x feature df
## the names of dat will become the suffixes of the features in a given data type ( to differentiate entries with the same name from different data types)
## output is a sample x feature data frame

collapse.data <- function(dat,intron="_"){

    if(length(dat)>1){
        same.names <- sapply(2:length(dat),function(x) identical(rownames(dat[[x]]),
                                                             rownames(dat[[x-1]])))

        stopifnot(sum(same.names)!=(length(same.names)-1))
    }
    dat.adj <- lapply(names(dat),function(x) {
        res <- dat[[x]]
        colnames(res) <- paste(colnames(res),x,sep=intron)
        res
    })

    out <- do.call(cbind,dat.adj)
    out <- data.frame(out)
    colnames(out) <- clean.names(colnames(out))
    return(out)
}


create.grp.suff <- function(dat.grp,intron="_") {
    paste0(intron,sapply(dat.grp,paste,collapse=intron))
    ##paste0(intron,paste0(dat.grp,collapse=intron))
}


split.set.suff <- function(combined,dat.grp,...){
    grp.suff <- create.grp.suff(dat.grp,...)

    fset.ind <- strsplit(combined,paste(grp.suff,collapse="|"))[[1]][1]
    feat.ind <- which(grp.suff==strsplit(combined,fset.ind,fixed=TRUE)[[1]][2])
    feat.suff <- if(length(feat.ind)>0) dat.grp[[feat.ind]] else NULL
    ## fset.ind <- strsplit(combined,grp.suff)[[1]][1]

    ## return(list(fset.ind=fset.ind,feat.suff=dat.grp))

    return(list(fset.ind=fset.ind,feat.suff=feat.suff))
}

extract.features<-function(allfeats,infeats,ids){
    allfeats[grepl(paste(infeats,collapse="|"),allfeats,ignore.case = TRUE)&
             grepl(paste(ids,collapse="|"),allfeats,ignore.case = TRUE)]
}


##the rows of metrics nees to be ranked from most to least desirable
which_keep <- function(metrics,comm.names,suffs){
    ind.keep <- c()
    for(i in comm.names){

        idx <- which(rownames(metrics)%in%paste0(i,suffs))
        if (length(idx)>0) ind.keep <- c(ind.keep,idx[1])
    }

    return(rownames(metrics)[sort(ind.keep)])
}

cv_grid <- function(nkern=500,len=250,lam.b=c(-5,5)){
        lam1 <- 2^runif(len, lam.b[1], lam.b[2])
        ##lam1 <- c(runif(len-2),0,1)
        lam2 <- 2^runif(len, lam.b[1], lam.b[2])

        nkern <- sample(1:nkern, size = len, replace = TRUE)
        newGrid <- data.frame(lam1=lam1,lam2=lam2,nkern=nkern)
        return(newGrid)
}


##pars - a data.frame with columns named lam1, lam2 & nkern giving the parameters to be tried with MKL
kernel.cv <- function(kerns,lbls,pars=cv_grid(nkern=dim(kerns)[3]),nfolds=5,type="binary",measure="bacc",wghts=NULL) {

    ############################
    kfolds <- createFolds(if(type=="regression") lbls else factor(lbls),
                          k=nfolds,
                          list=TRUE,
                          returnTrain = TRUE)

    .loss <- switch(type,
                    unsupervised=,
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
                spicer(kerns[kfolds[[j]],kfolds[[j]],1:pars[i,"nkern"] ,drop=FALSE],
                           lbls[kfolds[[j]]],
                           C=c(pars[i,"lam1"],pars[i,"lam2"]),
                       opt=list(loss=.loss,regname="elasticnet",
                           display=1,wghts=wghts[kfolds[[j]]]))
            },error=function(err){
                return(NULL)
            })
            .rtype <- switch(type,
                             unsupervised=,
                             multiclass=,
                             binary="probability",
                             regression="response")

            switch(type,
                   unsupervised=,
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

                metric <- select.metric(preds,lbls[-kfolds[[j]]],type,measure)

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
    ##o <- order(pars$lam1+pars$lam2)
    o <- order(res$metric,decreasing=TRUE)
    res <- res[o,]
    pars <- pars[o,]

    ##pars.id <- which.max(res$metric)
    pars.id <- round(0.1*nrow(res))
    return(list(res=res,pars=pars,best.id=pars.id))
}



is.valid.forest <- function(rf,minl=3) {
   medianl <- median(sapply(rf$forest$child.nodeIDs, function(x) length(x[[1]])))
   return(medianl>=minl)

}

correct.regression <- function(preds,lbls) {

      sq.err.f <- function(par) {
                    sum((lbls- par*preds )^2)
                }

       oo <- tryCatch({
                    optim(1,sq.err.f,method="L-BFGS-B",lower=0)
                }, error=function(cond) { return(list(par=0,value=Inf)) } )

      if(is.infinite(oo$value)) {
                    return (preds)
                } else {
                    return(oo$par*preds)
                }

}
## the binary/multiclass predictions need to be of type probability (i.e. dim=2)
select.metric <- function(preds,lbls,ttype,measure="roc") {

    metric <- switch(ttype,
           unsupervised=,
           binary={

               ##auroc
               rocr.pred <- ROCR::prediction(preds[,2],lbls)

               confM <- caret::confusionMatrix(
                   factor(apply(preds,1,function(x) colnames(preds)[which.max(x)]),
                          levels=levels(lbls)),lbls)

               msr <- switch(measure,
                             roc={
                                 ROCR::performance(rocr.pred,"auc")@y.values[[1]]
                                },
                             pr={
                                 perf <- performance(rocr.pred,"prec","rec")
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
                                       ebacc=,
                                       bacc={
                                           unname(mean(
                                               confM$byClass[,"Balanced Accuracy"]))
                                       },
                                       emar=,
                                       mar={
                                           unname(mean(
                                               confM$byClass[,"Sensitivity"]))
                                       })
                         msr
                     },
           error("Trying a method that is not implemented yet!"))

    return(metric)
}


compose.features <- function(dat,suff="composite",deg=3,topn=ncol(dat),idx=rownames(dat)) {

     comp.list <- foreach(i=iter(2:deg))%do%{
         combos <- combn(colnames(dat),i)
         out <- apply(combos,2,function(x) Emcdf::emcdf(dat[idx,x,drop=FALSE],dat[idx,x,drop=FALSE]))
         ##out <- apply(combos,2,function(x) copula::C.n(pobs(dat[idx,x,drop=FALSE]),dat[idx,x,drop=FALSE],ties.method="average"))
         rownames(out) <- idx
         colnames(out) <- paste0(apply(combos,2,function(x) paste0(x,collapse="_")),
                                    "_",suff)
         out
     }

     res <- do.call(cbind,comp.list)

     idx.test <- setdiff(rownames(dat),idx)
     res1 <- NULL

     if(length(idx.test)>0) {

     comp.list <- foreach(i=iter(2:deg))%do%{
         combos <- combn(colnames(dat),i)
         out <- apply(combos,2,function(x) Emcdf::emcdf(dat[idx,x,drop=FALSE],dat[idx.test,x,drop=FALSE]))
         ##out <- apply(combos,2,function(x) copula::C.n(pobs(dat[idx.test,x,drop=FALSE]),dat[idx,x,drop=FALSE],ties.method="average"))
         rownames(out) <- idx.test
         colnames(out) <- paste0(apply(combos,2,function(x) paste0(x,collapse="_")),
                                    "_",suff)
         out
     }
     res1 <- do.call(cbind,comp.list)
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

quantize.data <- function(dat,nbr=5,idx=rownames(dat)) {

    idx.test <- setdiff(rownames(dat),idx)
    res <- foreach(i=1:ncol(dat),.combine=cbind)%docomb%{
        breaks <- unique(quantile(dat[idx,i],seq(0,1,1/nbr)))
        breaks[length(breaks)] <- Inf
        breaks[1] <- -Inf
        as.numeric(cut(dat[idx,i],breaks=breaks))
        }

    rownames(res) <- idx

    if(length(idx.test)>0){

        res1 <- foreach(i=1:ncol(dat),.combine=cbind)%docomb%{
            breaks <- unique(quantile(dat[idx,i],seq(0,1,1/nbr)))
            breaks[length(breaks)] <- Inf
            breaks[1] <- -Inf
            as.numeric(cut(dat[idx.test,i],breaks=breaks))
        }
        rownames(res1) <- idx.test
        res <- rbind(res,res1)
    }


    colnames(res) <- colnames(dat)
    res <- res[rownames(dat),]

    return(res)
}


##default parameters for RF selected individually for each feature set
## you still need to keep the individual ones in case they are not included in oob.cv
rf.pars.default <- function(rf.pars=list()) {

    ##rf params defaults
    if(is.null(rf.pars$ntree)) rf.pars$ntree <- 1000
    if(is.null(rf.pars$min.node.prop)) rf.pars$min.node.prop <- 0.01
    if(is.null(rf.pars$min.nfeat)) rf.pars$min.nfeat <- 15
    if(is.null(rf.pars$mtry.prop)) rf.pars$mtry.prop <- 0.25
    if(is.null(rf.pars$regression.q)) rf.pars$regression.q <- 0.05
    rf.pars$replace <- if(is.null(rf.pars$replace)) FALSE else as.logical(match.arg(as.character(rf.pars$replace),c("TRUE","FALSE")))
    if(is.null(rf.pars$sample.frac)) rf.pars$sample.frac <- ifelse(rf.pars$replace,1,0.5)
    rf.pars$ttype <- match.arg(rf.pars$ttype,c("binary","regression","unsupervised","multiclass"))
    rf.pars$importance <- match.arg(rf.pars$importance,c("impurity_corrected","permutation","impurity"))

    ##this needs to be quoted out since we allow multiple metrics now
    ##rf.pars$bin.perf <- match.arg(rf.pars$bin.perf,c("roc","pr","acc","bacc","mar","rmse","rsq","mae","pearson","spearman"))
    ##oob.cv needs to have at least two entries - ntree + some parameter, otherwise errors occur
    if (is.null(rf.pars$oob.cv)) rf.pars$oob.cv <- data.frame(min.node.prop=0.05,mtry.prop=0.1,ntree=1000)
    return(rf.pars)
}

aklimate_pars_default<-function(akl_pars=list()) {

  if(is.null(akl_pars$topn)) akl_pars$topn <- 5
  if(is.null(akl_pars$cvlen)) akl_pars$cvlen <- 100
  if(is.null(akl_pars$nfold)) akl_pars$nfold <- 5
  if(is.null(akl_pars$lamb)) akl_pars$lamb <- c(-8,3)


  akl_pars$subsetCV <- if(is.null(akl_pars$subsetCV)) TRUE else as.logical(match.arg(as.character(akl_pars$subsetCV),c("TRUE","FALSE")))

  if(is.null(akl_pars$type)) akl_pars$type <- "response"
  return(akl_pars)

}
##Helper function that trains RF for each feature set

##lbls -a data.frame of 1 column
##dat -a data.frame of predictors
##rf.pars - list of parameters for RF
##always.split - vector of variables that will be added to the ones considered for splitting for each tree
## returns a trained RF model
rf.train <- function(dat,lbls,rf.pars,always.split=NULL) {

    colnames(lbls) <- c("labels")

    idx.train <- rownames(lbls)

    dat <- mlr::createDummyFeatures(dat)

        ## if(ncol(dat)>1) {
        ##     composite <- compose.features(as.matrix(dat),"composite",3,ncol(dat),idx.train)

        ##     dat <- cbind(dat,data.frame(composite))
        ## }

    dat <- cbind(lbls,dat[idx.train,,drop=FALSE])
    switch(rf.pars$ttype,
           unsupervised=,
           multiclass=,
           binary={


               rf <- ranger(data=dat[idx.train,,drop=FALSE],
                            dependent.variable.name="labels",
                            always.split.variables=always.split,
                            classification = TRUE,
                            sample.fraction = rf.pars$sample.frac,
                            num.trees=rf.pars$ntree,
                            mtry=ceiling((ncol(dat)-1)*rf.pars$mtry.prop),
                            min.node.size=ceiling(nrow(dat)*rf.pars$min.node.prop),
                            case.weights=rf.pars$weights,
                            class.weights = rf.pars$class.weights,
                            num.threads=1,
                            probability=TRUE,
                            respect.unordered.factors = FALSE,
                            importance=rf.pars$importance,
                            write.forest=TRUE,
                            keep.inbag=TRUE,
                            replace=rf.pars$replace)

           },
           regression={


               rf <- ranger(data=dat[idx.train,,drop=FALSE],
                            dependent.variable.name="labels",
                            always.split.variables=always.split,
                            sample.fraction = rf.pars$sample.frac,
                            num.trees=rf.pars$ntree,
                            mtry=ceiling((ncol(dat)-1)*rf.pars$mtry.prop),
                            min.node.size=ceiling(nrow(dat)*rf.pars$min.node.prop),
                            case.weights=rf.pars$weights,
                            num.threads=1,
                            importance=rf.pars$importance,
                            write.forest=TRUE,
                            respect.unordered.factors = FALSE,
                            keep.inbag=TRUE,
                            replace=rf.pars$replace)
           },
           error("Trying a method that is not implemented yet!"))


    return(rf)
}

##Helper function for selection of most relevant RFs from all RFs trained
##dat - a samples x featres df with all features together
## dat.grp - list of vectors, each vector containing a group of data type suffixes to be considered together
##fsets - a list of feature sets
##lbls -a df of one column - a factor - for classification and  muti-class analysis
## (status column should be 0-1 with 1's corresponding to deaths)
##always.add - vector of dat columns that should be added to each fset before RF trainins

## returns a ranking feature sets based on RF out-of-bag prediction performance on trainin set
train.forest.stats <- function(dat,dat.grp,fsets,lbls,rf.pars.global=rf.pars.default(list()),always.add=NULL,sep="_",verbose=FALSE){

    rf.pars.global <- rf.pars.default(rf.pars.global)

        idx.train <- rownames(lbls)

    grp.suff <- create.grp.suff(dat.grp,intron=sep)


    ##train models

    if(verbose) {
        print(date())
        print("Starting feature set forests")
    }

    zz <- foreach(k=1:nrow(rf.pars.global$oob.cv))%do%{

        rf.pars.local <- rf.pars.global
        rf.pars.local[colnames(rf.pars.global$oob.cv)] <- rf.pars.global$oob.cv[k,]

        vv <- foreach(j=iter(dat.grp)) %do% {
            tt <- foreach(i=iter(fsets)) %docomb% {

                feats<-extract.features(colnames(dat),i,j)
                ##feats <- colnames(dat)[colnames(dat)%in%expand.names(i,j,sep)]

                if(length(feats)<rf.pars.global$min.nfeat) return(NULL)

                curr.dat <- dat[idx.train,unique(c(feats,always.add)),drop=FALSE]

                rf <- rf.train(curr.dat,
                               lbls,
                               rf.pars.local,
                               NULL)

                imps <- importance(rf)
                imps <- imps[order(imps,decreasing=TRUE)]
                imps <- imps[imps>0]
                imps <- imps


                if(!is.valid.forest(rf,3)) return(list(metric=NA,imps=NA,preds=NA))
                ##metrics <- sapply(rf.pars.global$bin.perf,function(x) select.metric(rf$predictions,lbls[,1],rf.pars.global$ttype,x))
                ##metric <- abs(prod(metrics))*ifelse(any(sign(metrics)<0),-1,1)
                metric <- select.metric(rf$predictions,lbls[,1],rf.pars.global$ttype,rf.pars.global$bin.perf)


                ############
                preds <- if(length(dim(rf$predictions))==2) {
                    apply(rf$predictions,1,function(x) levels(lbls$labels)[which.max(x)])
                } else {
                    rf$predictions
                }
                names(preds) <- rownames(lbls)
                return(list(metric=metric,imps=imps,preds=preds))

            }

            names(tt) <- paste0(names(fsets),paste0(sep,paste(j,collapse=sep)))
            tt <- tt[!sapply(tt,is.null)]

            mms.in <- sapply(tt,function(x) x$metric)

            imps.in <- lapply(tt,function(x) x$imps)

            preds.in <- lapply(tt,function(x) x$preds)
            list(mms.in=mms.in,imps.in=imps.in,preds.in=preds.in)
        }

        mms <- unlist(lapply(vv,function(x) x$mms.in))

        imps <- do.call(c,lapply(vv, function(x) x$imps.in))

        preds <- do.call(c,lapply(vv, function(x) x$preds.in))
        list(mms=mms,imps=imps,preds=preds)
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
                                    unlist(rf.pars.global$oob.cv[x,colnames(rf.pars.global$oob.cv),drop=FALSE])
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


    msums <- rowMeans(ranks,na.rm=TRUE)


    o <- order(msums,decreasing = TRUE)

    msums <- msums[o]
    metrics <- metrics[o,,drop=FALSE]
    pars <- pars[o,,drop=FALSE]
    ranks <- ranks[o,,drop=FALSE]
    imps <- imps[o]
    preds <- preds[o]


    ##############################
    if(length(dat.grp)>1){
        ##keep only highest scoring RF corresponding to a given feature set
        kept <- which_keep(metrics,names(fsets),grp.suff)
        metrics <- metrics[kept,,drop=FALSE]
        msums <- msums[kept]
        pars <- pars[kept,,drop=FALSE]
        ranks <- ranks[kept,,drop=FALSE]
        imps <- imps[kept]
        preds <- preds[kept]
    }
    ##################################

    preds <- t(sapply(preds,function(x) x))

    ##needs to be extended for regression (refactored into a separate function)
    if(rf.pars.global$ttype=="regression"){
        preds.match <- t(apply(preds,1,function(x) x-lbls[idx.train,1]))
        preds.match <- apply(abs(preds.match),2,function(x) x<quantile(x,rf.pars.global$regression.q))

    } else {
        preds.match <- t(apply(preds,1,function(x) x==levels(lbls[idx.train,1])[lbls[idx.train,1]]))
    }
    return(list(pars.global=rf.pars.global,pars.local=pars,res=zz,stats=metrics,msums=msums,rankmax=rankmax,importance=imps,predictions=preds,predictions.match=preds.match))
}



##lbls - a data.frame of 1 or 2 columns
##rf.pars.local - the individual parameters for each forest (determined by oob error)- a df with colnames corresponding to parameters in rf.pars.default - rows correspond to each forest
train.forest.kernels <- function(dat,dat.grp,fsets,lbls,rf.pars.local,rf.pars.global=rf.pars.default(list()),always.add=NULL,sep="_",verbose=FALSE) {

    idx.train <- rownames(lbls)

    if(verbose) {
        print(date())
        print("Started training forest")
    }

    ii <- if(is.null(rf.pars.local)) names(fsets) else rownames(rf.pars.local)

    rf.out <- foreach(i=iter(ii),.options.multicore=list(preschedule=TRUE))%docomb%{

        ##fset.ind is the name of the set
        ##feat.suff is a vector of suffices that denominate data types
        expand(split.set.suff(i,dat.grp),nms=c("fset.ind","feat.suff"))

        feats<-extract.features(colnames(dat),fsets[[fset.ind]],feat.suff)

        ##feats <- colnames(dat)[colnames(dat)%in%expand.names(fsets[[fset.ind]],feat.suff,sep)]

        curr.dat <- dat[idx.train,unique(c(feats,always.add)),drop=FALSE]

        if (!is.null(rf.pars.local)) rf.pars.global[colnames(rf.pars.local)] <- rf.pars.local[i,]
        rf <- rf.train(curr.dat,
                       lbls,
                       rf.pars.global,
                       NULL)

        return(rf)

    }
    names(rf.out) <- ii

    ##gc()
    if(verbose) {
        print(date())
        print("Finished training forests")
    }

    return(rf.out)
}

#####
##dat - same input as junkle; by default idx.train are presumed to be all indices in rows of dat - if some of those are idx.test, specify accordingly
##dat rownames (samples) should be ordered accordingly

forest.to.kernel <- function(rf.models,dat,dat.grp,fsets,always.add=NULL,idx.train=rownames(dat),sep="_",verbose=FALSE){

    library(abind)

    ##check all idx.train indices are in the data provided
    stopifnot(length(idx.train)==length(intersect(rownames(dat),idx.train)))


    idx.test <- setdiff(rownames(dat),idx.train)
    ##if no idx.test , then we ar econstructing trianing kernels
    if (length(idx.test)==0) idx.test <- idx.train

    if(verbose) {
        print(date())
        print("Started kernel construction")
    }

    medianl <- sapply(1:length(rf.models),function(x) median(sapply(rf.models[[x]]$forest$child.nodeIDs, function(x) length(x[[1]]))))
    names(medianl) <- names(rf.models)


    kerns <- foreach(i=iter(names(rf.models)),.options.multicore=list(preschedule=FALSE))%docomb%{

        simn<-10
        if(match(i,names(rf.models))<2*simn){
            Sys.sleep(10*floor(match(i,names(rf.models))/simn))
        }


        expand(split.set.suff(i,dat.grp),nms=c("fset.ind","feat.suff"))

        feats<-extract.features(colnames(dat),fsets[[fset.ind]],feat.suff)


        curr.dat <- dat[,unique(c(feats,always.add)),drop=FALSE]
        curr.dat <- mlr::createDummyFeatures(curr.dat)

        ## if(ncol(curr.dat)>1) {
        ##     composite <- compose.features(as.matrix(curr.dat),"composite",3,ncol(curr.dat),idx.train)
        ##     curr.dat <- cbind(curr.dat,data.frame(composite))
        ## }

        q <- 1
        p <- 2

        if(medianl[i]>1){

            preds <- predict(rf.models[[i]],
                             curr.dat,
                             predict.all = TRUE)$predictions

            rownames(preds) <- rownames(curr.dat)

            multic<-length(dim(preds))>2 && dim(preds)[2]>2

            if(multic){
                extensions<-rf.models[[i]]$multic.rel
                counter<-which(rf.models[[i]]$forest$levels%in%rf.models[[i]]$multic.rel)
            } else {
                extensions<-"combined"
            }
            knames<-c(paste(i,extensions,sep=sep))

            proximity<-replicate(length(extensions),matrix(1,length(idx.train),length(idx.test)))

            dimnames(proximity)[[1]] <- idx.train
            dimnames(proximity)[[2]] <- idx.test
            dimnames(proximity)[[3]] <- knames

            if(multic && length(counter)>0) {
                for(j in 1:length(counter)){

                    pA <- Similarity::distance(preds[idx.train,counter[j],],
                                               preds[idx.test,counter[j],],
                                               method="euclidian",p=p)

                    if(identical(idx.train,idx.test)){
                        sA <- quantile(pA,probs=q)
                    } else {
                        sA <- quantile(Similarity::distance(preds[idx.train,counter[j],],
                                                            preds[idx.train,counter[j],],
                                                            method="euclidian",p=p),
                                       probs=q)
                    }

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

             prox <- Similarity::distance(preds[idx.train,],
                               preds[idx.test,],
                               method="euclidian",p=p)

             if(identical(idx.train,idx.test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma <- quantile(Similarity::distance(preds[idx.train,],
                                            preds[idx.train,],
                                            method="euclidian",p=p),probs=q)
             }

             proximity[,,1]<-1-(prox/sigma)


             ###############################
             preds1 <- predict(rf.models[[i]],
                               curr.dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr.dat)

             prox <- Similarity::distance(preds1[idx.train,],preds1[idx.test,],method="euclidian",p=p)

             if(identical(idx.train,idx.test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma <- quantile(Similarity::distance(preds1[idx.train,],
                                  preds1[idx.train,],
                                  method="euclidian",p=p),probs=q)

             }

             proximity[,,1]<-proximity[,,1]*(1-(prox/sigma))

             #####################################################

             prox <- Similarity::proximityMatrixRanger(curr.dat[idx.train,],
                                         curr.dat[idx.test,],
                                          rf.models[[i]])

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


        rm(curr.dat,preds,preds1,prox)
        ##gc()
        proximity
    }

    kerns <- abind(kerns,along=3)

    if(verbose) {
        print(date())
        print("Finished kernel construction")
    }

    rm(rf.models)

    return(kerns)
}

forest.to.kernel.oob <- function(rf.models,dat,dat.grp,fsets,always.add=NULL,idx.train=rownames(dat),sep="_",verbose=FALSE){

    library(abind)

    ##check all idx.train indices are in the data provided
    stopifnot(length(idx.train)==length(intersect(rownames(dat),idx.train)))


    idx.test <- setdiff(rownames(dat),idx.train)
    ##if no idx.test , then we ar econstructing trianing kernels
    if (length(idx.test)==0) idx.test <- idx.train

    if(verbose) {
        print(date())
        print("Started kernel construction")
    }

    medianl <- sapply(1:length(rf.models),function(x) median(sapply(rf.models[[x]]$forest$child.nodeIDs, function(x) length(x[[1]]))))
    names(medianl) <- names(rf.models)

    kerns <- foreach(i=iter(names(rf.models)),.options.multicore=list(preschedule=FALSE))%docomb%{

        simn<-10
        if(match(i,names(rf.models))<2*simn){
            Sys.sleep(10*floor(match(i,names(rf.models))/simn))
        }


        expand(split.set.suff(i,dat.grp),nms=c("fset.ind","feat.suff"))

        feats<-extract.features(colnames(dat),fsets[[fset.ind]],feat.suff)

        ##feats <- colnames(dat)[colnames(dat)%in%expand.names(fsets[[fset.ind]],feat.suff,sep)]

        curr.dat <- dat[,unique(c(feats,always.add)),drop=FALSE]
        curr.dat <- mlr::createDummyFeatures(curr.dat)

        q <- 1
        p <- 2

        if(medianl[i]>1){

            mask <- rf.models[[i]]$inbag.counts
            mask <- t(matrix(unlist(mask),ncol=length(mask[[1]]),byrow=TRUE))
            mask <- mask!=0


            preds <- predict(rf.models[[i]],
                             curr.dat,
                             predict.all = TRUE)$predictions

            rownames(preds) <- rownames(curr.dat)



            multic<-length(dim(preds))>2 && dim(preds)[2]>2

            if(multic){
                extensions<-rf.models[[i]]$multic.rel
                counter<-which(rf.models[[i]]$forest$levels%in%rf.models[[i]]$multic.rel)
            } else {
                extensions<-"combined"
            }

            knames<-c(paste(i,extensions,sep=sep))

            proximity<-replicate(length(extensions),matrix(1,length(idx.train),length(idx.test)))

            dimnames(proximity)[[1]] <- idx.train
            dimnames(proximity)[[2]] <- idx.test
            dimnames(proximity)[[3]] <- knames

            if(multic && length(counter)>0) {
                for(j in 1:length(counter)){
                    pl <- preds[,counter[j],]

                    pl[mask]<-NA

                    pl<-as.data.frame(t(pl))

                    pA<-get_cross_table(pl[,idx.train],
                                             pl[,idx.test],
                                             sample.sim.vect.eucl)

                    if(identical(idx.train,idx.test)){
                        sA <- quantile(pA,probs=q)
                    } else {

                        sA<-quantile(get_cross_table(pl[,idx.train],
                                                 pl[,idx.train],
                                                 sample.sim.vect.eucl),
                                     probs=q)
                    }

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
                        ##preds[mask] <- 0
                        preds[mask] <- NA
                    } else {
                        rps <- dim(preds)[2]
                    ##multiclass classification -
                    ##output is again 3D array but second dimension corresponds to more
                    ##classes and you can't just pick one to compute distance on
                    ##will have to concatenate them
                        preds <- t(apply(preds,1,c))
                    ##need to replicate mask matrix to match
                        mask.extend <- matrix(data=apply(mask,2,rep,rps),
                                   ncol=ncol(mask)*rps)
                        ##preds[mask.extend] <- 0
                        preds[mask.extend]<-NA
                    }
                } else {

                    preds[mask]<-NA
                    ##preds[mask]<-0
                }

                preds<-as.data.frame(t(preds))

                prox<-get_cross_table(preds[,idx.train],
                                         preds[,idx.test],
                                         sample.sim.vect.eucl)

             if(identical(idx.train,idx.test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma<-quantile(get_cross_table(preds[,idx.train],
                                          preds[,idx.train],
                                          sample.sim.vect.eucl),
                                 probs=q)
             }


                proximity[,,1]<-1-(prox/sigma)


             ###############################
             preds1 <- predict(rf.models[[i]],
                               curr.dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr.dat)

                ##preds1[mask]<-0
                preds1[mask] <- NA

                preds1<-as.data.frame(t(preds1))

                prox<-get_cross_table(preds1[,idx.train],
                                           preds1[,idx.test],
                                           sample.sim.vect.eucl)


             if(identical(idx.train,idx.test)){
                 sigma <- quantile(prox,probs=q)
             } else {
                 sigma<-quantile(get_cross_table(preds1[,idx.train],
                                                      preds1[,idx.train],
                                                      sample.sim.vect.eucl),
                                 probs=q)

             }


             proximity[,,1]<-proximity[,,1]*(1-(prox/sigma))

             #####################################################


             preds1 <- predict(rf.models[[i]],
                               curr.dat,
                               type="terminalNodes")$predictions

             rownames(preds1) <- rownames(curr.dat)
                preds1[mask] <- NA

              preds1 <- as.data.frame(t(preds1))

                prox  <-  get_cross_table(preds1[,idx.train],
                                               preds1[,idx.test],
                                               sample.sim.vect)


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


        rm(curr.dat,preds,preds1,prox,mask)
        proximity
    }

    kerns <- abind(kerns,along=3)

    if(verbose) {
        print(date())
        print("Finished kernel construction")
    }

    rm(rf.models)

    return(kerns)
}



rank.features.jklm <- function(jklmobj) {

    if(jklmobj$rf.pars.global$ttype=="multiclass"){
        imps.list <- foreach(k=iter(jklmobj$junkle.model))%do%{
                    wghts <- k$sorted_kern_weight
                    imps <- foreach(i=iter(names(wghts)))%do%{
                        ##you need to anchor the patterns at the start of the string
                        ## in case there are pathway names that are exact substrings of
                        ##other pathway names
                        ##this will become a problem if the matching substring is
                        ##at the beginning of the longer name
                        ind<-which(stringr:::str_detect(i,paste0("^",names(jklmobj$rf.models))))

                        res <- sort(ranger:::importance(jklmobj$rf.models[[ind]]),decreasing=TRUE)
                        res <- res[res>0]
                        res
                    }
                    names(imps) <- names(wghts)

                    scaled.imps <- lapply(names(wghts), function(x) imps[[x]]*wghts[x]/sum(imps[[x]]))
                    names(scaled.imps) <- names(wghts)

                    all.names <- unique(unlist(lapply(scaled.imps,names)))
                    comb.imps <- rep(0,length(all.names))
                    names(comb.imps) <- all.names

                    for(i in names(scaled.imps)) {
                        comb.imps[names(scaled.imps[[i]])] <- comb.imps[names(scaled.imps[[i]])] + scaled.imps[[i]]
                    }
                    comb.imps <- sort(comb.imps,decreasing=TRUE)
        }

        full.names <- unique(names(unlist(imps.list)))
        out <- rep(0,length(full.names))
        names(out) <- full.names

        for(l in 1:length(imps.list)){
            out[names(imps.list[[l]])] <-  out[names(imps.list[[l]])] + imps.list[[l]]
        }
        out <- sort(out,decreasing=TRUE)/sum(out)


    } else {
        wghts <- jklmobj$junkle.model$sorted_kern_weight
        imps <- foreach(i=iter(names(wghts)))%do%{
            ind<-which(stringr:::str_detect(i,names(jklmobj$rf.models)))
            res <- sort(ranger:::importance(jklmobj$rf.models[[ind]]),decreasing=TRUE)
            res <- res[res>0]
            res
        }
        names(imps) <- names(wghts)

        scaled.imps <- lapply(names(wghts), function(x) imps[[x]]*wghts[x]/sum(imps[[x]]))
        names(scaled.imps) <- names(wghts)

        all.names <- unique(unlist(lapply(scaled.imps,names)))
        comb.imps <- rep(0,length(all.names))
        names(comb.imps) <- all.names

        for(i in names(scaled.imps)) {
            comb.imps[names(scaled.imps[[i]])] <- comb.imps[names(scaled.imps[[i]])] + scaled.imps[[i]]
        }
        out <- sort(comb.imps,decreasing=TRUE)


    }

    return(out)
}

##dat - a data frame samples x features with all data types combined (can have factors)
##clean up names of pathways beforehand
##you should pre-filter set.wghts and feat.wghts to the # you want plotted!!
##the first entry in should.scale is the one whose column heatmap is used if cluster.col is TRUE
##the last column in lbls is the one that we use to group columns
plot.heatmap.jklm <- function(feat.wghts,dat.grp,dat,lbls,pdf.title,fset.wghts=NULL,fsets=NULL,should.scale=NULL,ordered=1,fsize=6,cluster_rows=TRUE,cluster_cols=TRUE,hsplit=NULL){


    library(ComplexHeatmap)
    library(circlize)
    library(foreach)
    ##require(colorRamps)
    require(RColorBrewer)

    ##color pallettes
    qual_col_pals  <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
    colors <-  unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))


    cont_pals <- lapply(c("RdYlBu","RdBu","RdYlGn","PuOr","RdGy","PRGn","BrBG"),
                        function(x) brewer.pal(n = 11,name = x))

    seq_col_pals <- lapply(rownames(brewer.pal.info)[brewer.pal.info$category == 'seq'],
                        function(x) brewer.pal(n = 9,name = x))

    ##prep data
    factor.feats <- names(feat.wghts)[sapply(dat[,names(feat.wghts)],is.factor)]
    dat[,factor.feats] <- apply(dat[,factor.feats,drop=FALSE],2,function(x) as.numeric(addNA(x)))

    include.fsets <- !is.null(fsets) && length(fset.wghts) >1

    stopifnot(length(intersect(rownames(dat.comb),rownames(lbls)))==nrow(lbls))

    dat <- dat[rownames(lbls),]
    ##annotation labels should be in a sample x annotation labels df
    oc <- do.call(order,lbls[,ordered,drop=FALSE])
    lbls <- lbls[oc,,drop=FALSE]
    dat <- dat[oc,,drop=FALSE]

    ##Data Preprocessing
    prsd.feats <- t(sapply(names(feat.wghts),function(x) {
        split.set.suff(x,unique(unlist(dat.grp)))
    }))

    if(!is.null(fsets)){

        nfset <- gsub(paste0(paste0("_",sapply(dat.grp,
                                               function(x) paste0(x,collapse ="_"))),
                             collapse="|"),
                      "",names(fset.wghts))
        names(fset.wghts) <- nfset

        fsets.trimmed <- lapply(fsets[names(fset.wghts)],function(x) {
            rownames(prsd.feats)[prsd.feats[,"fset.ind"]%in%x]
        })


    fset.wghts <- fset.wghts[names(fsets.trimmed)]
    fset.annot <- list.to.matrix.membership(fsets.trimmed)

    ##if just showing one feature set, trim features to match it
    if(length(fsets.trimmed)==1) prsd.feats <- prsd.feats[fsets.trimmed[[1]],,drop=FALSE]

    }




    #######################
    uniq.dtypes <- unlist(unique(prsd.feats[,"feat.suff"]))
    uniq.feats <-   unlist(unique(prsd.feats[,"fset.ind"]))

    if(is.null(should.scale)) should.scale <- setNames(rep(FALSE,length(uniq.dtypes)),
                                                       uniq.dtypes)
    num.feats <- table(as.character(prsd.feats[,2]))
    num.feats <- num.feats[uniq.dtypes]

    ##add pseudocounts for small groups
    num.feats <- num.feats+4
    ##account for annotations
    num.feats[1] <- 2*num.feats[1]

    if(include.fsets) {
        ##this is for the features that did not fall in any of the top pathways
        empty.annot <- matrix(0,nrow=nrow(prsd.feats)-nrow(fset.annot),ncol=ncol(fset.annot))
        rownames(empty.annot) <- setdiff(rownames(prsd.feats),rownames(fset.annot))
        fset.annot <- rbind(fset.annot,empty.annot)
        o <- order(fset.wghts,decreasing=TRUE)
        fset.annot <- fset.annot[,names(fset.wghts)[o],drop=FALSE]

        fset.annot.col <- sapply(1:ncol(fset.annot),function(x) {
            out <- rep(0,nrow(fset.annot))
            out[fset.annot[,x]>0] <- x
            out
        })
        rownames(fset.annot.col) <- rownames(fset.annot)
        colnames(fset.annot.col) <- colnames(fset.annot)

        fset.annot <- fset.annot.col
        rm(fset.annot.col)
    }



    ##############################
    ##Heatmap column annotations

    ##color mapping
    cc <- 0
    cc_seq <- 1
    lbl.col <- list()

    for (x in colnames(lbls)){
        if(!is.factor(lbls[,x])) {
                    br <- unname(quantile(lbls[,x],probs=seq(0,1,0.12),na.rm=TRUE))
                    br[1] <- floor(br[1])
                    br[length(br)] <- ceiling(br[length(br)])

                    lbl.col <- c(lbl.col,list(cols=colorRamp2(br,seq_col_pals[[cc_seq]])))
                    cc_seq <- cc_seq+1

         } else  {
            lbl.col <-  c(lbl.col,list(cols=setNames(colors[cc+1:length(levels(lbls[,x]))], levels(lbls[,x]))))
             cc <- cc+length(levels(lbls[,x]))
         }
    }
    names(lbl.col) <- colnames(lbls)


    ##unit(rep(3,ncol(lbls)),"mm")
    ##since the annotation height is included in the top heatmap height,
    ##best to keep the units relative since this is how we determine the
    ## relative size of viewports on the layout
    ha.top <- HeatmapAnnotation(df=lbls,
                                   col=lbl.col,
                                   na_col = "whitesmoke",
                                   show_annotation_name = TRUE,
                                   annotation_name_offset = unit(1,"mm"),
                                   annotation_name_gp=gpar(fontsize=fsize*1.5),
                                   annotation_name_side = "left",
                                   height=unit(0.05*ncol(lbls),"npc"),
                                annotation_height=if(ncol(lbls)<5) unit(rep(0.05,ncol(lbls)),"npc") else unit(rep(0.5/ncol(lbls),ncol(lbls)),"npc"),
                                   show_legend=FALSE)
    legend.list.top <- lapply(colnames(lbls),function(x) color_mapping_legend(ha.top@anno_list[[x]]@color_mapping,nrow=NULL,ncol=2,plot = FALSE))

    names(legend.list.top) <- colnames(lbls)

    legend.list.bottom <- list()
    ############################

    ##Plotting and Display

    pdf(pdf.title,paper="a4",
        width=unit(8,"inches"),
        height=unit(12*sum(num.feats)/100,"inches"))
    grid.newpage()
    main_vp <- viewport(name="main_vp",
                        layout = grid.layout(nr = length(uniq.dtypes),
                              nc = 1,
                              height=c(num.feats/sum(num.feats),0.05)))
    pushViewport(main_vp)

    if(include.fsets) {
        col.bplot.annot <- HeatmapAnnotation(fset_weights=column_anno_barplot(fset.wghts[colnames(fset.annot)],
                                                 baseline=floor(0.7*min(fset.wghts)),
                                             bar_width=0.5,
                                                ylim=c(floor(0.7*min(fset.wghts)),
                                                    max(fset.wghts)),
                                                 border=FALSE,axis=FALSE,
                                                 axis_side="left",
                                             gp=gpar(fill="#CCCCCC",col="#CCCCCC")),
                                             which="column",
                                         annotation_height = unit(15,"mm"))
    }


    for(i in 1:length(uniq.dtypes)) {
        ind <- uniq.dtypes[i]
        curr.feat <- rownames(prsd.feats)[prsd.feats[,2]==ind]
        names(curr.feat) <- as.character(prsd.feats[prsd.feats[,2]==ind,1])

        if(should.scale[ind]) {
            dat.curr <- scale(dat[rownames(lbls),curr.feat,drop=FALSE])
        } else {
            dat.curr <- dat[rownames(lbls),curr.feat,drop=FALSE]
        }

            dat.curr <- t(dat.curr)
            rownames(dat.curr) <- names(curr.feat)


        ##annotation_width=unit(2,"cm"),
            row.bplot.annot <-  HeatmapAnnotation(feat_weights=row_anno_barplot(feat.wghts[curr.feat],
                                             baseline=0,
                                             size=unit(min(fsize/4,2),"mm"),
                                                  ylim=c(0,1.1*max(feat.wghts)),
                                             border=FALSE,axis=FALSE,axis_side="top",
                                             gp=gpar(fill="#CCCCCC",col="#CCCCCC")),
                                              which="row",
                                              width=unit(0.067,"npc"))
            row.bplot.annot@anno_list[[1]]@name <- paste0(ind,"_feat_weights")
            ############

        ## if(i==1) {

        ##     chclust <- group_hclust(1-cor(dat.curr,method="spearman"),lbls)
        ## }


        ##viridis::inferno(11)),
                                     ##column_dend_reorder=TRUE,
                                     ##show_column_dend=if(i==1) TRUE else FALSE,

        breaks <- unname(quantile(dat.curr,probs=seq(0,1,0.1)))
        breaks[1] <- floor(breaks[1])
        breaks[length(breaks)] <- ceiling(breaks[length(breaks)])

        hm1 <- Heatmap(dat.curr,
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
                                     top_annotation= if(i==1) ha.top else NULL)

        legend.list.top <- c(legend.list.top,
                             list(color_mapping_legend(hm1@matrix_color_mapping,
                                                       ncol=2,
                                                   legend_direction="horizontal",
                                                   plot = FALSE)))


        layer <- hm1+row.bplot.annot

        if(include.fsets){


            hm2 <- Heatmap(name="fset_membership",
                            fset.annot[curr.feat,,drop=FALSE],
                           col=if((max(fset.annot[curr.feat,])-min(fset.annot[curr.feat,]))==0) "whitesmoke" else c("whitesmoke",colors[(length(colors)-ncol(fset.annot)+1):length(colors)]),
                            cluster_rows=FALSE,
                            cluster_columns=FALSE,
                            width=1,
                            top_annotation= if(i==1) col.bplot.annot else NULL,
                            show_row_names=FALSE,
                            show_column_names=FALSE,
                            show_heatmap_legend=FALSE,
                            heatmap_legend_param=list(
                                title="pathway",
                                at=1:ncol(fset.annot),
                                color_bar="discrete",
                                labels=colnames(fset.annot)))

            ## hm2@matrix_color_mapping
            ## you need to set the legend colors with the full heatmap set of colors!
            ## if you just extract hm2@matrix_color_mapping from the last dtype,
            ## if not all pathways have members in that dtype, weird things start
            ## happening with the legend colors
            if (i==length(uniq.dtypes) && include.fsets) legend.list.bottom <- c(legend.list.bottom,
                         list(color_mapping_legend(ColorMapping(colors=setNames(
                                                                    colors[(length(colors)-ncol(fset.annot)+1):length(colors)],
                                                               colnames(fset.annot))),
                                                   title="pathway",
                                                   at=colnames(fset.annot),
                                                   color_bar="discrete",

                                                   legend_direction="horizontal",
                                                   labels=sapply(colnames(fset.annot),function(x) substr(x,1,60)),
                                                   plot = FALSE)))

                    layer <- hm1 + hm2 + row.bplot.annot
        }

        ###################################################################

        pushViewport(viewport(layout.pos.row =i, layout.pos.col = 1))

        .spl <- if(nrow(dat.curr)>10) hsplit else NULL
        draw(layer,
             split=.spl,
             gap=unit(c(1.25,0.25),units="cm"),
             newpage=FALSE)


        ##separator lines
        if(i < length(uniq.dtypes))  grid.lines(x=unit(c(0.05,0.95),"npc"),
                                              y=unit(c(0,0),"npc"),
                                              gp=gpar(lty="dashed",lwd=2))



        if(include.fsets) {
            decorate_annotation("fset_weights", { grid.yaxis(at =seq(floor(0.7*min(fset.wghts)),
                                                             round(1.1*(max(fset.wghts)),digits=2),length.out=5),
                                                   label = seq(floor(0.7*min(fset.wghts)),
                                                       round(1.1*(max(fset.wghts)),digits=2),length.out=5),
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
                decorate_annotation(paste0(ind,"_feat_weights_1"), { grid.xaxis(at =seq(0,round(1.33*(max(feat.wghts)),digits=2),length.out=5),
                                                   label = seq(0,round(1.33*(max(feat.wghts)),digits=2),length.out=5),
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
    ##pushViewport(viewport(layout.pos.row =1, layout.pos.col = 1))
    ##necessary because you cannot plot the list of legends directly
    if(length(legend.list.top)>5)   {

        draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend.list.top[1:5],
                                          heatmap_legend_side="left",
             newpage=TRUE)



        draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend.list.top[6:length(legend.list.top)],
                                          heatmap_legend_side="left",
             newpage=TRUE)


   }  else   {
       draw(Heatmap(matrix(nrow=0,ncol=1),
                                                  show_row_names=FALSE,
                                                  show_column_names=FALSE,
                                                  show_heatmap_legend = FALSE),
                                          heatmap_legend_list=legend.list.top,
                                          heatmap_legend_side="left",
                                          newpage=TRUE)

   }
    ##upViewport()

    if(include.fsets) draw(Heatmap(matrix(nrow=0,ncol=1),
                 show_row_names=FALSE,
                 show_column_names=FALSE,
                 show_heatmap_legend = FALSE),
        heatmap_legend_list=legend.list.bottom,
       heatmap_legend_side="left",
             newpage=TRUE)

    dev.off()
}

################################
##need to fix labels for feature sets so that they appear as numbers in graph and have their separate legend

plot.network.jklm <- function(jklm,fsets,title,sep="_",topn=c(10,30)) {

    require(RColorBrewer)
    require(scales)
    qual_col_pals  <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
    colors <-  unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

    set.wghts <- jklm$junkle.model$sorted_kern_weight
    feat.wghts <- rank.features.jklm(jklm)

    set.wghts <- set.wghts[1:min(length(set.wghts),topn[1])]
    feat.wghts <- feat.wghts[1:min(length(feat.wghts),topn[2])]


    adj.df <- foreach(i=iter(names(set.wghts)))%docomb%{
        expand(split.set.suff(i,jklm$dat.grp),nms=c("fset.ind","feat.suff"))
        feats <- jklm$rf.models[[i]]$forest$independent.variable.names[jklm$rf.models[[i]]$forest$independent.variable.names%in%names(feat.wghts)]
        if(length(feats)>0){
            return(data.frame(from=toupper(fset.ind),to=toupper(feats),weight=set.wghts[i]*feat.wghts[feats],color="grey",stringsAsFactors = FALSE,check.names=FALSE) )
        } else {
            return (NULL)
        }
    }

    adj.df <- adj.df[!sapply(adj.df,is.null)]
    adj.df <- do.call(rbind,adj.df)


    nodes.set <- foreach(i=iter(names(set.wghts)),.combine=rbind)%docomb%{
        expand(split.set.suff(i,jklm$dat.grp),nms=c("fset.ind","feat.suff"))
        data.frame(nodeID=toupper(fset.ind),
                   nodeID.clean=gsub("genesigdb_|_up$|_down$|_dn$|^go_","",fset.ind,ignore.case = TRUE),
                   type=TRUE,
                   weight=set.wghts[i],
                   color=colors[1],stringsAsFactors = FALSE,check.names=FALSE)
        ##colors.sets[which(names(set.wghts)%in%i)]

    }

    unique.grps <- unique(unlist(jklm$dat.grp))
    colors.feats <- colors[2:(length(unique.grps)+1)]
    nodes.feat <- foreach(i=iter(names(feat.wghts)),.combine=rbind)%docomb%{
        expand(split.set.suff(i,unique.grps),nms=c("feat","feat.suff"))

        data.frame(nodeID=toupper(i),
                   nodeID.clean=feat,
                   type=FALSE,
                   weight=feat.wghts[i],
                   color=colors.feats[unique.grps%in%feat.suff],
                   stringsAsFactors = FALSE,check.names=FALSE)

    }



    nodes <- rbind(nodes.set,nodes.feat)
    ig <- graph_from_data_frame(d=adj.df,vertices=nodes,directed=FALSE)

    V(ig)$degree <- degree(ig,mode="all")
    ig <- delete_vertices(ig,V(ig)[V(ig)$degree==0])
    V(ig)$shape <- ifelse(V(ig)$type, "square", "circle")

   ## edge.coff <- mean(E(ig)$weight)
   ## ig <- delete_edges(ig,E(ig)[E(ig)$weight<edge.coff])


    E(ig)$width <- 0.1+0.5*log(E(ig)$weight/min(E(ig)$weight))
    V(ig)$size <- 3+log(V(ig)$weight/min(V(ig)$weight),1.2)

    l <- layout.davidson.harel(ig)
    l <- norm_coords(l, ymin=-2, ymax=2, xmin=-2, xmax=2)

    color.sets <- colors[(length(unique.grps)+2):(length(unique.grps)+length(set.wghts)+1)]
    ##this sets the transparency of the shades to 50%, seems reasonable for now
    color.sets <- scales::alpha(color.sets,alpha=0.5)
        ends <- ends(ig,E(ig))
    grp.mark <- lapply(unique(ends[,1]),function(x) {
        which(names(V(ig))%in%ends[ends[,1]==x,2])
    })
    names(grp.mark) <- unique(ends[,1])

    pdf(title)
    plot(ig,mark.groups=grp.mark,mark.color=color.sets,mark.border=NA,layout=l,vertex.label.font=10,vertex.label=V(ig)$nodeID.clean,vertex.label.cex=0.7,vertex.label.dist=1.1,vertex.label.degree=pi/2)
    dev.off()
}


