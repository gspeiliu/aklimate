**AKLIMATE: Algorithm for Kernel Learning with Approximating Tree Ensembles**

This is an R package of AKLIMATE as described in [1]. AKLIMATE is a stacked kernel learning algorithm that integrates multi-modal data and "prior knowledge" feature sets.\
It does so by : 
1. Creating a collection of Random Forest base learners, each trained on the features in an individual feature set.
2. Converting the most relevant RF base learners into integrative RF kernels
3. Combining the RF kernels into an optimal meta-kernel via Multiple Kernel Learning.

AKLIMATE can handle binary and multiclass classification as well as regression tasks.
It can operate on continuous, binary, categorical, ordinal and count data, requires minimal pre-processing and can query thousands of feature sets.

**Training**\
*Inputs*\
*dat* : samples x features data frame where columns might be of different type \
*dat_grp* : a list of vectors, each consisting of suffixes for data types that match the ones used in dat. Each vector corresponds to a particular combination of data types that will be tested for each component RF. Only the data type combination with the
best performance for a given feature set is retained. 
The data type suffixes should be distinct from one another so that none is a proper substring of another - i.e. c('CNV','CNV_GISTIC') is not OK, but c('MUTA:HOT','MUTA:NONSENSE') is. \
*fsets* : list of prior knowledge feature sets. \
*lbls* : vector of training data labels. \
*always_add* : vector of dat column names that are to be included with each fset. \
*rf_pars* : list of parameters for RF base kernels run
* ntree : Number of trees for RF kernel construction. Default is 1000.
* min_node_prop : Minimal size of leaf nodes (unit is proportion of training set size). Default is 0.01. 
* min_nfeat : Minimal size of feature set (across all data modalities) for an RF to be constructed. Default is 15.
* mtry_prop : Proportion of features to be considered for each splitting decision. Default is 0.25
* regression_q : For regression predictions only. Quantile of the per-sample empirical distribution of absolute differences between RF sample predictions and sample label. Used for binarization of sample predictions during best RF selection. Default 0.05.
* replace : TRUE/FALSE. Is subsampling to be done with replacement? Default is FALSE.
* sample_frac : Fraction of training data points to subsample for each tree. Default is 0.5 for sampling without replacement and 1 for bootstrapping.
* ttype : Type of learning task - choices are "binary","multiclass", and "regression". Default is "binary".
* split_rule : Type of splitting criteria- choices are "gini","hellinger","variance",and "beta". See ranger documentation for more details. Default is "gini".
* importance : Rule for calculating feature and feature set importance - choices are "impurity_corrected","permutation",and "impurity". Default is "impurity_corrected".
* metric : Metric for ranking RF base learner performance used in the selection of best RFs. Choices are "roc","pr","acc","bacc","mar","rmse","rsq","mae","pearson", and "spearman". Default is "roc".
* unordered_factors : How to treat unordered factors. Choices are "order","ignore", and "partition". See ranger for more details. Default is "order".
* oob_cv : A data frame of parameters to tune during trainings of all RF base learners, with OOB metric performance (from choices above) used to select best combination. Each row of the data frame includes a different combination of RF hyperparameters. The data frame has to contain at least two columns, with one column being "ntree". Having too many hyperparameter combinations can lead to significant slowdown in computation time. Default is a data frame of 1 row using the "min_node_prop","mtry_prop", and "ntree"/2 values of the rf_pars list. This argument is experimental - we recommend using the default setting. 
*akl_pars* list of parameters for RF best kernel selection and MKL meta-learner
* topn : number of RF kernels (ranked by metric specified in rf_pars) that correctly predict a given sample to be included in best RF list. Default is 5.
* cvlen : Number of random MKL hyperparameter combinations to be tested during MKL CV step. Default is 100.
* nfold : Number of folds to be used in MKL CV. Default is 5.
* lamb : Interval bounds from which random MKL hyperparameter combinations are drawn (log2 units). Default is (-20,0).
* celnet : Hyperparameters for MKL elastic net run. Should be a vector of length 2. Default is NULL - hyperparameters are tuned via internal cross-validation.
* subsetCV : TRUE/FALSE. When TRUE, the MKL CV step also randomly varies the number of RF kernels in addition to the MKL regularization hyperparameters. It does so by training on a subset of kernels of size K, randomly selected on the (0,number best RF kernels) interval. Once K is selected, the top K kernels (ranked by metric specified in rf_pars) are included in current CV run. Default is TRUE
* type : Type of predictions - possible choices are "response" and "probability". Default is "response".
*store_kernels* TRUE/FALSE. Should the model store the training RF kernels. Default is FALSE. \
*verbose* TRUE/FALSE. Should the model print verbose progress statements. Default is FALSE.

*Output*\
An AKLIMATE model with the following components:
* rf_stats : List of metrics and predictions from training run on all RF base learners.
* kernels : RF kernels used in MKL training step. NULL if store_kernels is set to FALSE. 
* kern_cv : if akl_pars$celnet is NULL, hyperparameter vectors examined during MKL cross-validation, along with matching metric scores.
* rf_models : Set of RF base learners used to produce RF kernels for stacked MKL.
* akl_model : Trained spicer MKL model, with either user-supplied elastic net hyperparameters, or the hyperparameters selected via CV tuning.
* rf_pars_global : rf_pars argument
* rf_pars_local : optimal RF parameters for each RF base learner. Those will be the same (with the exception of ntree) as the rf_pars_global parameters unless rf_pars$oob_cv was specified by the user.
* akl_pars : akl_pars argument
* dat_grp : dat_grp argument
* idx_train : Vector of training data instances.
* preds_train : AKLIMATE predictions on training set.

**Prediction**\
*Inputs* \
*akl_obj* : an AKLIMATE model \
*dat* : samples x features data frame where columns might be of different type \
*fsets* : list of prior knowledge feature sets \
*kernels* : Pre-computed kernels for the MKL prediction part. By default the value is NULL, i.e. test kernels are computed on the fly from the RF models stored in the aklimate object. \
*store_kernels* : Boolean variable indicating whether test_kernels should be stored in output. Default is FALSE. 

*Output* \
List with following elements:
* preds : Predictions.If akl_obj$akl_pars$type is "response", this is a vector class predictions or regression estimates. If akl_obj$akl_pars$type is "probability" (classification tasks only), it is a matrix with predicted probabilities of membership for each class.
* kernels : If store_kernels is TRUE, the test kernels. Default is FALSE.

**Citations** 
1.  V. Uzunangelov, C. K. Wong, and J. Stuart. *Highly Accurate Cancer Phenotype Prediction with AKLIMATE, a Stacked Kernel Learner Integrating Multimodal Genomic Data and Pathway Knowledge.* bioRxiv, July 2020.

(c) Vladislav Uzunangelov 2019  
uzunangelov@gmail.com
