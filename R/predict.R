#' @title Prediction function for the latent class CTMC model
#'
#' @description Input an "lctmc" model object and individuals' disease history. This function will predict the individuals' latent class.
#'
#' @name predict
#'
#' @param object an object of class "lctmc_2x2" or "lctmc_3x3". Obtained by fitting the latent class CTMC model using the `lctmc` functions
#' @param ... the following argument must be specified as part of `...`
#' \describe{
#'   \item{df_past}{a data.frame object containing the history of observed disease status. \cr
#'         This is used to perform latent clustering as well as predict disease state at a future time point}
#'   \item{pred_threshold}{a numeric value between 0 and 1. By specifying this threshold, the predicted class will be assigned as the class with probability \cr
#'         greater than this threshold **and** has the largest probability mass. Defaults to 0.}
#' }
#'
#' @return A data.frame object containing the predicted latent class.
#'
#' @note The `df_past` data frame should contain the disease history of each person that is to be predicted. \cr
#' All covariates (x0, x1, x2, w0, w1, w2) should be included as well in the `df_past` data frame.
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]
#'
#' @example inst/examples/ex_predict.R
NULL

#' @rdname predict
#' @exportS3Method
predict.lctmc_2x2 = function(object, ...) {
  ### unpack `...`
  supps = list(...)
  df_past = supps$df_past
  pred_threshold = ifelse(is.null(supps$pred_threshold), 0, supps$pred_threshold)

  ### check argument specification
  lctmc = object
  if (!is.data.frame(df_past)) {
    stop("`df_past` must be of class 'data.frame'")
  }
  if (!all(c("id", "obsTime") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'id' and 'obsTime'")
  } else if (all(c("id", "obsTime") %in% colnames(df_past))) {
    temp = unique(df_past[c("id", "obsTime")])
    temp = as.numeric(table(temp$id)) >= 2
    if (!all(temp)) {
      stop("there should be at minimal two observation (one transition) per ID in `df_past`")
    }
  }
  if (!all(c("state_at_obsTime", "x1", "x2", "w1", "w2") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'x1', 'x2', 'w1', 'w2', 'state_at_obsTime'")
  }
  if (!all(df_past$state_at_obsTime %in% c(1, 2))) {
    stop("for a 2x2 model, 'state_at_obsTime' should only take values 1 or 2")
  }
  if (pred_threshold < 0 || pred_threshold > 1) {
    stop("`pred_threshold` should be between 0 and 1 inclusive.")
  }

  ### sort
  df_past = df_past[order(df_past$id, df_past$obsTime), ]

  ### format data for computing `bik`
  my_df1 = fmt_rowwise_trans(
    data = df_past,
    type = "2x2",
    X_names = lctmc$X_names,
    W_names = lctmc$W_names,
    scaling = 1
  )
  my_df1.IDs = unique(my_df1$df_trans[["id"]])
  my_df1.Xmat = my_df1$Xmat
  my_df1.Wmat = my_df1$Wmat
  my_df1.dt = my_df1$dt
  my_df1 = my_df1$df_trans

  ### which parameter vector to use
  model_mle = lctmc$SE$SE$mle_theta
  names(model_mle) = lctmc$SE$SE$names

  ### compute `bik`
  bik = bik_all_2x2(
    theta = model_mle,
    data = my_df1,
    Xmat = my_df1.Xmat,
    Wmat = my_df1.Wmat,
    dt = my_df1.dt,
    P.rs = FALSE,
    K = lctmc$K,
    theta.names = gen_theta_names(K = lctmc$K, type = "2x2", purpose = "bik")
  )

  ### compute probability of each individual's latent class assignment probability
  PYi = Reduce(`+`, bik)
  P.z = Map(f = function(x) x/PYi, bik)
  names(P.z) = paste("P.z", 1:lctmc$K, sep = "")

  ### generate predicted class (i.e., the class with largest probability)
  max_probs1 = Reduce(f = function(u, v) ifelse(u >= v, u, v), P.z)
  pred_class = Map(f = function(x) x == max_probs1, P.z)
  for (i in seq_along(pred_class)) {
    pred_class[[i]] = i * as.numeric(pred_class[[i]])
  }
  pred_class = Reduce(`+`, pred_class) # a numeric vector

  ### output
  latent_class = cbind(
    data.frame(id = my_df1.IDs, pred_class = pred_class),
    as.data.frame(P.z)
  )

  ### subjects that do not meet probability threshold
  preds_not_meet_threshold = apply(
    X = as.matrix(as.data.frame(P.z)),
    MARGIN = 1,
    FUN = function(x) {
      max(x) < pred_threshold
    }
  )

  ### those that do not meet threshold have predicted class turned into NA
  latent_class$pred_class_strict = ifelse(preds_not_meet_threshold, NA, latent_class$pred_class)

  ### return
  return(latent_class)
}

#' @rdname predict
#' @exportS3Method
predict.lctmc_3x3 = function(object, ...) {
  ### unpack `...`
  supps = list(...)
  df_past = supps$df_past
  pred_threshold = ifelse(is.null(supps$pred_threshold), 0, supps$pred_threshold)

  ### check argument specification
  lctmc = object
  if (!is.data.frame(df_past)) {
    stop("`df_past` must be of class 'data.frame'")
  }
  if (!all(c("id", "obsTime") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'id' and 'obsTime'")
  } else if (all(c("id", "obsTime") %in% colnames(df_past))) {
    temp = unique(df_past[c("id", "obsTime")])
    temp = as.numeric(table(temp$id)) >= 2
    if (!all(temp)) {
      stop("there should be at minimal two observation (one transition) per ID in `df_past`")
    }
  }
  if (!all(c("state_at_obsTime", "x1", "x2", "w1", "w2") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'x1', 'x2', 'w1', 'w2', 'state_at_obsTime'")
  }
  if (!all(df_past$state_at_obsTime %in% c(1, 2, 3))) {
    stop("for a 3x3 model, 'state_at_obsTime' should only take values 1, 2, or 3")
  }
  if (pred_threshold < 0 || pred_threshold > 1) {
    stop("`pred_threshold` should be between 0 and 1 inclusive.")
  }

  ### sort
  df_past = df_past[order(df_past$id, df_past$obsTime), ]

  ### format data for computing `bik`
  my_df1 = fmt_rowwise_trans(
    data = df_past,
    type = "3x3",
    X_names = lctmc$X_names,
    W_names = lctmc$W_names,
    scaling = 1
  )
  my_df1.IDs = unique(my_df1$df_trans[["id"]])
  my_df1.Xmat = my_df1$Xmat
  my_df1.Wmat = my_df1$Wmat
  my_df1.dt = my_df1$dt
  my_df1 = my_df1$df_trans


  ### which parameter vector to use
  model_mle = lctmc$SE$SE$mle_theta
  names(model_mle) = lctmc$SE$SE$names

  ### compute `bik`
  bik = bik_all_3x3(
    theta = model_mle,
    data = my_df1,
    Xmat = my_df1.Xmat,
    Wmat = my_df1.Wmat,
    dt = my_df1.dt,
    P.rs = FALSE,
    K = lctmc$K,
    theta.names = gen_theta_names(K = lctmc$K, type = "3x3", purpose = "bik")
  )

  ### compute probability of each individual's latent class assignment probability
  PYi = Reduce(`+`, bik)
  P.z = Map(f = function(x) x/PYi, bik)
  names(P.z) = paste("P.z", 1:lctmc$K, sep = "")

  ### generate predicted class (i.e., the class with largest probability)
  max_probs1 = Reduce(f = function(u, v) ifelse(u >= v, u, v), P.z)
  pred_class = Map(f = function(x) x == max_probs1, P.z)
  for (i in seq_along(pred_class)) {
    pred_class[[i]] = i * as.numeric(pred_class[[i]])
  }
  pred_class = Reduce(`+`, pred_class) # a numeric vector

  ### output
  latent_class = cbind(
    data.frame(id = my_df1.IDs, pred_class = pred_class),
    as.data.frame(P.z)
  )

  ### subjects that do not meet probability threshold
  preds_not_meet_threshold = apply(
    X = as.matrix(as.data.frame(P.z)),
    MARGIN = 1,
    FUN = function(x) {
      max(x) < pred_threshold
    }
  )

  ### those that do not meet threshold have predicted class turned into NA
  latent_class$pred_class_strict = ifelse(preds_not_meet_threshold, NA, latent_class$pred_class)

  ### return
  return(latent_class)
}
