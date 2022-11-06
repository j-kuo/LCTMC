#' @title Prediction function for the latent class CTMC model
#'
#' @description Input a "lctmc" object and performs two types of prediction via the estimated model: \cr
#' \enumerate{
#'   \item predicts the latent cluster classification. \cr
#'   \item predicts the disease status at a user specified time
#' }
#'
#' @name predict
#'
#' @param object an object of class "lctmc_2x2" or "lctmc_3x3". Obtained by fitting the latent class CTMC model using the `lctmc` functions
#' @param ... the following argument must be specified as part of `...`
#' \describe{
#'   \item{df_pred}{a data.frame object containing two columns: 'id' and 'obsTime'. \cr
#'         The first is the identifier for each individual and the second is a numeric continuous variable, indicating the future time point that we wish to predict the disease status.
#'         This time variable should be greater than the individual-wise maximum time specified in the `df_past` data frame}
#'   \item{df_past}{a data.frame object containing the history of observed disease status. This is used to perform latent clustering as well as predict disease state at a future time point}
#'   \item{param.type}{a character scalar. Either "mle" or "kmeans". \cr
#'         If "mle" then the estimated MLE of the latent class CTMC is used as the parameter for prediction. If "kmeans" then Step 1 of the initial value generation process (K-Means algorithm) is used as the prediction parameters.}
#' }
#'
#' @return A list object containing two elements:
#' \itemize{
#'   \item **df_pred.lc** is the predicted latent class classification. The predicted class is the class with the largest probability mass.
#'   \item **df_pred.ds** is the predicted disease status. This should be the same size as the input argument `df_pred` with an additional column added for the predicted class at the respective observation time.
#' }
#'
#' @note In the `df_pred` data frame, there should only be one ID per person and "obsTime" must be a future time point relative to the `df_past` data frame.
#' The `df_past` data frame should contain the disease history of each person that is to be predicted. \cr
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
  df_pred = supps$df_pred
  df_past = supps$df_past
  param.type = supps$param.type

  ### check argument specification
  lctmc = object
  if (!(is.data.frame(df_pred) && is.data.frame(df_past))) {
    stop("`df_pred` and `df_past` both must be of class 'data.frame'")
  }
  if (length(unique(df_pred$id)) != nrow(df_pred)) {
    stop("Each ID number in `df_pred` must only appear at max one time")
  }
  if (!all(df_pred$id %in% df_past$id)) {
    stop("All IDs that appear in `df_pred` should also appear in the `df_past` data frame")
  } else {
    temp = merge(df_pred, df_past, by = "id", all.x = TRUE, sort = FALSE, suffixes = c(".pred", ".past"))
    temp = (temp$obsTime.pred - temp$obsTime.past) > 0
    if (!all(temp)) {
      stop("observation times in `df_pred` must be larger greater than those in `df_past`")
    }
  }
  if (!all(c("id", "obsTime") %in% colnames(df_pred)) || (ncol(df_pred) != 2)) {
    stop("`df_pred` should contain ONLY two columns: 'id' and 'obsTime'")
  }
  if (!all(c("id", "obsTime") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'id' and 'obsTime'")
  } else if (all(c("id", "obsTime") %in% colnames(df_past))) {
    temp = unique(df_past[c("id", "obsTime")])
    temp = as.numeric(table(temp$id)) >= 2
    if (!all(temp)) {
      stop("there should be at minimal two observation (one transition) per ID in `df_past`")
    }
  } else {
    if (!all(c("state_at_obsTime", "x1", "x2", "w1", "w2") %in% colnames(df_past))) {
      stop("`df_past` should contain columns 'x1', 'x2', 'w1', 'w2', 'state_at_obsTime'")
    }
  }
  if (!all(df_past$state_at_obsTime %in% c(1, 2))) {
    stop("for a 2x2 model, 'state_at_obsTime' should only take values 1 or 2")
  }


  ### only predict on IDs common in both input data frames
  common_ID = intersect(df_pred$id, df_past$id)
  df_past = df_past[df_past$id %in% common_ID, ]
  df_pred = df_pred[df_pred$id %in% common_ID, ]

  ### sort
  df_past = df_past[order(df_past$id, df_past$obsTime), ]
  df_pred = df_pred[order(df_pred$id, df_pred$obsTime), ]

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
  if (param.type == "mle") {
    model_param = lctmc$SE$SE$mle_theta
    names(model_param) = lctmc$SE$SE$names
  } else if (param.type == "kmeans") {
    model_param = lctmc$init01$step1_full$theta
  } else {
    stop("`param.type` should be either 'mle' or 'kmeans'")
  }

  ### compute `bik`
  bik = bik_all_2x2(
    theta = model_param,
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

  ### output (1)
  latent_class = cbind(
    data.frame(id = my_df1.IDs, pred_class = pred_class),
    as.data.frame(P.z)
  )


  ### data frame of latent class probability (long format id x z)
  pivot_vars = paste("P.z", 1:lctmc$K, sep = "")
  df.prob_lc = latent_class[!colnames(latent_class) %in% c("pred_class")]
  df.prob_lc = lapply(
    X = pivot_vars,
    FUN = function(v) {
      cbind(
        df.prob_lc[!colnames(df.prob_lc) %in% pivot_vars],
        data.frame(z = sub("P\\.z", "", v), P.z = df.prob_lc[[v]])
      )
    }
  )
  df.prob_lc = do.call(`rbind`, df.prob_lc)

  ### a df of last known observation for each person (from `df_past`)
  my_df2 = Map(f = function(d) unique(d[d$obsTime == max(d$obsTime), ]), split(x = df_past, f = df_past$id))
  my_df2 = do.call(`rbind`, my_df2)
  my_df2$x0 = my_df2$w0 = 1
  my_df2 = merge(my_df2, df_pred, by = 'id', sort = FALSE, all.x = TRUE, suffixes = c(".past", ".pred"))

  ### data objects needed to compute probabilities
  my_df2.IDs = my_df2$id
  my_df2.Xmat = as.matrix(my_df2[c("x0", "x1", "x2")])
  my_df2.Wmat = as.matrix(my_df2[c("w0", "w1", "w2")])
  my_df2.dt = my_df2$obsTime.pred - my_df2$obsTime.past

  ### previous state (CTMC, future depend on closest past)
  prev_state = my_df2$state_at_obsTime
  names(prev_state) = my_df2$id

  ### get transition probability
  P.rs_k = bik_all_2x2(
    theta = model_param,
    data = my_df2,
    Xmat = my_df2.Xmat,
    Wmat = my_df2.Wmat,
    dt = my_df2.dt,
    P.rs = TRUE,
    K = lctmc$K,
    theta.names = gen_theta_names(K = lctmc$K, type = "2x2", purpose = "bik")
  )

  ### transition probabilities convert into data frame format
  df.P.rs_k = NULL
  for (i in seq_along(P.rs_k)) {
    ## i^th class probability
    d = P.rs_k[[i]]

    ## ID col + transition prob in columns + previous state vector
    d = cbind(id = names(prev_state), prev_state = prev_state, d)
    rownames(d) = NULL

    ## depending on previous state, probably to future state is different
    for (future_state in 1:2) {
      # new variable name
      new_v = paste("P.ds", future_state, sep = "")
      # which transition probability column to use
      p_from1_v = paste("P1", future_state, sep = "")
      p_from2_v = paste("P2", future_state, sep = "")
      # new column
      d[[new_v]] = ifelse(d$prev_state == 1, d[[p_from1_v]], d[[p_from2_v]])
    }

    ## append + `z` is an indicator for current class in long format data frame
    d$z = as.character(i)
    df.P.rs_k = rbind(df.P.rs_k, d)
  }
  df.P.rs_k = df.P.rs_k[order(df.P.rs_k$id, df.P.rs_k$z), ]

  ### join with latent class probability (these act as prediction weights)
  df.P.rs_k = merge(df.P.rs_k, df.prob_lc, by = c("id", "z"), all.x = TRUE, sort = FALSE)
  df.P.rs_k$P.ds1.z = df.P.rs_k$P.ds1 * df.P.rs_k$P.z
  df.P.rs_k$P.ds2.z = df.P.rs_k$P.ds2 * df.P.rs_k$P.z


  ### prediction probability for the disease state
  df.P.rs = split(x = df.P.rs_k, f = df.P.rs_k$id)
  df.P.rs = Map(
    f = function(d) {
      data.frame(id = unique(d$id),
                 P.trans_to1 = sum(d$P.ds1.z),
                 P.trans_to2 = sum(d$P.ds2.z))
    },
    df.P.rs
  )
  df.P.rs = do.call(`rbind`, df.P.rs)


  ### out of the possible disease states, whats the largest probability mass?
  P.rs = list(P1 = df.P.rs$P.trans_to1, P2 = df.P.rs$P.trans_to2)
  max_probs2 = Reduce(f = function(u, v) ifelse(u >= v, u, v), P.rs)

  ### disease state prediction is based on the state with the highest probability mass
  pred_state = Map(f = function(x) x == max_probs2, P.rs)
  for (i in seq_along(pred_state)) {
    pred_state[[i]] = i * as.numeric(pred_state[[i]])
  }
  df.P.rs$pred_state = Reduce(`+`, pred_state)

  ### output (2)
  disease_states = merge(df_pred, df.P.rs, by = "id", all.x = TRUE, sort = FALSE)
  lead_cols = c("id", "obsTime", "pred_state")
  trail_cols = colnames(disease_states)[!colnames(disease_states) %in% lead_cols]
  disease_states = disease_states[c(lead_cols, trail_cols)] # re-order columns

  ### return
  return(list(df_pred.lc = latent_class, df_pred.ds = disease_states))
}

#' @rdname predict
#' @exportS3Method
predict.lctmc_3x3 = function(object, ...) {
  ### unpack `...`
  supps = list(...)
  df_pred = supps$df_pred
  df_past = supps$df_past
  param.type = supps$param.type

  ### check argument specification
  lctmc = object
  if (!(is.data.frame(df_pred) && is.data.frame(df_past))) {
    stop("`df_pred` and `df_past` both must be of class 'data.frame'")
  }
  if (length(unique(df_pred$id)) != nrow(df_pred)) {
    stop("Each ID number in `df_pred` must only appear at max one time")
  }
  if (!all(df_pred$id %in% df_past$id)) {
    stop("All IDs that appear in `df_pred` should also appear in the `df_past` data frame")
  } else {
    temp = merge(df_pred, df_past, by = "id", all.x = TRUE, sort = FALSE, suffixes = c(".pred", ".past"))
    temp = (temp$obsTime.pred - temp$obsTime.past) > 0
    if (!all(temp)) {
      stop("observation times in `df_pred` must be larger greater than those in `df_past`")
    }
  }
  if (!all(c("id", "obsTime") %in% colnames(df_pred)) || (ncol(df_pred) != 2)) {
    stop("`df_pred` should contain ONLY two columns: 'id' and 'obsTime'")
  }
  if (!all(c("id", "obsTime") %in% colnames(df_past))) {
    stop("`df_past` should contain columns 'id' and 'obsTime'")
  } else if (all(c("id", "obsTime") %in% colnames(df_past))) {
    temp = unique(df_past[c("id", "obsTime")])
    temp = as.numeric(table(temp$id)) >= 2
    if (!all(temp)) {
      stop("there should be at minimal two observation (one transition) per ID in `df_past`")
    }
  } else {
    if (!all(c("state_at_obsTime", "x1", "x2", "w1", "w2") %in% colnames(df_past))) {
      stop("`df_past` should contain columns 'x1', 'x2', 'w1', 'w2', 'state_at_obsTime'")
    }
  }
  if (!all(df_past$state_at_obsTime %in% c(1, 2, 3))) {
    stop("for a 3x3 model, 'state_at_obsTime' should only take values 1, 2, or 3")
  }


  ### only predict on IDs common in both input data frames
  common_ID = intersect(df_pred$id, df_past$id)
  df_past = df_past[df_past$id %in% common_ID, ]
  df_pred = df_pred[df_pred$id %in% common_ID, ]

  ### sort
  df_past = df_past[order(df_past$id, df_past$obsTime), ]
  df_pred = df_pred[order(df_pred$id, df_pred$obsTime), ]

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
  if (param.type == "mle") {
    model_param = lctmc$SE$SE$mle_theta
    names(model_param) = lctmc$SE$SE$names
  } else if (param.type == "kmeans") {
    model_param = lctmc$init01$step1_full$theta
  } else {
    stop("`param.type` should be either 'mle' or 'kmeans'")
  }

  ### compute `bik`
  bik = bik_all_3x3(
    theta = model_param,
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

  ### output (1)
  latent_class = cbind(
    data.frame(id = my_df1.IDs, pred_class = pred_class),
    as.data.frame(P.z)
  )


  ### data frame of latent class probability
  pivot_vars = paste("P.z", 1:lctmc$K, sep = "")
  df.prob_lc = latent_class[!colnames(latent_class) %in% c("pred_class")]
  df.prob_lc = lapply(
    X = pivot_vars,
    FUN = function(v) {
      cbind(
        df.prob_lc[!colnames(df.prob_lc) %in% pivot_vars],
        data.frame(z = sub("P\\.z", "", v), P.z = df.prob_lc[[v]])
      )
    }
  )
  df.prob_lc = do.call(`rbind`, df.prob_lc)

  ### a df of last known observation for each person (from `df_past`)
  my_df2 = Map(f = function(d) unique(d[d$obsTime == max(d$obsTime), ]), split(x = df_past, f = df_past$id))
  my_df2 = do.call(`rbind`, my_df2)
  my_df2$x0 = my_df2$w0 = 1
  my_df2 = merge(my_df2, df_pred, by = 'id', sort = FALSE, all.x = TRUE, suffixes = c(".past", ".pred"))

  ### data objects needed to compute probabilities
  my_df2.Xmat = as.matrix(my_df2[c("x0", "x1", "x2")])
  my_df2.Wmat = as.matrix(my_df2[c("w0", "w1", "w2")])
  my_df2.dt = my_df2$obsTime.pred - my_df2$obsTime.past

  ### previous state (CTMC, future depend on closest past)
  prev_state = my_df2$state_at_obsTime
  names(prev_state) = my_df2$id

  ### get transition probability
  P.rs_k = bik_all_3x3(
    theta = model_param,
    data = my_df2,
    Xmat = my_df2.Xmat,
    Wmat = my_df2.Wmat,
    dt = my_df2.dt,
    P.rs = TRUE,
    K = lctmc$K,
    theta.names = gen_theta_names(K = lctmc$K, type = "3x3", purpose = "bik")
  )

  ### transition probabilities convert into data frame format
  df.P.rs_k = NULL
  for (i in seq_along(P.rs_k)) {
    ## i^th class probability
    d = P.rs_k[[i]]

    ## ID col + transition prob in columns + previous state vector
    d = cbind(id = names(prev_state), prev_state = prev_state, d)
    rownames(d) = NULL

    ## depending on previous state, probably to future state is different
    for (future_state in 1:3) {
      # new variable name
      new_v = paste("P.ds", future_state, sep = "")
      # which transition probability column to use
      p_from1_v = paste("P1", future_state, sep = "")
      p_from2_v = paste("P2", future_state, sep = "")
      # new column
      d[[new_v]] = ifelse(d$prev_state == 1, d[[p_from1_v]], d[[p_from2_v]])
      d[[new_v]][d$prev_state == 3] = as.numeric(future_state == 3)
    }

    ## append + `z` is an indicator for current class in long format data frame
    d$z = as.character(i)
    df.P.rs_k = rbind(df.P.rs_k, d)
  }
  df.P.rs_k = df.P.rs_k[order(df.P.rs_k$id, df.P.rs_k$z), ]

  ### join with latent class probability (these act as prediction weights)
  df.P.rs_k = merge(df.P.rs_k, df.prob_lc, by = c("id", "z"), all.x = TRUE, sort = FALSE)
  df.P.rs_k$P.ds1.z = df.P.rs_k$P.ds1 * df.P.rs_k$P.z
  df.P.rs_k$P.ds2.z = df.P.rs_k$P.ds2 * df.P.rs_k$P.z
  df.P.rs_k$P.ds3.z = df.P.rs_k$P.ds3 * df.P.rs_k$P.z

  ### prediction probability for the disease state
  df.P.rs = split(x = df.P.rs_k, f = df.P.rs_k$id)
  df.P.rs = Map(
    f = function(d) {
      data.frame(id = unique(d$id),
                 P.trans_to1 = sum(d$P.ds1.z),
                 P.trans_to2 = sum(d$P.ds2.z),
                 P.trans_to3 = sum(d$P.ds3.z))
    },
    df.P.rs
  )
  df.P.rs = do.call(`rbind`, df.P.rs)

  ### out of the possible disease states, whats the largest probability mass?
  P.rs = list(P1 = df.P.rs$P.trans_to1, P2 = df.P.rs$P.trans_to2, P3 = df.P.rs$P.trans_to3)
  max_probs2 = Reduce(f = function(u, v) ifelse(u >= v, u, v), P.rs)

  ### disease state prediction is based on the state with the highest probability mass
  pred_state = Map(f = function(x) x == max_probs2, P.rs)
  for (i in seq_along(pred_state)) {
    pred_state[[i]] = i * as.numeric(pred_state[[i]])
  }
  df.P.rs$pred_state = Reduce(`+`, pred_state)

  ### output (2)
  disease_states = merge(df_pred, df.P.rs, by = "id", all.x = TRUE, sort = FALSE)
  lead_cols = c("id", "obsTime", "pred_state")
  trail_cols = colnames(disease_states)[!colnames(disease_states) %in% lead_cols]
  disease_states = disease_states[c(lead_cols, trail_cols)] # re-order columns

  ### return
  return(list(df_pred.lc = latent_class, df_pred.ds = disease_states))
}
