# load data ---------
data("df_to_predict", package = "LCTMC")
data("model_2x2", package = "LCTMC")

# data creation ---------
  ### disease history ~
   ## get df past
  my_df_past = df_to_predict$`2x2`
  my_df_past = split(x = my_df_past, f = my_df_past$id)

  ### disease history ~
   ## filter out the last observation for each ID number
  my_df_past = lapply(
    X = my_df_past,
    FUN = function(d) {
      d[d$obsTime < max(d$obsTime), ]
    }
  )
  my_df_past = do.call(`rbind`, my_df_past)
  rownames(my_df_past) = NULL


  ### future time point to predict on ~
   ## create data frame
  my_df_predict = df_to_predict$`2x2`[c("id", "obsTime", "state_at_obsTime")]
  my_df_predict = split(x = my_df_predict, f = my_df_predict$id)

  ### future time point to predict on ~
   ## filter out the last observation for each ID number
  my_df_predict = lapply(
    X = my_df_predict,
    FUN = function(d) {
      d[d$obsTime == max(d$obsTime), ]
    }
  )
  my_df_predict.answer = do.call(`rbind`, my_df_predict)
  my_df_predict = my_df_predict.answer[c("id", "obsTime")]
  rownames(my_df_predict.answer) = rownames(my_df_predict) = NULL


# predict --------
  ### run function
  my_predictions = predict(
    object = model_2x2,
    df_past = my_df_past, df_pred = my_df_predict,
    param.type = "mle"
  )


  ### check result ~
  ## disease state
  check_result2 = merge(
    my_predictions$df_pred.ds[c("id", "pred_state")],
    my_df_predict.answer[c("id", "state_at_obsTime")],
    by = "id"
  )
  table(check_result2$state_at_obsTime, check_result2$pred_state)
  mean(check_result2$state_at_obsTime == check_result2$pred_state)


  ### check result ~
   ## latent class
  check_result1 = merge(
    my_predictions$df_pred.lc[c("id", "pred_class")],
    unique(df_to_predict$`2x2`[c("id", "latent_class")]),
    by = "id"
  )
  table(check_result1$pred_class, check_result1$latent_class)
  mean(check_result1$pred_class == check_result1$latent_class)

# --------------------------------------------------------------------------- #
# ## NOTE: #### #### #### #
# ##  * Even though the cross-table shows 0% accuracy, it is actually 100%.
# ##  * As latent class's labels are meaningless,
# ##     it only matters if the corresponding transition rates are aligned
# ##     (see `LCTMC::align_MLE_2x2` for more info)
# ##  * the lctmc_2x2() and lctmc_3x3() functions will always
# ##     assign the class with greatest `alpha0` estimate as class 1
# ##     and the class with smallest `alpha0` as the last class, `K`.
# --------------------------------------------------------------------------- #
