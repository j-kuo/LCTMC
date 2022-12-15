# load data
data("df_to_predict", package = "LCTMC")
data("model_2x2", package = "LCTMC")

# create data sets (full data)
my_df_past.full = df_to_predict$`2x2`


# create data sets (truncated, 3 obs per person)
my_df_past.only3 = lapply(
  X = split(x = my_df_past.full, f = my_df_past.full$id),
  FUN = function(d) {
    d = d[1:3, ]
    d
  }
)
my_df_past.only3 = do.call(`rbind`, my_df_past.only3)


# create data sets (truncated, 6 obs per person)
my_df_past.only6 = lapply(
  X = split(x = my_df_past.full, f = my_df_past.full$id),
  FUN = function(d) {
    d = d[1:6, ]
    d
  }
)
my_df_past.only6 = do.call(`rbind`, my_df_past.only6)


# predict (1) ~ only 3 obs per person in disease history
my_predictions = predict(
  object = model_2x2,
  df_past = my_df_past.only3,
  pred_threshold = 0.7
)
table(my_predictions$pred_class, my_predictions$pred_class_strict)


# predict (2) ~ only 6 obs per person in disease history
my_predictions = predict(
  object = model_2x2,
  df_past = my_df_past.only6,
  pred_threshold = 0.7
)
table(my_predictions$pred_class, my_predictions$pred_class_strict)


# predict (3) ~ max obs per person in disease history
my_predictions = predict(
  object = model_2x2,
  df_past = my_df_past.full,
  pred_threshold = 0.7
)
table(my_predictions$pred_class, my_predictions$pred_class_strict)


# predict (4) ~ extreme threshold criteria, leading to all NA predictions
my_predictions = predict(
  object = model_2x2,
  df_past = my_df_past.only6,
  pred_threshold = 0.9999
)
table(my_predictions$pred_class, my_predictions$pred_class_strict)
