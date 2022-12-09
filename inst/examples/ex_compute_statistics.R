# load data
data("model_2x2", package = "LCTMC")
data("model_3x3", package = "LCTMC")


# Q matrix ---------
# 2x2
compute_Q(m = model_2x2, x1 = 1, x2 = 1)
# 3x3
compute_Q(m = model_3x3, x1 = 1, x2 = 1)


# Sojourn Times ---------
# 2x2
compute_sojourn(m = model_2x2, x1 = 1, x2 = 1)
# 3x3
compute_sojourn(m = model_3x3, x1 = 1, x2 = 1)


# Rate Ratios ---------
# 2x2
compute_RR(m = model_2x2, x1 = 1, x2 = 1, ci_method = 'delta')
# 3x3
compute_RR(m = model_3x3, x1 = 1, x2 = 1, ci_method = 'delta')


# Transition Probability ---------
# 2x2
compute_P(m = model_2x2, x1 = 1, x2 = 1, dt = 30)
# 3x3
compute_P(m = model_3x3, x1 = 1, x2 = 1, dt = 30)



# Transition Probability when t --> infty
lapply(
  X = compute_P(m = model_3x3, x1 = 1, x2 = 1, dt = 999999),
  FUN = function(x) round(x, 10)
)


# Transition Probability when t = 0
lapply(
  X = compute_P(m = model_3x3, x1 = 1, x2 = 1, dt = 0),
  FUN = function(x) round(x, 10)
)
