# 2x2 case ------------------
# transition rates
q12 = c(1, 0.5, 0.1)
q21 = c(1, 2, 3)

# time interval between observations
time_diffs = c(1, 5, 10)

# compute
LCTMC:::get_P_2x2(
  q12 = q12,
  q21 = q21,
  dt = time_diffs
)


# 3x3 case ------------------
# transition rates
q12 = c(1, 0.5, 0.1)
q21 = c(1, 2, 3)
q23 = c(1, 1, 1)

# time interval between observations
time_diffs = c(1, 5, 10)

# compute
LCTMC:::get_P_3x3(
  q12 = q12,
  q21 = q21,
  q23 = q23,
  dt = time_diffs
)
