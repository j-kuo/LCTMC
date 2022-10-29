# 2x2 case -------------
# generate P
q12 = q21 = c(1.1, 0.5, 0.1, 2.2, 0.3)
time_diffs = c(1, 2, 5.5, 8.1, 13.5)
my_P = LCTMC:::get_P_2x2(q12 = q12, q21 = q21, dt = time_diffs)

# generate data
my_ids = c(1,1,1, 2,2, 3,3,3)

c1 = c(0,1,1.5, 0,1, 0,2,5)
c2 = c(1,2,2, 1,1, 2,1,2)

x1 = c(rnorm(n = 3), rnorm(n = 2), rnorm(n = 3))
x2 = c(1,1,1, 0,0, 1,1,1)

w1 = c(rexp(n = 3, rate = 1), rexp(n = 2, rate = 1), rexp(n = 3, rate = 1))
w2 = c(rexp(n = 3, rate = 10), rexp(n = 2, rate = 10), rexp(n = 3, rate = 10))

my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)
my_df2 = LCTMC:::fmt_rowwise_trans(
  data = my_df,
  type = "2x2",
  X_names = c("x0", "x1", "x2"),
  W_names = c("w0", "w1", "w2"),
  scaling = 1
)$df_trans

# use LCTMC:::Li_2x2()
LCTMC:::Li_2x2(P = my_P, data = my_df2)


# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - #


# 3x3 case -------------
# generate P
q12 = q21 = q23 = c(1.1, 0.5, 0.1, 2.2, 0.3)
time_diffs = c(1, 2, 5.5, 8.1, 13.5)

my_P = LCTMC:::get_P_3x3(q12 = q12, q21 = q21, q23 = q23, dt = time_diffs)

# generate data
my_ids = c(1,1,1, 2,2, 3,3,3)

c1 = c(0,1,1.5, 0,1, 0,2,5)
c2 = c(1,2,3, 1,1, 1,2,3)

x1 = c(rnorm(n = 3), rnorm(n = 2), rnorm(n = 3))
x2 = c(1,1,1, 0,0, 1,1,1)

w1 = c(rexp(n = 3, rate = 1), rexp(n = 2, rate = 1), rexp(n = 3, rate = 1))
w2 = c(rexp(n = 3, rate = 10), rexp(n = 2, rate = 10), rexp(n = 3, rate = 10))

my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)
my_df2 = LCTMC:::fmt_rowwise_trans(
  data = my_df,
  type = "3x3",
  X_names = c("x0", "x1", "x2"),
  W_names = c("w0", "w1", "w2"),
  scaling = 1
)$df_trans

# use Li_3x3()
LCTMC:::Li_3x3(P = my_P, data = my_df2)
