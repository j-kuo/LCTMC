# id column
my_ids = c(1, 1, 1,
           2, 2,
           3, 3, 3)

# observation times
c1 = c(0, 1, 1.5,
       0, 1,
       0, 2, 5)

# state at observation time
c2 = c(1, 2, 2,
       1, 1,
       2, 1, 3)

# some covariate data
x1 = c(rnorm(n = 3), rnorm(n = 2), rnorm(n = 3))
x2 = c(rep(1, 3), rep(0, 2), rep(1, 3))
w1 = c(rexp(n = 3, rate = 1), rexp(n = 2, rate = 1), rexp(n = 3, rate = 1))
w2 = c(rexp(n = 3, rate = 10), rexp(n = 2, rate = 10), rexp(n = 3, rate = 10))

# create data frame
my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)

# scaling parameters
my_scaling = c(x0 = 1, x1 = 0.01, x2 = 5, w0 = 1, w1 = 1000, w2 = 1, dt = 100)

# use function ~ 2x2 variation ~ note the 1-3 transition is not captured
LCTMC:::fmt_rowwise_trans(
  data = my_df,
  type = "2x2",
  X_names = c("x0", "x1", "x2"),
  W_names = c("w0", "w1", "w2"),
  scaling = my_scaling
)

# use function ~ 3x3 variation
LCTMC:::fmt_rowwise_trans(
  data = my_df,
  type = "3x3",
  X_names = c("x0", "x1", "x2"),
  W_names = c("w0", "w1", "w2"),
  scaling = my_scaling
)
