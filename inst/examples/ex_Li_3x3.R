# generate P (see documentation for get_P_3x3())
q12 = c(1, 0.5, 0.1, 2.2, 3)
q21 = c(1, 2, 3, 2, 1)
q23 = c(1, 1, 1, 1, 1)
time_diffs = c(1, 2, 5.5, 8.1, 13.5)

my_P = LCTMC:::get_P_3x3(q12 = q12, q21 = q21, q23 = q23, dt = time_diffs)

# - # - # - # - # - # -

# generate data (see documentation for fmt_rowwise_3x3trans())
my_ids = c(1,1,1, 2,2, 3,3,3)

c1 = c(0,1,1.5, 0,1, 0,2,5)
c2 = c(1,2,3, 1,1, 1,2,3)
x1 = c(rnorm(n = 3), rnorm(n = 2), rnorm(n = 3))
x2 = c(1,1,1, 0,0, 1,1,1)
w1 = c(rexp(n = 3, rate = 1), rexp(n = 2, rate = 1), rexp(n = 3, rate = 1))
w2 = c(rexp(n = 3, rate = 10), rexp(n = 2, rate = 10), rexp(n = 3, rate = 10))
latent_class = c(1,1,1, 1,1 ,2,2)

my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)
my_df2 = fmt_rowwise_3x3trans(data = my_df)

# - # - # - # - # - # -

# use Li_3x3()
Li_3x3(P = my_P, data = my_df2)
