# id column
my_ids = c(1,1,1, 2,2, 3,3,3)

# observation times
c1 = c(0,1,1.5, 0,1, 0,2,5)

# state at observation time
c2 = c(1,2,2, 1,1, 2,1,2)

# some covariate data
x1 = c(rnorm(n = 3), rnorm(n = 2), rnorm(n = 3))
x2 = c(1,1,1, 0,0, 1,1,1)
w1 = c(rexp(n = 3, rate = 1), rexp(n = 2, rate = 1), rexp(n = 3, rate = 1))
w2 = c(rexp(n = 3, rate = 10), rexp(n = 2, rate = 10), rexp(n = 3, rate = 10))

# create data frame
my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)

# use function
fmt_rowwise_2x2trans(data = my_df)
