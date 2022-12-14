# generate data
my_ids = c(1, 1, 1,
           2, 2,
           3, 3, 3)
c1 = c(0, 1, 1.5,
       0, 1,
       0, 2, 5)
c2 = c(1, 2, 2,
       1, 1,
       2, 1, 3)
x1 = c(rep(rnorm(n = 1), 3), rep(rnorm(n = 1), 2), rep(rnorm(n = 1), 3))
x2 = c(rep(1, 3), rep(0, 2), rep(1, 3))
w1 = c(rep(runif(n = 1), 3), rep(runif(n = 1), 2), rep(runif(n = 1), 3))
w2 = c(rep(0, 3), rep(0, 2), rep(1, 3))

# create data frame and matrix objects
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
)
xmat = my_df2$Xmat
wmat = my_df2$Wmat
dt = my_df2$dt
my_df2 = my_df2$df_trans

# - # - # - # - # - # - # - # - # -

# create theta vector
my_theta = c("alpha0.1" = 1.0, "alpha1.1" = 1.0, "alpha2.1" = 1.0,
             "alpha0.2" = 0.5, "alpha1.2" = 0.2, "alpha2.2" = 0.2,
             "beta0.12_1" = 1.1, "beta1.12_1" = 1.1, "beta2.12_1" = 1.1,
             "beta0.21_1" = 1.5, "beta1.21_1" = 1.5, "beta2.21_1" = 1.5,
             "beta0.12_2" = 1.8, "beta1.12_2" = 1.8, "beta2.12_2" = 1.8,
             "beta0.21_2" = 0.3, "beta1.21_2" = 0.3, "beta2.21_2" = 0.3,
             "beta0.12_3" = 1.4, "beta1.12_3" = 1.4, "beta2.12_3" = 1.4,
             "beta0.21_3" = 0.4, "beta1.21_3" = 0.4, "beta2.21_3" = 0.4)

# use bik_all_2x2()
my_b = LCTMC:::bik_all_2x2(
  theta = my_theta,
  data = my_df2,
  Xmat = xmat,
  Wmat = wmat,
  dt = dt,
  P.rs = FALSE,
  K = 3L,
  theta.names = gen_theta_names(K = 3L, type = "2x2", purpose = "bik")
)

# manually create negative and zero values
my_b$bi1[2] = 0
my_b$bi1[3] = -1
my_b$bi3[1:3] = c(0, -1 , 0)

print(my_b)                          # before imputation
my_b = LCTMC:::impute_bik(x = my_b)  # imputes
print(my_b)                          # after imputation
