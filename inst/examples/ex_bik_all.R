# This is a 2x2 case example
# generate data
my_ids = c(1,1,1, 2,2, 3,3,3)
c1 = c(0,1,1.5, 0,1, 0,2,5)
c2 = c(1,2,2, 1,1, 2,1,2)
x1 = c(rep(rnorm(n = 1), 3), rep(rnorm(n = 1), 2), rep(rnorm(n = 1), 3))
x2 = c(rep(1, 3), rep(0, 2), rep(1, 3))
w1 = c(rep(runif(n = 1), 3), rep(runif(n = 1), 2), rep(runif(n = 1), 3))
w2 = c(rep(0, 3), rep(0, 2), rep(1, 3))
latent_class = c(1,1,1, 1,1 ,2,2)

# create data frame and matrix objects
my_df = data.frame(
  id = my_ids,
  obsTime = c1, state_at_obsTime = c2,
  x1 = x1, x2 = x2,
  w1 = w1, w2 = w2
)

my_df2 = LCTMC:::fmt_rowwise_trans(data = my_df, type = "2x2")
xmat = as.matrix(my_df2[c("x0", "x1", "x2")])
wmat = unique(my_df2[c("id", "w0", "w1", "w2")])
wmat = as.matrix(wmat[!colnames(wmat) %in% c("id")])

# - # - # - # - # - # - # - # - # -

# create theta vector ~ case of 2 latent classes
my_theta = c("alpha0.1" = 1.0, "alpha1.1" = 1.0, "alpha2.1" = 1.0,
             "beta0.12_1" = 1.1, "beta1.12_1" = 1.1, "beta2.12_1" = 1.1,
             "beta0.21_1" = 1.5, "beta1.21_1" = 1.5, "beta2.21_1" = 1.5,
             "beta0.12_2" = 1.8, "beta1.12_2" = 1.8, "beta2.12_2" = 1.8,
             "beta0.21_2" = 0.3, "beta1.21_2" = 0.3, "beta2.21_2" = 0.3)

# use bik_all_2x2() ~ case of 2 latent classes
LCTMC:::bik_all_2x2(
  theta = my_theta,
  data = my_df2, dt = my_df2$dt,
  Xmat = xmat,
  Wmat = wmat,
  P.rs = FALSE,
  K = 2L,
  theta.names = gen_theta_names(K = 2L, type = "2x2", purpose = "bik")
)

# - # - # - # - # - # - # - # - # -

# create theta vector ~ case of 3 latent classes
my_theta = c("alpha0.1" = 1.0, "alpha1.1" = 1.0, "alpha2.1" = 1.0,
             "alpha0.2" = 0.5, "alpha1.2" = 0.2, "alpha2.2" = 0.2,
             "beta0.12_1" = 1.1, "beta1.12_1" = 1.1, "beta2.12_1" = 1.1,
             "beta0.21_1" = 1.5, "beta1.21_1" = 1.5, "beta2.21_1" = 1.5,
             "beta0.12_2" = 1.8, "beta1.12_2" = 1.8, "beta2.12_2" = 1.8,
             "beta0.21_2" = 0.3, "beta1.21_2" = 0.3, "beta2.21_2" = 0.3,
             "beta0.12_3" = 1.4, "beta1.12_3" = 1.4, "beta2.12_3" = 1.4,
             "beta0.21_3" = 0.4, "beta1.21_3" = 0.4, "beta2.21_3" = 0.4)

# use bik_all_2x2() ~ case of 3 latent classes
LCTMC:::bik_all_2x2(
  theta = my_theta,
  data = my_df2, dt = my_df2$dt,
  Xmat = xmat,
  Wmat = wmat,
  P.rs = FALSE,
  K = 3L,
  theta.names = gen_theta_names(K = 3L, type = "2x2", purpose = "bik")
)

# - # - # - # - # - # - # - # - # -

# create theta vector ~ case of 4 latent classes
my_theta = c("alpha0.1" = 1.0, "alpha1.1" = 1.0, "alpha2.1" = 1.0,
             "alpha0.2" = 0.5, "alpha1.2" = 0.2, "alpha2.2" = 0.2,
             "alpha0.3" = -0.5, "alpha1.3" = -1.0, "alpha2.3" = -0.2,
             "beta0.12_1" = 1.1, "beta1.12_1" = 1.1, "beta2.12_1" = 1.1,
             "beta0.21_1" = 1.5, "beta1.21_1" = 1.5, "beta2.21_1" = 1.5,
             "beta0.12_2" = 1.8, "beta1.12_2" = 1.8, "beta2.12_2" = 1.8,
             "beta0.21_2" = 0.3, "beta1.21_2" = 0.3, "beta2.21_2" = 0.3,
             "beta0.12_3" = 1.4, "beta1.12_3" = 1.4, "beta2.12_3" = 1.4,
             "beta0.21_3" = 0.4, "beta1.21_3" = 0.4, "beta2.21_3" = 0.4,
             "beta0.12_4" = 0.9, "beta1.12_4" = 0.9, "beta2.12_4" = 0.9,
             "beta0.21_4" = 0.8, "beta1.21_4" = 0.8, "beta2.21_4" = 0.8)

# use bik_all_2x2() ~ case of 4 latent classes
LCTMC:::bik_all_2x2(
  theta = my_theta,
  data = my_df2, dt = my_df2$dt,
  Xmat = xmat,
  Wmat = wmat,
  P.rs = FALSE,
  K = 4L,
  theta.names = gen_theta_names(K = 4L, type = "2x2", purpose = "bik")
)
