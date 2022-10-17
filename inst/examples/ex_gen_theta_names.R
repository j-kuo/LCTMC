# 2x2 case with 2 latent classes
theta_names.2x2.K2 = list(
  # alphas, class 1
  c("alpha0.1", "alpha1.1", "alpha2.1"),
  # betas, class 1
  c("beta0.12_1", "beta1.12_1", "beta2.12_1"),
  c("beta0.21_1", "beta1.21_1", "beta2.21_1"),
  # betas, class 2
  c("beta0.12_2", "beta1.12_2", "beta2.12_2"),
  c("beta0.21_2", "beta1.21_2", "beta2.21_2")
)

# 2x2 case with 3 latent classes
theta_names.2x2.K3 = list(
  # alphas, class 1
  c("alpha0.1", "alpha1.1", "alpha2.1"),
  # alphas, class 1
  c("alpha0.2", "alpha1.2", "alpha2.2"),
  # betas, class 1
  c("beta0.12_1", "beta1.12_1", "beta2.12_1"),
  c("beta0.21_1", "beta1.21_1", "beta2.21_1"),
  # betas, class 2
  c("beta0.12_2", "beta1.12_2", "beta2.12_2"),
  c("beta0.21_2", "beta1.21_2", "beta2.21_2"),
  # betas, class 3
  c("beta0.12_3", "beta1.12_3", "beta2.12_3"),
  c("beta0.21_3", "beta1.21_3", "beta2.21_3")
)

# print TRUE
identical(
  theta_names.2x2.K2,
  gen_theta_names(K = 2L, type = "2x2", purpose = "em")
)
identical(
  theta_names.2x2.K3,
  gen_theta_names(K = 3L, type = "2x2", purpose = "em")
)

# 3x3 case with 2 latent classes
theta_names.3x3.K2 = list(
  # alphas, class 1
  c("alpha0.1", "alpha1.1", "alpha2.1"),
  # betas, class 1
  c("beta0.12_1", "beta1.12_1", "beta2.12_1"),
  c("beta0.21_1", "beta1.21_1", "beta2.21_1"),
  c("beta0.23_1", "beta1.23_1", "beta2.23_1"),
  # betas, class 2
  c("beta0.12_2", "beta1.12_2", "beta2.12_2"),
  c("beta0.21_2", "beta1.21_2", "beta2.21_2"),
  c("beta0.23_2", "beta1.23_2", "beta2.23_2")
)

# prints TRUE
identical(
  theta_names.3x3.K2,
  gen_theta_names(K = 2L, type = "3x3", purpose = "em")
)
