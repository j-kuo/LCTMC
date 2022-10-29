# The example below demonstrates how to specify the 'lctmc' functions correctly
#   to perform the latent class CTMC model.
# Both the 2x2 and 3x3 examples are listed below, however it should be noted
#   that the models may take a while to run (approx. 10-30min on a high spec laptop)
# In both examples, we set the number of latent classes, K, to be 2.
# In the 2x2 model, we constrained the model so that the covariate effects on
#   the latent class probability are equal to 0.
# In the 3x3 model, we left the model un-constrained.


\dontrun{
  ## this is a 2x2 example, with 3 latent classes
  data("example_df2x2", package = "LCTMC")
  model_2x2 = LCTMC::lctmc_2x2(
    # data
    data = example_df2x2,
    # any scaling transformation
    dt_scale = c(dt = 1/365),
    x_scale = c(x0 = 1, x1 = 1/7, x2 = 1),
    w_scale = c(w0 = 1, w1 = 1/9, w2 = 1),
    # general model specification
    K = 2L,
    X_names = c('x0', 'x1', 'x2'),
    W_names = c('w0', 'w1', 'w2'),
    par_constraint = c(alpha1.1 = 0, alpha2.1 = 0),
    # controls ~ gen_inits
    N_sub = Inf,
    pct_keep = seq(0.40, 1.00, 0.001),
    parallelize = TRUE,
    which_step1 = "best",
    # controls ~ EM
    theta.names = LCTMC::gen_theta_names(K = 2L, type = "2x2", purpose = "em"),
    EM_controls = list(maxit = 50, ELL_tol = 1e-1, LPY_tol = 1e-3, par_tol = 1e-3),
    optim_controls = list(fnscale = -nrow(example_df2x2), maxit = 1e8, factr = 1e-4),
    # misc.
    test_if_global_optim = list(test = FALSE, true_params = NA),
    parallel_optim = list(run = TRUE, cl = parallel::makeCluster(spec = parallel::detectCores()-1)),
    MyModelName = "My 2x2 (K2) Model"
  )

  # - # - # - # - # - # - # - # - #

  ## this is a 3x3 example, with 3 latent classes
  data("example_df3x3", package = "LCTMC")
  model_3x3 = LCTMC::lctmc_3x3(
    # data
    data = example_df3x3,
    # any scaling transformation
    dt_scale = c(dt = 1/365),
    x_scale = c(x0 = 1, x1 = 1/8, x2 = 1),
    w_scale = c(w0 = 1, w1 = 1/8, w2 = 1),
    # general model specification
    K = 2L,
    X_names = c('x0', 'x1', 'x2'),
    W_names = c('w0', 'w1', 'w2'),
    par_constraint = NULL,
    # controls ~ gen_inits
    N_sub = Inf,
    pct_keep = seq(0.40, 1.00, 0.001),
    parallelize = TRUE,
    which_step1 = "best",
    # controls ~ EM
    theta.names = LCTMC::gen_theta_names(K = 2L, type = "3x3", purpose = "em"),
    EM_controls = list(maxit = 50, ELL_tol = 1e-1, LPY_tol = 1e-3, par_tol = 1e-3),
    optim_controls = list(fnscale = -nrow(example_df3x3), maxit = 1e8, factr = 1e-4),
    # misc.
    test_if_global_optim = list(test = FALSE, true_params = NA),
    parallel_optim = list(run = TRUE, cl = parallel::makeCluster(spec = parallel::detectCores()-1)),
    MyModelName = "My 3x3 (K2) Model"
  )
}
