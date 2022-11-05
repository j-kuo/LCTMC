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
  ctrl_2x2 = LCTMC::create_controls(type = "2x2", data = example_df2x2)

  model_2x2 = LCTMC::lctmc_2x2(
    # data
    data = example_df2x2,
    # general model specification
    K = 2L,
    X_names = c('x0', 'x1', 'x2'),
    W_names = c('w0', 'w1', 'w2'),
    par_constraint = c(alpha1.1 = 0, alpha2.1 = 0),
    # misc.
    controls = ctrl_2x2,
    parallel_optim = list(run = TRUE, cl = parallel::makeCluster(spec = parallel::detectCores()-1)),
    MyModelName = "My 2x2 (K2) model"
  )

  # - # - # - # - # - # - # - # - #

  ## this is a 3x3 example, with 3 latent classes
  data("example_df3x3", package = "LCTMC")
  ctrl_3x3 = LCTMC::create_controls(type = "3x3", data = example_df3x3)

  model_3x3 = LCTMC::lctmc_3x3(
    # data
    data = example_df3x3,
    # general model specification
    K = 2L,
    X_names = c('x0', 'x1', 'x2'),
    W_names = c('w0', 'w1', 'w2'),
    par_constraint = NULL,
    # misc.
    controls = ctrl_3x3,
    parallel_optim = list(run = T, cl = parallel::makeCluster(spec = parallel::detectCores()-1)),
    MyModelName = "My 3x3 (K2) Model"
  )
}
