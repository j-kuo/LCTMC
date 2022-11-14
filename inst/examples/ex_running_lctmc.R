# The example below demonstrates how to specify the 'lctmc' functions correctly
#   to perform the latent class CTMC model.
# -------------------------------------------------------------------------------------
# Both the 2x2 and 3x3 examples are listed below, however it should be noted
#   that the models may take a while to run (approx. 10-30min on a high spec laptop)
# -------------------------------------------------------------------------------------
# In the 2x2 model, we constrained the model so that the covariate effects on
#   the latent class probability are equal to 0.
# -------------------------------------------------------------------------------------
# In the 3x3 model, we left the model un-constrained.
# -------------------------------------------------------------------------------------
# In this code, we also use a for-loop to fit models with 2, 3, and 4 latent classes
#   the actual number of latent classes within these data sets is 2.
# -------------------------------------------------------------------------------------
# The model with the smallest BIC is the "best" model.
# -------------------------------------------------------------------------------------

\dontrun{
  ## this is a 2x2 example, with 3 latent classes
  data("example_df2x2", package = "LCTMC")
  ctrl_2x2 = LCTMC::create_controls(type = "2x2", data = example_df2x2)

  m2x2_list = list()
  for (k in 2:4) {
    m2x2 = LCTMC::lctmc_2x2(
      # data
      data = example_df2x2,
      # general model specification
      K = k,
      X_names = c('x0', 'x1', 'x2'),
      W_names = c('w0', 'w1', 'w2'),
      par_constraint = c(alpha1.1 = 0, alpha2.1 = 0),
      # misc.
      controls = ctrl_2x2,
      parallel_optim = list(
        run = TRUE, cl = parallel::makeCluster(spec = parallel::detectCores()-1)
      ),
      MyModelName = paste("My 2x2 (K=", k, ") model", sep = "")
    )

    m2x2_list[[paste("k", k, sep = "")]] = m2x2
  }

  BIC = sapply(
    X = m2x2_list[-1],
    FUN = function(x) {
      log_like = LCTMC::test_global_optim(m = x, data = example_df2x2)$L_mle
      k = x$n_pars
      n = nrow(example_df2x2) - length(unique(example_df2x2$id))
      -2*log_like + k*log(n)
    }
  )
  BIC
  which.min(BIC)


  # - # - # - # - # - # - # - # - #


  ## this is a 3x3 example, with 3 latent classes
  data("example_df3x3", package = "LCTMC")
  ctrl_3x3 = LCTMC::create_controls(type = "3x3", data = example_df3x3)

  m3x3_list = list()
  for (k in 2:4) {
    m3x3 = LCTMC::lctmc_3x3(
      # data
      data = example_df3x3,
      # general model specification
      K = k,
      X_names = c('x0', 'x1', 'x2'),
      W_names = c('w0', 'w1', 'w2'),
      par_constraint = NULL,
      # misc.
      controls = ctrl_3x3,
      parallel_optim = list(
        run = T,
        cl = parallel::makeCluster(spec = parallel::detectCores()-1)
      ),
      MyModelName = paste("My 3x3 (K=", k, ") model", sep = "")
    )

    m3x3_list[[paste("k", k, sep = "")]] = m3x3
  }

  BIC = sapply(
    X = m3x3_list[-1],
    FUN = function(x) {
      log_like = LCTMC::test_global_optim(m = x, data = example_df3x3)$L_mle
      k = x$n_pars
      n = nrow(example_df3x3) - length(unique(example_df3x3$id))
      -2*log_like + k*log(n)
    }
  )
  BIC
  which.min(BIC)
}
