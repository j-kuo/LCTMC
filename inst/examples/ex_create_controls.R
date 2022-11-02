create_controls(
  # type & aspect
  type = '2x2',
  # ...
  scaling = c(x1 = 1),
  trace = FALSE
)


create_controls(
  # type & aspect
  type = '3x3',
  # ...
  scaling = c(x1 = 1),
  trace = FALSE,
  pct_keep = seq(0.4, 1, 0.1),
  parallelize.init01 = FALSE
)


create_controls(
  # type & aspect
  type = '2x2',
  # ...
  scaling = c(x1 = 1),
  trace = FALSE,
  pct_keep = seq(0.4, 1, 0.1),
  parallelize.init01 = TRUE,
  which_step1 = "best"
)


create_controls(
  # type & aspect
  type = '3x3',
  # ...
  scaling = c(x1 = 1),
  trace = FALSE,
  pct_keep = seq(0.4, 1, 0.1),
  parallelize.init01 = TRUE,
  which_step1 = "best",
  EM.maxit = 10,
  EM.ELL_tol = 0.1,
  EM.LPY_tol = 0.001,
  EM.par_tol = 0.001,
  LBFGSB.fnscale = -1,
  LBFGSB.maxit = 1e5,
  LBFGSB.factr = 0.003
)


create_controls(
  # type & aspect
  type = '2x2',
  # ...
  scaling = c(x1 = 1),
  trace = FALSE,
  pct_keep = seq(0.4, 1, 0.1),
  parallelize.init01 = TRUE,
  which_step1 = "best",
  EM.maxit = 10,
  EM.ELL_tol = 0.1,
  EM.LPY_tol = 0.001,
  EM.par_tol = 0.001,
  LBFGSB.fnscale = -1,
  LBFGSB.maxit = 1e5,
  LBFGSB.factr = 0.003,
  solve_tol = 1e-6,
  # symmetric_tol = 1e-10,
  eigen0_tol = 1e-10
)
