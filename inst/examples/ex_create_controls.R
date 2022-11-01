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
