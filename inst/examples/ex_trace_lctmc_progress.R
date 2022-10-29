# width = 120
{
  LCTMC:::trace_lctmc_progress(width = 120L,
                               section = 'header1', type = "format",
                               MyModelName = "ABCDEF")
  LCTMC:::trace_lctmc_progress(width = 120L,
                               section = 'header2', type = "format",
                               MyModelName = "ABCDEF")
}

# width = 100
{
  LCTMC:::trace_lctmc_progress(width = 100L,
                               section = 'header1', type = "format",
                               MyModelName = "Hello World !!")
  LCTMC:::trace_lctmc_progress(width = 100L,
                               section = 'header2', type = "init1",
                               MyModelName = "Hello World !!")
}

# an example with closing messages too
{
  xx = "Hellow World !!"
  tt = as.POSIXct(1, origin = "2022-10-15", tz = "UTC")

  LCTMC:::trace_lctmc_progress(width = 101L,
                               section = "header1", type = "format",
                               ref_t = tt,
                               MyModelName = xx)
  LCTMC:::trace_lctmc_progress(width = 101L,
                               section = "header2", type = "init2",
                               ref_t = tt,
                               MyModelName = xx)

  cat(" * some notes here \n * xyz \n * abc \n * wow~ \n * more notes \n")

  LCTMC:::trace_lctmc_progress(width = 101L,
                               section = "tail1", type = "rescale",
                               ref_t = tt,
                               MyModelName = xx)
  LCTMC:::trace_lctmc_progress(width = 101L,
                               section = "tail2", type = "rescale",
                               ref_t = tt,
                               MyModelName = xx)
}
