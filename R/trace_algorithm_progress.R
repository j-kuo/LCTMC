#' @title Trace algorithm progress
#'
#' @description Trace algorithm progress by printing section titles
#'
#' @param width an integer scalar. Use this to indicate output message width. Minimum width is 100. \cr
#' If any thing less than 100 is specified then function will automatically default it to `100L`.
#' @param section a character scalar. If 'header1' then `MyModelName` is printed. If equals to 'header2' then the output message depends on `type`. \cr
#' If 'tail1' then prints the closing message depending on `type`. If 'tail2' then prints the time stamp comparing to the argument `ref_t`
#' @param type a character scalar. when `section` is equal to 'header2' this argument is used to decide the type of message that is printed. \cr
#' Ignored if `section` is equal to 'header1'
#' @param ref_t a scalar time object obtained from `Sys.time()`.
#' @param MyModelName a character scalar. This will be used when `section` is equal to 'header1'
#'
#' @return NULL. \cr
#' This function prints messages and outputs NULL
#'
#' @note This function is usually called within the [lctmc_2x2()] or [lctmc_3x3()] function
#'
#' @seealso [lctmc_2x2()], [lctmc_3x3()]
#'
#' @example inst/examples/ex_trace_lctmc_progress.R

trace_lctmc_progress = function(width = 100L,
                                section = c("header1", "header2", "tail1", "tail2"),
                                type = c("format", "init1", "init2", "em", "se", 'rescale'),
                                ref_t = Sys.time(),
                                MyModelName = "") {
  ### check arguments
  if (!is.integer(width)) {
    stop("`width` should be an integer value. Specified it by using 'L'. For example, `100L` or `150L`")
  }
  if (!section %in% c("header1", "header2", "tail1", "tail2")) {
    stop("`section` must be either 'header1' or 'header2'")
  }
  if (!type %in% c("format", "init1", "init2", "em", "se", 'rescale')) {
    stop("`type` must be either 'format', 'init1', 'init2', 'em', 'se', or 'rescale'")
  }
  ### width is minimum 100 characters
  width = max(width, 100)

  ### create messages ~ headers
  if (section == "header1") {
    ## header line 1
    h1_lead = c("#####~{", MyModelName, "}~")
    h1_tail = rep("#", times = width - 9 - nchar(MyModelName))

    ## out msg
    cat(h1_lead, h1_tail,"\n", sep = "")
  } else if (section == "header2") {
    ## header line 2 ~ depends on `type` (i.e., middle section)
    if (type == "format") {
      h2_mid = " Format Data (my_df , X & W matrix , dt) "
    } else if (type == "init1") {
      h2_mid = " Generating Initial Values -- Step 1 "
    } else if (type == "init2") {
      h2_mid = " Generating Initial Values -- Step 2 "
    } else if (type == "em") {
      h2_mid = " Initiating Expectation-Maximization Algorithm "
    } else if (type == "se") {
      h2_mid = " Computing Hessian Matrix for SE approx. "
    } else if (type == "rescale") {
      h2_mid = " Re-scaling Parameters "
    }
    h2_lead = rep("#", times = 5)
    h2_tail = rep("#", times = width - length(h2_lead) - nchar(h2_mid))

    ## out msg
    cat(h2_lead, h2_mid, h2_tail, "\n", sep = "")
  }

  ### create messages ~ closing
  if (section == "tail1") {
    ## closing line 1 ~ depending on value of `type`
    if (type == "format") {
      t1_mid = " finished formatting "
    } else if (type == "init1") {
      t1_mid = " initial values 1/2 "
    } else if (type == "init2") {
      t1_mid = " initial values 2/2 "
    } else if (type == "em") {
      t1_mid = " EM converged/maxit reached "
    } else if (type == "se") {
      t1_mid = " obtained hessian "
    } else if (type == "rescale") {
      t1_mid = " rescaled parameters to original data's units "
    }
    t1_tail = "#####"
    t1_head = rep(".", times = width - nchar(t1_mid) - nchar(t1_tail))

    ## out msg
    cat(t1_head, t1_mid, t1_tail, "\n", sep = "")
  } else if (section == "tail2") {
    ## tail 2nd line ~ time stamp
    time_diff = difftime(time1 = Sys.time(), time2 = ref_t, units = 'hours')
    t2_tail = paste(" Timestamp: [", sprintf("%.2f", time_diff), "] hr #####", sep = "")
    t2_lead = rep(".", times = width - nchar(t2_tail))

    ## out msg
    cat(t2_lead, t2_tail, "\n\n", sep = "")
  }
}
