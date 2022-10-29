#' A simulated data set for binary outcome.
#'
#' The data set was simulated using the 'LCTMC.simulate' package. This data set
#' is meant to demonstrate how the input data should be formatted for a **2x2** model,
#' and it can be used in conjunction with the example code under [lctmc_2x2()].
#'
#' @usage data(example_df2x2)
#'
#' @format A data frame with dimension **60895 x 7**, columns have the following interpretation:
#' \describe{
#'   \item{id}{ID number as the individual/person identifier (character)}
#'   \item{obsTime}{time of observation/data collection (numeric)}
#'   \item{state_at_obsTime}{the observed outcome corresponding to the time of observation/data collection (numeric; 1 or 2)}
#'   \item{x1}{variable #1 that can affect the CTMC process (numeric)}
#'   \item{x2}{variable #2 that can affect the CTMC process (numeric)}
#'   \item{w1}{variable #1 that can affect the latent class probabilities (numeric)}
#'   \item{w2}{variable #2 that can affect the latent class probabilities (numeric)}
#' }
"example_df2x2"

#' A simulated data set for 3-category outcome with one absorbing state.
#'
#' The data set was simulated using the 'LCTMC.simulate' package. This data set
#' is meant to demonstrate how the input data should be formatted for a **3x3** model,
#' and it can be used in conjunction with the example code under [lctmc_3x3()].
#'
#' @usage data(example_df3x3)
#'
#' @format A data frame with dimension **52243 x 7**, columns have the following interpretation:
#' \describe{
#'   \item{id}{ID number as the individual/person identifier (character)}
#'   \item{obsTime}{time of observation/data collection (numeric)}
#'   \item{state_at_obsTime}{the observed outcome corresponding to the time of observation/data collection (numeric; 1, 2, or 3). \cr
#'   in this case, stage 3 is set to be an absorbing state (e.g., death)}
#'   \item{x1}{variable #1 that can affect the CTMC process (numeric)}
#'   \item{x2}{variable #2 that can affect the CTMC process (numeric)}
#'   \item{w1}{variable #1 that can affect the latent class probabilities (numeric)}
#'   \item{w2}{variable #2 that can affect the latent class probabilities (numeric)}
#' }
"example_df3x3"
