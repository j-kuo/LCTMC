% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_controls.R
\name{create_controls}
\alias{create_controls}
\title{Creates the control object for fitting LCTMC}
\usage{
create_controls(type, ...)
}
\arguments{
\item{type}{a character scalar. Equals either "2x2" or "3x3" to indicate the number of category in the outcome variable.}

\item{...}{can be the following values:
\describe{
\item{data}{certain parameters will require an input data frame to compute the default values. \cr
Use \code{data} to specify this data frame.}
\item{scaling}{controls what factor covariate should be scaled to ensure algorithm stability. \cr
See \code{\link[=fmt_rowwise_trans]{fmt_rowwise_trans()}} for more info. Default value depends on \code{data}}
\item{trace}{controls whether data formatting function should print function details. \cr
See \code{\link[=fmt_rowwise_trans]{fmt_rowwise_trans()}} for more info}
\item{N_sub}{controls how many subjects to use for step 1 of initial value generation. \cr
See \code{\link[=gen_inits01_lctmc_2x2]{gen_inits01_lctmc_2x2()}} for more info}
\item{pct_keep}{controls how much of left/right tail should be truncated for the k-means algorithm. \cr
See \code{\link[=gen_inits01_lctmc_2x2]{gen_inits01_lctmc_2x2()}} for more info}
\item{init01.parallelize}{controls whether step 1 of initial value generation should be parallelized. \cr
See \code{\link[=gen_inits01_lctmc_2x2]{gen_inits01_lctmc_2x2()}} for more info}
\item{which_step1}{a character scalar. Controls which kmeans result to use for step 2 of initial value generation. \cr
Default is "best" which uses the parameter vector that yields the highest \eqn{\log(P(Y))}.
The alternative is "all" which uses all of the kmeans estimate without truncating either tail.}
\item{init02.maxit}{a numeric scalar. Controls the maximum number of iteration for step 2 of initial value generation. \cr
See \code{\link[=gen_inits02_lctmc_2x2]{gen_inits02_lctmc_2x2()}} for more info.}
\item{init02.factr}{a numeric scalar. Controls the convergence tolerance number for step 2 of initial value generation. \cr
See \code{\link[=gen_inits02_lctmc_2x2]{gen_inits02_lctmc_2x2()}} for more info.}
\item{EM.maxit}{controls how many EM iterations the algorithm will perform. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} for more info}
\item{EM.ELL_tol}{controls the convergence tolerance on the expected conditional log-likelihood value. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} for more info}
\item{EM.LPY_tol}{controls the convergence tolerance on the observed log-likelihood value. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} for more info}
\item{EM.par_tol}{controls the convergence tolerance on the magnitude change in parameter values. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} for more info}
\item{LBFGSB.fnscale}{controls the \code{fnscale} argument for \code{optim()}. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} & \code{\link[=optim]{optim()}} for more info. Default value depends on \code{data}.}
\item{LBFGSB.maxit}{controls the \code{maxit} argument for \code{optim()}. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} & \code{\link[=optim]{optim()}} for more info}
\item{LBFGSB.factr}{controls the \code{factr} argument for \code{optim()}. \cr
See \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}} & \code{\link[=optim]{optim()}} for more info}
\item{solve_tol}{controls the tolerance value for detecting linear dependency when \code{solve()} is called. \cr
See \code{\link[=get_SE_lctmc_2x2]{get_SE_lctmc_2x2()}} & \code{\link[=solve]{solve()}} for more info}
\item{symmetric_tol}{controls the tolerance value when checking if the hessian matrix is symmetric. \cr
See \code{\link[=get_SE_lctmc_2x2]{get_SE_lctmc_2x2()}} for more info}
\item{eigen0_tol}{controls the tolerance value for treating eigen values as essentially 0. \cr
See \code{\link[=get_SE_lctmc_2x2]{get_SE_lctmc_2x2()}} & \code{\link[=eigen]{eigen()}} for more info}
}}
}
\value{
a custom class object that acts like a list. The output contains the following:
\describe{
\item{fmt_data}{contains elements: \code{scaling} and \code{trace}.}
\item{init01}{contains elements: \code{N_sub}, \code{pct_keep}, and \code{init01.parallelize}}
\item{init02}{contains elements: \code{which_step1}, \code{maxit}, and \code{factr}}
\item{EM}{contains elements: \code{EM.maxit}, \code{EM.ELL_tol}, \code{EM.LPY_tol}, \code{EM.par_tol},
\code{LBFGS.fnscale}, \code{LBFGS.maxit}, and \code{LBFGSB.factr}.}
\item{SE}{contains elements: \code{solve_tol}, \code{symmetric_tol}, and \code{eigen0_tol}.}
\item{rescale}{currently there are no control options for rescaling parameters.}
\item{type}{controls whether contorls is for a '2x2' model or a '3x3' model.}
}
}
\description{
Create a custom object that acts like a list which serve as object that holds all control
parameters for fitting the LCTMC model.
}
\examples{
# use only default
create_controls(type = "2x2")


# use only default ~ but use specify input data to generate some arguments
data("example_df2x2", package = "LCTMC")
create_controls(type = "3x3", data = example_df2x2)
}
\seealso{
\code{\link[=lctmc_2x2]{lctmc_2x2()}}; \code{\link[=lctmc_3x3]{lctmc_3x3()}}
}
