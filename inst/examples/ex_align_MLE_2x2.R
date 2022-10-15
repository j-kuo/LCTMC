# ~ ####################### ~ NOTE ~ ################################ ~ #
# note that the latent classes are re-arranged
# such that class 1 always have the largest alpha0 estimate
# transformation are done following multinomial logistic regression model
# ~ ################################################################# ~ #

## how many latent classes?
my_K = 3

## create names of alpha parameter for the MLE
names.a = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(my_K-1))
names.a = paste(names.a$Var1, ".", names.a$Var2, sep = "")

## create names of beta parameter for the MLE
names.b = expand.grid(c("beta0", "beta1", "beta2"), c("12", "21"), 1:my_K)
names.b = paste(names.b$Var1, ".", names.b$Var2, "_", names.b$Var3, sep = "")

## create a vector of numeric value to serve as the MLE
my_mle = 1:(length(names.a) + length(names.b))
names(my_mle) = c(names.a, names.b)

## use function in the case when true parameters are NOT supplied
## works for any `my_K` >= 2
## warning message is normal
align_MLE_2x2(mle = my_mle, K = my_K)

\dontrun{
  # use the LCTMC.simulate package to generate some true values
  my_true = LCTMC.simulate::gen_true_param(K_class = 3, M_state = 2)

  # use function in the case when true parameters are supplied
  # only works when `my_K` is equal to 3 & requires `my_true` to be defined
  align_MLE_2x2(true = my_true, mle = my_mle, K = my_K)
}

