## 2x2 example
data("example_df2x2", package = "LCTMC")
data("model_2x2", package = "LCTMC")

LCTMC::test_global_optim(
  m = model_2x2,
  data = example_df2x2
)


## 3x3 example
data("example_df3x3", package = "LCTMC")
data("model_3x3", package = "LCTMC")

LCTMC::test_global_optim(
  m = model_3x3,
  data = example_df3x3
)
