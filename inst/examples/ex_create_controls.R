# use only default
create_controls(type = "2x2")


# use only default ~ but use specify input data to generate some arguments
data("example_df2x2", package = "LCTMC")
create_controls(type = "3x3", data = example_df2x2)
