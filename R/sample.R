library(tidyverse)

iris |> 
	dplyr::select(Sepal.Length) |>
  head()