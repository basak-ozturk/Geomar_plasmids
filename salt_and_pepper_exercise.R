
library(tibble)
library(dplyr)

saltandpepper <- tibble(
  salt = paste("salt", 1:100),
  pepper = paste("pepper", 1:100)
)

print(saltandpepper)
saltandpepper[,1]
saltandpepper |> select(salt)
saltandpepper |> pull(salt)
saltandpepper |> select(salt) |> slice(1:50)
saltandpepper |> slice(2:25) |> pull(salt)
saltandpepper |> slice_sample(n=20) |> pull(salt)
saltandpepper |> slice_max(order_by = salt, n=20) |> pull(salt)
