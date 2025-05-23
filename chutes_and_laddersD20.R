pacman::p_load(wrappedtools, tidyverse, ggplot2)

# board parameters
n_squares <- 49
ladder_positions <- c(3, 7, 23, 27, 33)
ladder_destinations <- c(18, 21, 34, 42, 46)
chute_positions <- c(10, 16, 32, 45)
chute_destinations <- c(6, 14, 20, 41)

# helper function to get the destination of a ladder or chute
get_destination <- function(square) {
  #check for ladders or chutes. If the next square us a ladder or chute, we take the index of the ladder or chute position and reassign next square the the corresponding index in the destination square using which()
  if (square %in% ladder_positions) {
    return(ladder_destinations[which(ladder_positions == square)])
  } else if (square %in% chute_positions) {
    return(chute_destinations[which(chute_positions == square)])
  }
  return(square)
}

# define a function to simulate a single game
chutes_and_ladders <- function() {
  visit_counter <- rep(0, n_squares) #The rep() function creates an empty vector which can be filled with the number of visits
  names(visit_counter) <- 1:n_squares
  current_square <- 1  
  die_roll_count <- 0
  while (current_square < n_squares) {
    roll <- sample(1:20, 1)  # die roll
    next_square <- current_square + roll
    
    # ensure the player doesn't go past the last square
    if (next_square > n_squares) {
      next_square <- n_squares
    }
    
    # count the visit to the destination square (accounting for ladders and chutes)
    destination_square <- get_destination(next_square)
    visit_counter[as.character(destination_square)] <- visit_counter[as.character(destination_square)] + 1 #you need to use as.character because the vector names are strings. Otherwise it doesn't work.
    
    # update the current square
    current_square <- destination_square
    die_roll_count <- die_roll_count+1
    # check if the game has ended (reached or passed square 49)
    if (current_square == n_squares) {
      break  # End the game
    }
  }
  
  return(list(visits = visit_counter, rolls = die_roll_count)) #R doesn't allow multiple return statements
}

# Run the simulation multiple times and accumulate results. Counting vector to keep track of visits to squares (excluding ladder/chute squares). 

runs <- 10000
total_visits <- rep(0, n_squares)
names(total_visits) <- 1:n_squares
die_rolls_per_game <- numeric(runs) #create numeric vector to store die roll counts

for (i in 1:runs) {
  result <- chutes_and_ladders()
  total_visits <- total_visits + result$visits
  die_rolls_per_game[i] <- result$rolls
}

print(total_visits)
print(summary(die_rolls_per_game)) #tried to print the whole thing, don't do it bad idea

#ggplot dataframe preparation

# create the tibble for plotting
results_tibble <- tibble(
  square = names(total_visits[-49]),
  visits = as.vector(total_visits[-49])
)


ggplot(results_tibble, aes(x = square, y = visits)) +
  geom_col(fill = "salmon") +
  theme_minimal() +
  labs(
    title = "Number of Visits per Square in Chutes and Ladders After 10000 Runs",
    x = "Square Number",
    y = "Number of Visits"
  ) +
  scale_x_discrete(breaks = seq(1, 48)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

