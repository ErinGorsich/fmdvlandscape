setwd("~/My Drive/Warwick/students/Diana/Buffalo behavior_Dylan")
library(sf)

# Population features from Aug 2025 survey
population_data <- read.csv("buffalo_census_2025.csv")
pop_df <- data.frame(
    pop_id = seq(1, length(population_data[,1]) ),
    total = population_data$TOTAL,
    bulls = population_data$Bulls,
    latitude = population_data$Latitude,
    longitude = population_data$Longitude)
pop_df$herdtype <- NA
pop_df$land_type <- NA
for (i in 1:length(pop_df[, 1])){
    ifelse((pop_df$total[[i]] - pop_df$bulls[[i]]) > 0, 
    pop_df$herdtype[[i]] <- "mixed",
    pop_df$herdtype[[i]] <- "bulls")
}
n_groups <- length(pop_df[,1])
n_mixed_groups <- length(pop_df[pop_df$herdtype == "mixed", ])
n_buffalo <- sum(pop_df$total)
par(mfrow = c(1, 2))
hist(pop_df$total[pop_df$herdtype == "mixed"], xlab = "Size, Mixed herd", main = "")
hist(pop_df$total[pop_df$herdtype == "bulls"], xlab = "Size, Bachelor herds", main = "")

# pop_sf <- st_as_sf(pop_df, coords = c("longitude","latitude"))
# ecozone <- 
# pop_df$ecozone <- overlay(r1, r2, fun=function(x,y){return(x/y)})


# Set initial conditions. Put one infected in the second herd... (UPDATE ME LATER)
make_initial_states <- function(pop_df, n_groups) {
    state <- matrix(0L, nrow = n_groups, ncol = 5L,
        dimnames = list(NULL, c("S", "E", "I", "R", "C")) )
    state[, "S"] <- pop_df$total - c(0, 1, rep(0, n_groups - 2))
    state[, "I"] <- c(0, 1, rep(0, n_groups - 2))
    return(state)
}
initial_state <- make_initial_states(pop_df, n_groups)


# Set parameters for stochastic disease dynamics, non-interacting populations
params <- list(beta = 2.8,
    epsilon = 1 / 0.5,     # E -> I rate, per day
    gamma = 1/ 5.7,        # I -> R rate, per days
    p = 0.9,               # fraction of I transitions going to C
    xi = 1/243.0,          # C -> S rate, per day
    tau = 1.0              # tau leap
)


# Run one stochastic ODE for any population
run_one_stochastic_dz = function(state, days, params){
    n_steps <- ceiling(days / params$tau)
    p_EI <- 1 - exp(-params$epsilon * params$tau)
    p_I_exit <- 1 - exp(-params$gamma * params$tau)
    p_CS <- 1 - exp(-params$xi * params$tau)
    
    for (step in seq_len(n_steps)) {
        S <- state[, "S"]
        E <- state[, "E"]
        I <- state[, "I"]
        R <- state[, "R"]
        C <- state[, "C"]
        N <- S + E + I + R + C
        
        force_of_infection <- ifelse(N > 0, params$beta * I / N, 0)
        p_SE <- 1 - exp(-force_of_infection * params$tau)
        p_SE <- pmin(pmax(p_SE, 0), 1)
        
        new_exposed <- rbinom(length(S), size = S, prob = p_SE)
        new_infectious <- rbinom(length(E), size = E, prob = p_EI)
        infectious_exits <- rbinom(length(I), size = I, prob = p_I_exit)
        to_C <- rbinom(length(I), size = infectious_exits, prob = params$p)
        to_R <- infectious_exits - to_C
        C_to_S <- rbinom(length(C), size = C, prob = p_CS)
        
        state[, "S"] <- S - new_exposed + C_to_S
        state[, "E"] <- E + new_exposed - new_infectious
        state[, "I"] <- I + new_infectious - infectious_exits
        state[, "R"] <- R + to_R
        state[, "C"] <- C + to_C - C_to_S
    }
    return(state)
}

# Test run, only one seed infection, one takes off)
state <- run_one_stochastic_dz(state = initial_state, days = 14, params)


# Set up herd parameters
# Randomly divide one population into two populations while preserving all
# compartment totals. A 50:50 expected split is used here.
split_population <- function(population) {
    daughter_1 <- as.integer(mapply(
        FUN = function(total) rbinom(1L, size = total, prob = 0.5),
        total = population
    ))
    daughter_2 <- as.integer(population - daughter_1)
    rbind(daughter_1, daughter_2)
}

test <- split_population(state[2,])

# Apply merge/split/stay events to the current population state.
apply_population_events <- function(state, params) {
    n <- nrow(state)
    
    p_merge <- rbeta(
        1L, params$merge_beta_shape1, params$merge_beta_shape2
    )
    p_split <- rbeta(
        1L, params$split_beta_shape1, params$split_beta_shape2
    )
    
    # Ensure the three event probabilities form a valid categorical model.
    if (p_merge + p_split > 1) {
        scale <- 1 / (p_merge + p_split)
        p_merge <- p_merge * scale
        p_split <- p_split * scale
    }
    
    event <- sample(
        c("merge", "split", "stay"),
        size = n,
        replace = TRUE,
        prob = c(p_merge, p_split, 1 - p_merge - p_split)
    )
    
    merge_candidates <- which(event == "merge")
    split_candidates <- which(event == "split")
    merge_pairs <- split(merge_candidates, ceiling(seq_along(merge_candidates) / 2))
    
    output <- list()
    output_event <- character(0)
    used <- rep(FALSE, n)
    
    # Merge candidates in random pairs. An unpaired candidate stays unchanged.
    if (length(merge_pairs) > 0L) {
        for (pair in merge_pairs) {
            if (length(pair) == 2L) {
                merged <- colSums(state[pair, , drop = FALSE])
                output[[length(output) + 1L]] <- merged
                output_event <- c(output_event, "merge")
                used[pair] <- TRUE
            }
        }
    }
    
    # Unpaired merge candidates and stay candidates remain as single populations.
    unchanged <- which(!used & event != "split")
    for (index in unchanged) {
        output[[length(output) + 1L]] <- state[index, ]
        output_event <- c(output_event, ifelse(event[index] == "merge", "merge_unpaired", "stay"))
        used[index] <- TRUE
    }
    
    # Split selected populations into two daughter populations.
    for (index in split_candidates) {
        daughters <- split_population(state[index, ])
        output[[length(output) + 1L]] <- daughters[1, ]
        output[[length(output) + 1L]] <- daughters[2, ]
        output_event <- c(output_event, "split", "split")
        used[index] <- TRUE
    }
    
    new_state <- do.call(rbind, output)
    colnames(new_state) <- colnames(state)
    rownames(new_state) <- NULL
    
    list(
        state = matrix(as.integer(new_state), ncol = 5L,
                       dimnames = list(NULL, colnames(state))),
        event = output_event,
        p_merge = p_merge,
        p_split = p_split,
        n_before = n,
        n_after = nrow(new_state)
    )
}

run_model <- function(params, seed = 2026) {
    set.seed(seed)
    state <- make_initial_state(params)
    history <- vector("list", params$cycles)
    event_log <- vector("list", params$cycles)
    
    for (cycle in seq_len(params$cycles)) {
        state_before_event <- simulate_block(state, params$block_days, params)
        event_result <- apply_population_events(state_before_event, params)
        state <- event_result$state
        
        history[[cycle]] <- list(
            cycle = cycle,
            end_day = cycle * params$block_days,
            state_before_event = state_before_event,
            state_after_event = state,
            p_merge = event_result$p_merge,
            p_split = event_result$p_split
        )
        
        event_log[[cycle]] <- data.frame(
            cycle = cycle,
            end_day = cycle * params$block_days,
            p_merge = event_result$p_merge,
            p_split = event_result$p_split,
            populations_before = event_result$n_before,
            populations_after = event_result$n_after,
            merge_events = sum(event_result$event == "merge"),
            split_events = sum(event_result$event == "split"),
            stay_events = sum(event_result$event == "stay"),
            unpaired_merge_candidates = sum(event_result$event == "merge_unpaired"),
            row.names = NULL
        )
    }
    
    list(
        history = history,
        event_log = do.call(rbind, event_log),
        final_state = state,
        parameters = params
    )
}

# Run 26 cycles of 14 days: total simulated time = 364 days.
result <- run_model(parameters, seed = 2026)

# Event summary.
print(result$event_log)

# Plot total compartment counts across all currently existing populations.
times <- vapply(result$history, function(x) x$end_day, numeric(1))
total_state <- do.call(rbind, lapply(
    result$history,
    function(x) colSums(x$state_before_event)
))
matplot(
    times, total_state,
    type = "l", lty = 1, lwd = 2,
    col = c("steelblue", "darkorange", "firebrick", "darkgreen", "purple"),
    xlab = "Time (days)", ylab = "Total number of animals",
    main = "Blockwise stochastic disease dynamics"
)
legend(
    "topright", legend = colnames(total_state),
    col = c("steelblue", "darkorange", "firebrick", "darkgreen", "purple"),
    lty = 1, lwd = 2, bty = "n"
)

# Plot the number of populations after each event cycle.
plot(
    result$event_log$end_day,
    result$event_log$populations_after,
    type = "o", pch = 16,
    xlab = "Time (days)", ylab = "Number of populations",
    main = "Population count after merge/split events"
)



