
#ifndef FUNCTIONS_H
#define FUNCTIONS_H


//////////////////////////////
// Header Files
//////////////////////////////


#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include<algorithm>


//////////////////////////////
// Function Prototypes
//////////////////////////////


// Function to populate action space
std::vector<double> populate_action_space_fn(
        int m,
        double A_min,
        double A_step_size
);


// Cartesian product of action space
std::vector<std::vector<double>> action_space_cartesian_product_fn(
        const std::vector<double>& action_space,
        int n,
        int m
);

// Function to determine quantity
std::vector<std::vector<double>> quantity_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<int>& a,
        int n,
        int a_0,
        double mu,
        int S_cardinality
);

// Function to determine total market revenue for each state
std::vector<double> revenue_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& quantity,
        int S_cardinality,
        int n
);

// Function to determine consumer surplus for each state
std::vector<double> cs_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<int>& a,
        double mu,
        int a_0,
        int n,
        int S_cardinality
);

// Function to determine profits
std::vector<std::vector<double>> profits_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& quantity,
        const std::vector<int>& marginal_cost,
        int S_cardinality,
        int n
);

// Function to determine average profits for both firms from setting each possible price
std::vector<std::vector<double>> avg_profits_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& profits,
        const std::vector<double>& action_space,
        int n,
        int m,
        double delta
);

// Function to populate Q-matrix
std::vector<std::vector<double>> initial_Q_matrix_fn(
        const std::vector<std::vector<double>>& avg_profits_for_states_matrix,
        int n,
        int m,
        int S_cardinality
);

// Function to get logit competitive equilibrium outcome
std::vector<std::vector<double>> logit_competitive_outcome_fn(
        const std::vector<int>& a,
        const std::vector<int>& mc,
        double mu,
        int n,
        int a_0,
        int iteration_limit
);

// Function to print main results after all episodes
void print_results(
        const std::vector<double>& avg_profits,
        const std::vector<double>& avg_prices,
        const std::vector<double>& avg_quantities,
        const std::vector<double>& avg_revenues,
        const std::vector<int>& a,
        const std::vector<int>& mc,
        double avg_cs,
        double average_seconds,
        double mu,
        int average_time_steps,
        int n,
        int a_0,
        int iteration_limit
);



#endif //FUNCTIONS_H
