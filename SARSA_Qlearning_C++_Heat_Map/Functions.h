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


//////////////////////////////
// Function Prototypes
//////////////////////////////


// Function to populate action space
std::vector<double> populate_A(
    int m,
    double A_min,
    double step_size
);


// Cartesian product of action space
std::vector<std::vector<double>> generate_combinations(
    const std::vector<double>& A,
    int n,
    int m
);

// Function to determine quantity for each state for each agent
std::vector<std::vector<double>> q_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<int>& a,
    int n,
    int a_0,
    double mu,
    int S_cardinality
);

// Function to determine total market revenue for each state for each agent
std::vector<double> rvn_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<std::vector<double>>& quantity,
    int S_cardinality,
    int n
);

// Function to determine profit for each state for each agent
std::vector<std::vector<double>> r_fn(
        const std::vector<std::vector<double>>& price,
        const std::vector<std::vector<double>>& quantity,
        const std::vector<int>& marginal_cost,
        int S_cardinality,
        int n
);

// Function to determine consumer surplus for each state
std::vector<double> cs_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<int>& a,
    double mu,
    int a_0,
    int n,
    int S_cardinality
);

// Function to determine average profits for each agent when setting each of the m possible prices
std::vector<std::vector<double>> Qi0_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<std::vector<double>>& profit,
    const std::vector<double>& A,
    int n,
    int m,
    double delta
);

// Function to populate Q-matrix for each agent
std::vector<std::vector<std::vector<double>>> Q_fn(
        const std::vector<std::vector<double>>& Qi0,
        int n,
        int m,
        int S_cardinality
);

// Function to generate all possible values of alpha
std::vector<double> alpha_fn(
        int alpha_values,
        double alpha_min,
        double alpha_step_size
);

// Function to generate all possible values of beta
std::vector<double> beta_fn(
        int beta_values,
        double beta_min,
        double beta_step_size
);


#endif //FUNCTIONS_H