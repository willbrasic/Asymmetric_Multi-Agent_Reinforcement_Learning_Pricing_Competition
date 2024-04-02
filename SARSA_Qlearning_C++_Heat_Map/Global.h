
#include <vector>
#include <cmath>

#ifndef GLOBAL_H
#define GLOBAL_H


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Parameters demand model
    extern const int n;                     // Number of agents                          
    extern const int a_0;                   // Inverse index of aggregate demand
    extern const double mu;                 // Horizontal differentiation index
    extern const std::vector<int> mc;       // Marginal cost
    extern const std::vector<int> a;        // Product quality index

    // Learning parameters
    extern const int k;                     // Memory (number of periods)
    extern const int m;                     // Number of equally spaced price points
    extern const int E;                     // Number of episodes
    extern const double delta;              // Discount factor

    // Logit equilibrium competitive and collusive profits for each firm
    extern const std::vector<double> competitive_profits;
    extern const std::vector<double> collusive_profits;

    // Minimum and maximum values in action space
    extern const double A_min;
    extern const double A_max;

    // Step size for action space
    extern const double step_size;

    // How many values of alpha and beta for the heat map
    extern const int alpha_values;
    extern const int beta_values;

    // Minimum and maximum values for alpha and beta space
    extern const double alpha_min;
    extern const double alpha_max;
    extern const double beta_min;
    extern const double beta_max;

    // Step size for alpha and beta space
    extern const double alpha_step_size;
    extern const double beta_step_size;

    // Size of state space
    extern const int S_cardinality;

    // Calculate argmax_a(Q(a,s)) for each s once every convergence_check time steps to check for convergence
    extern const int convergence_check;

    // Calculate results over last results_check time steps prior to convergence
    extern const int results_check;

    // Number of time steps needed for argmax_a(Q(a,s)) to be constant for each s to achieve convergence
    extern const int norm;

    // Number of time steps allowed per episode
    extern const int maxt;

    // Each agent's action at time step t
    extern std::vector<int> actions_t;

    // Each agent's price at time step t
    extern std::vector<double> prices_t;

    // Each agent's quantity at time step t
    extern std::vector<double> quantities_t;

    // Each agent's price at time step t
    extern std::vector<double> profits_t;

    // Each agent's price at time step t
    extern std::vector<double> profits_t;

    // Each agent's Delta at time step t
    extern std::vector<double> Delta_t;

    // SARSA agent's Delta averaged over last results_check time steps prior to convergence for all episodes
    extern std::vector<double> avg_Delta_0;

    // Q-learning agent's Delta averaged over last results_check time steps prior to convergence for all episodes
    extern std::vector<double> avg_Delta_1;

    // Each agent's Delta for each (alpha, beta) pair averaged across all episodes
    extern std::vector<std::vector<std::vector<double>>> Delta;

    // Average convergence time for each (alpha, beta) combination
    extern std::vector<std::vector<std::vector<int>>> time_steps_avg;


} // namespace Global

#endif //GLOBAL_H
