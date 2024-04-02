
#include <vector>
#include <cmath>

#ifndef GLOBAL_H
#define GLOBAL_H


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Parameters demand model
    extern const int n;                     // Number of firms
    extern const int a_0;                   // Inverse index of aggregate demand
    extern const double mu;                 // Horizontal differentiation index
    extern const std::vector<int> mc;       // Marginal cost
    extern const std::vector<int> a;        // Product quality index

    // Learning parameters
    extern const int k;                     // Memory (number of periods)
    extern const int m;                     // Number of equally spaced price points
    extern const int E;                     // Number of episodes
    extern const double delta;              // Discount factor
    extern const double alpha;              // Learning rate
    extern const double beta;               // epsilon-greedy experimentation parameter

    // Minimum and maximum values in action space
    extern const double A_min;
    extern const double A_max;

    // Step size for action space
    extern const double A_step_size;

    // Size of state space
    extern const int S_cardinality;

    // Calculate argmax(Q) once every *convergence_check* time steps to check for convergence
    extern const int convergence_check;

    // Number of time steps needed for argmax(Q) to be constant for convergence
    extern const int norm;

    // Number of time steps allowed per episode
    extern const int max_t;

    // Display results for episode e every *results_check* time steps
    extern const int results_check;

    // How many points for learning curve plot
    extern const int learning_curve_points;

    // Store results for learning curve every learning_curve_check time steps
    extern const int learning_curve_check;

    // Vector to store average profits for learning curve
    extern std::vector<std::vector<double>> learning_curve_profits;

    // Vector store how much time it took episode loop to run
    extern std::vector<double> timings;

    // Vector to store how many time steps needed for convergence for each episode
    extern std::vector<int> time_steps;

    // Vector to store prices in episode e averaged over the last results_check time steps prior to convergence
    extern std::vector<std::vector<double>> prices_e;

    // Vector to store quantities in episode e averaged over the last results_check time steps prior to convergence
    extern std::vector<std::vector<double>> quantities_e;

    // Vector to store profits in episode e averaged over the last results_check time steps prior to convergence
    extern std::vector<std::vector<double>> profits_e;

    // Vector to store revenues in episode e averaged over the last results_check time steps prior to convergence
    extern std::vector<std::vector<double>> revenues_e;

    // Vector to store consumer surplus in episode e averaged over the last results_check time steps prior to convergence
    extern std::vector<double> cs_e;

    // Vector to store consumer surplus averaged over the last results_check time steps prior to convergence in episode e
    extern std::vector<double> avg_cs_e;


} // namespace Global


#endif //GLOBAL_H
