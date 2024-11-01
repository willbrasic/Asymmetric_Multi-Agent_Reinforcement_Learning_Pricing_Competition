
#include "Global.h"


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Parameters demand model
    const int n{ 2 };                              // Number of firms
    const int a_0{ 0 };                            // Inverse index of aggregate demand
    const double mu{ 1.0 / 4.0 };                  // Horizontal differentiation index
    const std::vector<int> mc(n, 1);         // Marginal cost
    const std::vector<int> a(n, 2);          // Product quality index

    // Learning parameters
    const int k{ 1 };                              // Memory (number of periods)
    const int m{ 15 };                             // Number of equally spaced price points
    const int E{ 100 };                              // Number of episodes
    const double delta{ 0.95 };                    // Discount factor
    const double alpha{ 0.15 };                    // Learning rate
    const double beta{ 1e-5 };                     // epsilon-greedy experimentation parameter

    // Minimum and maximum values in action space
    const double A_min{ 1.0 };
    const double A_max{ 2.1 };

    // Step size for action space
    const double A_step_size = (A_max - A_min) / (m - 1);

    // Size of state space
    const int S_cardinality = static_cast<int>(pow(m, n * k));

    // Calculate argmax(Q) once every *convergence_check* time steps to check for convergence
    const int convergence_check{ 100 };

    // Number of time steps needed for argmax(Q) to be constant for convergence
    const int norm{ 100000 / convergence_check };

    // Number of time steps allowed per episode
    const int max_t{ 10000000 };

    // Display results for episode e every results_check time steps
    const int results_check{ 100000 };

    // How many points for learning curve plot
    const int learning_curve_points{ 70 };

    // Store results for learning curve every learning_curve_check time steps
    const int learning_curve_check{ 10000 };

    // Vector to store average profits for learning curve
    std::vector<std::vector<double>> learning_curve_profits(learning_curve_points, std::vector<double>(E, 0.0));

    // Vector store how much time it took episode loop to run
    std::vector<double> timings(E, 0.0);

    // Vector to store how many time steps needed for convergence for each episode
    std::vector<int> time_steps(E, 0);

    // Vector to store prices in episode e averaged over the last results_check time steps prior to convergence
    std::vector<std::vector<double>> prices_e(n, std::vector<double>(E, 0.0));

    // Vector to store quantities in episode e  averaged over the last results_check time steps prior to convergence
    std::vector<std::vector<double>> quantities_e(n, std::vector<double>(E, 0.0));

    // Vector to store profits in episode e averaged over the last results_check time steps prior to convergence
    std::vector<std::vector<double>> profits_e(n, std::vector<double>(E, 0.0));

    // Vector to store revenues in episode e averaged over the last results_check time steps prior to convergence
    std::vector<std::vector<double>> revenues_e(n, std::vector<double>(E, 0.0));

    // Vector to store consumer surplus in episode e
    std::vector<double> cs_e;

    // Vector to store consumer surplus averaged over the last results_check time steps prior to convergence in episode e
    std::vector<double> avg_cs_e(E, 0.0);

} //Global
