#include "Global.h"


//////////////////////////////
// Global Variables Namespace
//////////////////////////////


namespace Global {

    // Parameters demand model
    const int n{ 2 };                             // Number of agents
    const int a_0{ 0 };                           // Inverse index of aggregate demand
    const double mu{ 1.0 / 4.0 };                 // Horizontal differentiation index
    const std::vector<int> mc(n, 1);        // Marginal cost
    const std::vector<int> a(n, 2);         // Product quality index

    // Learning parameters
    const int k{ 1 };                             // Memory (number of periods)
    const int m{ 15 };                            // Number of equally spaced price points
    const int E{ 100 };                           // Number of episodes
    const double delta{ 0.95 };                   // Discount factor

    // Logit equilibrium competitive and collusive profits for each firm (use MATLAB or C++ script to get these values; make sure they align with n)
    const std::vector<double> competitive_profits(n, 0.2229);
    const std::vector<double> collusive_profits(n, 0.3375);

    // Minimum and maximum values in action space
    const double A_min{ 1.0 };
    const double A_max{ 2.1 };

    // Step size for action space
    const double step_size = (A_max - A_min) / (m - 1);

    // How many values of alpha and beta for the heat map
    const int alpha_values{ 100 };
    const int beta_values{ 100 };

    // Minimum and maximum values for alpha and beta space
    const double alpha_min{ 0.01 };
    const double alpha_max{ 0.3 };
    const double beta_min{ 1e-6 };
    const double beta_max{ 2 * 1e-5 };

    // Step size for alpha and beta space
    const double alpha_step_size = (alpha_max - alpha_min) / (alpha_values - 1);
    const double beta_step_size = (beta_max - beta_min) / (beta_values - 1);

    // Size of state space
    const int S_cardinality = static_cast<int>(pow(m, n * k));

    // Calculate argmax_a(Q(a,s)) for each s once every convergence_check time steps to check for convergence
    const int convergence_check{ 100 };

    // Calculate results over last results_check time steps prior to convergence
    const int results_check{ 100000 };

    // Number of time steps needed for argmax_a(Q(a,s)) to be constant for each s to achieve convergence
    const int norm{ 100000 / convergence_check };

    // Number of time steps allowed per episode
    const int maxt{ 10000000 };

    // Each agent's action at time step t
    std::vector<int> actions_t(n, 0);

    // Each agent's price at time step t
    std::vector<double> prices_t(n, 0.0);

    // Each agent's quantity at time step t
    std::vector<double> quantities_t(n, 0.0);

    // Each agent's profits at time step t
    std::vector<double> profits_t(n, 0.0);

    // Each agent's price at time step t
    std::vector<double> Delta_t(n, 0.0);

    // SARSA agent's Delta averaged over last results_check time steps prior to convergence for all episodes
    std::vector<double> avg_Delta_0(E, 0.0);

    // Q-learning agent's Delta averaged over last results_check time steps prior to convergence for all episodes
    std::vector<double> avg_Delta_1(E, 0.0);

    // Each agent's Delta for each (alpha, beta) pair averaged across all episodes
    std::vector<std::vector<std::vector<double>>> Delta(alpha_values, std::vector<std::vector<double>>(beta_values, std::vector<double>(n, 0.0)));

    // Average convergence time for each (alpha, beta) combination
    std::vector<std::vector<std::vector<int>>> time_steps_avg(alpha_values, std::vector<std::vector<int>>(beta_values, std::vector<int>(1, 0.0)));

} //Global


