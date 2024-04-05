/////////////////////////////////////////////////////
// Q-learning and SARSA Pricing Competition
// William B. Brasic
// The University of Arizona
// wbrasic@arizona.edu
// Website:
// December 2023, Last revision: 5 April 2024
//
// This project allows SARSA and Q-learning agents
// to engage in Bertrand-Markov pricing game to
// compute Δ values for 10,000
// (α, β) pairs averaged over 100,000 time steps
// prior to convergence for E = 100 episodes.
/////////////////////////////////////////////////////


//////////////////////////////
// Header Files
//////////////////////////////


#include "Functions.h"
#include "Global.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace Global;


//////////////////////////////
// Main Function
//////////////////////////////


int main()
{

    ////////////////////////////////////
    // Preliminaries
    ////////////////////////////////////


    // Set seed
    std::mt19937 gen(1024);

    // Create discrete uniform distribution over {0, 1, ..., m-1} to draw initial action from
    std::uniform_int_distribution<> dist_1(0, m - 1);

    // Create uniform continuous uniform distribution over [0, 1] for epsilon-greedy exploration
    std::uniform_real_distribution<> dist_2(0, 1);

    // Fixed surpresses scientific notation and setprecision(4) rounds to four decimal places
    std::cout << std::fixed << std::setprecision(6);


    ////////////////////////////////////
    // main Initialization
    ////////////////////////////////////


    // Create actions space
    std::vector<double> A = populate_A(m, A_min, step_size);

    // Variables for each state
    std::vector<std::vector<double>> p = generate_combinations(A, n, m);
    std::vector<std::vector<double>> q = q_fn(p, a, n, a_0, mu, S_cardinality);
    std::vector<double> rvn = rvn_fn(p, q, S_cardinality, n);
    std::vector<double> cs = cs_fn(p, a, mu, a_0, n, S_cardinality);
    std::vector<std::vector<double>> r = r_fn(p, q, mc, S_cardinality, n);
    std::vector<std::vector<double>> Qi0 = Qi0_fn(p, r, A, n, m, delta);

    // Populate alpha vector
    std::vector<double> alpha_vector = alpha_fn(alpha_values, alpha_min, alpha_step_size);

    // Populate beta vector
    std::vector<double> beta_vector = beta_fn(beta_values, beta_min, beta_step_size);

    // Vector to store each agents' preliminary actions (these actions are only used to get the next state)
    std::vector<int> actions_preliminary(n, 0);

    // Vector to store value for exploration for all possible time steps
    std::vector<double> minus_beta_times_t(maxt, 0.0);

    // Value used to determine next state
    const int pow_value = static_cast<int>(pow(m, n - 1));

    // Vector storing all possible quantity numerators for each firm (I pull this out of the below loop to increase speed)
    std::vector<std::vector<double>> exp_quantity_vector(m, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < m; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            exp_quantity_vector[j][i] = exp( (a[i] - A[j]) / mu);
        }
    }
    double exp_outside_option = exp(a_0 / mu);


    ////////////////////////////////////
    // (alpha, beta) Loop
    ////////////////////////////////////


    // Simulation header
    std::cout << "\n***********************\n";
    std::cout << "* Simulation Results *\n";
    std::cout << "***********************\n";


    // Loop over values for learning rate
    for (int alpha{ 0 }; alpha < alpha_values; ++alpha)
    {
        // Current value of alpha
        double current_alpha = alpha_vector[alpha];

        // Loop over values for experimentation rate
        for (int beta{ 0 }; beta < beta_values; ++beta)
        {
            // Current value of beta
            double current_beta = beta_vector[beta];

            // Get exploration number for all possible time steps
            for (int t{ 0 } ; t <= maxt; ++t)
                minus_beta_times_t[t] = exp(-current_beta * t);

            // Print the current value of alpha
            std::cout << "Alpha " << alpha + 1 << " of " << alpha_values << " is: " << current_alpha << "\n";

            // Print the current value of beta
            std::cout << "Beta " << beta + 1 << " of " << beta_values << " is: " << current_beta << "\n\n";


            ////////////////////////////////////
            // Episode Loop
            ////////////////////////////////////


            // Vector to store how many time steps needed for convergence for each episode
            std::vector<int> time_steps(E, 0);

            // Start time measurement
            auto start_time = std::chrono::steady_clock::now();

            // Loop over all possible episodes
            for (int e{ 0 }; e < E; ++e)
            {
                // Initialize convergence counter to 0
                int convergence_counter{ 0 };

                // Initialize time step to 0
                int t{ 0 };

                // Vectors to store each agent Q-value maximizing action (used to test for convergence)
                std::vector<std::vector<int>> argmax_1(S_cardinality, std::vector<int>(n, 0));
                std::vector<std::vector<int>> argmax_2(S_cardinality, std::vector<int>(n, 0));

                // Populate each agent's Q-matrix
                std::vector<std::vector<std::vector<double>>> Q = Q_fn(Qi0, n, m, S_cardinality);

                // Initialize vector to store each agent's Delta values over the current episode
                std::vector<double> Delta_0;
                std::vector<double> Delta_1;

                // Draw each agents action from discrete uniform distribution over {0, 1, ..., m-1} to get initial states
                actions_preliminary[0] = dist_1(gen);
                actions_preliminary[1] = dist_1(gen);

                // Use preliminary actions to get the initial state
                int current_state = actions_preliminary[0] + (actions_preliminary[1]) * pow_value;


                ////////////////////////////////////
                // Time Step 0
                ////////////////////////////////////


                // I PULL TIME STEP 0 OUT OF THE BELOW WHILE LOOP SO IT DOESN'T HAVE TO CHECK FOR TIME STEP 0 AT EACH ITERATION
                // THE PURPOSE OF THIS IS TO INCREASE RUNTIME SPEED ALTHOUGH IT INDUCES CODE REDUNDANCY


                // Get actions for time step 0
                actions_t[0] = dist_1(gen);
                actions_t[1] = dist_1(gen);

                // Determine subsequent state
                int next_state = actions_t[0] + (actions_t[1]) * pow_value;

                // Get each agent's price at time step 0
                prices_t[0] = A[actions_t[0]];
                prices_t[1] = A[actions_t[1]];

                // Compute quantity
                double exp_quantity_0 = exp_quantity_vector[actions_t[0]][0];
                double exp_quantity_1 = exp_quantity_vector[actions_t[1]][1];
                double quantity_denominator = exp_outside_option + exp_quantity_0  + exp_quantity_1;
                quantities_t[0] = exp_quantity_0 / quantity_denominator;
                quantities_t[1] = exp_quantity_1 / quantity_denominator;

                // Get each agent's profit at time step t
                profits_t[0] = (prices_t[0] - mc[0]) * quantities_t[0];
                profits_t[1] = (prices_t[1] - mc[1]) * quantities_t[1];

                // Take uniform random action for next time step (since t = 0, SARSA would explore w/ probability of roughly 1)
                int SARSA_next_action = dist_1(gen);

                // SARSA update
                Q[current_state][actions_t[0]][0] = (1 - current_alpha) * Q[current_state][actions_t[0]][0] + current_alpha * (profits_t[0] + delta * Q[next_state][SARSA_next_action][0]);

                // Find maximum Q-value for Q-learning in next state (can do j+=2 since m is even)
                double max_Q = std::numeric_limits<double>::lowest();
                for (int j{ 0 }; j + 1 < m; j += 2)
                {
                    double q1 = Q[next_state][j][1];
                    double q2 = Q[next_state][j + 1][1];
                    double local_max = std::max(q1, q2);
                    if (local_max > max_Q)
                        max_Q = local_max;
                }
                // Q-learning update
                Q[current_state][actions_t[1]][1] = (1 - current_alpha) * Q[current_state][actions_t[1]][1] + current_alpha * (profits_t[1] + delta * max_Q);

                // Transition to next state
                current_state = next_state;

                // t is now 1
                ++t;


                ////////////////////////////////////
                // While Loop
                ////////////////////////////////////


                // Start while loop
                while ((t < maxt) && (convergence_counter < norm))
                {
                    // Exploration condition for getting current action
                    bool explore_condition_1 = dist_2(gen) < minus_beta_times_t[t];

                    // Exploration condition for getting next action
                    bool explore_condition_2 = dist_2(gen) < minus_beta_times_t[t + 1];

                    // For t > 0, action for i = 0 (SARSA agent) is already determined (SARSA_next_action)
                    actions_t[0] = SARSA_next_action;

                    // Explore for Q-learning agent
                    if (explore_condition_1)
                        // Take uniform random action
                        actions_t[1] = dist_1(gen);
                    // Exploit for Q-learning agent
                    else
                    {
                        max_Q = std::numeric_limits<double>::lowest();
                        for (int j{ 0 }; j < m; ++j)
                        {
                            if (Q[current_state][j][1] > max_Q)
                            {
                                max_Q = Q[current_state][j][1];
                                actions_t[1] = j;
                            }
                        }
                    }

                    // Determine subsequent state
                    next_state = actions_t[0] + (actions_t[1]) * pow_value;

                    // Get each agent's price at time step t
                    prices_t[0] = A[actions_t[0]];
                    prices_t[1] = A[actions_t[1]];

                    // Compute quantities
                    exp_quantity_0 = exp_quantity_vector[actions_t[0]][0];
                    exp_quantity_1 = exp_quantity_vector[actions_t[1]][1];
                    quantity_denominator = exp_outside_option + exp_quantity_0  + exp_quantity_1;
                    quantities_t[0] = exp_quantity_0 / quantity_denominator;
                    quantities_t[1] = exp_quantity_1 / quantity_denominator;

                    // Get each agent's profit at time step t
                    profits_t[0] = (prices_t[0] - mc[0]) * quantities_t[0];
                    profits_t[1] = (prices_t[1] - mc[1]) * quantities_t[1];

                    // Each agent's Delta at time step t
                    for (int i{ 0 }; i < n; ++i)
                        Delta_t[i] = (profits_t[i] - competitive_profits[i]) / (collusive_profits[i] - competitive_profits[i]);
                    Delta_0.emplace_back(Delta_t[0]);
                    Delta_1.emplace_back(Delta_t[1]);

                    // Explore for SARSA's action in the next state
                    if (explore_condition_2)
                        // Take uniform random action for next time step
                        SARSA_next_action = dist_1(gen);
                    // Exploit for SARSA's action in the next state
                    else
                    {
                        // Find index of greedy action
                        int greedy = 0;
                        max_Q = std::numeric_limits<double>::lowest();
                        for (int j{ 0 }; j < m; ++j)
                        {
                            if (Q[next_state][j][0] > max_Q)
                            {
                                max_Q = Q[next_state][j][0];
                                greedy = j;
                            }
                        }
                        // Take greedy action at next time step
                        SARSA_next_action = greedy;
                    }

                    // SARSA update
                    Q[current_state][actions_t[0]][0] = (1 - current_alpha) * Q[current_state][actions_t[0]][0] + current_alpha * (profits_t[0] + delta * Q[next_state][SARSA_next_action][0]);

                    // Find maximum Q-value for Q-learning in next state (can do j+=2 since m is even)
                    max_Q = std::numeric_limits<double>::lowest();
                    for (int j{ 0 }; j + 1 < m; j += 2)
                    {
                        double q_1 = Q[next_state][j][1];
                        double q_2 = Q[next_state][j+1][1];
                        double local_max = std::max(q_1, q_2);
                        if (local_max > max_Q)
                        {
                            max_Q = local_max;
                        }
                    }

                    // Q-learning update
                    Q[current_state][actions_t[1]][1] = (1 - current_alpha) * Q[current_state][actions_t[1]][1] + current_alpha * (profits_t[1] + delta * max_Q);

                    // Transition to the next state
                    current_state = next_state;

                    // Check for convergence every convergence_check time steps
                    if (t % convergence_check == 0)
                    {
                        // If we have already checked for convergence once before, compare old argmax_a(Q(s,a)) to current argmax_a(Q(s,a))
                        if (t > convergence_check)
                        {
                            // Get actions that give maximum Q-values for each agent for each state
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_cardinality; ++s)
                                {
                                    int argmax{ 0 };
                                    double max_Q = std::numeric_limits<double>::lowest();
                                    for (int j{ 0 }; j + 1 < m; j+=2)
                                    {
                                        if (Q[s][j][i] > max_Q)
                                        {
                                            max_Q = Q[s][j][i];
                                            argmax = j;
                                        }
                                        if (Q[s][j + 1][i] > max_Q)
                                        {
                                            max_Q = Q[s][j + 1][i];
                                            argmax = j + 1;
                                        }
                                    }
                                    argmax_2[s][i] = argmax;
                                }
                            }
                            // Calculate how many times previous argmax calculation convergence_check time steps ago is the same as new argmax calculation
                            int sum_argmax{ 0 };
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_cardinality; ++s)
                                {
                                    if (argmax_2[s][i] == argmax_1[s][i])
                                        ++sum_argmax;
                                }
                            }
                            // If both agents' optimal actions have not changed, increase the convergence counter
                            if (sum_argmax == n * S_cardinality)
                                ++convergence_counter;
                            else
                                convergence_counter = 0;

                            // Set the old argmaxes to the new ones for next convergence check
                            argmax_1 = argmax_2;
                        }
                        // Calculate initial agent optimal actions at time step convergence_check
                        else
                        {
                            for (int i{ 0 }; i < n; ++i)
                            {
                                for (int s{ 0 }; s < S_cardinality; ++s)
                                {
                                    int argmax{ 0 };
                                    double max_Q = std::numeric_limits<double>::lowest();
                                    for (int j{ 0 }; j < m; ++j)
                                    {
                                        if (Q[s][j][i] > max_Q)
                                        {
                                            max_Q = Q[s][j][i];
                                            argmax = j;
                                        }
                                    }
                                    argmax_1[s][i] = argmax;
                                }
                            }
                        }
                    }

                    // Increment t
                    ++t;

                // End while loop
                }

                // Add time step at convergence to time steps storage vector (t - 1 b/c we increment t after convergence check)
                time_steps[e] = t - 1;

                // Delta for each agent averaged over the last results_check time steps prior to convergence for current episode
                avg_Delta_0[e] = std::accumulate(Delta_0.begin() + t - results_check, Delta_0.begin() + t, 0.0) / results_check;
                avg_Delta_1[e] = std::accumulate(Delta_1.begin() + t - results_check, Delta_1.begin() + t, 0.0) / results_check;


            // End episode loop
            }

            // Average time steps across each episode for the current (alpha, beta) pair
            time_steps_avg[alpha][beta][0] = static_cast<int>(std::round(std::accumulate(time_steps.begin(), time_steps.end(), 0.0) / E));

            // Store Delta for each agent averaged across each episode
            Delta[alpha][beta][0] = std::accumulate(avg_Delta_0.begin(), avg_Delta_0.end(), 0.0) / E;
            Delta[alpha][beta][1] = std::accumulate(avg_Delta_1.begin(), avg_Delta_1.end(), 0.0) / E;

            // End time measurement for beta loop
            auto end_time = std::chrono::steady_clock::now();

            // Calculate the duration
            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            // Print episode loop duration
            std::cout << "Episode loop duration in seconds: " << duration_seconds << "\n\n";

        // End beta loop
        }

    // End alpha loop
    }

    // Write Delta results to a .csv file
    std::ofstream Delta_File(R"(C:\Users\wbras\OneDrive\Desktop\UA\2nd_Year_Paper\2nd_Year_Paper_Code_C++\SARSA_Qlearning_Heat_Map\Delta_Results.csv)");
    for (int alpha{ 0 }; alpha < alpha_values; ++alpha)
    {
        for (int beta{ 0 }; beta < beta_values; ++beta)
        {
            Delta_File << alpha_vector[alpha] << "," << beta_vector[beta] << "," << Delta[alpha][beta][0] << "," << Delta[alpha][beta][1] << "\n";
        }
    }

    // Write time step results to a .csv file
    std::ofstream Time_Step_File(R"(C:\Users\wbras\OneDrive\Desktop\UA\2nd_Year_Paper\2nd_Year_Paper_Code_C++\SARSA_Qlearning_Heat_Map\Time_Step_Results.csv)");
    for (int alpha{ 0 }; alpha < alpha_values; ++alpha)
    {
        for (int beta{ 0 }; beta < beta_values; ++beta)
        {
            Time_Step_File << alpha_vector[alpha] << "," << beta_vector[beta] << "," << time_steps_avg[alpha][beta][0] << "\n";
        }
    }

    return 0;
}
