
/////////////////////////////////////////////////////
// SARSA and Q-learning Bertrand-Markov Game
//
// William B. Brasic
// The University of Arizona
// wbrasic@arizona.edu
// Website:
// December 2023, Last revision: 25 March 2024
//
// This project allows SARSA and Q-learning agents
// to engage in a Bertrand-Markov pricing game.
// I leverage OOP principles to create a
// cleaner code. THIS IS STILL A WORK IN PROGRESS.
// THERE IS A BUG SOMEWHERE THAT I HAVE TO FIND.
// I GET GET SIGNIFICANTLY DIFFERENT
// RESULTS ACROSS EPISODES.
/////////////////////////////////////////////////////


////////////////////////////////////
// Header Files
////////////////////////////////////


#include "Agent.h"
#include "SARSA_Agent.h"
#include "Qlearning_Agent.h"
#include "Global.h"
#include "Functions.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>


using namespace Global;


int main()
{


    ////////////////////////////////////
    // Preliminaries
    ////////////////////////////////////

    // Discrete uniform distribution over {0, 1, ..., m-1} to draw actions from
    std::uniform_int_distribution<> dist_1(0, m - 1);

    // Continuous uniform distribution over [0, 1] for epsilon-greedy exploration
    std::uniform_real_distribution<> dist_2(0, 1);

    // Suppress scientific notation and round output to four decimal places
    std::cout << std::fixed << std::setprecision(4);


    ////////////////////////////////////
    // main Initialization
    ////////////////////////////////////


    // Create action space
    std::vector<double> action_space = populate_action_space_fn(m, A_min, A_step_size);

    // Variables for each state
    std::vector<std::vector<double>> prices = action_space_cartesian_product_fn(action_space, n, m);
    std::vector<std::vector<double>> quantities = quantity_for_states_fn(prices, a, n, a_0, mu, S_cardinality);
    std::vector<double> revenues = revenue_for_states_fn(prices, quantities, S_cardinality, n);
    std::vector<double> cs = cs_for_states_fn(prices, a, mu, a_0, n, S_cardinality);
    std::vector<std::vector<double>> profits = profits_for_states_fn(prices, quantities, mc, S_cardinality, n);

    // Average profit across states when firm sets one of the m prices
    std::vector<std::vector<double>> avg_profits_for_states = avg_profits_for_states_fn(prices, profits, action_space, n, m, delta);

    // Exploration value for each possible time step
    std::vector<double> minus_beta_times_t(max_t, 0.0);
    for (int t{ 0 }; t < max_t + 1; ++t)
        minus_beta_times_t[t] = exp(-beta * t);

    // Power vector used to get next state for all agents
    std::vector<int> power_vector(n, 0.0);
    for (int i{ 0 }; i < n; ++i)
        power_vector[0] = static_cast<int>(pow(m, i));


    ////////////////////////////////////
    // Episode Loop
    ////////////////////////////////////


    // Simulation header
    std::cout << "\n***********************\n";
    std::cout << "* Simulation Results *\n";
    std::cout << "***********************\n";


    // Start episode loop
    for (int e{ 0 }; e < E; ++e)
    {

        ////////////////////////////////////
        // Episode Loop Initialization
        ////////////////////////////////////


        // Start time measurement
        auto start_time = std::chrono::steady_clock::now();

        // Print episode e of E total
        std::cout << "Episode " << e + 1 << " of " << E << " total" << "\n\n";

        // Generate initial Q-matrix for each agent
        std::vector<std::vector<double>> initial_Q_matrix = initial_Q_matrix_fn(avg_profits_for_states, n, m, S_cardinality);

        // Create Agent objects from SARSA and Qlearning Agent subclasses
        SARSA_Agent SARSA_1(initial_Q_matrix, e);
        SARSA_Agent SARSA_2(initial_Q_matrix, e + 1);
        Qlearning_Agent Qlearning_1(initial_Q_matrix, e + 2);

        // Vector to store pointers to subclass instances
        std::vector<Agent*> Agents;
        Agents.push_back(&SARSA_1);
        Agents.push_back(&SARSA_2);
        Agents.push_back(&Qlearning_1);

        // Ensure number of agent objects is equal to Global::n
        int n_agents{ 0 };
        for (const auto& agent : Agents)
            ++n_agents;
        if (n_agents != n)
        {
            std::cerr << "Error: The number of agent objects created is not equal to n.\n";
            std::exit(EXIT_FAILURE);
        }

        // Ensure each agent object has correct value of Global::mc and Global::a
        for (int i{ 0 }; i < n ; ++i)
        {
            int mc = Agents[i]->get_mc();
            int a = Agents[i]->get_a();
            if (a != Global::a[i] || mc != Global::mc[i])
            {
                std::cerr << "Error: Agent " << i << " does not have correct value of a or mc as defined in Global.cpp.";
                std::exit(EXIT_FAILURE);
            }
        }

        // Vector to store preliminary actions only to get initial state
        std::vector<int> actions_preliminary(n, 0);
        for (int i{ 0 }; i < n; ++i)
        {
            std::mt19937& gen = Agents[i]->get_gen();
            actions_preliminary[i] = dist_1(gen);
        }

        // Use preliminary actions to get the initial state
        int current_state{ 0 };
        for (int i{ 0 }; i < n; ++i)
        {
            int action_preliminary = actions_preliminary[i];
            current_state += action_preliminary * power_vector[i];
        }

        // Initialize time step
        int t{ 0 };

        // Initialize convergence counter to 0
        int convergence_counter{ 0 };

        // Counter for learning curve results storage
        int learning_curve_counter{ 0 };


        ////////////////////////////////////
        // While Loop
        ////////////////////////////////////


        // Start while loop
        while ( (t < max_t) && (convergence_counter < norm) )
        {
            // Set actions for time step t
            for(auto& agent : Agents)
            {
                // Determine if agent explores or not
                std::mt19937& gen = agent->get_gen();
                bool explore_condition = dist_2(gen) < minus_beta_times_t[t];

                // Set agent's action at time step t
                agent->set_action_t( t, explore_condition, current_state);
                agent->add_action_t();
            }

            // Determine subsequent state
            int next_state{ 0 };
            for (int i{ 0 }; i < n; ++i)
            {
                int action_t = Agents[i]->get_action_t();
                next_state += action_t * power_vector[i];
            }

            // Set each agent's price at time step t
            for (auto& agent : Agents)
            {
                int action_t = agent->get_action_t();
                double price_t = action_space[action_t];
                agent->set_price_t(price_t);
                agent->add_price_t();
            }

            // Compute quantity denominator for time step t
            std::vector<double> exp_quantities(n, 0.0);
            double quantity_denominator = exp(a_0 / mu);
            for (int i{ 0 }; i < n; ++i)
            {
                int a = Agents[i]->get_a();
                double price_t = Agents[i]->get_price_t();
                exp_quantities[i] = exp((a - price_t) / mu);
                quantity_denominator += exp_quantities[i];
            }

            // Set each agent's quantity at time step t
            for (int i{ 0 }; i < n; ++i)
            {
                double quantity_t = exp_quantities[i] / quantity_denominator;
                Agents[i]->set_quantity_t(quantity_t);
                Agents[i]->add_quantity_t();
            }

            // Set each agent's revenue at time step t
            for (auto& agent : Agents)
            {
                double price = agent->get_price_t();
                double quantity = agent->get_quantity_t();
                double revenue = price * quantity;
                agent->set_revenue_t(revenue);
                agent->add_revenue_t();
            }

            // Set each agent's profit at time step t
            for (auto& agent : Agents)
            {
                int mc = agent->get_mc();
                double price_t = agent->get_price_t();
                double quantity_t = agent->get_quantity_t();
                double profit_t = (price_t - mc) * quantity_t;
                agent->set_profit_t(profit_t);
                agent->add_profit_t();
            }

            // Get each agent's next action for the Q-matrix update (only SARSA ends up actually taking it; Qlearning just finds the greedy action to use in update)
            for (const auto& agent : Agents)
            {
                std::mt19937& gen = agent->get_gen();
                bool explore_condition = dist_2(gen) < minus_beta_times_t[t + 1];
                agent->set_next_action(next_state, explore_condition);
            }

            // Q-matrix update for each agent
            for (const auto& agent : Agents)
            {
                int action_t = agent->get_action_t();
                int next_action = agent->get_next_action();
                double profit_t = agent->get_profit_t();
                agent->Q_matrix_update(current_state, next_state, action_t, next_action, profit_t, delta, alpha);
            }

            // Log of summation term in consumer surplus formula
            double exp_cs_sum = exp( a_0 / mu );
            for (const auto& agent : Agents)
            {
                int a = agent->get_a();
                double price_t = agent->get_price_t();
                exp_cs_sum += exp((a - price_t) / mu);
            }
            double log_exp_cs_sum = log(exp_cs_sum);

            // Consumer surplus at time step t in episode e
            double cs = mu * log_exp_cs_sum;
            cs_e.emplace_back(cs);

            // Add results for learning curve
            if ( (t % learning_curve_check == 0) && (t < 700001) && (t > 0) )
            {
                // Add learning curve profits which are averaged over agents and the last learning_curve_check time steps
                double agent_avg_profits{ 0.0 };
                for (const auto& agent : Agents)
                {
                    std::vector<double> profits = agent->get_profits_e();
                    agent_avg_profits += std::accumulate(profits.begin() + t - learning_curve_check, profits.begin() + t, 0.0);
                }
                double agents_avg_profits = agent_avg_profits / (n * learning_curve_check);
                learning_curve_profits[learning_curve_counter][e] = agents_avg_profits;

                ++learning_curve_counter;
            }

            // Check for convergence
            if ( (t % convergence_check == 0) && (t != 0) )
            {
                // Counter for when all agents' previous argmaxes are the same this time as last time they were checked (convergence_check time steps ago)
                int total_argmax_sum{ 0 };
                bool time_step_above_convergence_check = t > convergence_check;
                for (auto& agent : Agents)
                {
                    // Determine if agent has same argmax as that of convergence_check time steps ago
                    bool same_argmax = agent->convergence_check_fn(time_step_above_convergence_check);

                    // If agent does have the same argmax action for each state as last check, then update the counter
                    if (same_argmax)
                        ++total_argmax_sum;
                }
                // If each agents' argmaxes have not changed since last convergence check, then update the convergence counter
                if (total_argmax_sum == n)
                    ++convergence_counter;
                else
                    convergence_counter = 0;
            }

            // Transition to the next state
            current_state = next_state;

            // Increase time step counter
            ++t;

        // End of while loop
        }

        // End time measurement
        auto end_time = std::chrono::steady_clock::now();

        // Calculate the duration
        auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

        // Add seconds it took for episode e to converge
        timings[e] = duration_seconds;

        // Add number of time steps needed for episode e to converge (t - 1 b/c t gets incremented after convergence is checked for)
        time_steps[e] = t - 1;

        // Print the time it took for episode e to end
        std::cout << "\n\nLoop execution time: " << duration_seconds << " seconds" << "\n";

        // Check if converged in episode e
        if (t == max_t)
            std::cout << "Did not converge.";
        else
            // Print the time step at convergence (t - 1 b/c t gets incremented after convergence is checked for)
            std::cout << "Time steps until convergence " << t - 1 << "\n\n";

        // Key variables for each agent averaged over the last results_check time steps prior to convergence
        for (int i{ 0 }; i < n; ++i)
        {
            // Prices
            std::vector<double> price_e = Agents[i]->get_prices_e();
            double avg_price_e = std::accumulate(price_e.begin() + t - results_check, price_e.begin() + t, 0.0) / results_check;
            prices_e[i][e] = avg_price_e;

            // Quantities
            std::vector<double> quantity_e = Agents[i]->get_quantities_e();
            double avg_quantity_e = std::accumulate(quantity_e.begin() + t - results_check, quantity_e.begin() + t, 0.0) / results_check;
            quantities_e[i][e] = avg_quantity_e;

            // Revenues
            std::vector<double> revenue_e = Agents[i]->get_revenues_e();
            double avg_revenue_e = std::accumulate(revenue_e.begin() + t - results_check, revenue_e.begin() + t, 0.0) / results_check;
            revenues_e[i][e] = avg_revenue_e;

            // Profits
            std::vector<double> profit_e = Agents[i]->get_profits_e();
            double avg_profit_e = std::accumulate(profit_e.begin() + t - results_check, profit_e.begin() + t, 0.0) / results_check;
            profits_e[i][e] = avg_profit_e;
        }

        // Average of consumer surplus over the last results_check time steps prior to convergence
        avg_cs_e[e] = std::accumulate(cs_e.begin() + t - results_check, cs_e.begin() + t, 0.0) / results_check;

    // End of episode loop
    }

    // Key variables averaged over each episode e for the last results_check time steps prior to convergence
    int total_elements = E;
    std::vector<double> avg_prices(n, 0.0), avg_quantities(n, 0.0), avg_revenues(n, 0.0), avg_profits(n, 0.0);
    for (int i{ 0 }; i < n; ++i)
    {
        avg_prices[i] = std::accumulate(prices_e[i].begin(), prices_e[i].end(), 0.0) / total_elements;
        avg_quantities[i] = std::accumulate(quantities_e[i].begin(), quantities_e[i].end(), 0.0) / total_elements;
        avg_revenues[i] = std::accumulate(revenues_e[i].begin(), revenues_e[i].end(), 0.0) / total_elements;
        avg_profits[i] = std::accumulate(profits_e[i].begin(), profits_e[i].end(), 0.0) / total_elements;
    }

    // Consumer surplus averaged over each episode e for the last results_check time steps prior to convergence
    double avg_cs = std::accumulate(avg_cs_e.begin(), avg_cs_e.end(), 0.0) / total_elements;

    // Average seconds and time steps until convergence across episodes
    double average_seconds = std::accumulate(timings.begin(), timings.end(), 0.0) / timings.size();
    int average_time_steps = static_cast<int>(std::accumulate(time_steps.begin(), time_steps.end(), 0.0) / time_steps.size());


    ////////////////////////////////////
    // Print Results
    ////////////////////////////////////


    std::vector<std::vector<double>> competitive_outcome = logit_competitive_outcome_fn(a, mc, mu, n, a_0, 150);
    for (const auto& vec : competitive_outcome)
    {
        for (const auto& j : vec)
        {
            std::cout << j << " ";
        }
        std::cout << "\n\n";
    }

    // Use print results function to print results averaged across all episodes
    print_results(avg_profits, avg_prices, avg_quantities, avg_revenues,
                  a, mc, avg_cs, average_seconds, mu, average_time_steps, n,
                  a_0, 150);

    // Write learning_curve_profits to a CSV file
    std::ofstream learning_curve_profits_file(R"(C:\Users\wbras\OneDrive\Desktop\UA\2nd_Year_Paper\2nd_Year_Paper_Code_C++\SARSA_Qlearning\SARSA_Qlearning_Results\Learning_Curve_Profits.csv)");
    for (const auto &vec : learning_curve_profits)
    {
        for (size_t v{ 0 }; v < vec.size(); ++v)
        {
            // Place element v of vec in the file
            learning_curve_profits_file << vec[v];

            // If v is not the last element of vec, then place a comma as we still have more of the row to write to the file
            if (v != vec.size() - 1)
                learning_curve_profits_file << ",";
        }
        // Start a new row in the file
        learning_curve_profits_file << "\n";
    }



    return 0;
}
