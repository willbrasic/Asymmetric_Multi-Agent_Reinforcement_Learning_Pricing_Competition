
#ifndef AGENT_H
#define AGENT_H

#include "Global.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

using namespace Global;


class Agent
{

protected:
    // Agent's Q matrix
    std::vector<std::vector<double>> Q;

    // Each agent's random seed
    unsigned seed;

    // Agent's action at time step t
    int action_t;

    // Random number generator
    std::mt19937 gen;

    // Discrete uniform distribution over {0, 1, ..., m-1} to draw actions from
    std::uniform_int_distribution<> dist_1;

    // Vectors to store each agent's Q-value maximizing action for all states (used to test for convergence)
    std::vector<int> argmax_1 = std::vector<int>(S_cardinality, 0);
    std::vector<int> argmax_2 = std::vector<int>(S_cardinality, 0);

    // Agent's marginal cost (ensure this aligns with Global.cpp)
    int mc { 1 };

    // Agent's product quality index (ensure this aligns with Global.cpp)
    int a { 2 };

    // Agents' actions for each time step in episode e
    std::vector<int> actions_e;

    // Agent's price at time step t
    double price_t;

    // Agents' actions for each time step in episode e
    std::vector<double> prices_e;

    // Agent's quantity at time step t
    double quantity_t;

    // Agents' quantities for each time step in episode e
    std::vector<double> quantities_e;

    // Agent's revenue at time step t
    double revenue_t;

    // Agents' revenues for each time step in episode e
    std::vector<double> revenues_e;

    // Agent's quantity at time step t
    double profit_t;

    // Agents' quantities for each time step in episode e
    std::vector<double> profits_e;

public:
    // Constructor
    Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed);

    // Destructor
    virtual ~Agent() {  };

    // Get Q-matrix
    const std::vector<std::vector<double>>& get_Q() const { return Q; }

    // Get seed
    const unsigned get_seed() const
    {
        return seed;
    }

    // Return Agent's random number generator
    std::mt19937& get_gen()
    {
        return gen;
    }

    // Get marginal cost
    const int& get_mc() const { return mc; }

    // Get product quality index
    const int& get_a() const { return a; }

    // Set action at time step t (overwritten in subclasses)
    virtual void set_action_t(int time_step, bool explore_condition, int current_state) = 0;

    // Get action at time step t
    const int& get_action_t() const { return action_t; }

    // Add action at time step t to actions_e
    void add_action_t() { actions_e.emplace_back(get_action_t()); }

    // Get actions_e
    const std::vector<int>& get_actions_e() const { return actions_e; }

    // Set price at time step t
    void set_price_t(double price) { price_t = price; }

    // Get price at time step t
    const double& get_price_t() const { return price_t; }

    // Add price at time step t to actions_e
    void add_price_t() { prices_e.emplace_back(get_price_t()); }

    // Get prices_e
    const std::vector<double>& get_prices_e() const { return prices_e; }

    // Set quantity at time step t
    void set_quantity_t(double quantity) { quantity_t = quantity; }

    // Get quantity at time step t
    const double& get_quantity_t() const { return quantity_t; }

    // Add quantity at time step t to actions_e
    void add_quantity_t() { quantities_e.emplace_back(get_quantity_t()); }

    // Get quantities_e
    const std::vector<double>& get_quantities_e() const { return quantities_e; }

    // Set revenue at time step t
    void set_revenue_t(double revenue) { revenue_t = revenue; }

    // Get revenue at time step t
    const double& get_revenue_t() const { return revenue_t; }

    // Add revenue at time step t to actions_e
    void add_revenue_t() { revenues_e.emplace_back(get_revenue_t()); }

    // Get revenues_e
    const std::vector<double>& get_revenues_e() const { return revenues_e; }

    // Set profit at time step t
    void set_profit_t(double profit) { profit_t = profit; }

    // Get profit at time step t
    const double& get_profit_t() const { return profit_t; }

    // Add profit at time step t to profits_e
    void add_profit_t() { profits_e.emplace_back(get_profit_t()); }

    // Get profits_e
    const std::vector<double>& get_profits_e() const { return profits_e; }

    // Set next_action (to be overwritten by derived classes)
    virtual void set_next_action(int next_state, bool explore_condition) = 0;

    // Get next_action
    virtual const int& get_next_action() const = 0;

    // Q-matrix update
    void Q_matrix_update(int current_state,
                         int next_state,
                         int current_action,
                         int next_action,
                         double current_profit,
                         double delta,
                         double alpha);

    bool convergence_check_fn(bool time_step_above_convergence_check);

};


#endif //AGENT_H
