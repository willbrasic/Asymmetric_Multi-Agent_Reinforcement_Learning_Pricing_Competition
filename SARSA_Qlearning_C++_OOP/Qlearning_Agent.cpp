
#include "Qlearning_Agent.h"

// Using superclass constructor
Qlearning_Agent::Qlearning_Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed)
        : Agent(Q_matrix, seed)
{
    std::cout << "Qlearning Agent Created\n";
}

void Qlearning_Agent::set_action_t(int time_step, bool explore_condition, int current_state)
{
    // Explore
    if (time_step == 0 || explore_condition)
    {
        int action = dist_1(gen);
        action_t = action;
    }
    // Exploit
    else
    {
        double max_Q_val = std::numeric_limits<double>::lowest();
        for (int j{ 0 }; j < m; ++j)
        {
            if (Q[current_state][j] > max_Q_val)
            {
                max_Q_val = Q[current_state][j];
                action_t = j;
            }
        }
    }
}

void Qlearning_Agent::set_next_action(int next_state, bool explore_condition)
{
    // Get action yielding highest Q-value in the next state
    int next_action{ 0 };
    double max_Q_for_next_state = std::numeric_limits<double>::lowest();
    for (int j{ 0 }; j < m; ++j)
    {
        if (Q[next_state][j] > max_Q_for_next_state)
        {
            max_Q_for_next_state = Q[next_state][j];
            next_action = j;
        }
    }
    Qlearning_next_action = next_action;
}


