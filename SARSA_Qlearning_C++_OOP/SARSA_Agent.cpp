

/////////////////////////////
// Methods for SARSA_Agent.h
//////////////////////////////


#include "SARSA_Agent.h"

// Using superclass constructor
SARSA_Agent::SARSA_Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed)
    : Agent(Q_matrix, seed)
{
    std::cout << "SARSA Agent Created\n";
}

// Set SARSA agent's action at time step t
void SARSA_Agent::set_action_t(int time_step, bool explore_condition, int current_state)
{
    if (time_step == 0)
    {
        int action = dist_1(gen);
        action_t = action;
    }
    else
        action_t = SARSA_next_action;
}

// Set SARSA agent's action at time step t + 1
void SARSA_Agent::set_next_action(int next_state, bool explore_condition)
{
    // Explore
    if (explore_condition)
    {
        int action = dist_1(gen);
        SARSA_next_action = action;
    }
    // Exploit
    else
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
        SARSA_next_action = next_action;
    }
}

