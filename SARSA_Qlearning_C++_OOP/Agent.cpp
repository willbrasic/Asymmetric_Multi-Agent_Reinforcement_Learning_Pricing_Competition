

/////////////////////////////
// Methods for Agent.h
//////////////////////////////


#include "Agent.h"


Agent::Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed)
    : Q{Q_matrix}, gen(seed), dist_1(0, m - 1)
{
    this->seed = seed;
}

// Q-matrix update
void Agent::Q_matrix_update(int current_state,
                            int next_state,
                            int current_action,
                            int next_action,
                            double current_profit,
                            double delta,
                            double alpha)
{
    Q[current_state][current_action] = (1 - alpha) * Q[current_state][current_action] + alpha * (current_profit + delta * Q[next_state][next_action]);
}

bool Agent::convergence_check_fn(bool time_step_above_convergence_check)
{
    // Is argmax_2[s] == argmax_1[s] for each state s? Initialize to false
    bool same_argmax{ false };

    // Initialize sum for when argmax_2[s] == argmax_1[s]
    int argmax_sum{ 0 };

    // If time step is greater than convergence_check, then update argmax_2 for each agent
    if (time_step_above_convergence_check)
    {
        // Look for maximum Q-values in each state
        for (int s{ 0 }; s < S_cardinality; ++s)
        {
            int argmax{ 0 };
            double max_Q = std::numeric_limits<double>::lowest();
            for (int j{ 0 }; j + 1 < m; j += 2)
            {
                if (Q[s][j] > max_Q) {
                    max_Q = Q[s][j];
                    argmax = j;
                }
                if (Q[s][j + 1] > max_Q) {
                    max_Q = Q[s][j + 1];
                    argmax = j + 1;
                }
            }
            argmax_2[s] = argmax;
            if (argmax_2[s] == argmax_1[s])
                ++argmax_sum;
        }
    }
    // If time step is not greater than convergence_check, then update argmax_1 (this chunk only gets executed once)
    else
    {
        // Look for maximum Q-values in each state
        for (int s{ 0 }; s < S_cardinality; ++s)
        {
            int argmax{ 0 };
            double max_Q = std::numeric_limits<double>::lowest();
            for (int j{ 0 }; j + 1 < m; j += 2)
            {
                if (Q[s][j] > max_Q) {
                    max_Q = Q[s][j];
                    argmax = j;
                }

                if (Q[s][j + 1] > max_Q) {
                    max_Q = Q[s][j + 1];
                    argmax = j + 1;
                }
            }
            argmax_1[s] = argmax;
        }
    }

    // If argmax is the same as last check, then update same_argmax to be true
    if (argmax_sum == S_cardinality)
        same_argmax = true;

    // Move argmax_2 to the old argmax collection
    argmax_1 = argmax_2;

    return same_argmax;
}
