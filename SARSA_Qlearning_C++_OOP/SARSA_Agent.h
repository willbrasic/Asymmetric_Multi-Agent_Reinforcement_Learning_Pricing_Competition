
#ifndef SARSA_AGENT_H
#define SARSA_AGENT_H

#include "Agent.h"


class SARSA_Agent: public Agent
{
private:
    // SARSA's action in next time step
    int SARSA_next_action;

public:
    // Constructor
    SARSA_Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed);

    // Set SARSA agent's action at time step t
    void set_action_t(int time_step, bool explore_condition, int current_state) override;

    // Set SARSA_next_action at time step t
    void set_next_action(int next_state, bool explore_condition) override;

    const int& get_next_action() const override { return SARSA_next_action; }

};


#endif //SARSA_AGENT_H
