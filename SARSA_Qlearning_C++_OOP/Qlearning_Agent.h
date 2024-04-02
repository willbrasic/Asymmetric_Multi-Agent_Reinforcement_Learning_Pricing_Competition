
#ifndef QLEARNING_AGENT_H
#define QLEARNING_AGENT_H

#include "Agent.h"

class Qlearning_Agent: public Agent
{
private:
    int Qlearning_next_action;
public:
    // Constructor
    Qlearning_Agent(std::vector<std::vector<double>>& Q_matrix, unsigned seed);

    // Set SARSA action
    void set_action_t(int time_step, bool explore_condition, int current_state) override;

    // Set SARSA_next_action at time step t
    void set_next_action(int next_state, bool explore_condition) override;

    const int& get_next_action() const override { return Qlearning_next_action; }

 };


#endif //QLEARNING_AGENT_H
