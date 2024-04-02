# Heterogeneous Multi-Agent Reinforcement Learning Pricing Competition

This project allows different algorithms to engage in Bertrand pricing competition to determine if they can achieve levels above the competitive Bertrand-Nash outcome and sustain them using reward-punishment (trigger) strategies.

I first ran experiments using the well-known SARSA and Q-learning reinforcement learning algorithms. This is coded in MATLAB and subsequently re-coded in C++ to gain more computing power.

I then moved beyond the tabular based algorithms to more sophisticated deep reinforcement learning algorithms: Twin Delayed Deep Deterministic Policy Gradient (TD3) and Soft-Actor Crtic (SAC). This allows for for continuous action space, a feature likely employed by firms in practice. This is coded using Python to make use of the popular deep learning library PyTorch.
