# Heterogeneous Multi-Agent Reinforcement Learning Pricing Competition

Hi there! Thank you for checking out my repository! This README.md file gives
details on the Heterogeneous Multi-Agent Reinforcement Learning Pricing Competition project
of the repository.

## Overview

This repository contains code for allowing two distinct tabular based reinforcement learning methods,
SARSA and Q-learning, in a Bertrand-Markov pricing game. The purpose of this project is to determine if such
algorithms engaging in an unknown environment can learn collusive outcomes and sustain them in equilibrium
using trigger strategies. The environment is coded in MATLAB and subsequently re-coded in C++
to gain more computing power and make use of object-oriented programming (OOP) design.
This project is largely inspired by the seminal Calvano et al. (2020) paper.

My results indicate that, indeed, simple asymmetric reinforcement learning algorithms
can learn anti-competitive outcomes by interacting with each other over time and such
outcomes can be sustained via using reward-punishment schemes.

## Table of Contents

- [Heterogeneous Multi-Agent Reinforcement Learning Pricing Competition](#project-name)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [Getting Started](#getting-started)
    - [Description of Repository](#description-of-repository)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Training Curves](#training_curves)
  - [Results](#results)

  ## Getting Started

  Below are some instructions on how to get the project up and running.

  ### Description of Repository

  The Repository contains three main folders: SARSA_Qlearning_Matlab,
  SARSA_Qlearning_C++_OOP, and SARSA_Qlearning_C++_Heat_Map.

  SARSA_Qlearning_Matlab contains the code SARSA_Qlearning_Base_Case.m  to
  replicate the baseline results in MATLAB as well as the plotting script
  SARSA_Qlearning_Results.m to generate the plots presented in the paper.

  SARSA_Qlearning_C++_OOP contains a much neater implementation of
  SARSA_Qlearning_Base_Case.m used in MATLAB. In this version, I leverage
  object-oriented programming (OOP) design in C++ to create a concise
  version of the environment. Given the increased computing power of C++,
  this code runs roughly 50 times faster than the original MATLAB code when
  optimizing the compiler for speed by using -0fast (see CMakeLists.txt).

  Lastly, SARSA_Qlearning_C++_Heat_Map constructs the environment in C++
  using procedural programming with a focus on optimizing for speed. The purpose
  of the code in this folder is to run the experiment over 10,000 different
  combinations of the learning rate α and the epsilon-greedy
  experimentation parameter β. The results for each parameter combination are
  subsequently averaged over one-hundred episode runs and these results are
  used to construct heat maps.

  ### Prerequisites

  The main dependencies are MATLAB and C++.

  ### Installation

  ```bash
  # Clone the repository
  git clone https://github.com/willbrasic/Heterogeneous_MARL_Pricing_Competition.git

  # Navigate to the project directory
  cd Heterogeneous_MARL_Pricing_Competition
  ```

  ## Training Curves

  Here are training curves for the profit measure Δ (defined in the paper) and
  consumer surplus:

  ![Picture 1](https://github.com/willbrasic/Heterogeneous_MARL_Pricing_Competition/blob/main/Heterogeneous_MARL_Pricing_Competition_Pictures/SARSA_Qlearning_Base_Case_Learning_Curves.png)

  
