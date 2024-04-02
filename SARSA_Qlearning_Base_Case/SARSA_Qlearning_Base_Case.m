%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SARSA and Q-learning Base Case
% William B. Brasic
% The University of Arizona
% wbrasic@arizona.edu
% Website:
% October 2023; Last Revision: 16 March 2024
%
% This project simulates SARSA and Q-learning price competition
% for the base case model.
%
% Before executing script:
% 1. Ensure R is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear workspace
clear;

% Do not show warnings
warning off all;

% Numbers are rounded without scientific notation
format longG;

% Reset random number generator
rng(0,'twister');

% Number of episodes
E = 100;

% Run (for saving results)
R = 1;

% Version of SARSA and Q-learning simulation
version = 'Base_Case';

% File name for storing results averaged over *E* episodes
results_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Results_', num2str(R), '.csv');

% File name for storing final results averaged over *E* episodes
results_final_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Results_Final_', num2str(R), '.csv');

% File name for storing learning curve results
lc_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Learning_Curves_', num2str(R), '.csv');

% File name for storing action distribution results
action_distribution_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Action_Distribution_', num2str(R), '.csv');

% File name for storing reward-punishment (RP) scheme results
rp_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_RP_Firm_');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Demand Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add path to logit equilibrium solver scripts
addpath('Logit_Equilibrium');

% Solve for logit competitive equilibrium using fixed point iteration
run('Logit_Equilibrium\Logit_Competitive_Equilibrium.m');

% Solve for logit collusive equilibrium using fixed point iteration
run('Logit_Equilibrium\Logit_Collusive_Equilibrium.m');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Primitives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear console
clc;

% Model Parameters
n = 2;                 % Number of products/firms
a = 2 * ones(n, 1);    % Product quality index
a_0 = 0;               % Inverse index of aggregate demand
mu = 1/4;              % Horizontal differentiation index
c = 1 * ones(n, 1);    % Marginal costs

% Learning parameters
k = 1;                 % Memory (number of periods)
m = 15;                % Number of equally spaced price points
delta = 0.95;          % Discount factor
alpha = 0.15;          % Learning rate
beta = 1e-5;           % Experimentation rate

% Minimum and maximum prices
Amin = 1.0;
Amax = 2.1;

% Discretization of the action space A (m equally spaced points)
A = linspace(Amin, Amax, m);

% Size of state space
S_cardinality = m^(n * k);

% All possible combinations of actions for firms
As = A;
for i = 1:n-1
    As = combvec(As, A);
end

% Prices set by each firm in each state
p = As;

% Quantity of each firm in each state
q = exp((a - p) ./ mu) ./ (sum(exp((a - p) ./ mu)) + exp(a_0 ./ mu));

% Total revenue in the market
rvn = sum(p .* q);

% Consumer surplus
cs = mu .* log(sum(exp((a - p) ./ mu)) + exp(a_0 ./ mu));

% Profits for each firm in each state
r = (p - c) .* q;

% Find average profit of the states when firm i sets price j
for i = 1:n
    for j = 1:m
        Qi0(i, j) = mean(r(i, p(i, :) == A(j))) ./ (1 - delta);
    end
end

% Simulation header
fprintf(1,'\n**********************\n');
fprintf(1,'* SIMULATION RESULTS *\n');
fprintf(1,'**********************\n');

% Delete old file storing results of each episode
delete(results_file_name);

% Start new file storing results for each episode
results = fopen(results_file_name,'at');

% If file opened successfully, write file headers
if results ~= -1
    fprintf(results, 'Session\t');
    fprintf(results, 'Repetitions\t');
    fprintf(results, 'Profit\t');
    fprintf(results, 'Price\t');
    fprintf(results, 'Quantity\t');
    fprintf(results, 'Price_Market\t');
    fprintf(results, 'Revenue\t');
    fprintf(results, 'CS\t');
    fprintf(results, '\n');
    fclose(results);
end

% Calculate argmax(Q) once every *convergence_check* repetitions to check for convergence
convergence_check = 100;

% Number of time steps needed for argmax(Q) to be constant for convergence
norm = 100000 / convergence_check;

% Number of time steps allowed per episode
maxt = 10000001;

% Compute learning curve data every *lc_check* time steps
lc_check = 10000;

% Display results for episode *e* every *results_check* time steps
results_check = 100000;

% Time steps used for reward punishment (RP) scheme test
t_rp = 19;

% Initialize matrix to store taken actions
action_store = zeros(n, m, E);

% Initialize matrix to store profit learning curve results
r_lc = zeros(700000 / lc_check, E);

% Initialize matrix to store consumer surplus learning curve results
cs_lc = zeros(700000 / lc_check, E);

% Price vector to store evolution before and after RP test
p_rp = zeros(n, t_rp, n, E);

% Initialize Q-matrix
Q = zeros(S_cardinality, m, n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Episode Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Start looping over episodes
for e = 1:E

    % Start timer
    tic

    % Clear variables at the start of each episode
    clear action cs expl p r q state

    % Each firm's Q-matrix is initialized with average profits when setting each of the m possible prices
    for i = 1:n
        Q(:, :, i, e) = ones(S_cardinality, 1) * Qi0(i, :);
    end

    % Randomly determine intitial actions
    for i = 1:n
        a_t0(i, 1) = randperm(m, 1);
    end

    % Determine initial state
    temp = a_t0(1, 1);
    for i = 2:n
        % Modify temp in some way that it lies in {1,2,...,S_cardinality}
        temp = temp + (a_t0(i, 1) - 1) * (m^(i - 1));
    end

    % Get initial state which is a function of the intitial action
    state(1) = temp;

    % Initalize convergence counter to 0
    convergence_count = 0;

    % Initialize time step to 1
    t = 1;

    % Initialize a counter for plotting training curves
    lc_count = 1;

    % Print episode e of E total
    fprintf('\nEpisode %1.0f of %1.0f total. Time Step: ', [e;E]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin Time Step Loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % While max time steps not reached and convergence not achieved
    while ((t < maxt) && (convergence_count < norm))

        % Get actions for each firm
        for i = 1:n
            % Only pick action for i = 1 (SARSA) in t = 1 b/c it's done below
            % for the remainder of time steps
            if ((i == 1 && t == 1) || (i == 2))
                % Explore
                if rand(1, 1) < exp(-beta * t)
                    % Take uniform random action
                    action(i, t) = randperm(m, 1);
                % Exploit
                else
                    % Take greedy action
                    action(i, t) = find(Q(state(t), :, i, e) == max(Q(state(t), :, i, e)), 1);
                end
            end
            % Increment counter for action taken
            action_store(i, action(i, t), e) = action_store(i, action(i, t), e) + 1;
        end

        % Determine subsequent state
        temp = action(1, t);
        for i = 2:n
            % Modify temp in some way that it lies in {1,2,...,S_cardinality}
            temp = temp + (action(i, t) - 1) * (m^(i - 1));
        end

        % Subsequent state which is a function of prior period action
        state(t + 1) = temp;

        % Loop over firms to get their price
        for i = 1:n
            % Set firm i's price for time step t
            p(i, t) = A(action(i, t));
        end

        % Quantity for each firm at time step t
        q(:, t) = exp((a - p(:, t)) ./ mu) ./ (sum(exp((a - p(:, t)) ./ mu)) + exp(a_0 ./ mu));

        % Consumer surplus at time step t
        cs(t) = mu .* log(sum(exp((a - p(:, t)) ./ mu)) + exp(a_0 ./ mu));

        % Profits for each firm at time step t
        r(:, t) = (p(:, t) - c) .* q(:, t);

        % Store data for learning curves every lc_check time steps
        if ((mod(t, lc_check) == 0) && (t <= 700001))
            % Store mean profits across firms
            r_lc(lc_count, e) = mean(mean(r(:, t-9999:t), 2));
            % Store mean consumer surplus
            cs_lc(lc_count, e) = mean(cs(t - 9999:t));
            % Increment the learning curve counter
            lc_count = lc_count + 1;
        end

        % Update Q-matrix
        for i = 1:n
            if i == 1
                % Get i's action in next state to use for SARSA update
                if rand(1, 1) < exp(-beta * (t + 1))
                    action(i, t + 1) = randperm(m, 1);
                else
                    action(i, t + 1) = find(Q(state(t + 1), :, i, e) == max(Q(state(t + 1), :, i, e)), 1);
                end
                % SARSA Update
                Q(state(t), action(i, t), i, e) = (1 - alpha) .* Q(state(t), action(i,t), i, e) + ...
                    alpha .* (r(i, t) + delta .* Q(state(t + 1), action(i, t + 1), i, e));
            else
                % Q-learning Update
                Q(state(t), action(i, t), i, e) = (1 - alpha) .* Q(state(t), action(i, t), i, e) + ...
                    alpha .* (r(i, t) + delta .* max(Q(state(t + 1), :, i, e)));
            end
        end

        % Check for convergence every *convergence_check* time steps
        if mod(t, convergence_check) == 0
            if (t / convergence_check) > 1
                % Find only the indices of the maximal actions for each state for each firm
                [~, amax2] = max(Q(:, :, :, e), [], 2);
                % Check if prior period's optimal actions for each firm
                % in each state is the same as this period
                if sum(sum(amax2 == amax1)) == n * S_cardinality
                    convergence_count = convergence_count + 1;
                else
                    convergence_count = 0;
                end
                % Set the old optimal action indices to the new ones
                amax1 = amax2;
            else
                % Find only the indices of the maximal actions for each state for each firm
                [~, amax1] = max(Q(:, :, :, e), [], 2);
            end
        end

        % Display counter every *results_check* time steps
        if mod(t, results_check) == 0
            if t > results_check
                % Delete previous counter display
                for j = 0:log10(t - 1)
                    fprintf('\b');
                end
            end
            % Print the time
            fprintf('%d', t);
            % Allows time for display to update
            pause(.05);
        end

        % Update time step
        t = t + 1;

    % End while loop
    end

    % Print new line
    fprintf('\n')

    % Stop timer
    tt = toc;

    % Averge profit over last *results_check* time steps of episode *e*
    r_e(:, e) = mean(r(1:n, t - results_check:t - 1), 2);

    % Average Delta over last *results_check* time steps of episode *e*
    Delta_e(:, e) = (r_e(:, e) - comp_pi) ./ (coll_pi - comp_pi);

    % Average price over last *results_check* time steps of episode *e*
    p_e(:, e) = mean(p(:, t - results_check:t - 1), 2);

    % Average quantity over last *results_check* time steps of episode *e*
    q_e(:, e) = mean(q(:, t - results_check:t - 1), 2);

    % Average consumer surplus over last *results_check* time steps of episode *e*
    cs_e(:, e) = mean(cs(:, t - results_check:t - 1), 2);

    % Average profit in market over last *results_check* time steps of episode *e*
    p_market_e(:, e) = mean(sum(p(:, t - results_check:t - 1) .* q(:, t - results_check:t - 1) ...
        ./ sum(q(:, t - results_check:t - 1))), 2);

    % Average revenue over last *results_check* time steps of episode *e*
    rvn_e(:, e) = mean((p(:, t - results_check:t - 1) .* q(:, t - results_check:t - 1)), 2);

    % Track iterations until convergence for episode *e*
    converge(e) = t - 1;

    % If t is the maximum time step (so the algorithm didn't converge)
    if t == maxt
        fprintf(1, 'Did not converge.\n');
    else
        fprintf(1, '# of time steps until convergance: %1.0f\n', t - 1);
    end

    % Results for episode *e* of *E* averaged over last *results_check* time steps
    fprintf(1,'It took about %1.0f minutes and %1.0f seconds until convergence\n', ...
        floor(tt / 60), round(tt - floor(tt / 60) * 60))
    fprintf(1,'\nAveraged across last 100,000 time steps:\n');
    fprintf(1,'                                 Firms                 \n');
    fprintf(1,'                       ----------------------------------\n');
    fprintf(1,'             Tot/Avg');
    fprintf(1,'         %1.0f', [1:n]')
    fprintf(1,'\n---------------------------------------------------------\n');
    fprintf(1,'\nDelta           %1.4f', mean(Delta_e(:, e)))
    fprintf(1,'    %1.4f', Delta_e(:, e))
    fprintf(1,'\nProfits         %1.4f', sum(r_e(:, e)))
    fprintf(1,'    %1.4f', r_e(:, e))
    fprintf(1,'\nDemand          %1.4f', sum(q_e(:, e)))
    fprintf(1,'    %1.4f', q_e(:, e))
    fprintf(1,'\nPrices          %1.4f', p_market_e(:, e))
    fprintf(1,'    %1.4f', p_e(:, e))
    fprintf(1,'\nRevenue         %1.4f', sum(rvn_e(:, e)))
    fprintf(1,'    %1.4f', rvn_e(:, e))
    fprintf(1,'\nCS              %1.4f', cs_e(:, e))
    fprintf(1,'\n---------------------------------------------------------\n');

    % Open file to write results to
    results = fopen(results_file_name, 'at');

    % If file opened successfully, write results of episode *e*
    if results ~= -1
        fprintf(results, '%1.0f\t', e);
        fprintf(results, '%1.0f\t', converge(e));
        fprintf(results, '%1.14f\t', mean(r_e(:, e)));
        fprintf(results, '%1.14f\t', mean(p_e(:, e)));
        fprintf(results, '%1.14f\t', mean(q_e(:, e)));
        fprintf(results, '%1.14f\t', p_market_e(:, e));
        fprintf(results, '%1.14f\t', sum(rvn_e(:, e)));
        fprintf(results, '%1.14f\t', cs_e(:, e));
        fprintf(results, '\n');
        fclose(results);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reward-Punishment Test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Agent j corresponds to the deviator
    for j = 1:n

        % Initialization for RP test
        Q_rp = Q;
        state_rp = state;
        action_rp = action;
        q_rp = q;
        r_rp = r;

        % Re-initialize t to convergence achieved for episode *e*
        t = converge(e);

        % Second index for p_rp
        p_rp_2 = 1 + t - converge(e);

        % While t is less than or equal to t_rp
        while t <= converge(e) + t_rp

            % Get the actions for each firm
            for i = 1:n
                % Ensure both firms have action at first time step
                if ((i == 1 && t == converge(e)) || (i == 2))
                    % Explore
                    if rand(1, 1) < exp(-beta * (t))
                        % Take uniform random action
                        action_rp(i, t) = randperm(m, 1);
                    % Exploit
                    else
                        % Take greedy action
                        action_rp(i, t) = find(Q_rp(state_rp(t), :, i, e) ...
                            == max(Q_rp(state_rp(t), :, i, e)), 1);
                    end
                end
            end

            % If SARSA agent is deviator and it's the second time step
            if ((j == 1) && (t == converge(e) + 1))
                action_rp(j, t + 1) = 7; % A(7) is roughly the Bertrand-Nash price
            end
            % If Q-learner is deviator and it's the third time step
            if ((j == 2) && (t == converge(e) + 2))
                action_rp(j, t) = 7;
            end

            % Determine subsequent state
            temp = action_rp(1, t);
            for i = 2:n
                % Modify temp in some way that it lies in {1,...,S_cardinality}
                temp = temp + (action_rp(i, t) - 1) * (m^(i - 1));
            end

            % Subsequent state which is a function of prior period action
            state_rp(t + 1) = temp;

            % Loop over firms to get their price
            for i = 1:n
                % Set firm i's price
                p_rp(i, p_rp_2, j, e) = A(action_rp(i, t));
            end

            % Quantity for each firm for time step t
            q_rp(:, t) = exp((a - p_rp(:, p_rp_2, j, e)) ./ mu) ...
                ./ (sum(exp((a - p_rp(:, p_rp_2, j, e)) ./ mu)) + exp(a_0 ./ mu));

            % Profits for each firm for time step t
            r_rp(:, t) = (p_rp(:, p_rp_2, j, e) - c) .* q_rp(:, t);

            % Ensure SARSA agent chooses actions properly
            for i = 1:n
                % If SARSA agent's turn and not second time step so above deviation
                % doesn't get overwritten or when j = 2 and i = 1 so SARSA agent
                % has an action selection at third time step
                if ((i == 1 && t ~= converge(e) + 1) || (i == 1 && j == 2))
                    % Explore
                    if rand(1, 1) < exp(-beta * (t + 1))
                        % Take uniform random action
                        action_rp(i, t + 1) = randperm(m,1);
                    % Exploit
                    else
                        % Take greedy action
                        action_rp(i, t + 1) = find(Q_rp(state_rp(t + 1), :, i, e) ...
                            == max(Q_rp(state_rp(t + 1), :, i, e)), 1);
                    end
                end
            end

            % Update Q-matrix
            for i = 1:n
                if i == 1
                    % SARSA Update
                    Q_rp(state_rp(t), action_rp(i, t), i, e) = (1 - alpha) ...
                        .* Q_rp(state_rp(t), action_rp(i, t), i, e) ...
                         + alpha .* (r_rp(i, t) + delta .* Q_rp(state_rp(t + 1), action_rp(i, t + 1), i, e));
                else
                    % Q-learning Update
                    Q_rp(state_rp(t), action_rp(i, t), i, e) = (1 - alpha) ...
                        *Q_rp(state_rp(t), action_rp(i, t), i, e) ...
                        + alpha.*(r_rp(i, t) + delta .* max(Q_rp(state_rp(t + 1), :, i, e)));
                end
            end

            % Update time step
            t = t + 1;

            % Update second index of p_rp
            p_rp_2 = p_rp_2 + 1;

        % End while loop
        end

    % End reward-punishment test
    end

% End episode loop
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Average results across *E* episodes
fprintf(1,'\n*********************************************************\n');
fprintf(1,'* OVERAL RESULTS                  ***********************\n');
fprintf(1,'* (Averaged across %4.0f episodes) ***********************\n',E);
fprintf(1,'*********************************************************\n');
fprintf(1,'\nAverage number of time steps until convergence: %1.0f\n', mean(converge));
fprintf(1,'\n                                 Firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'             Tot/Avg');
fprintf(1,'         %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nDelta           %1.4f', mean(mean(Delta_e), 2))
fprintf(1,'    %1.4f', mean(Delta_e, 2))
fprintf(1,'\nProfits         %1.4f', mean(sum(r_e), 2))
fprintf(1, '    %1.4f', mean(r_e, 2));
fprintf(1,'\nDemand          %1.4f', mean(sum(q_e), 2))
fprintf(1,'    %1.4f', mean(q_e, 2))
fprintf(1,'\nPrices          %1.4f', mean(p_market_e, 2))
fprintf(1,'    %1.4f', mean(p_e, 2))
fprintf(1,'\nRevenue         %1.4f', mean(sum(rvn_e), 2))
fprintf(1,'    %1.4f', mean(rvn_e, 2))
fprintf(1,'\nCS              %1.4f', mean(cs_e, 2))
fprintf(1,'\n---------------------------------------------------------\n');

% Average percentage change from competitive outcome across *E* episodes
fprintf(1,'\n*********************************************************\n');
fprintf(1,'* PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ************\n');
fprintf(1,'* (Averaged across %4.0f episodes)            ************\n',E);
fprintf(1,'*********************************************************\n');
fprintf(1,'\nAverage number of time steps until convergence: %1.0f\n', mean(converge));
fprintf(1,'\n                                 Firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'               Tot/Avg');
fprintf(1,'        %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nProfits         %1.2f%%', 100 * (mean(sum(r_e), 2) - sum(comp_pi))/sum(comp_pi))
fprintf(1,'    %1.2f%%', 100 * (mean(r_e, 2) - mean(comp_pi))/mean(comp_pi))
fprintf(1,'\nDemand          %1.2f%%', 100 * (mean(sum(q_e), 2) - sum(comp_q))/sum(comp_q))
fprintf(1,'    %1.2f%%', 100 * (mean(q_e, 2) - mean(comp_q))/mean(comp_q))
fprintf(1,'\nPrices          %1.2f%%', 100 * (mean(sum(p_e), 2) - sum(comp_p))/sum(comp_p))
fprintf(1,'    %1.2f%%', 100 * (mean(p_e, 2) - mean(comp_p))/mean(comp_p))
fprintf(1,'\nRevenue          %1.2f%%', 100 * (mean(sum(rvn_e), 2) - sum(comp_rvn)) / sum(comp_rvn))
fprintf(1,'     %1.2f%%', 100 * (mean(rvn_e, 2) - mean(comp_rvn))/mean(comp_rvn))
fprintf(1,'\nCS             %1.2f%%', 100*(mean(cs_e, 2) - comp_cs)/comp_cs)
fprintf(1,'\n---------------------------------------------------------\n');

% Delete old final results file storing results averaged over *E* episodes
delete(results_final_file_name);

% Start new final results file storing results averaged over *E* episodes
results_final = fopen(results_final_file_name, 'at');

% If file opened successfully, write file headers
if results_final ~= -1
    fprintf(results_final, 'Mean Convergence\t');
    fprintf(results_final, 'Mean_Market_Delta\t');
    for i = 1:n
        fprintf(results_final, strcat('Mean_Firm_', num2str(i), '_Delta\t'));
    end
    fprintf(results_final, 'Mean_Market_Profit\t');
    for i = 1:n
        fprintf(results_final, strcat('Mean_Firm_', num2str(i), '_Profit\t'));
    end
    fprintf(results_final, 'Mean_Market_Quantity\t');
    for i = 1:n
        fprintf(results_final, strcat('Mean_Firm_', num2str(i), '_Quantity\t'));
    end
    fprintf(results_final, 'Mean_Market_Price\t');
    for i = 1:n
        fprintf(results_final, strcat('Mean_Firm_', num2str(i), '_Price\t'));
    end
    fprintf(results_final, 'Mean_Market_Revenue\t');
    for i = 1:n
        fprintf(results_final, strcat('Mean_Firm_', num2str(i), '_Revenue\t'));
    end
    fprintf(results_final, 'Mean_CS\t');
    fprintf(results_final, '\n');
    fclose(results_final);
end

% Open file to write final results to
results_final = fopen(results_final_file_name, 'at');

% If file opened successfully, write final results averaged over *E* episodes
if results_final ~= -1
    fprintf(results_final, '%1.0f\t', mean(converge));
    fprintf(results_final, '%1.14f\t', mean(mean(Delta_e), 2));
    fprintf(results_final, '%1.14f\t', mean(Delta_e, 2));
    fprintf(results_final, '%1.14f\t', mean(sum(r_e), 2));
    fprintf(results_final, '%1.14f\t', mean(r_e, 2));
    fprintf(results_final, '%1.14f\t', mean(sum(q_e), 2));
    fprintf(results_final, '%1.14f\t', mean(q_e, 2));
    fprintf(results_final, '%1.14f\t', mean(p_market_e, 2));
    fprintf(results_final, '%1.14f\t', mean(p_e, 2));
    fprintf(results_final, '%1.14f\t', mean(sum(rvn_e), 2));
    fprintf(results_final, '%1.14f\t', mean(rvn_e, 2));
    fprintf(results_final, '%1.14f\t', mean(cs_e, 2));
    fprintf(results_final, '\n');
    fclose(results_final);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning Curves Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a table to store learning curve results
lc_data = table(r_lc, cs_lc);

% Write the table to a CSV file
writetable(lc_data, lc_file_name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Action Distribution Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Average taken actions across *E* episodes
action_store_avg = mean(action_store, 3);

% Create a table to store taken action distribution results
actions_distribution_data = table(action_store_avg, Qi0);

% Write the table to a CSV file
writetable(actions_distribution_data, action_distribution_file_name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reward-Punishment Test Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Average p_rp across *E* episodes
p_rp_avg = mean(p_rp, 4);

% Create a table to store RP results for each deviating firm j
for j = 1:n
    p_rp_data = table(p_rp_avg(:, :, j));
    writetable(p_rp_data, ...
        strcat(rp_file_name, num2str(j), '_', num2str(R), '.csv'));
end
