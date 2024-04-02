%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SARSA and Q-learning Results
% William B. Brasic
% The University of Arizona
% wbrasic@arizona.edu
% Website:
% October 2023; Last Revision: 2 March 2024
%
% This project obtains results for SARSA and Q-learning
% agents engaging in price competition.
%
% Before executing script:
% 1. Ensure R is correct
% 2. Ensure version is correct
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

% Colors for plots
color_1 = [0 0 0];
color_2 = [0.3010 0.7450 0.9330];

% Number of episodes
E = 100;

% Run to compute results for
R = 1;

% Number of firms
n = 2;

% Version of SARSA and Q-learning simulation
version = 'Seed_Robustness';

% File name for storing final results averaged over *E* episodes
results_final_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Results_Final_', num2str(R), '.csv');

% File name for storing learning curve results
lc_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_', ...
    version, '_Learning_Curves_', num2str(R), '.csv');

% File name for storing action distribution results
action_distribution_file_name = strcat('SARSA_Qlearning_', version, '\SARSA_Qlearning_',...
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
% Simulation Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear console
clc;

% Read in final results
results_final = readtable(results_final_file_name);

% Mean convergence
mean_converge = results_final.MeanConvergence;

% Mean market and firm Delta
mean_market_Delta = results_final.Mean_Market_Delta;
mean_firm_Delta = [results_final.Mean_Firm_1_Delta; results_final.Mean_Firm_2_Delta];

% Mean market and firm profit
mean_market_r = results_final.Mean_Market_Profit;
mean_firm_r = [results_final.Mean_Firm_1_Profit; results_final.Mean_Firm_2_Profit];

% Mean market and firm quantity
mean_market_q = results_final.Mean_Market_Quantity;
mean_firm_q = [results_final.Mean_Firm_1_Quantity; results_final.Mean_Firm_2_Quantity];

% Mean market and firm price
mean_market_p = results_final.Mean_Market_Price;
mean_firm_p = [results_final.Mean_Firm_1_Price; results_final.Mean_Firm_2_Price];

% Mean market and firm revenue
mean_market_rvn = results_final.Mean_Market_Revenue;
mean_firm_rvn = [results_final.Mean_Firm_1_Revenue; results_final.Mean_Firm_2_Revenue];

% Mean consumer surplus
mean_cs = results_final.Mean_CS;

% Average results across *E* episodes
fprintf(1,'\n*********************************************************\n');
fprintf(1,'* OVERAL RESULTS                  ***********************\n');
fprintf(1,'* (Averaged across %4.0f episodes) ***********************\n',E);
fprintf(1,'*********************************************************\n');
fprintf(1,'\nAverage number of time steps until convergence: %1.0f\n', mean_converge);
fprintf(1,'\n                                 Firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'             Tot/Avg');
fprintf(1,'         %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nDelta           %1.4f', mean_market_Delta)
fprintf(1,'    %1.4f', mean_firm_Delta)
fprintf(1,'\nProfits         %1.4f', mean_market_r)
fprintf(1,'    %1.4f', mean_firm_r)
fprintf(1,'\nDemand          %1.4f', mean_market_q)
fprintf(1,'    %1.4f', mean_firm_q)
fprintf(1,'\nPrices          %1.4f', mean_market_p)
fprintf(1,'    %1.4f', mean_firm_p)
fprintf(1,'\nRevenue         %1.4f', mean_market_rvn)
fprintf(1,'    %1.4f', mean_firm_rvn)
fprintf(1,'\nCS              %1.4f', mean_cs)
fprintf(1,'\n---------------------------------------------------------\n');

% Average percentage change from competitive outcome across *E* episodes
fprintf(1,'\n*********************************************************\n');
fprintf(1,'* PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ************\n');
fprintf(1,'* (Averaged across %4.0f episodes)            ************\n',E);
fprintf(1,'*********************************************************\n');
fprintf(1,'\nAverage number of time steps until convergence: %1.0f\n', mean(mean_converge));
fprintf(1,'\n                                 firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'               Tot/Avg');
fprintf(1,'         %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nProfits         %1.2f%%', 100 * (mean_market_r - sum(comp_pi))/sum(comp_pi))
fprintf(1,'    %1.2f%%', 100 * (mean_firm_r - mean(comp_pi))/mean(comp_pi))
fprintf(1,'\nDemand          %1.2f%%', 100 * (mean_market_q - sum(comp_q))/sum(comp_q))
fprintf(1,'    %1.2f%%', 100 * (mean_firm_q - mean(comp_q))/mean(comp_q))
fprintf(1,'\nPrices          %1.2f%%', 100 * (sum(mean_firm_p) - sum(comp_p))/sum(comp_p))
fprintf(1,'    %1.2f%%', 100 * (mean_firm_p - mean(comp_p))/mean(comp_p))
fprintf(1,'\nRevenue          %1.2f%%', 100 * (mean_market_rvn - sum(comp_rvn)) / sum(comp_rvn))
fprintf(1,'     %1.2f%%', 100 * (mean_firm_rvn - mean(comp_rvn))/mean(comp_rvn))
fprintf(1,'\nCS             %1.2f%%', 100 * (mean_cs - comp_cs)/comp_cs)
fprintf(1,'\n---------------------------------------------------------\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read in learning curve results
lc = readtable(lc_file_name);

% Learning curve profit averaged across firms for each episode
r_lc = table2array(lc(:, 1:100));

% Learning curve consumer surplus for each episode
cs_lc = table2array(lc(:, 101:200));

% Replace zero values for profit learning curve using
% last non-zero element in epsisode *e* for firm i
for e = 1:size(r_lc, 2)
    % Replace zero values for r_lc using last non_zero element for
    % episode *e*
    last_non_zero_r_lc = find(r_lc(:, e) ~= 0, 1, 'last');
    r_lc(r_lc(:, e) == 0, e) = r_lc(last_non_zero_r_lc, e);
end

% Replace zero values for consumer surplus learning curve using
% last non-zero element in epsisode *e* for firm i
for e = 1:size(cs_lc, 2)
    % Replace zero values for cs_lc using last non_zero element for
    % episode e
    last_non_zero_cs_lc = find(cs_lc(:, e) ~= 0, 1, 'last');
    cs_lc(cs_lc(:, e) == 0, e) = cs_lc(last_non_zero_cs_lc, e);
end

% Learning curve Delta averaged across firms for each episode
Delta_lc = (r_lc - mean(comp_pi)) ./ (mean(coll_pi) - mean(comp_pi));

% Learning curve Delta averaged across episodes and then across firms
Delta_lc_avg = mean(Delta_lc, 2);

% Learning curve consumer surplus averaged across epsiodes
cs_lc_avg = mean(cs_lc, 2);

% Significance level for confidence intervals
alpha = 0.05;

% Lower and upper bounds for average Delta across firms
lower_bound_Delta_lc_avg = quantile(Delta_lc, alpha/2, 2);
upper_bound_Delta_lc_avg = quantile(Delta_lc, 1 - alpha/2, 2);

% Lower and upper bounds for consumer surplus for each episode
lower_bound_cs_lc_avg = quantile(cs_lc, alpha/2, 2);
upper_bound_cs_lc_avg = quantile(cs_lc, 1 - alpha/2, 2);

% Plot learning curves
figure;

% Plot Delta_lc_avg
subplot(2, 1, 1);
fill([1:length(Delta_lc_avg), fliplr(1:length(Delta_lc_avg))], ...
    [upper_bound_Delta_lc_avg', fliplr(lower_bound_Delta_lc_avg')], ...
    color_1, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
hold on;
plot(Delta_lc_avg, '-o', 'Color', color_1, 'LineWidth', 1.5, 'HandleVisibility','off');
yline(1, '--red', 'LineWidth', 2, 'DisplayName', 'Perfectly Collusive');
yline(0, '--blue', 'LineWidth', 2, 'DisplayName', 'Perfectly Competitive');
hold off;
title('Learning Curves', 'FontSize', 12);
subtitle('SARSA vs. Q-learning', 'FontSize', 12);
ylabel('\Delta');
ylim([min(Delta_lc_avg) - 0.1, 1 + 0.1]);
xticklabels([]);
legend('Location', 'Northwest');
grid on;

% Plot cs_lc_avg
subplot(2, 1, 2);
fill([1:length(cs_lc_avg), fliplr(1:length(cs_lc_avg))], ...
    [upper_bound_cs_lc_avg', fliplr(lower_bound_cs_lc_avg')], color_2, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;
plot(cs_lc_avg, '-o', 'Color', color_2, 'LineWidth', 1.5);
yline(coll_cs, '--red', 'LineWidth', 2);
yline(comp_cs, '--blue', 'LineWidth', 2);
hold off;
ylabel('Consumer Surplus');
ylim([coll_cs - 0.1, comp_cs + 0.1]);
xlabel('Time Step');
xtickangle(45);
xticklabels({'10K', '100K', '200K', '300K', '400K', '500K', '600K', '700K'});
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Action_Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read in action distribution data
action_distribution = readtable(action_distribution_file_name);

% Action distribution averaged over *E* episodes
action_store_avg = table2array(action_distribution(:,1:15));

% Initial Q-matrix
Qi0 = table2array(action_distribution(:,16:30));

% Highlight the bar that is the profit maximizing action for firm 1
highlighted_bar_1 = find(Qi0(1, :) == max(Qi0(1, :)));

% % Highlight the bar that is the profit maximizing action for firm 2
highlighted_bar_2 = find(Qi0(2, :) == max(Qi0(2, :)));

% Plot action distribution
figure;

% Plot action_store_avg for firm 1
subplot(2, 1, 1);
h1 = bar(action_store_avg(1, :), 'FaceColor', color_1, 'EdgeColor', 'black');
hold on;
h2 = bar(highlighted_bar_1 , action_store_avg(1, highlighted_bar_1 ), 'FaceColor', color_2);
hold off;
ax1 = gca;
ax1.YAxis.Exponent = 0;
title('Distribution of Executed Actions', 'FontSize', 12);
subtitle('SARSA', 'FontSize', 12);
ylabel('Count');
ylim([0, max(action_store_avg(1, :)) + 50000]);
yticklabels({'0', '100K', '200K', '300K', '400K'});
xticklabels([]);
legend([h2], {'Mean Profit Maximizing Action'}, 'Location', 'Northwest');
grid on;

% Plot action_store_avg for firm 2
subplot(2, 1, 2);
h1 = bar(action_store_avg(2, :), 'FaceColor', color_1, 'EdgeColor', 'black');
hold on;
h2 = bar(highlighted_bar_2, action_store_avg(2, highlighted_bar_2), 'FaceColor', color_2);
hold off;
ax2 = gca;
ax2.YAxis.Exponent = 0;
subtitle('Q-learning', 'FontSize', 12);
ylabel('Count');
ylim([0, max(action_store_avg(2, :)) + 50000]);
yticks(ax2, [0 100000 200000 300000 400000, 500000]); % Set y-ticks positions
yticklabels(ax2, {'0', '100K', '200K', '300K', '400K', '500K'}); % Set y-tick labels
xlabel('Action');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RP Test Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read in reward-punishment data for each deviator firm j
for j = 1:n
    p_rp_table = readtable(strcat(rp_file_name, ...
        num2str(j), '_', num2str(R), '.csv'));
    p_rp_array = table2array(p_rp_table);
    p_rp_avg(:, :, j) = p_rp_array;
end

% Plot reward-punishment scheme test
figure;

% Plot p_rp_avg when agent 1 deviates
subplot(2, 1, 1);
hold on;
plot(p_rp_avg(1, :, 1), '--o', 'Color', color_1, 'DisplayName', 'Deviating Agent (SARSA)');
plot(p_rp_avg(2, :, 1), '--o', 'Color', color_2, 'DisplayName', 'Non-Deviating Agent (Q-learning)');
yline(max(coll_p), '--red', 'LineWidth', 2, 'DisplayName', 'Perfectly Collusive');
yline(min(comp_p), '--blue', 'LineWidth', 2, 'DisplayName', 'Perfectly Competitive');
hold off;
title('Reward-Punishment Scheme Test', 'FontSize', 12);
subtitle('SARSA vs. Q-learning', 'FontSize', 12);
ylabel('Price');
ylim([min(comp_p) - 0.1, max(coll_p) + 0.1]);
xticklabels([]);
legend('Location', 'Northwest');
grid on;

% Plot p_rp_avg when agent 2 deviates
subplot(2, 1, 2);
hold on;
plot(p_rp_avg(2, :, 2), '--o', 'Color', color_2, 'DisplayName', 'Deviating Agent (Q-learning)');
plot(p_rp_avg(1, :, 2), '--o', 'Color', color_1, 'DisplayName', 'Non-Deviating Agent (SARSA)');
yline(max(coll_p), '--red', 'LineWidth', 2, 'DisplayName', 'Perfectly Collusive');
yline(min(comp_p), '--blue', 'LineWidth', 2, 'DisplayName', 'Perfectly Competitive');
hold off;
ylabel('Price');
ylim([min(comp_p) - 0.1, max(coll_p) + 0.1]);
xlabel('Time Step');
xtickangle(45);
legend('Location', 'Northwest');
grid on;
