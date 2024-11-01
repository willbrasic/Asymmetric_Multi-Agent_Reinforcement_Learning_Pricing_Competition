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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logit Demand Equilibrium 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Solve for logit competitive equilibrium using fixed point iteration
run('Logit_Competitive_Equilibrium.m');

% Solve for logit collusive equilibrium using fixed point iteration
run('Logit_Collusive_Equilibrium.m');

% Read in learning curve results
learning_curve_profits = table2array(readtable("Learning_Curve_Profits.csv"));

% Replace zero values for profit learning curve using 
% last non-zero element in epsisode *e* for firm i
for e = 1:size(learning_curve_profits, 2)
    % Replace zero values for learning_curve_profits using last non_zero element for
    % episode *e*
    last_non_zero_learning_curve_profit = find(learning_curve_profits(:, e) ~= 0, 1, 'last');
    learning_curve_profits(learning_curve_profits(:, e) == 0, e) = learning_curve_profits(last_non_zero_learning_curve_profit, e); 
end

% Learning curve Delta averaged across firms for each episode 
learning_curve_Delta = (learning_curve_profits - mean(comp_pi)) ./ (mean(coll_pi) - mean(comp_pi));

% Average results across episodes
avg_learning_curve_Delta = mean(learning_curve_Delta, 2);

% Significance level for confidence intervals
alpha = 0.05;

% Lower and upper bounds for average Delta across firms
lower_bound_avg_learning_curve_Delta = quantile(learning_curve_Delta, alpha/2, 2);
upper_bound_avg_learning_curve_Delta = quantile(learning_curve_Delta, 1 - alpha/2, 2);

% Plot learning curves
figure;

% Plot Delta_lc_avg
fill([1:length(avg_learning_curve_Delta), fliplr(1:length(avg_learning_curve_Delta))], ...
    [upper_bound_avg_learning_curve_Delta', fliplr(lower_bound_avg_learning_curve_Delta')], ...
    color_1, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
hold on;
plot(avg_learning_curve_Delta, '-o', 'Color', color_1, 'LineWidth', 1.5, 'HandleVisibility','off');
yline(1, '--red', 'LineWidth', 2, 'DisplayName', 'Perfectly Collusive');
yline(0, '--blue', 'LineWidth', 2, 'DisplayName', 'Perfectly Competitive');
hold off;
title('Learning Curves', 'FontSize', 12);
subtitle('SARSA vs. Q-learning', 'FontSize', 12);
ylabel('\Delta');
ylim([min(avg_learning_curve_Delta) - 0.1, 1 + 0.1]);  
xticklabels([]);
legend('Location', 'Northwest');
grid on;


















