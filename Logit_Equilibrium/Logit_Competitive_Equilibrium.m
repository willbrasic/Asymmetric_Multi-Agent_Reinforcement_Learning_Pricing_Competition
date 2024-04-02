% Title: Logit Competitive Equilibrium Solver Using Fixed Point Iteration
% Author: William B. Brasic and Matthijs Wildenbeest
% Work Address: The University of Arizona
% Email: wbrasic@arizona.edu and wildenbeest@arizona.edu
% Website: 
% October 2023; Last revision: 11 December 2023


%------------- BEGIN CODE --------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not show warnings
warning off all;   

% Numbers are rounded
format longG;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Primitives   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of products/firms
n = 2;

% Product quality index (Demand characteristic influencing utility)
a = 2*ones(n, 1);

% Marginal costs for each firm
c = ones(n, 1);

% Horizontal differentiation index (Scale parameter of T1EV distributed)
mu = 1/4;

% Inverse index of aggregate demand (Mean utility of outside option)
a0 = 0;

% Initial prices for each firm
p0 = 2.5*ones(n, 1);

% Assign pj0 to the initialized prices
pj0 = p0;

% Initialize the max distance of the two prices
norm = 1;

% Initialize the mean distance of the two prices
avgnorm = 1;

% When max distance is smaller than this number we converged
tol = 1e-8;

% "Dampening" parameter
eta = 0.5;

% Counter for number of convergence iterations
i = 0;

% While both the max and average distances are greater than tol
while norm > tol && avgnorm > tol

    % Set pj equal to the initialized prices
    pj = pj0;
    
    % Compute demand for each firm using pj
    comp_q = exp((a-pj)./mu)/(exp(a0./mu)+sum(exp((a-pj)./mu)));

    % Compute prices (unique symmetric eq. price) with eta parameter 
    pj0 = eta*pj+(1-eta)*(c+mu./(1-comp_q));
    
    % Obtain the distance between old prices pj and new prices pj0
    t = abs(pj0-pj);

    % Get the max distance out of the two pricing distances
    norm = max(t);

    % Get the average distance out of the two pricing distances
    avgnorm = mean(t);

    % Increase the counter
    i = i + 1;

    % Print the counter and the average distance out of the two prices
    %fprintf(1,'iteration %1.0f---avgnorm %1.12f.\n',i,avgnorm);

    % If we get to 50 loops, we probably are not convering so quit the
    % program
    if i == 1000 
        fprintf(1,'Probably not converging---quiting fixed point routine (avgnorm = %1.6f).\n',avgnorm)
        norm = tol; avgnorm = tol; 
    end
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Competitive prices across the firms 
comp_p = pj0;

% Competitive demand across the firms
comp_q = exp((a-comp_p)./mu)/(exp(a0./mu)+sum(exp((a-comp_p)./mu)));

% Competitive revenue across the firms
comp_rvn = comp_p.*comp_q;

% Competitive average profits across the firms
comp_pi = (comp_p-c).*comp_q;

% Competitive consumer surplus
comp_cs = mu.*log(sum(exp((a-comp_p)./mu))+exp(a0./mu));

fprintf(1,'\n****************************************************\n');
fprintf(1,'********** OVERALL RESULTS (Competitve) *************\n');
fprintf(1,'****************************************************\n');
fprintf(1,'\nNumber of repetitions until convergence: %1.0f\n',i);
fprintf(1,'\n                                 Firms                 \n');
fprintf(1,'                       ----------------------------------\n');
fprintf(1,'             Tot/Avg');
fprintf(1,'         %1.0f', [1:n]')
fprintf(1,'\n---------------------------------------------------------\n');
fprintf(1,'\nprofits         %1.4f', sum(comp_pi))
fprintf(1,'    %1.4f', comp_pi)
fprintf(1,'\ndemand          %1.4f', sum(comp_q))
fprintf(1,'    %1.4f', comp_q)
fprintf(1,'\nprices          %1.4f', mean(comp_p))
fprintf(1,'    %1.4f', comp_p)
fprintf(1,'\nrevenue         %1.4f', sum(comp_rvn))
fprintf(1,'    %1.4f', comp_rvn)
fprintf(1,'\nCS              %1.4f', comp_cs)
fprintf(1,'\n---------------------------------------------------------\n');


%------------- END OF CODE --------------




