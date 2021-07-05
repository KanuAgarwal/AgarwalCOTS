function [t_vec, x, y] = simulate_reefs(num_reefs, t_end, params, initial_state, control_effort)

% This function simulates starfish and coral populations at n reefs for t 
% years. 
%
% The function inputs are:
%   num_reefs      = number of reefs to simulate
%   t_end          = final timestep i.e. number of years to run simlation
%   params         = a struct with parameter values for model equations
%   intial_state   = a struct with the initial state of the system i.e. the 
%                    initial coral cover and starfish populations
%   control_effort = an matrix (size n x t-1) with control effort at each
%                    reef over time, or a constant for constant or no
%                    control effort 
%
% The function outputs are:
%   t_vec = a vector of time (in years) that the simulation is run over
%   x     = a matrix of coral cover for each reef over time
%   y     = a matrix of starfish poopulations for each reef over time


% PARAMETERS
% Parameters in model equations
alpha_c = params.alpha_c;       % growth rate of coral 
beta_sc = params.beta_sc;       % mortality rate of coral from starfish
beta_cs = params.beta_cs;       % starfish response to coral
alpha_s = params.alpha_s;       % natural mortality rate of starfish
rho = params.rho;               % proliferation rate of starfish larvae 
K = params.K;                   % starfish larvae carrying capacity
r_c = params.r_c;               % coral larvae production rate
r_s = params.r_s;               % starfish larvae production rate
kappa_c = params.kappa_c;       % coral larvae connectivity matrix
kappa_s = params.kappa_s;       % starfish larvae connectivity matrix
A = params.A;                   % reef area i.e. coral carrying capacity

% Initial state of the system
x_0 = initial_state.x_0;        % initial coral cover
y_0 = initial_state.y_0;        % initial starfish population


% SETUP
% Initialise time vector
t_0 = 0;
t_vec = t_0:1:t_end;

% Initialise matrices for coral and starfish
x = zeros(num_reefs, length(t_vec));
y = zeros(num_reefs, length(t_vec));
x(:, 1) = x_0;
y(:, 1) = y_0;

% Matrix for control effort
% If there is no control effort, create a zero  matrix
if control_effort == 0
    k = zeros(num_reefs, length(t_vec)-1);
% Otherwise, leave the control effort matrix as is
else
    k = control_effort;
end

% Larval recruitment
% Keep track of percentage of larvae that die at each reef
kappa_c_dead = 1 - sum(kappa_c, 2);
kappa_s_dead = 1 - sum(kappa_s, 2);
% Matrices for larval recruitment before mortality
sigma = zeros(num_reefs, length(t_vec)-1);              % coral
tau = zeros(num_reefs, length(t_vec)-1);                % starfish
% Matrices for larval recruitment after mortality
phi = zeros(num_reefs, length(t_vec)-1);                % coral
gamma = zeros(num_reefs, length(t_vec)-1);              % starfish


% SOLVE
% Loop over and calculate population size through time
for t = t_0+1:t_end
    % Loop over and calculate population at each reef
    for i = 1:num_reefs
        % Calculate starfish larval recruitment
        for j = 1:num_reefs
            tau(i, t) = tau(i, t) + kappa_s(j, i) * y(j, t) * r_s;
        end

        % Use Beverton-Holt model for starfish survival
        gamma(i, t) = (rho * tau(i, t)) / (1 + (rho - 1)/K(i) * tau(i, t));
        
        % Check gamma isn't negative
        if gamma(i, t) < 0
            gamma(i, t) = 0;
        end
       
        % Calculate coral larval recuitment
        for j = 1:num_reefs
            sigma(i, t) = sigma(i, t) + kappa_c(j, i) * x(j, t) * r_c;
        end
        
        % Calculate coral settlement based on area left on reef
        if sigma(i, t) <= (A(i) - x(i, t))
            phi(i, t) = sigma(i, t);
        else
            phi(i, t) = A(i) - x(i, t);
        end
        
        % Calculate population sizes 
        y(i, t+1) = y(i, t) * (beta_cs * x(i, t) - alpha_s) ...
                        - k(i, t) * y(i, t) + gamma(i, t);            
        x(i, t+1) = x(i, t) * (alpha_c - beta_sc * y(i, t)) + phi(i, t);
        
        % Round no. of starfish to nearest integer
        y(i, t+1) = round(y(i, t+1));
        
        % Check no. of starfish isn't less than zero
        if y(i, t+1) < 0
            y(i, t+1) = 0;
        end

        % Check coral cover isn't less than zero
        if x(i, t+1) < 0
            x(i, t+1) = 0;
        end
        
        % Check coral cover isn't larger than the reef size
        if x(i, t+1) > A(i)
            x(i, t+1) = A(i);
        end
    end
end
