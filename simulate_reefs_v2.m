function [t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort)

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


% PARAMETERS ==============================================================
% Parameters in model equations
p_tilde = params.p_tilde;       % effect of fast-growing coral on COTS     
M_cots = params.M_cots;         % natural mortality of COTS 
p_1_f = params.p_1_f;           % effect of COTS on fast-growing coral
% p_1_m = params.p_1_m;           % effect of COTS on slow-growing coral
r_f = params.r_f;               % intrinsic growth rate of fast-growing coral
% r_m = params.r_m;               % intrinsic growth rate of slow-growing coral
K_f = params.K_f;               % carrying capacity of fast-growing coral
% K_m = params.K_m;               % carrying capacity of slow-growing coral
p_2_f = params.p_2_f;           % effect of COTS on fast-growing coral
% p_2_m = params.p_2_m;           % effect of COTS on slow-growing coral
% rho = params.rho;               % proliferation rate of starfish larvae 
% K = params.K;                   % starfish larvae carrying capacity
r_c = params.r_c;               % coral larvae production rate
r_s = params.r_s;               % starfish larvae production rate
omega_c = params.omega_c;       % coral larvae connectivity matrix
omega_s = params.omega_s;       % starfish larvae connectivity matrix
% A = params.A;                   % reef area i.e. coral carrying capacity

% Initial state of the system
C_0_f = initial_state.C_0_f;            % fast-growing coral
% C_0_m = initial_state.C_0_m;            % slow-growing coral
N_0_2 = initial_state.N_0_2;            % age 2+ COTS
N_0_1 = initial_state.N_0_1;            % age 1 COTS
N_0_0 = initial_state.N_0_0;            % age 0 COTS


% SETUP ===================================================================
% Initialise time vector
t_0 = 0;
t_vec = t_0:1:t_end;

% Initialise matrices for coral and starfish
C_y_f = zeros(num_reefs, length(t_vec));
% C_y_m = zeros(num_reefs, length(t_vec));
N_y_2 = zeros(num_reefs, length(t_vec));
N_y_1 = zeros(num_reefs, length(t_vec));
N_y_0 = zeros(num_reefs, length(t_vec));
C_y_f(:, 1) = C_0_f;
% C_y_m(:, 1) = C_0_m;
N_y_2(:, 1) = N_0_2;
N_y_1(:, 1) = N_0_1;
N_y_0(:, 1) = N_0_0;

% Initialise matrix for control effort
% If there is no control effort, create a zero  matrix
if control_effort == 0
    k = zeros(num_reefs, length(t_vec)-1);
% Otherwise, leave the control effort matrix as is
else
    k = control_effort;
end

% Initialise matrices for larval dispersal
% Matrices for larval recruitment before mortality
sigma = zeros(num_reefs, length(t_vec)-1);              % coral
tau = zeros(num_reefs, length(t_vec)-1);                % starfish
% Matrices for larval recruitment after mortality
phi = zeros(num_reefs, length(t_vec)-1);                % coral
gamma = zeros(num_reefs, length(t_vec)-1);              % starfish


% SOLVE ===================================================================
% Loop over and calculate population size through time
for t = t_0+1:t_end
    % Loop over and calculate population at each reef
    for i = 1:num_reefs
        % STARFISH ========================================================
        % Calculate starfish larval recruitment
        for j = 1:num_reefs
            tau(i, t) = tau(i, t) + omega_s(j, i) * N_y_2(j, t) * r_s;
        end

%         % Use Beverton-Holt model for starfish survival
%         gamma(i, t) = (rho * tau(i, t)) / (1 + K * tau(i, t));
%         gamma(i, t) = (rho * tau(i, t)) / (1 + (rho - 1)/K * tau(i, t));
        
%         % Check gamma isn't negative
%         if gamma(i, t) < 0
%             gamma(i, t) = 0;
%         end

%         % Calculate population sizes for age 0 COTS depending on year -
%         % based on Morello paper - testing base case model
%         if t == 1
%             N_y_0(i, t+1) = 1 + exp(4.292);
%         elseif t == 3
%             N_y_0(i, t+1) = exp(4.307) + 1;
%         else
%             N_y_0(i, t+1) = 2;
%         end
        
        % Calculate population sizes for age 0 COTS
        N_y_0(i, t+1) = tau(i, t);
        
        % Functions to help calculate age 1 and 2+ COTS population sizes
        f_of_C = 1 - p_tilde * (C_y_f(i, t) / (1 + C_y_f(i, t)));

        % Calculate population sizes for age 1 and 2+ COTS
        N_y_1(i, t+1) = N_y_0(i, t) * exp(-f_of_C * M_cots);
        N_y_2(i, t+1) = (N_y_1(i, t) + N_y_2(i, t)) * exp(-f_of_C * M_cots) ...
                            - k(i, t) * N_y_2(i, t);
        
        % CORAL ===========================================================
        % Calculate coral larval recuitment
        for j = 1:num_reefs
            sigma(i, t) = sigma(i, t) + omega_c(j, i) * C_y_f(j, t) * r_c;
        end
        
%         % Calculate coral settlement based on area left on reef
%         if sigma(i, t) <= (1 - C_y_f(i, t))
%             phi(i, t) = sigma(i, t);
%         else
%             phi(i, t) = 1 - C_y_f(i, t);
%         end
                        
        % Function for calculating different coral population sizes
        rho_y = exp(-5 * C_y_f(i, t) / K_f);
%         rho_y = 1 / (1 + exp(-70 * (C_y_f(i, t) / (K_f - 0.1))));
        
        % Fast-growing coral mortality from COTS 
        Q_y_f = (1-rho_y) * (p_1_f * (N_y_1(i, t) + N_y_2(i, t)) * C_y_f(i, t)) ...
                    / (1 + exp( -(N_y_1(i, t) + N_y_2(i, t)) / p_2_f));

%         % Slow-growing coral mortality from COTS
%         Q_y_m = rho_y * (p_1_m * (N_y_1(i, t) + N_y_2(i, t)) * C_y_m(i, t)) ...
%                     / (1 + exp( -(N_y_1(i, t) + N_y_2(i, t))/p_2_m) );

        % Calculate population size for fast-growing coral
        C_y_f(i, t+1) = C_y_f(i, t) + r_f * C_y_f(i, t) * (1 - C_y_f(i, t)/K_f) ...
                            - Q_y_f + sigma(i, t);
        
%         % Calculate population size for slow-growing coral
%         C_y_m(i, t+1) = C_y_m(i, t) + r_m * C_y_m(i, t) * (1 - C_y_m(i, t)/K_m) - Q_y_m;

        % CHECKS ==========================================================
        % Ensure populations aren't less than zero
        if C_y_f(i, t+1) < 0
            C_y_f(i, t+1) = 0;          % fast-growing coral
        end

%         if C_y_m(i, t+1) < 0
%             C_y_m(i, t+1) = 0;          % slow-growing coral
%         end
        
        if N_y_0(i, t+1) < 0
            N_y_0(i, t+1) = 0;          % age 0 COTS
        end
        
        if N_y_1(i, t+1) < 0
            N_y_1(i, t+1) = 0;          % age 1 COTS
        end
        
        if N_y_2(i, t+1) < 0            % age 2+ COTS
            N_y_2(i, t+1) = 0;
        end
        % =================================================================
    end
end

% END =====================================================================