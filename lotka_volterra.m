clear all;

% x = coral cover
% y = COTS

% Define problem parameters
t_0 = 0;
t_end = 15;
xy_0 = [20; 20];            % x and y at time 0

% Solve using ode45
[t, xy] = ode45(@lotka, [t_0 t_end], xy_0);

% Plotting
figure(1), clf, hold on
plot(t, xy(:, 1), 'Linewidth', 2)
plot(t, xy(:, 2), 'Linewidth', 2)
title('Lotka-Volterra Model')
xlabel('Time (t)')
ylabel('Population Size')
legend('Coral Cover', 'COTS', 'Location', 'NorthEastOutside')


%--------------------------------------------------------------------------

function xy_dash = lotka(t, xy)
% Lotka-Volterra model

% Define parameters
alpha = 1;
beta = 0.01;
delta = 0.02;
gamma = 1;
k = 0.2;

% System of equations - two first order ODEs
% xy = [x; y]; 
xy_dash = [xy(1) * (alpha - beta*xy(2)); 
    xy(2) * (delta*xy(1) - gamma - k)];

end