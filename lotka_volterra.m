clear all;

% x = coral cover
% y = COTS

% Define problem parameters
t_0 = 0;
t_end = 15;
x_0 = 20;
y_0 = 20;
xy_0 = [x_0; y_0]; 

% Solve using ode45
[t_vals, xy] = ode45(@lotka, [t_0 t_end], xy_0);
x_vals = xy(:, 1);
y_vals = xy(:, 2);

% Plotting
figure(1), clf, hold on
plot(t_vals, x_vals, 'Linewidth', 2)
plot(t_vals, y_vals, 'Linewidth', 2)
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
x = xy(1);
y = xy(2);
xy_dash = [x * (alpha - beta*y); 
           y * (delta*x - gamma - k)];

end