function mu_s = calculate_cots_age1_reproduction(plot_gonad_weight)

% This function calculates mu_s or the percentage of age 1 COTS that can 
% reproduce, using Lucas (1984) and Babcock et al. (2016)

if nargin == 0
    plot_gonad_weight = 0;
end

% Initialise a time array - time is in months (age of cots)
t = 0:1:100;

% From Lucas (1984) where D = diameter in mm and t = time in months
D = 323 * (1 + 393 * exp(-0.294 * t)) .^ (-1);

% From Babcock et al. (2016) where W = gonad weight in grams i.e.
% reproductive output of cots
W = 3.384 * exp(0.0115 * D);

% Plot the gonad weight against age of cots (time in months)
if plot_gonad_weight == 1
    figure(1), clf, hold on
    plot(t, W, 'LineWidth', 2)
    xlabel('Age of COTS (months)', 'Interpreter', 'Latex', 'FontSize', 13)
    ylabel('Gonad weight (grams)', 'Interpreter', 'Latex', 'FontSize', 13)
end

% Calculate average gonad weight between 12 and 24 months
age_1_reproductive_output = mean(W(12:24));

% Calculate average gonad weight between 24 and 36 months
age_2_reproductive_output = mean(W(24:36));

% Ratio of age 1 reproductive output to age 2 reproductive output
mu_s = age_1_reproductive_output / age_2_reproductive_output;
