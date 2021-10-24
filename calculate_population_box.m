function [coral_box, starfish_age2_box, starfish_age1_box, starfish_age0_box] ...
    = calculate_population_box(t_end, C_y_f, N_y_2, N_y_1, N_y_0, num_reefs, lat, lon)

% Function which calculates the total population of coral and starfish in
% the outbreak initiation box.


% Initialise arrays for storing total populations
coral_box = zeros(1, t_end+1);
starfish_age2_box = zeros(1, t_end+1);
starfish_age1_box = zeros(1, t_end+1);
starfish_age0_box = zeros(1, t_end+1);

% Loop over each timestep
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_2_sum = 0;
    starfish_1_sum = 0;
    starfish_0_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f(i, t);
            starfish_2_sum = starfish_2_sum + N_y_2(i, t);
            starfish_1_sum = starfish_1_sum + N_y_1(i, t);
            starfish_0_sum = starfish_0_sum + N_y_0(i, t);
        end
    end
    % Assign the sum to the coral cover array
    coral_box(t) = coral_sum;
    starfish_age2_box(t) = starfish_2_sum;
    starfish_age1_box(t) = starfish_1_sum;
    starfish_age0_box(t) = starfish_0_sum;
end
