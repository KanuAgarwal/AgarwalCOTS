function [num_reefs, lon, lat, reef_area, omega] = match_reefs_to_areas(lg, lt, con_matrix)

% This code matches the reef areas in DataSources/ReefGazette.csv to the
% reef coordinates in IdentifyKeySources/original_centroids and keeps the
% intersection of the two datasets. It also updates the larval dispersal
% matrices in IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7
% accordingly. 


% LOAD DATA ---------------------------------------------------------------
% =========================================================================
% Put lat and long data into matrix
reef_coord_data = [lg, lt];

% =========================================================================
% Load reef areas 
% View import options in csv 
opts = detectImportOptions('DataSources/ReefGazette.csv');

% Only select the variables that we want to import 
opts.SelectedVariableNames = [15:16, 4];

% Load in the data with the right variables
M = readmatrix('DataSources/ReefGazette.csv', opts);

% Replace all entries that have NaN with 0
M(isnan(M)) = 0;

% Remove rows that have all zeros
reef_area_data = M(all(M, 2), :);

% =========================================================================


% MATCH REEFS AND AREAS ---------------------------------------------------
% =========================================================================
% Find the smallest Euclidean distance i.e. for each reef with an area,
% find the closest reef in my model
[D, I_coord] = pdist2(reef_coord_data, reef_area_data(:, 1:2), ...
    'Euclidean', 'Smallest', 1);

% Sort from smallest to largest distance
[D_sorted, I_area_sorted] = sort(D);
I_coord_sorted = I_coord(I_area_sorted);

% =========================================================================
% Remove reefs that we couldn't find a match for from my model
unique_reefs = unique(I_coord_sorted);
reef_coord_data_v2 = reef_coord_data(unique_reefs, :);

% Remove data from connectivity matrices for reefs with no area data
omega = con_matrix(unique_reefs, unique_reefs);

% =========================================================================
% Find the smallest three distances i.e. for each reef in my model, find 
% the closest 3 reefs with an area
[D_v2, I_v2_area] = pdist2(reef_area_data(:, 1:2), reef_coord_data_v2, ...
    'Euclidean', 'Smallest', 3);

% =========================================================================
% Create array to store reefs matches
num_matched = 0;
I_v2_area_matched = zeros(1, 1);
I_v2_coord_matched = zeros(1, 1);
num_repeats = 0;
I_v2_area_repeats = zeros(1, 1);
I_v2_coord_repeats = zeros(1, 1);

% Loop over each observation from euclidean distance algorithm
for i = 1:length(I_v2_area)
    current_reef = I_v2_area(1, i);
    % If current reef hasn't already been matched and it doesn't repeat, save it
    if ismember(current_reef, I_v2_area_matched) == 0 && ismember(current_reef, I_v2_area_repeats) == 0
        num_matched = num_matched + 1; 
        I_v2_area_matched(num_matched) = current_reef;
        I_v2_coord_matched(num_matched) = i;
    % Otherwise, if current reef has already been matched or repeated
    else
        % If the reef hasn't repeated before 
        if ismember(current_reef, I_v2_area_repeats) == 0
            % Find previous entry in array
            index_remove = find(I_v2_area_matched == current_reef);
            % Add previous entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = I_v2_area_matched(index_remove);
            I_v2_coord_repeats(num_repeats) = I_v2_coord_matched(index_remove);
            % Remove entry from matched reefs
            I_v2_area_matched(index_remove) = [];
            I_v2_coord_matched(index_remove) = [];
            num_matched = num_matched - 1;
            % Add new entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = current_reef;
            I_v2_coord_repeats(num_repeats) = i;
        % Otherwise if the reef has repeated before
        else 
            % Add new entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = current_reef;
            I_v2_coord_repeats(num_repeats) = i;
        end
    end
end

% =========================================================================
% Get all the distances and indexes that repeat
D_v2_repeats = D_v2(:, I_v2_coord_repeats);
I_v2_area_repeats_all = I_v2_area(:, I_v2_coord_repeats);
I_v2_area_repeats_unique = unique(I_v2_area_repeats);

% Create new array for area data
I_v2_area_actual = I_v2_area(1, :);

% Find the closest reef for any repeats
for i = 1:length(I_v2_area_repeats_unique)
    % Find all repeats of each unique element
    current_repeats = find(I_v2_area_repeats == I_v2_area_repeats_unique(i));
    % Make a new array with all for index and distance
    D_current = D_v2_repeats(:, current_repeats);
    I_area_current = I_v2_area_repeats_all(:, current_repeats);
    I_coord_current = I_v2_coord_repeats(current_repeats);
    % Find the smallest distance in the first row i.e. the closest reef
    [~, ind_min] = min(D_current(1, :));
    [~, ind_max] = max(D_current(1, :));
    % Update the index so that the closest reef is used 
    I_v2_area_actual(I_coord_current(ind_min)) = I_area_current(1, ind_min);
    % Figure out the closest reef for the other repeats
    % If there is only one other reef
    if size(D_current, 2) == 2
        % Make sure the second closest reef isn't already used
        if ~any(I_v2_area_actual == I_area_current(2, ind_max))
            I_v2_area_actual(I_coord_current(ind_max)) = I_area_current(2, ind_max);
        % Make sure the third closest reef isn't already used
        elseif ~any(I_v2_area_actual == I_area_current(3, ind_max))
            I_v2_area_actual(I_coord_current(ind_max)) = I_area_current(3, ind_max);
        % Otherwise place a zero and delete the entry later
        else 
            I_v2_area_actual(I_coord_current(ind_max)) = 0;
        end
    % If there are two other reefs with the same 
    elseif size(D_current, 2) == 3
        % Get the subset of indexes and values
        D_subset = D_current(2:3, 1:end ~= ind_min);
        I_area_subset = I_area_current(2:3, 1:end ~= ind_min);
        I_coord_subset = I_coord_current(1:end ~= ind_min);
        % Check that the second closest reefs aren't the same
        if I_area_subset(1, 1) ~= I_area_subset(1, 2)
            % Make sure the first reef isn't already used
            if ~any(I_v2_area_actual == I_area_subset(1, 1))
                I_v2_area_actual(I_coord_current(1)) = I_area_subset(1, 1);
            else
                fprintf('help reef 1 already used')
            end
            % Make sure the second reef isn't already used
            if ~any(I_v2_area_actual == I_area_subset(1, 2))
                I_v2_area_actual(I_coord_current(2)) = I_area_subset(1, 2);
            else
                fprintf('help reef 2 already used')
            end
        else
            fprintf('help same reefs again')
        end
    end
end

% =========================================================================
% Delete any entries from indexes that are zero - these are repeat reefs
index_delete = find(I_v2_area_actual == 0); 
I_v2_area_actual(index_delete) = [];

% Delete the corresponding data from everything else i.e. connectivity 
% matrix, lat/long data
reef_coord_data_v2(index_delete, :) = []; 
omega(index_delete, :) = [];
omega(:, index_delete) = [];

% Keep only the reefs that match
reef_area_data_actual = reef_area_data(I_v2_area_actual, :);

% =========================================================================


% RETURN VARIABLES --------------------------------------------------------
% =========================================================================
% Variables to return
num_reefs = length(reef_area_data_actual);
lon = reef_coord_data_v2(:, 1);
lat = reef_coord_data_v2(:, 2);
reef_area = reef_area_data_actual(:, 3);
omega = omega;

% =========================================================================
