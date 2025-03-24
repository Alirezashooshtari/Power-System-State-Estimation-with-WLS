%% Plot Voltage Magnitudes and Angles from WLS State Estimation
% This script generates two subplots to compare the estimated voltage magnitudes
% and angles from WLS against their true values at a specific time sample.

%% Clear Existing Figures
close all; % Close all open figure windows to start fresh

%% Define Parameters
sample_time = 30; % Specific time sample to plot (adjust as needed)
B_N = size(V_WLS_noisy, 1); % Number of buses (assumed from V_WLS_noisy dimensions)

%% Create Subplot for Voltage Magnitudes
subplot(1, 2, 1); % 1 row, 2 columns, first subplot
v_se = V_WLS_noisy(:, sample_time);   % Estimated voltage magnitudes at sample_time (pu)
v_true = V_true(:, sample_time);      % True voltage magnitudes at sample_time (pu)
t = linspace(1, B_N, B_N);            % X-axis: Bus indices from 1 to B_N

% Plot estimated and true voltage magnitudes
plot(t, v_se', 'g--*', ...            % Estimated values: Green dashed line with stars
     t, v_true', 'b--o', ...          % True values: Blue dashed line with circles
     'LineWidth', 1);                 % Set line width for clarity
legend('SE-WLS', 'True-value', ...    % Legend labels
       'Location', 'NorthWest');      % Position legend in northwest corner
ylabel('Voltage Magnitudes (p.u.)');  % Y-axis label with units (per-unit)
xlabel('Bus Numbers');                % X-axis label
grid on;                              % Enable grid lines for readability

%% Create Subplot for Voltage Angles
subplot(1, 2, 2); % 1 row, 2 columns, second subplot
theta_se = Theta_WLS_noisy(:, sample_time); % Estimated voltage angles at sample_time (radians)
theta_true = Theta_true(:, sample_time);    % True voltage angles at sample_time (radians)
t = linspace(1, B_N, B_N);                  % X-axis: Bus indices from 1 to B_N

% Plot estimated and true voltage angles
plot(t, theta_se', 'g--*', ...        % Estimated values: Green dashed line with stars
     t, theta_true', 'b--o', ...      % True values: Blue dashed line with circles
     'LineWidth', 1);                 % Set line width for clarity
legend('SE-WLS', 'True-value', ...    % Legend labels
       'Location', 'NorthWest');      % Position legend in northwest corner
ylabel('Voltage Angles (Radian)');    % Y-axis label with units (radians)
xlabel('Bus Numbers');                % X-axis label
grid on;                              % Enable grid lines for readability