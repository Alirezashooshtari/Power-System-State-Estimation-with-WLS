%% Chi-Squared Test for Bad Data Detection
% This script performs a chi-squared test to detect bad data in power system
% state estimation by comparing the objective function J_BDD against a
% critical threshold.

%% Time Vector Definition
% Define the time samples for analysis
T = 100; % Total number of time samples (replace with actual T value)
time_samples = 1:T; % Vector of time indices from 1 to T

%% Parameter Initialization
% Define key parameters for the chi-squared test
m = num_meas;           % Total number of measurements (e.g., from extractZData)
N = Var_N;              % Number of estimated parameters (e.g., state variables)
p = 0.01;               % Significance level (0.01 for 99% confidence, adjust as needed)

% Calculate degrees of freedom
df = m - N;             % Degrees of freedom: measurements minus parameters

%% Chi-Squared Critical Value Calculation
% Compute the critical value from the chi-squared distribution
chi_squared_critical = chi2inv(1 - p, df); % Threshold for bad data detection
% Note: 1 - p gives the upper tail critical value (e.g., 99th percentile)

%% Plotting Chi-Squared Test Results
% Create a figure to visualize J_BDD against the chi-squared threshold
figure; % Open a new figure window

% Plot J_BDD over time
plot(time_samples(1:T), Jbdd(1:T), 'b', 'LineWidth', 1, ...
    'DisplayName', 'J_{BDD}');
% Arguments:
% - time_samples(1:T): X-axis (time indices)
% - Jbdd(1:T): Y-axis (objective function values, assumed precomputed)
% - 'b': Blue solid line
% - 'LineWidth', 1: Line thickness

hold on; % Allow multiple plots on the same figure

% Plot the chi-squared critical threshold as a horizontal line
yline(chi_squared_critical, '--r', 'LineWidth', 1, ...
    'DisplayName', 'Threshold');
% Arguments:
% - chi_squared_critical: Y-value of the threshold
% - '--r': Red dashed line
% - 'LineWidth', 1: Line thickness

%% Figure Customization
% Add labels and formatting
xlabel('Time Sample');         % X-axis label
ylabel('J_{BDD}');             % Y-axis label (BDD = Bad Data Detection)
% ylim([0, 150]);              % Uncomment to set Y-axis limits (adjust as needed)
% title('χ²-test');            % Uncomment to add a title
legend('show');                % Display legend with specified DisplayNames
grid on;                       % Enable grid lines for better readability

hold off; % Release the hold on the figure