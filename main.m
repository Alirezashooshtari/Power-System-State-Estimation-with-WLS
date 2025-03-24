clc;
clearvars;
close all;
format short g;

disp(' ');
disp(' ============== START AC SE ================ ');
disp(' =========================================== ');
disp(' ');
rng(42); % Set random seed for reproducibility

%% ========================== Load Profile & Case Setup ==========================
load_profile = csvread("load.csv"); % Simulated load profile
T = 100;%numel(load_profile); % Number of time steps
casestudy = 'case14'; % MATPOWER case

% Load MATPOWER case for IEEE 14-bus system
mpc = loadcase(casestudy); 
mpopt = mpoption('out.all', 0, 'verbose', 0); % MATPOWER settings
unsuccess_count = 0; % Counter for unsuccessful OPFs
slack_bus = find(mpc.bus(:, 2)==3);
sb_angle = mpc.bus(slack_bus, 9) * pi/180;
sb_volatge = mpc.bus(slack_bus, 8);


%% ========================== Define Measurement Locations ==========================
% Full redundancy for IEEE 14-bus system
voltage_buses = 1:numel(mpc.bus(:,1));
power_buses = setdiff(1:numel(mpc.bus(:,1)),slack_bus);
power_branches = 1:numel(mpc.branch(:,1));

%power injection full redundancy for IEEE 14
% voltage_buses = 1:1;
% power_buses = setdiff(1:14,1);
% power_branches = setdiff(1:20,[20,18,17,16, 14, 12, 6, 2]);

%less redundancy for IEEE 14
% voltage_buses = 1:1;
% power_buses = setdiff(1:14,[1,4,5,6,9,13]);
% power_branches = setdiff(1:20,[20,18,17,16, 14, 12, 6, 2]);

% Extract measurement data
z = extractZData(casestudy, voltage_buses, power_buses, power_branches);

%% ========================== Preallocate Memory ==========================
B_N = str2double(regexp(casestudy, '\d+', 'match')); % Number of buses
num_meas = numel(z(:,1)); % Total number of measurements
Var_N = 2 * B_N - 1; % Number of state variables (excluding slack bus)

% Initialize matrices for storing results
V_WLS_noisy = zeros(B_N, T);
Theta_WLS_noisy = zeros(B_N, T);
meas_T_true = zeros(num_meas, T);
V_true = zeros(B_N, T);
Theta_true = zeros(B_N, T);

%% ========================== Time Loop ==========================
disp('Running Power Flow Simulation...');

for i = 1:T
    % Create a temporary copy of the case to modify load values
    mpc_temp = mpc; 
    mpc_temp.bus(:,3) = mpc_temp.bus(:,3) * load_profile(i);
    mpc_temp.bus(:,4) = mpc_temp.bus(:,4) * load_profile(i);

    % Run AC OPF
    res = runopf(mpc_temp, mpopt);
    
    % Track unsuccessful OPFs
    if res.success == 0
        unsuccess_count = unsuccess_count + 1;
    end

    % Display progress every 100 iterations
    if mod(i, 100) == 0
        fprintf('AC OPF at time step: %d | Unsuccessful OPFs: %d\n', i, unsuccess_count);
    end

    % Extract true voltage magnitudes and angles
    V_true(:, i) = res.bus(:,8);
    Theta_true(:, i) = res.bus(:,9) * pi / 180;

    % Compute true measurements using the power flow function
    meas_T_true(:, i) = powerflow(casestudy, V_true(:, i), Theta_true(:, i), voltage_buses, power_buses, power_branches);
end


%% ========================== Add Noise to Measurements ==========================
disp('Adding noise to measurements...');

accRT = sqrt(z(:, 6)); % Still relative error scaling
noise = accRT .* meas_T_true .* randn(num_meas, T); % Additive, scaled by magnitude
meas_T_noisy = meas_T_true + noise;
sig = accRT .* meas_T_true; % Same sig

%% ========================== Preallocate WLS Storage ==========================
Residuals_WLS_noisy = zeros(num_meas, T);
L2_norm = zeros(1, T);
Norm_Residuals_WLS_noisy = zeros(num_meas, T);
Obj_valu_WLS_noisy = zeros(num_meas, T);

%% ========================== Weighted Least Squares (WLS) ==========================
disp(' ');
disp(' ============== Running WLS Estimation ================ ');
disp(' ');

tic;
for i = 1:T
    %%% Check for pseudo-measurements or SLC  
    % When load is zero, `sigma = 0` → `sigma inverse would be inf`
    % This causes singularity, so we need to replace zeros with a small value
    % if any(abs(sig(:, i)) < 1e-8)
    %     sig(abs(sig(:, i)) < 1e-8, i) = 1e-4;
    % end
    if sum(abs(sig(:,i)) < 1e-8)~=0
        sig(abs(sig(:,i))< 1e-8,i) = 10^-4;
    end

    % Display progress every 100 iterations
    if mod(i, 100) == 0
        fprintf('WLS estimation at time step: %d\n', i);
    end

    % Weight matrix W = diag(1/sigma^2)
    W = diag(sig(:, i) .^ -2);

    % Run WLS state estimation
    [J, r, V, del] = wls(casestudy, meas_T_noisy(:, i), sig(:, i), voltage_buses, power_buses, power_branches, slack_bus, sb_angle ,sb_volatge);

    % Store estimated values
    V_WLS_noisy(:, i) = V;
    Theta_WLS_noisy(:, i) = del;
    Residuals_WLS_noisy(:, i) = r;
    Norm_Residuals_WLS_noisy(:, i) = r ./ sig(:, i);
    L2_norm(i) = norm(r, 2); % L2 norm of residuals
    Obj_valu_WLS_noisy(:, i) = W * (r .^ 2);
    Jbdd = sum(Obj_valu_WLS_noisy);
end
elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);


