function zdt = extractZData(casestudy, voltage_buses, power_buses, power_branches)
    % extractZData: Extracts selected measurement data from a MATPOWER case file.
    %
    % This function allows the user to select specific buses for voltage, power injections,
    % and branches for power flows.
    %
    % Inputs:
    % - casestudy: MATPOWER case file (e.g., 'case14')
    % - voltage_buses: List of bus indices to include for voltage magnitude measurements
    % - power_buses: List of bus indices to include for power injections
    % - power_branches: List of branch indices to include for power flow measurements
    %
    % Output:
    % - zdt: Measurement data matrix in the format:
    %   [Measurement ID | Type | Value | From Bus | To Bus | Rii]
    %

    % Load the MATPOWER case file
    try
        mpc = loadcase(casestudy);
    catch
        error('Invalid MATPOWER case file. Ensure the file exists and is properly formatted.');
    end

    % Extract data sizes
    n_buses = size(mpc.bus, 1);       % Total number of buses
    n_branches = size(mpc.branch, 1); % Total number of branches

    % Validate Inputs
    voltage_buses = unique(voltage_buses); % Ensure unique bus indices
    power_buses = unique(power_buses);
    power_branches = unique(power_branches);

    % Ensure selected buses and branches exist
    if any(voltage_buses > n_buses) || any(power_buses > n_buses)
        error('Selected voltage or power buses exceed total available buses.');
    end
    if any(power_branches > n_branches)
        error('Selected power branches exceed total available branches.');
    end

    % Define total measurements count
    total_measurements = length(voltage_buses) + 2 * length(power_buses) + 2 * length(power_branches);

    % Preallocate zdt for efficiency
    zdt = zeros(total_measurements, 6); 

    % Define measurement IDs
    row_idx = (1:total_measurements)';

    %% 1. Voltage Magnitudes (Type 1) - Custom Selection
    num_voltage = length(voltage_buses);
    zdt(1:num_voltage, :) = [row_idx(1:num_voltage), ...
                             ones(num_voltage, 1), ...
                             ones(num_voltage, 1), ...
                             voltage_buses(:), ...
                             zeros(num_voltage, 1), ...
                             repmat(4e-6, num_voltage, 1)];

    %% 2. Active Power Injections (Type 2) - Custom Selection
    row_start = num_voltage + 1;
    row_end = row_start + length(power_buses) - 1;
    zdt(row_start:row_end, :) = [row_idx(row_start:row_end), ...
                                 repmat(2, length(power_buses), 1), ...
                                 ones(length(power_buses), 1), ...
                                 power_buses(:), ...
                                 zeros(length(power_buses), 1), ...
                                 repmat(1e-4, length(power_buses), 1)];

    %% 3. Reactive Power Injections (Type 3) - Custom Selection
    row_start = row_end + 1;
    row_end = row_start + length(power_buses) - 1;
    zdt(row_start:row_end, :) = [row_idx(row_start:row_end), ...
                                 repmat(3, length(power_buses), 1), ...
                                 ones(length(power_buses), 1), ...
                                 power_buses(:), ...
                                 zeros(length(power_buses), 1), ...
                                 repmat(1e-4, length(power_buses), 1)];

    %% 4. Real Power Flows (Type 4) - Custom Selection
    row_start = row_end + 1;
    row_end = row_start + length(power_branches) - 1;
    zdt(row_start:row_end, :) = [row_idx(row_start:row_end), ...
                                 repmat(4, length(power_branches), 1), ...
                                 ones(length(power_branches), 1), ...
                                 mpc.branch(power_branches, 1), ...
                                 mpc.branch(power_branches, 2), ...
                                 repmat(2.5e-6, length(power_branches), 1)];

    %% 5. Reactive Power Flows (Type 5) - Custom Selection
    row_start = row_end + 1;
    row_end = row_start + length(power_branches) - 1;
    zdt(row_start:row_end, :) = [row_idx(row_start:row_end), ...
                                 repmat(5, length(power_branches), 1), ...
                                 ones(length(power_branches), 1), ...
                                 mpc.branch(power_branches, 1), ...
                                 mpc.branch(power_branches, 2), ...
                                 repmat(2.5e-6, length(power_branches), 1)];

end