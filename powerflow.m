% This has been used to calculate the true measurements

function h = powerflow(casestudy, V, del, voltage_buses, power_buses, power_branches)
    % wls_z: Computes the true measurements for the WLS estimation
    %
    % Inputs:
    % - casestudy: MATPOWER case file
    % - V: Bus voltage magnitudes
    % - del: Bus voltage angles (radians)
    % - voltage_buses: List of buses for voltage magnitude measurements
    % - power_buses: List of buses for power injection measurements
    % - power_branches: List of branches for power flow measurements
    %
    % Output:
    % - h: True measurement values for the selected measurements

    % Load system data
    ybus = computeYBus(casestudy); % Bus admittance matrix (YBus)
    zdata = extractZData(casestudy, voltage_buses, power_buses, power_branches); % Measurement data
    bpq = computeBBus(casestudy); % B-bus data

    % Extract system properties
    nbus = max(max(zdata(:, 4)), max(zdata(:, 5))); % Determine number of buses
    type = zdata(:, 2); % Measurement types
    fbus = zdata(:, 4); % From bus indices
    tbus = zdata(:, 5); % To bus indices

    % Extract real and imaginary parts of YBus
    G = real(ybus);
    B = imag(ybus);

    % Identify measurement types
    vi = find(type == 1); % Voltage magnitude measurements
    ppi = find(type == 2); % Real power injections
    qi = find(type == 3); % Reactive power injections
    pf = find(type == 4); % Real power flows
    qf = find(type == 5); % Reactive power flows

    % Get number of measurements per type
    nvi = length(vi);
    npi = length(ppi);
    nqi = length(qi);
    npf = length(pf);
    nqf = length(qf);

    % Preallocate measurement function h
    h1 = V(fbus(vi), 1); % Voltage magnitudes
    h2 = zeros(npi, 1);  % Active power injections
    h3 = zeros(nqi, 1);  % Reactive power injections
    h4 = zeros(npf, 1);  % Active power flows
    h5 = zeros(nqf, 1);  % Reactive power flows

    %% **1. Compute Power Injection Measurements
    if npi > 0
        for i = 1:npi
            m = fbus(ppi(i));
            h2(i) = V(m) * sum(V(:) .* (G(m, :)' .* cos(del(m) - del(:)) + B(m, :)' .* sin(del(m) - del(:))));
        end
    end

    if nqi > 0
        for i = 1:nqi
            m = fbus(qi(i));
            h3(i) = V(m) * sum(V(:) .* (G(m, :)' .* sin(del(m) - del(:)) - B(m, :)' .* cos(del(m) - del(:))));
        end
    end

    %% **2. Compute Power Flow Measurements
    if npf > 0
        for i = 1:npf
            m = fbus(pf(i));
            n = tbus(pf(i));
            h4(i) = -V(m)^2 * G(m, n) + V(m) * V(n) * (G(m, n) * cos(del(m) - del(n)) + B(m, n) * sin(del(m) - del(n)));
        end
    end

    if nqf > 0
        for i = 1:nqf
            m = fbus(qf(i));
            n = tbus(qf(i));
            h5(i) = V(m)^2 * (B(m, n) - bpq(m, n)) + V(m) * V(n) * (G(m, n) * sin(del(m) - del(n)) - B(m, n) * cos(del(m) - del(n)));
        end
    end

    % Combine all measurement values
    h = [h1; h2; h3; h4; h5];

end
