function [J, r, V, del] = wls(casestudy, z, sig, voltage_buses, power_buses, power_branches, slack_bus, sb_angle, sb_voltage)
    % Weighted Least Squares State Estimation for Power Systems
    % Inputs:
    %   casestudy      - System data structure
    %   z             - Measurement vector
    %   sig           - Measurement standard deviations
    %   voltage_buses - Voltage measurement bus indices
    %   power_buses   - Power injection measurement bus indices
    %   power_branches- Power flow measurement branch indices
    %   slack_bus     - Slack bus index
    %   sb_angle      - Slack bus angle (in radians)
    %   sb_voltage    - Slack bus voltage magnitude
    % Outputs:
    %   J             - Final objective function value
    %   r             - Measurement residuals
    %   V             - Estimated voltage magnitudes
    %   del           - Estimated voltage angles

    %% System Data Initialization
    % Load and compute system matrices
    ybus = computeYBus(casestudy);          % Bus admittance matrix (YBus)
    zdata = extractZData(casestudy, voltage_buses, power_buses, power_branches); % Measurement data
    bpq = computeBBus(casestudy);           % B-bus data for reactive power

    % Extract real and imaginary parts of YBus
    G = real(ybus);                         % Conductance matrix
    B = imag(ybus);                         % Susceptance matrix

    %% Measurement Data Processing
    % Extract measurement properties
    nbus = max(max(zdata(:, 4)), max(zdata(:, 5))); % Total number of buses
    type = zdata(:, 2);                     % Measurement types
    fbus = zdata(:, 4);                     % From bus indices
    tbus = zdata(:, 5);                     % To bus indices

    % Identify measurement type indices
    vi = find(type == 1);                   % Voltage magnitude measurements
    ppi = find(type == 2);                  % Real power injections
    qi = find(type == 3);                   % Reactive power injections
    pf = find(type == 4);                   % Real power flows
    qf = find(type == 5);                   % Reactive power flows

    % Count number of each measurement type
    nvi = length(vi);                       % Number of voltage measurements
    npi = length(ppi);                      % Number of real power injections
    nqi = length(qi);                       % Number of reactive power injections
    npf = length(pf);                       % Number of real power flows
    nqf = length(qf);                       % Number of reactive power flows

    %% State Vector Initialization
    V = sb_voltage * ones(nbus, 1);        % Initialize voltage magnitudes
    del = sb_angle * ones(nbus, 1);        % Initialize voltage angles
    del(slack_bus) = sb_angle;              % Set slack bus angle
    angle_indices = 1:nbus;                 % All bus indices
    angle_indices(slack_bus) = [];          % Exclude slack bus from estimation
    E = [del(angle_indices); V];            % State vector: [angles; voltages]

    % Measurement error covariance
    Ri_inv = 1 ./ (sig.^2);                 % Inverse of covariance matrix

    %% Iterative Solution
    iter = 1;                               % Iteration counter
    tol = 5;                                % Initial tolerance
    max_tol = 1e-6;                         % Convergence threshold

    while tol > max_tol
        %% Measurement Function Calculation (h)
        h1 = V(fbus(vi), 1);                % Voltage magnitude measurements
        h2 = zeros(npi, 1);                 % Real power injections
        h3 = zeros(nqi, 1);                 % Reactive power injections
        h4 = zeros(npf, 1);                 % Real power flows
        h5 = zeros(nqf, 1);                 % Reactive power flows
        
        % Real power injection calculations
        for i = 1:npi
            m = fbus(ppi(i));
            for k = 1:nbus
                h2(i) = h2(i) + V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
            end
        end
    
        % Reactive power injection calculations
        for i = 1:nqi
            m = fbus(qi(i));
            for k = 1:nbus
                h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end

        % Real power flow calculations
        for i = 1:npf
            m = fbus(pf(i));
            n = tbus(pf(i));
            h4(i) = -V(m)^2*G(m,n) - V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
        end
    
        % Reactive power flow calculations
        for i = 1:nqf
            m = fbus(qf(i));
            n = tbus(qf(i));
            h5(i) = -V(m)^2*(-B(m,n)+bpq(m,n)) - V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
        end

        h = [h1; h2; h3; h4; h5];           % Complete measurement function vector
        r = z - h;                          % Measurement residuals

        %% Jacobian Matrix (H) Calculation
        % H11: Voltage magnitude w.r.t. angles (always zero)
        H11 = zeros(nvi, nbus-1);

        % H12: Voltage magnitude w.r.t. voltages
        H12 = sparse(1:nvi, fbus(vi), 1, nvi, nbus);

    
        % H21: Real power injections w.r.t. angles
        H21 = zeros(npi, nbus-1);
        if npi > 0
            m = fbus(ppi);
            for i = 1:npi
                for k = 1:nbus-1
                    bus_k = angle_indices(k); % Map k to the actual bus number
                    if bus_k == m(i)
                        temp = V(m(i)) * V(:)' .* (-G(m(i),:) .* sin(del(m(i)) - del(:)') + ...
                               B(m(i),:) .* cos(del(m(i)) - del(:)'));
                        H21(i,k) = sum(temp) - V(m(i))^2 * B(m(i),m(i));
                    else
                        H21(i,k) = V(m(i)) * V(bus_k) * (G(m(i),bus_k) * sin(del(m(i)) - del(bus_k)) - ...
                                  B(m(i),bus_k) * cos(del(m(i)) - del(bus_k)));
                    end
                end
            end
        end
        
        % H22: Real power injections w.r.t. voltages
        H22 = zeros(npi, nbus);
        if npi > 0
            m = fbus(ppi);
            for i = 1:npi
                for k = 1:nbus
                    if k == m(i)
                        temp = V(:)' .* (G(m(i),:) .* cos(del(m(i)) - del(:)') + ...
                               B(m(i),:) .* sin(del(m(i)) - del(:)'));
                        H22(i,k) = sum(temp) + V(m(i)) * G(m(i),m(i));
                    else
                        H22(i,k) = V(m(i)) * (G(m(i),k) * cos(del(m(i)) - del(k)) + ...
                                  B(m(i),k) * sin(del(m(i)) - del(k)));
                    end
                end
            end
        end
        
        % H31: Reactive power injections w.r.t. angles
        H31 = zeros(nqi, nbus-1);
        if nqi > 0
            m = fbus(qi);
            for i = 1:nqi
                for k = 1:nbus-1
                    bus_k = angle_indices(k); % Map k to the actual bus number
                    if bus_k == m(i)
                        temp = V(m(i)) * V(:)' .* (G(m(i),:) .* cos(del(m(i)) - del(:)') + ...
                               B(m(i),:) .* sin(del(m(i)) - del(:)'));
                        H31(i,k) = sum(temp) - V(m(i))^2 * G(m(i),m(i));
                    else
                        H31(i,k) = V(m(i)) * V(bus_k) * (-G(m(i),bus_k) * cos(del(m(i)) - del(bus_k)) - ...
                                  B(m(i),bus_k) * sin(del(m(i)) - del(bus_k)));
                    end
                end
            end
        end
        
        % H32: Reactive power injections w.r.t. voltages
        H32 = zeros(nqi, nbus);
        if nqi > 0
            m = fbus(qi);
            for i = 1:nqi
                for k = 1:nbus
                    if k == m(i)
                        temp = V(:)' .* (G(m(i),:) .* sin(del(m(i)) - del(:)') - ...
                               B(m(i),:) .* cos(del(m(i)) - del(:)'));
                        H32(i,k) = sum(temp) - V(m(i)) * B(m(i),m(i));
                    else
                        H32(i,k) = V(m(i)) * (G(m(i),k) * sin(del(m(i)) - del(k)) - ...
                                  B(m(i),k) * cos(del(m(i)) - del(k)));
                    end
                end
            end
        end
        
        % H41: Real power flows w.r.t. angles
        H41 = zeros(npf, nbus-1);
        if npf > 0
            m = fbus(pf);
            n = tbus(pf);
            for i = 1:npf
                for k = 1:nbus-1
                    bus_k = angle_indices(k); % Map k to the actual bus number
                    if bus_k == m(i)
                        H41(i,k) = V(m(i)) * V(n(i)) * (-G(m(i),n(i)) * sin(del(m(i)) - del(n(i))) + ...
                                  B(m(i),n(i)) * cos(del(m(i)) - del(n(i))));
                    elseif bus_k == n(i)
                        H41(i,k) = -V(m(i)) * V(n(i)) * (-G(m(i),n(i)) * sin(del(m(i)) - del(n(i))) + ...
                                   B(m(i),n(i)) * cos(del(m(i)) - del(n(i))));
                    end
                end
            end
        end
        
        % H42: Real power flows w.r.t. voltages
        H42 = zeros(npf, nbus);
        if npf > 0
            m = fbus(pf);
            n = tbus(pf);
            for i = 1:npf
                for k = 1:nbus
                    if k == m(i)
                        H42(i,k) = -V(n(i)) * (-G(m(i),n(i)) * cos(del(m(i)) - del(n(i))) - ...
                                  B(m(i),n(i)) * sin(del(m(i)) - del(n(i)))) - 2 * G(m(i),n(i)) * V(m(i));
                    elseif k == n(i)
                        H42(i,k) = -V(m(i)) * (-G(m(i),n(i)) * cos(del(m(i)) - del(n(i))) - ...
                                  B(m(i),n(i)) * sin(del(m(i)) - del(n(i))));
                    end
                end
            end
        end
        
        % H51: Reactive power flows w.r.t. angles
        H51 = zeros(nqf, nbus-1);
        if nqf > 0
            m = fbus(qf);
            n = tbus(qf);
            for i = 1:nqf
                for k = 1:nbus-1
                    bus_k = angle_indices(k); % Map k to the actual bus number
                    if bus_k == m(i)
                        H51(i,k) = -V(m(i)) * V(n(i)) * (-G(m(i),n(i)) * cos(del(m(i)) - del(n(i))) - ...
                                  B(m(i),n(i)) * sin(del(m(i)) - del(n(i))));
                    elseif bus_k == n(i)
                        H51(i,k) = V(m(i)) * V(n(i)) * (-G(m(i),n(i)) * cos(del(m(i)) - del(n(i))) - ...
                                  B(m(i),n(i)) * sin(del(m(i)) - del(n(i))));
                    end
                end
            end
        end
        
        % H52: Reactive power flows w.r.t. voltages
        H52 = zeros(nqf, nbus);
        if nqf > 0
            m = fbus(qf);
            n = tbus(qf);
            for i = 1:nqf
                for k = 1:nbus
                    if k == m(i)
                        H52(i,k) = -V(n(i)) * (-G(m(i),n(i)) * sin(del(m(i)) - del(n(i))) + ...
                                  B(m(i),n(i)) * cos(del(m(i)) - del(n(i)))) - ...
                                  2 * V(m(i)) * (-B(m(i),n(i)) + bpq(m(i),n(i)));
                    elseif k == n(i)
                        H52(i,k) = -V(m(i)) * (-G(m(i),n(i)) * sin(del(m(i)) - del(n(i))) + ...
                                  B(m(i),n(i)) * cos(del(m(i)) - del(n(i))));
                    end
                end
            end
        end

        % Assemble complete Jacobian
        H = [H11 H12; H21 H22; H31 H32; H41 H42; H51 H52];
        
        %% State Update
        Gm = H' * (Ri_inv .* H);            % Gain matrix
        J = r' * (Ri_inv .* r);             % Objective function
        dE = Gm \ (H' * (Ri_inv .* r));     % State correction
        E = E + dE;                         % Update state vector

        % Update voltage angles and magnitudes
        del(angle_indices) = E(1:nbus-1);   % Update non-slack angles
        V = E(nbus:end);                    % Update all voltages

        % Check convergence
        iter = iter + 1;
        tol = max(abs(dE));
end
end