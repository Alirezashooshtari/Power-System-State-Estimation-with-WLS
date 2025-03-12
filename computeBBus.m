function bbus = computeBBus(casestudy)
    % computeBBus: Computes the B-bus (Shunt Admittance Matrix) from a MATPOWER case file.
    %
    % This function extracts line data from the MATPOWER case file and constructs the
    % B-bus matrix using shunt admittance (B/2 values).
    %
    % Input:
    % - casestudy: MATPOWER case file (e.g., 'case14')
    %
    % Output:
    % - bbus: Shunt admittance matrix (B-bus).
    %

    % Extract line data using the optimized linedatas function
    try
        linedata = extractLineData(casestudy);  % Ensure you use the optimized function
    catch
        error('Error in extracting line data. Ensure the MATPOWER case file is valid.');
    end

    % Extract required data columns
    fb = linedata(:, 1);  % From Bus
    tb = linedata(:, 2);  % To Bus
    b = linedata(:, 5);   % B/2 (Shunt Susceptance)

    % Determine number of buses and branches
    nbus = max(max(fb), max(tb));  % Total number of buses
    nbranch = length(fb);          % Total number of branches

    % Initialize B-bus matrix
    bbus = zeros(nbus, nbus);

    % Construct B-bus matrix using branch susceptance values
    for k = 1:nbranch
        bbus(fb(k), tb(k)) = b(k);
        bbus(tb(k), fb(k)) = b(k); % Ensure symmetry
    end
end
