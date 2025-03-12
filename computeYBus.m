function ybus = computeYBus(casestudy)
    % computeYBus: Computes the bus admittance matrix (Ybus) for a given MATPOWER case.
    %
    % This function loads a MATPOWER case file and calculates the bus admittance matrix (Ybus).
    % The function ensures the input file is valid and returns Ybus as a full matrix.
    %
    % Input:
    % - casestudy: MATPOWER case file (e.g., 'case14')
    %
    % Output:
    % - ybus: Full bus admittance matrix (Ybus) in complex form.
    %


    % Load the MATPOWER case file
    try
        mpc = loadcase(casestudy);
    catch
        error('Invalid MATPOWER case file. Ensure the file exists and is properly formatted.');
    end

    % Ensure required data fields exist
    if ~isfield(mpc, 'bus') || ~isfield(mpc, 'branch')
        error('MATPOWER case file must contain both bus and branch data.');
    end

    % Compute the bus admittance matrix (Ybus) and branch admittances (Yf, Yt)
    [Ybus, ~, ~] = makeYbus(mpc);

    % Convert Ybus to a full matrix for better readability
    ybus = full(Ybus);

end
