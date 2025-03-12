function linedt = extractLineData(casestudy)
    % extractLineData: Extracts line data from a MATPOWER case file.
    % 
    % This function loads a MATPOWER case file and extracts transmission line data.
    % If tap ratios are missing or zero, they are set to 1 by default.
    % The function also calculates B/2 values for each line.
    %
    % Input:
    % - casestudy: MATPOWER case file (e.g., 'case14')
    %
    % Output:
    % - linedt: Matrix containing line data in the format:
    %   [From Bus | To Bus | R (pu) | X (pu) | B/2 (pu) | Tap Ratio]
    %
    % Example Usage:
    %   linedata = extractLineData('case14');
    
    % Load MATPOWER case file
    try
        mpc = loadcase(casestudy);
    catch
        error('Invalid MATPOWER case file. Ensure the file exists and is properly formatted.');
    end

    % Ensure the required 'branch' field exists
    if ~isfield(mpc, 'branch')
        error('MATPOWER case file does not contain branch data.');
    end

    % Extract branch data
    branchData = mpc.branch;
    
    % Ensure the branch matrix has at least 9 columns (MATPOWER format)
    numCols = size(branchData, 2);
    if numCols < 9
        % Expand the branch matrix with zeros and assign default tap values
        branchData(:, 9) = 1;  
    end

    % Ensure all zero tap ratios are set to 1
    branchData(branchData(:, 9) == 0, 9) = 1;

    % Calculate B/2 from total branch susceptance (column 5)
    B_half = branchData(:, 5) / 2;

    % Construct the output matrix with required fields
    linedt = [branchData(:, 1), ...   % From Bus
              branchData(:, 2), ...   % To Bus
              branchData(:, 3), ...   % R (pu)
              branchData(:, 4), ...   % X (pu)
              B_half, ...             % B/2 (pu)
              branchData(:, 9)];      % Tap Ratio

end
