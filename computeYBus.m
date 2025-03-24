function Ybus = computeYBus(casestudy)
% Function to compute the bus admittance matrix (Ybus) for a given power system case.
% Input: 
%   caseNumber - Identifier for the power system case (used to load line data)
% Output:
%   Ybus - Bus admittance matrix (nbus x nbus complex matrix)

% Load line data using a separate function or file (e.g., linedata6.m)
lineData = extractLineData(casestudy);  

% Extracting line parameters from line data
fromBus = lineData(:,1);       % From bus number
toBus = lineData(:,2);         % To bus number
resistance = lineData(:,3);    % Line resistance R (in p.u.)
reactance = lineData(:,4);     % Line reactance X (in p.u.)
halfLineCharging = lineData(:,5); % Half-line charging susceptance B/2 (in p.u.)
tapRatio = lineData(:,6);      % Transformer tap ratio (if any), otherwise 1

% Calculate complex impedance and admittance
impedance = resistance + 1i * reactance;     % Complex line impedance Z = R + jX
admittance = 1 ./ impedance;                 % Line admittance Y = 1/Z
lineCharging = 1i * halfLineCharging;        % Line charging susceptance (imaginary)

% Determine system size
numBuses = max(max(fromBus), max(toBus));    % Total number of buses
numBranches = length(fromBus);               % Total number of branches

% Initialize the bus admittance matrix
Ybus = zeros(numBuses, numBuses);            % Complex Ybus matrix

%% Form off-diagonal elements of Ybus
for k = 1:numBranches
    Ybus(fromBus(k), toBus(k)) = Ybus(fromBus(k), toBus(k)) - admittance(k) / tapRatio(k);
    Ybus(toBus(k), fromBus(k)) = Ybus(fromBus(k), toBus(k));  % Ybus is symmetric
end

%% Form diagonal elements of Ybus
for m = 1:numBuses
    for n = 1:numBranches
        if fromBus(n) == m
            % If the current bus is the "from" bus
            Ybus(m, m) = Ybus(m, m) + admittance(n) / (tapRatio(n)^2) + lineCharging(n);
        elseif toBus(n) == m
            % If the current bus is the "to" bus
            Ybus(m, m) = Ybus(m, m) + admittance(n) + lineCharging(n);
        end
    end
end
end