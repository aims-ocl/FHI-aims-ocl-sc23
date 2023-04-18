%% Gives a new guess for Qsheet based on history. Based on a combination
% of the secant method and the double false position method, where the
% preferred choice is the secant result unless it is already known to be
% outside the limits bounding the position of the root, as discovered in
% previous iterations.
function next_Qsheet = mixQsheet(Qsheet_history, output_Qsheet_history, mixingParameter)

% Need at least 2 iterations to make an educated guess at Qsheet
if numel(Qsheet_history) < 2 || mixingParameter < 0
    % Mix next Qsheet simply
    next_Qsheet =  abs(mixingParameter) * output_Qsheet_history(end) + ...
        (1 - abs(mixingParameter)) * Qsheet_history(end);
    return
end

% Determine bounds for next guess
Qin_bounds = [-Inf; Inf];
% Find differences between input and output Qsheets
dQ_all = Qsheet_history - output_Qsheet_history;
% Sort dQ and output according to input Qsheets
[Qsheet_history_sorted, i] = sort(Qsheet_history);
dQ_sorted = dQ_all(i);
output_Qsheet_history_sorted = output_Qsheet_history(i);
% Find switches from negative to positive and vice versa
dQ_switches = diff(dQ_sorted > 0) ~= 0;
% Count switches from positive to negative
num_shifts = sum(dQ_switches);
% If one switch, narrow bounds for next guess
if num_shifts == 1
    % Only one crossing of (Qsheet, dQ). Find values closest to 0
    switch_index = find(dQ_switches ~= 0, 1);
    % Use as bounds
    Qin_bounds = Qsheet_history_sorted(switch_index : switch_index+1);
end

% Use secant method on last 2 iterations
Qin = Qsheet_history(end-1 : end);
Qout = output_Qsheet_history(end-1 : end);
intercept_Qin = getSecantGuess(Qin, Qout);

% Check that guess is within bounds (uses assumption of only one crossing)
if intercept_Qin < min(Qin_bounds) || intercept_Qin > max(Qin_bounds)
    % Simple secant guess is off. Use bounds to get next guess instead
    Qin = Qsheet_history_sorted(switch_index : switch_index+1);
    Qout = output_Qsheet_history_sorted(switch_index : switch_index+1);
    intercept_Qin = getSecantGuess(Qin, Qout);
end

% Mix next Qsheet simply
next_Qsheet =  mixingParameter * intercept_Qin + ...
    (1 - mixingParameter) * Qsheet_history(end);

end

%% Local function to calculate secant guess
function intercept_Qin = getSecantGuess(Qin, Qout)

% Find differences between input and output Qsheets
dQ = Qout - Qin;
% Extrapolate dQ(Qin) to dQ = 0
slope = diff(dQ) ./ diff(Qin);
dQ0 = dQ(1) - slope * Qin(1);
intercept_Qin = -dQ0 / slope;

end

