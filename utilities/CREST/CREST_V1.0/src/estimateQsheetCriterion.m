%% Estimate the required convergence criterion for the sheet charge from
% the convergence criterion on the band bending. Calculates Qsc for a range
% of dphi at distances = conv_crit_dphi, and takes minimal difference.
function cc_Qsheet_integrated = estimateQsheetCriterion(cc_dphi_b, cc_Ef, cc_Qsheet, ...
    dphi_b_bounds, Ef_bounds, T, Ndop, Nv, Nc, epsilon, Eg, ...
    useDOS, enLev_DOS, E_s, surfaceCellArea)

% CONSTANTS
physConsts = initPhysicalConsts;
k = physConsts.k;
q = physConsts.q;
e = physConsts.e;

ccs_Qsheet = [Inf Inf cc_Qsheet];
for i = [1 2]
    if i == 1
        % Bulk. Take smaller cc of cc_dphi_b, cc_Ef since both affect it
        cc_energy = min(e*cc_dphi_b, cc_Ef);
        % Unpack
        LB = min(e*dphi_b_bounds);
        UB = max(e*dphi_b_bounds);
    else
        % Slab. Can only calculate this if DOS was given.
        if ~useDOS; continue; end
        cc_energy = cc_Ef;
        % Unpack
        LB = min(Ef_bounds);
        UB = max(Ef_bounds);
        energyLevels = enLev_DOS(:, 1);
        DOS = enLev_DOS(:, 2);
    end
    
    if cc_energy == Inf; continue; end
    
    % Use an even spread of energy values which are the conv. crit. apart
    paramValues = (LB : cc_energy : UB)';
    
    % Calculate dQ for all values of the energy
    if i == 1
        % Calculate dQbulk = -Qsc values based on nondegenerate
        dQ = -solvePoisson_nondegenerate(paramValues, T, Ndop, Nv, Nc, epsilon, Eg);
    else
        % Calculate dQslab from integration of DOS
        dQ = zeros(size(paramValues));
        for ii = 1 : numel(paramValues)
            dQ(ii) = -sumStatesInInterval(q, k , T, energyLevels, DOS, [E_s paramValues(ii)]) / surfaceCellArea;     % e/cm^2
        end
    end
    
    % Take minimum difference between points
    minDiff = min(abs(diff(dQ)));
    
    % For too-high values of the conv. crit. the result may be []. Fix
    if isempty(minDiff); continue; end
    
    ccs_Qsheet(i) = minDiff;
end

% Return the smaller of the two derived convergence criteria
cc_Qsheet_integrated = min(ccs_Qsheet);

end

