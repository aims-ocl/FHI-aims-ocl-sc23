%% Given an initial interval for the Fermi level and a slab DOS profile,
% find the equilibrium value for Ef and the charging of the slab in light
% of band-bending inside the bulk.
% Input:
% E_s - Fermi level of the slab if no charge transfer
% E_b - Fermi level of the bulk if no charge transfer
% T - temperature (K)
% Ndop - Per-volume density of dopant-derived fixed charge in the bulk (e/cm^3)
% Eg - Gap energy in the bulk (eV)
% convergence_criterion_Qsheet - Convergence criterion for the sheet
% charge (e/cm^2)
% Output:
% equil_deltaQ - Charge transfer where charge conservation holds
% equil_dphi_b - Bulk band bending where charge conservation holds
% equil_Ef - Ef where charge conservation holds
function [equil_deltaQ, equil_dphi_b, equil_Ef] = calcCTBinary(E_s, E_b, ... 
        T, Ndop, Nv, Nc, epsilon, Eg, ...
        energyLevels, DOS, surfaceCellArea, ...     % cm^2
        cc_Qsheet_minimal, cc_dphi_b, cc_Ef, cc_Qsheet, ...     % e/cm^2,  e/cm^2, eV
        degeneracyLimit_kT, ignoreNumericFail, ...
        outputFile)

% CONSTANTS
physConsts = initPhysicalConsts;
k = physConsts.k;
q = physConsts.q;
e = physConsts.e;

% BINARY SEARCH
% Search for a band bending potential drop where Ef where the charge from the slab DOS matches the analytically
% calculated surface charge
Ef_interval = [E_s, E_b];

equilibrated = false;
while ~equilibrated
    Ef_guess = mean(Ef_interval);
    
    % Find slab charging caused by populating DOS from E_s to Ef_guess
    dQslab_guess = -sumStatesInInterval(q, k , T, energyLevels, DOS, [E_s Ef_guess]) / surfaceCellArea;     % e/cm^2
    if abs(dQslab_guess) < cc_Qsheet_minimal; dQslab_guess = 0; end
    
    % Find B.B. potential drop Ef_guess corresponds to. Convention is that
    % a negative potential drop bends bands UP.
    dphi_b_guess = -(E_b - Ef_guess) / e;     % eV / e = V
    
    % Use EM model to find the slab charge that matches the guess band
    % bending potential drop
    [Qsc, numericSol, degenerateCond] = findQsc(dphi_b_guess, ...
        T, Ndop, Nv, Nc, epsilon, Eg, ...
        cc_Qsheet_minimal, degeneracyLimit_kT, ...
        ignoreNumericFail, outputFile);
    dQbulk_guess = -Qsc;     % e/cm^2
    if abs(dQbulk_guess) < cc_Qsheet_minimal
        dQbulk_guess = 0;
        dphi_b_guess = 0;
    end
    if sign(dQbulk_guess) * sign(dQslab_guess) == -1
        % Something is wrong - these should not have opposing signs
        err = MException('MATLAB:calcCharging:IncompatibleSign', ...
            'Calculated eff. mass and atomistic slab charges have opposing signs');
        throw(err)
    end
    
    % Find overall conv. criterion on Qsheet for the given values.
    dphi_b_bounds = dphi_b_guess + [1 -1] * min(cc_dphi_b, cc_Ef);
    Ef_bounds = Ef_guess + [1 -1] * cc_Ef;
    cc_Qsheet_integrated = estimateQsheetCriterion(cc_dphi_b, cc_Ef, cc_Qsheet, ...
        dphi_b_bounds, Ef_bounds, T, Ndop, Nv, Nc, epsilon, Eg, ...
        true, [energyLevels, DOS], E_s, surfaceCellArea);
    
    % Compare the CT from the bulk to the CT to the slab. They should be
    % equal in magnitude. If the CT to the slab is larger (in absolute value)
    % than the CT from the bulk corresponding to the band bending, then the
    % band bending is too small (and the corresponding dQbulk is too small).
    % So the Ef should be closer to the E_s limit.
    % E_s limit ==> Maximal CT from bulk, dphi_b = -(E_b - E_s)/e.
    % E_b limit ==> CT from bulk --> 0, dphi_b --> 0.
    if abs(dQslab_guess - dQbulk_guess) < cc_Qsheet_integrated
        equilibrated = true;
    else
        if abs(dQslab_guess) > abs(dQbulk_guess)
            Ef_interval = [Ef_interval(1), Ef_guess];
        else
            Ef_interval = [Ef_guess, Ef_interval(2)];
        end
    end
    
    if numericSol ~= degenerateCond
        % This can actually happen only if the numerical solution failed
        % and findQsc fell back on the analytical solution. In this case,
        % change the degeneracy limit to ensure that only the analytical
        % solution is used, and restart.
        equilibrated = false;
        degeneracyLimit_kT = -Inf;
        Ef_interval = [E_s, E_b];
    end
end

% Take equilibrium value for deltaQ from the bulk (arbitrary choice).
equil_deltaQ = dQbulk_guess;
equil_Ef = Ef_guess;
equil_dphi_b = dphi_b_guess;

end