%% Returns total space charge in the semi-infinite bulk given the
% band-bending potential drop dphi_b and bulk parameters.
function [Qsc, numericSol, degenerateCond] = findQsc(dphi_b, ...
    T, Ndop, Nv, Nc, epsilon, Eg, ...
    cc_Qsheet_minimal, degeneracyLimit_kT, ...
    ignoreNumericFail, outputFile)

% Get surface charge using analytical method (even if the numeric solution
% is ultimately used the analytical result forms the initial guess)
Qsc_analytical = solvePoisson_nondegenerate(dphi_b, T, Ndop, Nv, Nc, epsilon, Eg);

% Handle zero situation
if abs(Qsc_analytical) < cc_Qsheet_minimal
    Qsc = 0;
    numericSol = false;
    degenerateCond = false;
    return
end
% If the user desires never to use the numeric solution, take analytical
if degeneracyLimit_kT == -Inf
    Qsc = Qsc_analytical;
    numericSol = false;
    degenerateCond = false;
    return
end

% Constants
physConsts = initPhysicalConsts;
q = physConsts.q;     % C/e
e = physConsts.e;     % e
k = physConsts.k;     % J/K

kT = k*T / q;     % J / C/e = CV / C/e = eV

% Bulk and surface electronic structure, all referenced to bulk Ev (VBM)
Ev_bulk = 0;
Ec_bulk = Eg;
Ev_surf = Ev_bulk - e * dphi_b;
Ec_surf = Ec_bulk - e * dphi_b;

% Find Ef as determined by neutrality in the bulk
Ef = findEMBulkEf(T, Ndop, Nv, Nc, Ev_bulk, Ec_bulk);

% Distance in kT of Ef from either band edge in the bulk and surface
relEfPositions = [Ef - [Ev_bulk; Ev_surf], [Ec_bulk; Ec_surf] - Ef]./ kT;

degenerateCond = any(any(relEfPositions < degeneracyLimit_kT));
numericSol = degenerateCond;
if degenerateCond
    % Degenerate doping conditions found at some extreme of the system.
    % Apply numerical solution.
    outputMessage(sprintf(['Identified degenerate conditions for delta_phi_b = %0.10f V. ' ...
        'Solving Poisson-Boltzmann equation numerically'], dphi_b), outputFile);
    
    [Qsc, ~, ~, ~, ~, ~, errmsg] = solvePoisson_numeric(dphi_b, T, ...
        Ndop, Nv, Nc, epsilon, Eg, degeneracyLimit_kT, Qsc_analytical);
    if ~strcmp(errmsg, '')
        % Numeric solution failed. Collect information about energies
        numericFailureMsg = sprintf(['Numeric Poisson solver may have failed. ' ...
            'System info:\n' ...
            'delta_phi_b = %0.10f V\n' ...
            'Ev_bulk - Ef = %0.10f kT\n' ...
            'Ec_bulk - Ef = %0.10f kT\n' ...
            'Ev_surf - Ef = %0.10f kT\n' ...
            'Ec_surf - Ef = %0.10f kT\n' ...
            'Error/warning message:\n%s'], ...
            dphi_b, (Ev_bulk - Ef) ./ kT, (Ec_bulk - Ef) ./ kT, ...
            (Ev_surf - Ef) ./ kT, (Ec_surf - Ef) ./ kT, errmsg);
        if ignoreNumericFail
            % Take analytical solution even though it may be wrong
            outputMessage(['Warning: ' numericFailureMsg], outputFile);
            outputMessage('Using analytical solution. Take care with results', outputFile);
            Qsc = Qsc_analytical;
            numericSol = false;
        else
            % Error
            err = MException('MATLAB:findQsurf:NumericPBSolverFailed', ...
                numericFailureMsg);
            throw(err)
        end
    end
else
    % Nondegenerate conditions. Trust analytical solution
    Qsc = Qsc_analytical;
end

end